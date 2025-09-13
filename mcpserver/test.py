import os
import sys
import json
import subprocess
from pathlib import Path

# ---- JSON-RPC over stdio helpers (MCP transport) ----

def _send_msg(proc, msg: dict):
    data = json.dumps(msg).encode("utf-8")
    header = f"Content-Length: {len(data)}\r\n\r\n".encode("utf-8")
    proc.stdin.write(header + data)
    proc.stdin.flush()

def _read_msg(proc) -> dict:
    # Read headers
    headers = {}
    while True:
        line = proc.stdout.readline()
        if not line:
            raise RuntimeError("Server closed stdout")
        line = line.decode("utf-8")
        if line in ("\r\n", "\n"):
            break
        if ":" in line:
            k, v = line.split(":", 1)
            headers[k.strip().lower()] = v.strip()
    content_length = int(headers.get("content-length", "0"))
    if content_length <= 0:
        raise RuntimeError("No Content-Length in response")
    body = proc.stdout.read(content_length)
    return json.loads(body.decode("utf-8"))

def _call(proc, method: str, params: dict, msg_id: int) -> dict:
    _send_msg(proc, {"jsonrpc": "2.0", "id": msg_id, "method": method, "params": params})
    resp = _read_msg(proc)
    if "error" in resp:
        raise RuntimeError(f"RPC error: {resp['error']}")
    return resp["result"]

def _extract_json_content(tool_result: dict):
    """
    MCP tool result shape:
    { "content": [ { "type": "json", "json": ... } ] }
    """
    content = tool_result.get("content", [])
    if not content:
        return None
    item = content[0]
    if item.get("type") == "json":
        return item.get("json")
    return None

# ---- Orchestrator to test your tools in one go ----

def run_end_to_end_workflow(
    symptoms: str | None = None,
    max_targets: int = 3,
    max_ligands: int = 3,
    include_admet: bool = True,
) -> dict:
    """
    Spawns mcpserver/main.py, then:
    - Calls diagnose_or_precautions (with your symptoms)
    - Calls analyze_targets (same symptoms)
    - For each target, calls get_ligands_for_protein (ADMET on ligands)
    Returns combined JSON dict.
    """

    # Ensure OPENROUTER_API_KEY is set if you plan to use LLM
    

    # Resolve path to mcpserver/main.py from project root
    root = Path(__file__).resolve().parents[1]  # tests/.. -> project root
    server_path = root / "mcpserver" / "main.py"
    if not server_path.exists():
        raise FileNotFoundError(f"Cannot find MCP server at {server_path}")

    # Start MCP server process
    proc = subprocess.Popen(
        [sys.executable, str(server_path)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        bufsize=0,
        cwd=str(root),  # run from project root
        env={**os.environ},  # inherit env (OPENROUTER_API_KEY, etc.)
    )

    try:
        # 0) initialize
        _call(proc, "initialize", {
            "clientInfo": {"name": "py-test", "version": "1.0"},
            "capabilities": {"tools": {}},
        }, 1)

        # 1) tools/list (optional sanity log)
        tools = _call(proc, "tools/list", {}, 2)
        tool_names = [t["name"] for t in tools.get("tools", [])]
        print("Tools available:", tool_names)

        # Ask user for symptoms if not provided
        if not symptoms:
            symptoms = input("Enter symptoms (comma/semicolon separated): ").strip()

        result = {
            "input_symptoms": symptoms,
            "diagnosis": None,
            "targets": {},
            "ligands_by_protein": {},
        }

        # 2) diagnose_or_precautions
        diag = _call(proc, "tools/call", {
            "name": "diagnose_or_precautions",
            "arguments": {"symptoms_input": symptoms, "use_llm": True},
        }, 3)
        diagnosis_json = _extract_json_content(diag)
        result["diagnosis"] = diagnosis_json
        print("\n--- diagnose_or_precautions ---")
        print(json.dumps(diagnosis_json, indent=2))

        # 3) analyze_targets
        at = _call(proc, "tools/call", {
            "name": "analyze_targets",
            "arguments": {"symptoms": symptoms, "max_targets": max_targets},
        }, 4)
        targets_json = _extract_json_content(at) or {}
        result["targets"] = targets_json
        print("\n--- analyze_targets ---")
        print(json.dumps(targets_json, indent=2))

        # Parse target_protein_1..N into a list
        targets_list = []
        for k in sorted(targets_json.keys()):
            if not k.startswith("target_protein_"):
                continue
            tp = targets_json[k]
            if not tp:
                continue
            targets_list.append(tp)

        # 4) For each target, call get_ligands_for_protein
        ligands_by_protein = {}
        msg_id = 5
        for tp in targets_list:
            pid = tp.get("protein_id")
            pname = tp.get("protein_name", "") or ""
            if not pid:
                continue
            lig = _call(proc, "tools/call", {
                "name": "get_ligands_for_protein",
                "arguments": {
                    "protein_id": pid,
                    "protein_name": pname,
                    "max_ligands": max_ligands,
                    "include_admet": include_admet,
                },
            }, msg_id)
            msg_id += 1
            lig_json = _extract_json_content(lig) or {}
            ligands_by_protein[pid] = lig_json.get("ligands", [])
            print(f"\n--- get_ligands_for_protein ({pid}) ---")
            print(json.dumps(lig_json, indent=2))

        result["ligands_by_protein"] = ligands_by_protein

        print("\n=== Combined Result ===")
        print(json.dumps(result, indent=2))
        return result

    finally:
        proc.kill()

# Run directly (no integration)
if __name__ == "__main__":
    # Optionally pass symptoms via CLI arg, else it will prompt.
    cli_symptoms = " ".join(sys.argv[1:]) if len(sys.argv) > 1 else None
    run_end_to_end_workflow(symptoms=cli_symptoms, max_targets=3, max_ligands=3, include_admet=True)