# mcpserver/main.py
import os
import sys
import json
import asyncio
from pathlib import Path
from typing import List, Dict, Any

import anyio

# Make src importable (for your RAG pipeline)
BASE_DIR = Path(__file__).resolve().parents[1]  # project root
SRC_DIR = BASE_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# RAG exact-match pipeline
from src.rag.main import run_pipeline

# Shared bio/LLM/services
from src.targetprotiensandtargetligands.shared_main import (
    OpenRouterLLMService,
    UniProtService,
    StringDBService,
    ChEMBLService,
    PubChemService,
    RDKitService,
    format_admet_for_spec,
)
from src.targetprotiensandtargetligands.schema import (_shape_target, _shape_ligand)

# Instantiate services once (shared by all tools)
llm = OpenRouterLLMService()
uniprot = UniProtService()
stringdb = StringDBService()
chembl = ChEMBLService()
pubchem = PubChemService()
rdkit = RDKitService()

class SimpleMCPServer:
    def __init__(self, name: str):
        self.name = name
        self.tools = {}

    def tool(self, name: str, description: str):
        def decorator(func):
            self.tools[name] = {
                "function": func,
                "description": description
            }
            return func
        return decorator

    async def handle_request(self, request):
        method = request.get("method")
        params = request.get("params", {})
        request_id = request.get("id")
       
        if method == "initialize":
            return {
                "id": request_id,
                "result": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {
                        "tools": {}
                    },
                    "serverInfo": {
                        "name": self.name,
                        "version": "1.0.0"
                    }
                }
            }
       
        elif method == "tools/list":
            tools_list = []
            for name, info in self.tools.items():
                tools_list.append({
                    "name": name,
                    "description": info["description"],
                    "inputSchema": {
                        "type": "object",
                        "properties": {}
                    }
                })
            return {
                "id": request_id,
                "result": {
                    "tools": tools_list
                }
            }
       
        elif method == "tools/call":
            tool_name = params.get("name")
            arguments = params.get("arguments", {})
           
            if tool_name in self.tools:
                try:
                    result = await self.tools[tool_name]["function"](**arguments)
                    return {
                        "id": request_id,
                        "result": result
                    }
                except Exception as e:
                    return {
                        "id": request_id,
                        "error": {
                            "code": -32603,
                            "message": str(e)
                        }
                    }
            else:
                return {
                    "id": request_id,
                    "error": {
                        "code": -32601,
                        "message": f"Tool {tool_name} not found"
                    }
                }
       
        return {
            "id": request_id,
            "error": {
                "code": -32601,
                "message": f"Method {method} not found"
            }
        }

    async def run(self):
        # LSP-style JSON-RPC over stdio: use Content-Length framed messages
        loop = asyncio.get_event_loop()

        def _read_message():
            # Read headers
            headers = {}
            while True:
                line = sys.stdin.buffer.readline()
                if not line:
                    return None
                line = line.decode("utf-8")
                if line in ("\r\n", "\n"):
                    break
                if ":" in line:
                    k, v = line.split(":", 1)
                    headers[k.strip().lower()] = v.strip()
            content_length = int(headers.get("content-length", "0"))
            if content_length <= 0:
                return None
            body = sys.stdin.buffer.read(content_length)
            return json.loads(body.decode("utf-8"))

        def _write_message(payload: dict):
            data = json.dumps(payload).encode("utf-8")
            header = f"Content-Length: {len(data)}\r\n\r\n".encode("utf-8")
            sys.stdout.buffer.write(header + data)
            sys.stdout.buffer.flush()

        while True:
            try:
                request = await loop.run_in_executor(None, _read_message)
                if request is None:
                    break

                response = await self.handle_request(request)
                await loop.run_in_executor(None, _write_message, response)

            except json.JSONDecodeError:
                # Malformed JSON body; report parse error
                await loop.run_in_executor(None, _write_message, {
                    "error": {"code": -32700, "message": "Parse error"}
                })
            except Exception as e:
                # Generic error
                await loop.run_in_executor(None, _write_message, {
                    "error": {"code": -32603, "message": str(e)}
                })

# Create server instance
server = SimpleMCPServer("drug_mcpserver")

@server.tool(
    "diagnose_or_precautions",
    "Given comma/semicolon-separated symptoms, return disease precautions if exact combo exists; else flag as new."
)
async def diagnose_or_precautions(symptoms_input: str = "", use_llm: bool = True):
    def _run():
        return run_pipeline(symptoms_input, use_llm=use_llm)
    data = await anyio.to_thread.run_sync(_run)
    return {"content": [{"type": "json", "json": data}]}


@server.tool(
    "analyze_targets",
    "Given user symptoms, return up to N related target proteins (no ligands, no ADMET)."
)
async def analyze_targets(symptoms: str = "", max_targets: int = 3):
    if not symptoms or not symptoms.strip():
        return {"content": [{"type": "json", "json": {"error": "symptoms required"}}]}

    max_targets = max(1, min(int(max_targets), 10))

    # 1) Normalize + extract keywords (robust fallbacks if LLM rate-limited)
    try:
        normalized = await llm.normalize_symptoms(symptoms)
    except Exception:
        normalized = symptoms
    try:
        keywords = await llm.extract_protein_keywords(normalized)
        if not keywords:
            raise ValueError("no keywords")
    except Exception:
        # Fallback: derive keywords from raw symptoms
        raw = [s.strip() for s in symptoms.replace(";", ",").split(",") if s.strip()]
        # take up to 5 words from each fragment
        kw_set = []
        for frag in raw:
            parts = [p for p in frag.split() if p.isalpha() and len(p) > 2]
            kw_set.extend(parts[:3])
        keywords = list(dict.fromkeys(kw_set))[:5]

    # 2) Query UniProt + STRING for candidates
    candidates: List[Dict[str, Any]] = []
    for kw in keywords[:5]:
        up = await uniprot.search(kw, limit=5)
        st = await stringdb.proteins_for_keyword(kw)
        candidates.extend(up)
        candidates.extend(st)

    # 3) Dedupe by protein_id and merge evidence/confidence
    uniq: Dict[str, Dict[str, Any]] = {}
    for c in candidates:
        pid = c.get("protein_id")
        if not pid:
            continue
        if pid not in uniq:
            uniq[pid] = c
        else:
            a = uniq[pid]
            a["confidence_score"] = max(a.get("confidence_score", 0.0), c.get("confidence_score", 0.0))
            a["evidence_sources"] = list({*(a.get("evidence_sources", [])), *(c.get("evidence_sources", []))})

    deduped = list(uniq.values())

    # 4) Relevance scoring (fallback if LLM unavailable)
    try:
        scored = await llm.score_protein_relevance(deduped, normalized)
        scored.sort(key=lambda x: x.get("confidence_score", 0.0), reverse=True)
    except Exception:
        # Fallback: sort by existing confidence_score if present
        scored = sorted(deduped, key=lambda x: x.get("confidence_score", 0.0), reverse=True)

    # If no candidates at all, synthesize pseudo-targets from keywords so downstream can proceed
    if not scored:
        pseudo = []
        for kw in (keywords or [])[:max_targets]:
            pseudo.append({
                "protein_id": None,
                "protein_name": kw,
                "gene_name": None,
                "organism": "Homo sapiens",
                "confidence_score": 0.2,
                "evidence_sources": ["keywords-fallback"],
                "function": None,
                "location": None,
            })
        scored = pseudo
    selected = scored[:max_targets]

    # 5) Optional: LLM description if function missing
    for p in selected:
        if not p.get("function"):
            try:
                desc = await llm.summarize_protein(p.get("protein_name", ""), p.get("function", ""))
                if desc:
                    p["function"] = desc
            except Exception:
                # leave function as-is
                pass

    # 6) Shape output and map to target_protein_1..N
    shaped = [_shape_target(p) for p in selected]
    resp: Dict[str, Any] = {}
    for i, item in enumerate(shaped, start=1):
        resp[f"target_protein_{i}"] = item

    return {"content": [{"type": "json", "json": resp}]}


@server.tool(
    "get_ligands_for_protein",
    "Fetch ligands for a target protein (ChEMBL/PubChem) and compute full ADMET (only for ligands)."
)
async def get_ligands_for_protein(protein_id: str = "", protein_name: str = "", max_ligands: int = 5, include_admet: bool = True):
    # Allow fallback by protein_name if no protein_id (to remain useful without UniProt)
    if not protein_id and not protein_name:
        return {"content": [{"type": "json", "json": {"error": "protein_id or protein_name required"}}]}

    max_ligands = max(1, min(int(max_ligands), 20))

    # 1) Prefer ChEMBL by UniProt accession when available
    ligands = []
    if protein_id:
        try:
            ligands = await chembl.get_protein_ligands(protein_id, limit=max_ligands)
        except Exception:
            ligands = []

    # 2) If not enough or no id path, fallback via PubChem name lookup
    if len(ligands) < max_ligands:
        fallback_kw = protein_name or protein_id
        if fallback_kw:
            more = await pubchem.search_compounds_by_name(fallback_kw, limit=max_ligands - len(ligands))
            ligands.extend(more)

    # 3) If still not enough, use LLM to suggest ligand names and query PubChem
    if len(ligands) < max_ligands and protein_name:
        try:
            suggested_names = await llm.suggest_ligands(protein_name, max_ligands - len(ligands))
            for name in suggested_names:
                more = await pubchem.search_compounds_by_name(name, limit=1)
                if more:
                    more[0]["source"] = "LLM-Suggested + PubChem"
                    ligands.extend(more)
        except Exception as e:
            pass  # Best effort

    # Dedupe by ligand_id
    seen = set()
    uniq: List[Dict[str, Any]] = []
    for l in ligands:
        lid = l.get("ligand_id") or l.get("chembl_id") or l.get("pubchem_cid")
        if not lid or lid in seen:
            continue
        seen.add(lid)
        uniq.append(l)

    # 3) Mechanism-of-action via LLM (best effort, robust to rate limits)
    for l in uniq:
        if not l.get("mechanism_of_action"):
            try:
                moa = await llm.summarize_ligand_moa(
                    l.get("ligand_name", ""), protein_name or protein_id, l.get("activity_type")
                )
                if moa:
                    l["mechanism_of_action"] = moa
            except Exception:
                pass

    # 4) ADMET only for ligands
    if include_admet:
        for l in uniq:
            smiles = l.get("smiles")
            if not smiles:
                continue
            raw = rdkit.admet(smiles)
            if raw:
                l["admet"] = format_admet_for_spec(raw)

    shaped = [_shape_ligand(l) for l in uniq[:max_ligands]]
    return {"content": [{"type": "json", "json": {"ligands": shaped}}]}


@server.tool(
    "react",
    "Perform chemical reaction between two SMILES compounds"
)
async def react(smiles1: str = "", smiles2: str = "", include_admet: bool = True):
    if not smiles1 or not smiles2:
        return {"content": [{"type": "json", "json": {"error": "Both smiles1 and smiles2 are required"}}]}

    # Import and use the actual reaction pipeline
    from src.reaction_mcpserver.main import run_reaction_pipeline

    def _run():
        return run_reaction_pipeline(smiles1, smiles2, include_admet=bool(include_admet))

    data = await anyio.to_thread.run_sync(_run)
    return {"content": [{"type": "json", "json": data}]}

@server.tool("health", "Health check")
async def health():
    return {"content": [{"type": "json", "json": {"status": "ok"}}]}

if __name__ == "__main__":
    # Do not print to stdout; MCP uses stdout. Use stderr for logs if needed.
    asyncio.run(server.run())