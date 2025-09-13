# mcpserver/reaction_mcpserver.py
import os
import sys
import time
import base64
import logging
import json
import asyncio
from io import BytesIO
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

import anyio

# Ensure project root is importable for "modules.*"
ROOT = Path(__file__).resolve().parents[1]  # project root
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, MolStandardize, rdChemReactions, Draw
RDLogger.DisableLog('rdApp.warning')

try:
    from rxn4chemistry import RXN4ChemistryWrapper
except Exception:
    RXN4ChemistryWrapper = None

# Your reaction knowledge bases
from modules.reaction_data import reaction_smarts_map
from modules.fun_grp_pattern import functional_group_patterns
from modules.cal_rxn_com import calculate_reaction_compatibility

logger = logging.getLogger("reaction-mcp")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Optional RXN4Chemistry config
RXN_API_KEY = os.getenv("RXN_API_KEY", "apk-c01d9c918e326ad080245e040134343481fe4a7513d859ae56743d11740cf93e")
RXN_PROJECT = os.getenv("RXN_PROJECT", "rxn4chemistry_tour")
rxn = None
if RXN4ChemistryWrapper and RXN_API_KEY:
    try:
        rxn = RXN4ChemistryWrapper(api_key=RXN_API_KEY)
        projects = rxn.get_projects()
        if not any(p.get("name") == RXN_PROJECT for p in projects):
            rxn.create_project(RXN_PROJECT)
        rxn.set_project(RXN_PROJECT)
        logger.info(f"RXN4Chemistry project set to {rxn.project_id}")
    except Exception as e:
        logger.warning(f"RXN init failed; using RDKit fallback. {e}")
        rxn = None
else:
    logger.info("RXN_API_KEY not set or rxn4chemistry not installed. Using RDKit rules only.")

# Simple MCP Server class (same as in main.py)
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
        while True:
            try:
                line = await asyncio.get_event_loop().run_in_executor(None, sys.stdin.readline)
                if not line:
                    break
               
                request = json.loads(line.strip())
                response = await self.handle_request(request)
               
                print(json.dumps(response), flush=True)
               
            except json.JSONDecodeError:
                continue
            except Exception as e:
                error_response = {
                    "error": {
                        "code": -32700,
                        "message": f"Parse error: {str(e)}"
                    }
                }
                print(json.dumps(error_response), flush=True)

# --------------------- Core molecule helpers ---------------------
def detect_functional_groups(mol):
    if not mol:
        return []
    groups = []
    for group_name, pattern in functional_group_patterns.items():
        if pattern and mol.HasSubstructMatch(pattern):
            groups.append(group_name)
    return groups

def standardize_molecule(mol):
    if not mol:
        return None
    try:
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
        standardizer = MolStandardize.Standardizer()
        mol = standardizer.standardize(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        mol = Chem.AddHs(mol)
        return mol
    except Exception:
        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return None

def compute_product_properties(mol):
    try:
        Chem.SanitizeMol(mol)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotbonds = Descriptors.NumRotatableBonds(mol)
        has_stereo = any(atom.HasProp('_ChiralityPossible') for atom in mol.GetAtoms()) or \
            any(bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds())
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True, allHsExplicit=False)
        return {
            "smiles": smiles,
            "molecular_weight": round(mw, 2),
            "num_atoms": mol.GetNumAtoms(),
            "logP": round(logp, 2),
            "tpsa": round(tpsa, 2),
            "num_h_donors": hbd,
            "num_h_acceptors": hba,
            "num_rotatable_bonds": rotbonds,
            "has_stereochemistry": has_stereo,
            "functional_groups": detect_functional_groups(mol),
        }
    except Exception as e:
        return {"error": f"Property computation failed: {e}", "functional_groups": []}

def compute_admet_properties(mol, reaction_type=None, reactant_groups=None, product_groups=None):
    Chem.SanitizeMol(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    arom_rings = Descriptors.NumAromaticRings(mol)
    fsp3 = Descriptors.FractionCSP3(mol)
    mol_refr = Descriptors.MolMR(mol)

    if reaction_type and reactant_groups and product_groups:
        ames_risk = "Yes" if arom_rings > 1 else "No"
        non_toxic_score = 20 if ames_risk == "Yes" else 80
        if "aromatic_addition" in (reaction_type or "").lower() and arom_rings > 0:
            ames_risk = "Yes"; non_toxic_score = 10
        herg_risk = "Yes" if (logp > 3 and mw > 400) else "No"
        herg_safe_score = 30 if herg_risk == "Yes" else 90
        if "lipophilic_addition" in (reaction_type or "").lower() and logp > 2:
            herg_risk = "Yes"; herg_safe_score = 20
        cyp3a4 = "Yes" if (mw > 400 and logp > 3) else "No"
        bioavailable_score = 40 if cyp3a4 == "Yes" else 80
        if "amide_formation" in (reaction_type or "").lower() and "amide" in (product_groups or []):
            cyp3a4 = "Yes"; bioavailable_score = 30
    else:
        ames_risk = "Yes" if arom_rings > 1 else "No"
        non_toxic_score = 20 if ames_risk == "Yes" else 80
        herg_risk = "Yes" if (logp > 3 and mw > 400) else "No"
        herg_safe_score = 30 if herg_risk == "Yes" else 90
        cyp3a4 = "Yes" if (mw > 400 and logp > 3) else "No"
        bioavailable_score = 40 if cyp3a4 == "Yes" else 80

    log_s = -logp + 0.5
    soluble_score = 90 if log_s > -3 else 50 if log_s > -5 else 20
    lipinski_violations = sum([mw > 500, logp > 5, h_donors > 5, h_acceptors > 10])
    if lipinski_violations > 1:
        bioavailable_score = max(10, bioavailable_score - (lipinski_violations * 20))

    gi_absorption = "High" if (tpsa < 140 and logp < 5) else "Low"
    caco2 = "High" if (mw < 500 and tpsa < 100 and logp < 3) else "Low"
    ppb = "High" if logp > 3 else "Low"
    vdss = "High" if logp > 3 and mw > 500 else "Low"
    bbb = "Yes" if (mw < 500 and 2 < logp < 5 and tpsa < 90 and h_donors < 3) else "No"
    clearance = "Low" if (logp > 3 and mw > 500) else "High"
    oct2 = "No"
    ld50 = "Moderate"

    return [
        {"section": "Physicochemical", "properties": [
            {"name": "Molecular Weight", "prediction": round(mw, 2), "units": "g/mol"},
            {"name": "LogP", "prediction": round(logp, 2), "units": "-"},
            {"name": "TPSA", "prediction": round(tpsa, 2), "units": "Å²"},
            {"name": "H-Bond Donors", "prediction": Descriptors.NumHDonors(mol), "units": "-"},
            {"name": "H-Bond Acceptors", "prediction": Descriptors.NumHAcceptors(mol), "units": "-"},
            {"name": "Rotatable Bonds", "prediction": Descriptors.NumRotatableBonds(mol), "units": "-"},
            {"name": "Aromatic Rings", "prediction": arom_rings, "units": "-"},
            {"name": "Fraction sp³ Carbons", "prediction": round(fsp3, 2), "units": "-"},
            {"name": "Molar Refractivity", "prediction": round(mol_refr, 2), "units": "-"},
        ]},
        {"section": "Absorption", "properties": [
            {"name": "Lipinski Violations", "prediction": lipinski_violations, "units": "-"},
            {"name": "GI Absorption", "prediction": gi_absorption, "units": "-"},
            {"name": "Aqueous Solubility (LogS)", "prediction": round(log_s, 2), "units": "-"},
            {"name": "Caco-2 Permeability", "prediction": caco2, "units": "-"},
        ]},
        {"section": "Distribution", "properties": [
            {"name": "Plasma Protein Binding", "prediction": ppb, "units": "-"},
            {"name": "Volume of Distribution", "prediction": vdss, "units": "-"},
            {"name": "BBB Permeability", "prediction": bbb, "units": "-"},
        ]},
        {"section": "Metabolism", "properties": [
            {"name": "CYP3A4 Inhibition", "prediction": "Yes" if (mw > 400 and logp > 3) else "No", "units": "-"},
            {"name": "CYP2D6 Inhibition", "prediction": "Yes" if (logp > 2 and Descriptors.NumHDonors(mol) > 0) else "No", "units": "-"},
            {"name": "CYP2C9 Inhibition", "prediction": "Yes" if (logp > 3 and arom_rings > 0) else "No", "units": "-"},
        ]},
        {"section": "Excretion", "properties": [
            {"name": "Total Clearance", "prediction": clearance, "units": "-"},
            {"name": "Renal OCT2 Substrate", "prediction": oct2, "units": "-"},
        ]},
        {"section": "Toxicity", "properties": [
            {"name": "hERG Inhibition", "prediction": "Yes" if (logp > 3 and mw > 400) else "No", "units": "-"},
            {"name": "AMES Mutagenicity", "prediction": "Yes" if arom_rings > 1 else "No", "units": "-"},
            {"name": "LD50 (Oral Rat)", "prediction": ld50, "units": "-"},
            {"name": "Skin Sensitization", "prediction": "No", "units": "-"},
        ]},
        {"section": "ADMET", "properties": [
            {"name": "Non-Toxic", "prediction": 80 if arom_rings <= 1 else 20, "units": "%"},
            {"name": "Soluble", "prediction": 90 if (-logp + 0.5) > -3 else 50 if (-logp + 0.5) > -5 else 20, "units": "%"},
            {"name": "Bioavailable", "prediction": 40 if (mw > 400 and logp > 3) else 80, "units": "%"},
            {"name": "hERG Safe", "prediction": 30 if (logp > 3 and mw > 400) else 90, "units": "%"},
        ]},
    ]

def predict_reaction_with_rxn(smiles1: str, smiles2: str) -> Tuple[Optional[str], float]:
    if not rxn:
        return None, 0.0
    try:
        reactants = f"{smiles1}.{smiles2}>>"
        resp = rxn.predict_reaction(reactants)
        if "prediction_id" not in resp:
            return None, 0.0
        time.sleep(2)
        pred_id = resp["prediction_id"]
        results = rxn.get_predict_reaction_results(pred_id)
        for _ in range(8):
            if "results" in results and results["results"]:
                break
            time.sleep(2)
            results = rxn.get_predict_reaction_results(pred_id)
        preds = results.get("results", [])
        if not preds:
            return None, 0.0
        top = preds[0]
        s = top.get("smiles", "")
        conf = float(top.get("confidence", 0.0))
        if conf < 0.1:
            return None, 0.0
        return s, conf
    except Exception:
        return None, 0.0

def predict_reaction_with_rdkit(m1, m2, s1: str, s2: str):
    g1 = detect_functional_groups(m1)
    g2 = detect_functional_groups(m2)
    applicable = []
    for rxn_type, data in reaction_smarts_map.items():
        comp = calculate_reaction_compatibility(g1, g2, rxn_type)
        if comp > 0.3:
            priority = data.get("priority", 5)
            applicable.append((rxn_type, comp, priority))
    applicable.sort(key=lambda x: (x[1], x[2]), reverse=True)
    if not applicable:
        return None, 0.0, None
    for rxn_type, comp, prio in applicable:
        try:
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map[rxn_type]["smarts"])
            rxn_smarts.Initialize()
            reactants = (m1, m2) if rxn_smarts.GetNumReactantTemplates() > 1 else ((m1,) if len(g1) >= len(g2) else (m2,))
            products = rxn_smarts.RunReactants(reactants)
            if not products:
                continue
            product_mols, seen, valid, total = [], set(), 0, 0
            for pset in products:
                for prod in pset:
                    if not prod:
                        continue
                    total += 1
                    prod = standardize_molecule(prod)
                    if not prod:
                        continue
                    smiles = Chem.MolToSmiles(prod, isomericSmiles=True, canonical=True)
                    m = Chem.MolFromSmiles(smiles)
                    if not m:
                        continue
                    try:
                        Chem.SanitizeMol(m)
                        if len(Chem.GetMolFrags(m)) == 1 and smiles not in seen:
                            seen.add(smiles)
                            product_mols.append(prod)
                            valid += 1
                    except Exception:
                        continue
            if product_mols:
                validity = valid / total if total else 0.0
                conf = min(1.0, (comp * 0.5 + validity * 0.4 + prio / 50.0))
                return product_mols, conf, rxn_type
        except Exception:
            continue
    return None, 0.0, None

def process_product_smiles(product_smiles: str) -> List[Chem.Mol]:
    if not product_smiles or product_smiles.strip() == "":
        return []
    products, seen = [], set()
    try:
        mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
        if mol:
            Chem.SanitizeMol(mol)
            mol = standardize_molecule(mol)
            if mol:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                if smi not in seen:
                    seen.add(smi)
                    products.append(mol)
        return products
    except Exception:
        pass
    for sep in ['.', '>>>', '>>', '>', '|']:
        product_smiles = "|".join(product_smiles.split(sep))
    for part in product_smiles.split("|"):
        part = part.strip()
        if not part or len(part) < 2:
            continue
        try:
            pm = Chem.MolFromSmiles(part, sanitize=False)
            if pm:
                Chem.SanitizeMol(pm)
                pm = standardize_molecule(pm)
                if pm and pm.GetNumAtoms() >= 3:
                    s = Chem.MolToSmiles(pm, isomericSmiles=True, canonical=True)
                    if s not in seen:
                        seen.add(s)
                        products.append(pm)
        except Exception:
            continue
    return products

def run_reaction_pipeline(smiles1: str, smiles2: str, include_admet: bool = True) -> Dict[str, Any]:
    if not smiles1 or not smiles2:
        return {"error": "Both smiles1 and smiles2 are required"}
    m1 = Chem.MolFromSmiles(smiles1)
    m2 = Chem.MolFromSmiles(smiles2)
    if not m1 or not m2:
        return {"error": "Invalid SMILES strings"}

    s1 = Chem.MolToSmiles(m1, isomericSmiles=True, canonical=True)
    s2 = Chem.MolToSmiles(m2, isomericSmiles=True, canonical=True)

    m1 = standardize_molecule(m1); m2 = standardize_molecule(m2)
    if not m1 or not m2:
        return {"error": "Molecule standardization failed"}
    Chem.SanitizeMol(m1); Chem.SanitizeMol(m2)

    g1 = detect_functional_groups(m1)
    g2 = detect_functional_groups(m2)

    reactants = []
    for smi, mol in [(s1, m1), (s2, m2)]:
        props = compute_product_properties(mol)
        admet = compute_admet_properties(mol) if include_admet else []
        reactants.append({"smiles": smi, "properties": props, "admet_properties": admet})

    # Try RXN, fallback RDKit
    product_smiles, confidence = predict_reaction_with_rxn(s1, s2)
    reaction_type = "RXN4Chemistry_Prediction"
    product_mols = []
    if product_smiles and confidence >= 0.15:
        product_mols = process_product_smiles(product_smiles)
    else:
        pmols, rd_conf, rd_type = predict_reaction_with_rdkit(m1, m2, s1, s2)
        if pmols and rd_conf > 0.3:
            product_mols = pmols
            confidence = rd_conf
            reaction_type = rd_type
            product_smiles = ".".join(Chem.MolToSmiles(x, isomericSmiles=True) for x in product_mols if x)
        else:
            return {
                "reactants": reactants,
                "reactionResults": [],
                "failedReactions": [{
                    "reactants": [s1, s2],
                    "reason": f"No viable reaction predicted (RXN: {confidence:.2f}, RDKit failed)",
                    "reactionType": "Unknown"
                }],
                "statistics": {"mean_mw": 0, "std_mw": 0, "min_mw": 0, "max_mw": 0},
                "createdAt": datetime.utcnow().isoformat() + "Z",
                "processedReactionTypes": {},
                "product_smiles": []
            }

    product_info = []
    for pm in product_mols:
        pm = standardize_molecule(pm)
        if not pm:
            continue
        Chem.SanitizeMol(pm)
        props = compute_product_properties(pm)
        if "error" in props:
            continue
        admet_props = compute_admet_properties(pm, reaction_type=reaction_type,
                                               reactant_groups=[g1, g2],
                                               product_groups=detect_functional_groups(pm)) if include_admet else []
        props["admet_properties"] = admet_props
        product_info.append(props)

    if not product_info:
        return {
            "reactants": reactants,
            "reactionResults": [],
            "failedReactions": [{
                "reactants": [s1, s2],
                "reason": "All products failed property computation",
                "reactionType": reaction_type
            }],
            "statistics": {"mean_mw": 0, "std_mw": 0, "min_mw": 0, "max_mw": 0},
            "createdAt": datetime.utcnow().isoformat() + "Z",
            "processedReactionTypes": {},
            "product_smiles": []
        }

    product_info.sort(key=lambda x: x["molecular_weight"], reverse=True)
    selected = product_info[:2]
    mws = [p["molecular_weight"] for p in selected]

    reaction_result = {
        "reactionType": reaction_type,
        "description": reaction_smarts_map.get(reaction_type, {}).get("description", "Predicted reaction"),
        "reactants": [s1, s2],
        "reactantGroups": [g1, g2],
        "products": selected,
        "confidence": float(round(confidence, 3)),
        "productSmiles": product_smiles or ""
    }

    return {
        "reactants": reactants,
        "reactionResults": [reaction_result],
        "failedReactions": [],
        "statistics": {
            "mean_mw": float(sum(mws)/len(mws)) if mws else 0.0,
            "std_mw": 0.0,
            "min_mw": min(mws) if mws else 0.0,
            "max_mw": max(mws) if mws else 0.0
        },
        "createdAt": datetime.utcnow().isoformat() + "Z",
        "processedReactionTypes": {f"{min(s1, s2)}:{max(s1, s2)}": [reaction_type]},
        "product_smiles": [reaction_result["productSmiles"]] if reaction_result["productSmiles"] else []
    }

def smiles_png_base64(smiles: str, size: int = 300) -> str:
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if not mol:
        raise ValueError("Invalid SMILES")
    mol = standardize_molecule(mol)
    if not mol:
        raise RuntimeError("Failed to standardize molecule")
    Chem.SanitizeMol(mol)
    img = Draw.MolToImage(mol, size=(size, size))
    buf = BytesIO(); img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("utf-8")

# --------------------- Registration for a parent server ---------------------
def register_reaction_tools(server):
    """Register reaction tools with any server that has a @server.tool decorator"""

    @server.tool("react", "Predict reaction between two SMILES (RXN4Chemistry + RDKit fallback).")
    async def react(smiles1: str = "", smiles2: str = "", include_admet: bool = True):
        def _run():
            return run_reaction_pipeline(smiles1, smiles2, include_admet=bool(include_admet))
        data = await anyio.to_thread.run_sync(_run)
        return {"content": [{"type": "json", "json": data}]}

    @server.tool("smiles_to_image", "Render SMILES to PNG as base64.")
    async def smiles_to_image(smiles: str = "", size: int = 300):
        def _run():
            return smiles_png_base64(smiles, size=int(size))
        b64png = await anyio.to_thread.run_sync(_run)
        return {"content": [{"type": "image", "data": b64png, "mimeType": "image/png"}]}

    @server.tool("reactions", "List available reaction SMARTS and metadata.")
    async def reactions():
        info = []
        for rtype, data in reaction_smarts_map.items():
            info.append({
                "type": rtype,
                "description": data.get("description", ""),
                "priority": data.get("priority", 5),
                "smarts": data.get("smarts", "")
            })
        info.sort(key=lambda x: x["priority"], reverse=True)
        return {"content": [{"type": "json", "json": {"total_reactions": len(info), "reactions": info}}]}

# --------------------- Standalone server mode -------------------
if __name__ == "__main__":
    # Create standalone server
    server = SimpleMCPServer("chemical-reaction")

    # Register reaction tools
    register_reaction_tools(server)

    # Add health check
    @server.tool("health", "Health check")
    async def health():
        return {"content": [{"type": "json", "json": {
            "status": "ok",
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "reactions_available": len(reaction_smarts_map),
            "functional_groups_available": len(functional_group_patterns),
            "rxn4chemistry_available": rxn is not None
        }}]}

    # Run the server
    asyncio.run(server.run())