from typing import List, Dict, Any

def _shape_target(p: Dict[str, Any]) -> Dict[str, Any]:
    # Aligns with TargetProtein schema (no ligands/ADMET here)
    return {
        "protein_id": p.get("protein_id"),
        "protein_name": p.get("protein_name"),
        "gene_name": p.get("gene_name"),
        "organism": p.get("organism", "Homo sapiens"),
        "confidence_score": float(p.get("confidence_score", 0.5)),
        "evidence_sources": list(p.get("evidence_sources", [])),
        "biological_function": p.get("function") or None,
        "subcellular_location": p.get("location") or None,
        "ligands": [],  # intentionally empty in this server
    }

def _shape_ligand(l: Dict[str, Any]) -> Dict[str, Any]:
    # Aligns with LigandInfo schema plus "admet" for ligands
    return {
        "ligand_id": l.get("ligand_id"),
        "ligand_name": l.get("ligand_name"),
        "smiles": l.get("smiles"),
        "molecular_weight": l.get("molecular_weight"),
        "binding_affinity": l.get("binding_affinity"),
        "drugbank_id": l.get("drugbank_id"),
        "chembl_id": l.get("chembl_id"),
        "pubchem_cid": l.get("pubchem_cid"),
        "mechanism_of_action": l.get("mechanism_of_action"),
        "clinical_status": l.get("clinical_status"),
        "admet": l.get("admet"),
    }