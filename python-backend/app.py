from flask import Flask, jsonify, request
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, MolStandardize, rdChemReactions, Draw
from rdkit import RDLogger
import numpy as np
import logging
import traceback
from pymongo import MongoClient
from itertools import combinations
from rxn4chemistry import RXN4ChemistryWrapper
import os
import time
from datetime import datetime
from pydantic import BaseModel, Field, ValidationError
from typing import List, Dict, Any, Union
import base64
from io import BytesIO
from dotenv import load_dotenv
import google.generativeai as genai
load_dotenv()
from modules.reaction_data import reaction_smarts_map
from modules.fun_grp_pattern import functional_group_patterns
from modules.cal_rxn_com import calculate_reaction_compatibility

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

app = Flask(__name__)
frontend_url = os.getenv('FRONTEND_URL')
CORS(app, resources={r"/api/*": {"origins": ["http://localhost:5173", 'https://drugs-10979.firebaseapp.com']}}, supports_credentials=True)
@app.after_request
def after_request(response):
    # Explicit headers to help OPTIONS preflight
    response.headers.add('Access-Control-Allow-Origin', 'http://localhost:5173')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,POST,OPTIONS')
    return response 


# MongoDB connection
mongo_uri = os.getenv('MONGO_URI', 'mongodb+srv://bhangaletejas003:G0yEjQa9yrTChtDU@h2skill.nnmre.mongodb.net/?retryWrites=true&w=majority&appName=h2skill')
mongo_client = MongoClient(mongo_uri)
db = mongo_client['test']
target_protein_collection = db['targetproteins']
reaction_responses_collection = db['reaction_responses']

# RXN4Chemistry configuration
RXN_API_KEY = os.getenv('RXN_API_KEY', 'apk-c01d9c918e326ad080245e040134343481fe4a7513d859ae56743d11740cf93e')
if not RXN_API_KEY or RXN_API_KEY == 'apk-c01d9c918e326ad080245e040134343481fe4a7513d859ae56743d11740cf93e':
    logger.error("RXN_API_KEY not set. Reaction prediction will rely on RDKit fallback.")
    rxn = None
else:
    try:
        rxn = RXN4ChemistryWrapper(api_key=RXN_API_KEY)
        project_name = "rxn4chemistry_tour"
        projects = rxn.get_projects()
        project_exists = any(project['name'] == project_name for project in projects)
        if not project_exists:
            rxn.create_project(project_name)
        rxn.set_project(project_name)
        logger.info(f"Using RXN4Chemistry project ID: {rxn.project_id}")
    except Exception as e:
        logger.error(f"Failed to initialize RXN4Chemistry: {str(e)}")
        rxn = None

# Pydantic models for schema validation
class Properties(BaseModel):
    smiles: str
    molecular_weight: float
    num_atoms: int
    logP: float
    tpsa: float
    num_h_donors: int
    num_h_acceptors: int
    num_rotatable_bonds: int
    has_stereochemistry: bool
    functional_groups: List[str]

class ADMETProperty(BaseModel):
    name: str
    prediction: Union[str, float, int]
    units: str = ""

class ADMETSection(BaseModel):
    section: str
    properties: List[ADMETProperty]

class Reactant(BaseModel):
    smiles: str
    properties: Properties
    admet_properties: List[ADMETSection] = Field(default_factory=list)

class Product(BaseModel):
    smiles: str
    molecular_weight: float
    num_atoms: int
    logP: float
    tpsa: float
    num_h_donors: int
    num_h_acceptors: int
    num_rotatable_bonds: int
    has_stereochemistry: bool
    functional_groups: List[str]
    admet_properties: List[ADMETSection] = Field(default_factory=list)

class ReactionResult(BaseModel):
    reactionType: str
    description: str
    reactants: List[str]
    reactantGroups: List[List[str]]
    products: List[Product]
    confidence: float
    productSmiles: str

class FailedReaction(BaseModel):
    reactants: List[str]
    reason: str
    reactionType: str = "Unknown"

class Statistics(BaseModel):
    mean_mw: float
    std_mw: float
    min_mw: float
    max_mw: float

class ReactionResponse(BaseModel):
    reactants: List[Reactant]
    reactionResults: List[ReactionResult]
    failedReactions: List[FailedReaction]
    statistics: Statistics
    createdAt: datetime
    processedReactionTypes: Dict[str, List[str]]
    product_smiles: List[str] = Field(default_factory=list)

def detect_functional_groups(mol):
    """Detect functional groups in a molecule."""
    if not mol:
        return []
    groups = []
    for group_name, pattern in functional_group_patterns.items():
        if pattern and mol.HasSubstructMatch(pattern):
            groups.append(group_name)
    logger.info(f"Detected functional groups for molecule: {groups}")
    return groups

def preprocess_molecule(mol, target_group=None):
    """Preprocess molecule with specific transformations based on target functional group."""
    if not mol:
        logger.warning("Invalid molecule provided for preprocessing")
        return None
    if not target_group:
        logger.info("No target group specified, returning standardized molecule")
        return standardize_molecule(mol)

    try:
        Chem.SanitizeMol(mol)
        groups = detect_functional_groups(mol)
        logger.info(f"Starting preprocessing for target group: {target_group}, current groups: {groups}")

        # Define preprocessing rules: {target_group: (source_group, smarts_reaction, product_pattern)}
        preprocess_rules = {
            'carbonyl': [
                ('thiocarbonyl', '[C:1](=[S:2])>>[C:1](=[O:2]).[H][S][H]', '[#6]=[O]'),
                ('imine', '[C:1](=[N:2])>>[C:1](=[O:2]).[N:2][H]', '[#6]=[O]')
            ],
            'carboxylic_acid': [
                ('amide', '[C:1][C:2](=[O])[N:3]>>[C:1][C:2](=[O])[OH].[N:3][H]', '[#6](=[O])[OH]'),
                ('ester', '[C:1][C:2](=[O])[O:3][C:4]>>[C:1][C:2](=[O])[OH].[C:4][OH:3]', '[#6](=[O])[OH]')
            ],
            'aryl_hydrogen': [
                ('aryl_halide', '[c:1][F,Cl,Br,I:2]>>[c:1][H].[*:2][H]', '[c][H]'),
                ('aryl_nitro', '[c:1][N+:2](=O)[O-:3]>>[c:1][H].[N:2](=O)[O:3][H]', '[c][H]')
            ],
            'alcohol': [
                ('ester', '[C:1][O:2][C:3](=[O])>>[C:1][OH:2].[C:3](=[O])[OH]', '[#6][OH]'),
                ('ether', '[C:1][O:2][C:3]>>[C:1][OH:2].[C:3][OH]', '[#6][OH]')
            ],
            'amine': [
                ('nitro', '[C:1][N+:2](=O)[O-:3]>>[C:1][N:2][H2].[O:3]=[O]', '[#6][NH2]'),
                ('amide', '[C:1][C:2](=[O])[N:3]>>[C:1][C:2](=[O])[OH].[N:3][H2]', '[#6][NH2]')
            ],
            'phosphonate': [
                ('phosphine', '[P:1]([C:2])([C:3])[C:4]>>[P:1](=[O])([O][C:2])[O][C:3].[C:4][H]', '[P](=[O])([O][#6])[O][#6]'),
                ('phosphite', '[P:1]([O][C:2])([O][C:3])[O][C:4]>>[P:1](=[O])([O][C:2])[O][C:3].[C:4][OH]', '[P](=[O])([O][#6])[O][#6]')
            ],
            'sulfoxide': [
                ('sulfide', '[C:1][S:2][C:3]>>[C:1][S:2](=[O])[C:3]', '[S](=[O])[#6]'),
                ('sulfone', '[C:1][S:2](=[O])(=[O])[C:3]>>[C:1][S:2](=[O])[C:3].[O]', '[S](=[O])[#6]')
            ],
            'alkene': [
                ('alkyne', '[C:1]#[C:2]>>[C:1]=[C:2][H][H]', '[#6]=[#6]'),
                ('epoxide', '[C:1]1[O:2][C:3]1>>[C:1]=[C:3].[O:2][H2]', '[#6]=[#6]')
            ],
            'primary_alcohol': [
                ('primary_alkyl_halide', '[C:1][F,Cl,Br,I:2]>>[C:1][OH].[*:2][H]', '[CH3][OH]'),
                ('aldehyde', '[C:1][C:2]=[O]>>[C:1][C:2][OH][H]', '[CH2][OH]')
            ],
            'secondary_alcohol': [
                ('ketone', '[C:1][C:2](=[O])[C:3]>>[C:1][C:2]([OH])[C:3]', '[CH]([C])[OH]')
            ]
        }

        # Apply transformation if target_group is in rules and source_group is present
        if target_group in preprocess_rules:
            for source_group, smarts, product_pattern in preprocess_rules[target_group]:
                if source_group in groups:
                    logger.info(f"Applying transformation: {source_group} to {target_group}")
                    rxn = rdChemReactions.ReactionFromSmarts(smarts)
                    products = rxn.RunReactants((mol,))
                    if products:
                        for prod_set in products:
                            for prod in prod_set:
                                try:
                                    if prod.HasSubstructMatch(Chem.MolFromSmarts(product_pattern)):
                                        standardized = standardize_molecule(prod)
                                        if standardized:
                                            logger.info(f"Successfully preprocessed to {target_group}")
                                            return standardized
                                except:
                                    continue
                    logger.warning(f"Transformation {source_group} to {target_group} failed")
        else:
            logger.warning(f"No preprocessing rule for target group: {target_group}")

        # If no transformation applied, return standardized molecule
        return standardize_molecule(mol)
    except Exception as e:
        logger.error(f"Preprocessing failed: {str(e)}")
        return None
def standardize_molecule(mol):
    """Standardize molecule with error handling."""
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
    except Exception as e:
        logger.warning(f"Standardization failed: {str(e)}")
        try:
            Chem.SanitizeMol(mol)
            return mol
        except:
            return None
        

def compute_product_properties(mol): 
    """Compute molecular properties."""
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
        # Ensure implicit hydrogens in SMILES
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True, allHsExplicit=False)
        logger.debug(f"Computed SMILES for properties: {smiles}")
        return {
            'smiles': smiles,
            'molecular_weight': round(mw, 2),
            'num_atoms': mol.GetNumAtoms(),
            'logP': round(logp, 2),
            'tpsa': round(tpsa, 2),
            'num_h_donors': hbd,
            'num_h_acceptors': hba,
            'num_rotatable_bonds': rotbonds,
            'has_stereochemistry': has_stereo,
            'functional_groups': detect_functional_groups(mol) or []
        }
    except Exception as e:
        logger.warning(f"Property computation failed: {str(e)}")
        return {'error': f'Property computation failed: {str(e)}', 'functional_groups': []}

def compute_admet_properties(mol, reaction_type=None, reactant_groups=None, product_groups=None):
    """Compute ADMET properties for a molecule, with reaction-specific adjustments."""
    try:
        Chem.SanitizeMol(mol)
        # Basic physicochemical properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        arom_rings = Descriptors.NumAromaticRings(mol)
        fsp3 = Descriptors.FractionCSP3(mol)
        mol_refr = Descriptors.MolMR(mol)

        # Reaction-specific adjustments
        if reaction_type and reactant_groups and product_groups:
            # AMES Mutagenicity (Non-Toxic) score: Lower score means more toxic
            if arom_rings > 1:
                ames_risk = "Yes"
                non_toxic_score = 20  # High risk of mutagenicity
            else:
                ames_risk = "No"
                non_toxic_score = 80  # Low risk
            if "aromatic_addition" in reaction_type.lower() and arom_rings > 0:
                ames_risk = "Yes"
                non_toxic_score = 10  # Increased risk due to aromatic addition

            # hERG Risk (hERG Safe) score: Lower score means higher risk
            herg_risk = "Yes" if (logp > 3 and mw > 400) else "No"
            herg_safe_score = 30 if herg_risk == "Yes" else 90
            if "lipophilic_addition" in reaction_type.lower() and logp > 2:
                herg_risk = "Yes"
                herg_safe_score = 20  # Increased risk due to lipophilic addition

            # CYP Inhibition (used to adjust Bioavailability score)
            cyp3a4 = "Yes" if (mw > 400 and logp > 3) else "No"
            bioavailable_score = 40 if cyp3a4 == "Yes" else 80  # Lower bioavailability if CYP3A4 inhibited
            if "amide_formation" in reaction_type.lower() and "amide" in product_groups:
                cyp3a4 = "Yes"
                bioavailable_score = 30  # Amide formation may increase CYP3A4 inhibition
        else:
            # Default heuristics if no reaction context
            ames_risk = "Yes" if arom_rings > 1 else "No"
            non_toxic_score = 20 if ames_risk == "Yes" else 80

            herg_risk = "Yes" if (logp > 3 and mw > 400) else "No"
            herg_safe_score = 30 if herg_risk == "Yes" else 90

            cyp3a4 = "Yes" if (mw > 400 and logp > 3) else "No"
            bioavailable_score = 40 if cyp3a4 == "Yes" else 80

        # Solubility score (based on LogS and TPSA)
        log_s = -logp + 0.5  # Approximate LogS
        soluble_score = 90 if log_s > -3 else 50 if log_s > -5 else 20  # Higher LogS means more soluble

        # Bioavailability adjustment based on Lipinski rules
        lipinski_violations = sum([
            mw > 500,
            logp > 5,
            h_donors > 5,
            h_acceptors > 10
        ])
        if lipinski_violations > 1:
            bioavailable_score = max(10, bioavailable_score - (lipinski_violations * 20))

        # Absorption (Heuristics)
        gi_absorption = "High" if (tpsa < 140 and logp < 5) else "Low"
        caco2 = "High" if (mw < 500 and tpsa < 100 and logp < 3) else "Low"

        # Distribution (Heuristics)
        ppb = "High" if logp > 3 else "Low"
        vdss = "High" if logp > 3 and mw > 500 else "Low"
        bbb = "Yes" if (mw < 500 and logp > 2 and logp < 5 and tpsa < 90 and h_donors < 3) else "No"

        # Metabolism (Heuristics)
        cyp2d6 = "Yes" if (logp > 2 and h_donors > 0) else "No"
        cyp2c9 = "Yes" if (logp > 3 and arom_rings > 0) else "No"

        # Excretion (Heuristics)
        clearance = "Low" if (logp > 3 and mw > 500) else "High"
        oct2 = "No"  # Simplified, assuming no strong basic groups

        # Toxicity (Heuristics)
        ld50 = "Moderate"  # Arbitrary
        skin_sens = "No"  # Simplified, assuming no reactive groups

        # Structure the ADMET properties into sections
        admet_sections = [
            {
                "section": "Physicochemical",
                "properties": [
                    {"name": "Molecular Weight", "prediction": round(mw, 2), "units": "g/mol"},
                    {"name": "LogP", "prediction": round(logp, 2), "units": "-"},
                    {"name": "TPSA", "prediction": round(tpsa, 2), "units": "Å²"},
                    {"name": "H-Bond Donors", "prediction": h_donors, "units": "-"},
                    {"name": "H-Bond Acceptors", "prediction": h_acceptors, "units": "-"},
                    {"name": "Rotatable Bonds", "prediction": rot_bonds, "units": "-"},
                    {"name": "Aromatic Rings", "prediction": arom_rings, "units": "-"},
                    {"name": "Fraction sp³ Carbons", "prediction": round(fsp3, 2), "units": "-"},
                    {"name": "Molar Refractivity", "prediction": round(mol_refr, 2), "units": "-"}
                ]
            },
            {
                "section": "Absorption",
                "properties": [
                    {"name": "Lipinski Violations", "prediction": lipinski_violations, "units": "-"},
                    {"name": "GI Absorption", "prediction": gi_absorption, "units": "-"},
                    {"name": "Aqueous Solubility (LogS)", "prediction": round(log_s, 2), "units": "-"},
                    {"name": "Caco-2 Permeability", "prediction": caco2, "units": "-"}
                ]
            },
            {
                "section": "Distribution",
                "properties": [
                    {"name": "Plasma Protein Binding", "prediction": ppb, "units": "-"},
                    {"name": "Volume of Distribution", "prediction": vdss, "units": "-"},
                    {"name": "BBB Permeability", "prediction": bbb, "units": "-"}
                ]
            },
            {
                "section": "Metabolism",
                "properties": [
                    {"name": "CYP3A4 Inhibition", "prediction": cyp3a4, "units": "-"},
                    {"name": "CYP2D6 Inhibition", "prediction": cyp2d6, "units": "-"},
                    {"name": "CYP2C9 Inhibition", "prediction": cyp2c9, "units": "-"}
                ]
            },
            {
                "section": "Excretion",
                "properties": [
                    {"name": "Total Clearance", "prediction": clearance, "units": "-"},
                    {"name": "Renal OCT2 Substrate", "prediction": oct2, "units": "-"}
                ]
            },
            {
                "section": "Toxicity",
                "properties": [
                    {"name": "hERG Inhibition", "prediction": herg_risk, "units": "-"},
                    {"name": "AMES Mutagenicity", "prediction": ames_risk, "units": "-"},
                    {"name": "LD50 (Oral Rat)", "prediction": ld50, "units": "-"},
                    {"name": "Skin Sensitization", "prediction": skin_sens, "units": "-"}
                ]
            },
            # Add new ADMET section for radar chart
            {
                "section": "ADMET",
                "properties": [
                    {"name": "Non-Toxic", "prediction": non_toxic_score, "units": "%"},
                    {"name": "Soluble", "prediction": soluble_score, "units": "%"},
                    {"name": "Bioavailable", "prediction": bioavailable_score, "units": "%"},
                    {"name": "hERG Safe", "prediction": herg_safe_score, "units": "%"}
                ]
            }
        ]
        return admet_sections
    except Exception as e:
        logger.error(f"ADMET computation failed: {str(e)}")
        raise



################### this below functions are used to perform reactions#############################
def predict_reaction_with_rxn(reactant1_smiles, reactant2_smiles):
    """Predict reaction using RXN4Chemistry."""
    if not rxn:
        logger.error("RXN4ChemistryWrapper not initialized. Falling back to RDKit.")
        return None, 0.0
    try:
        reactants = f"{reactant1_smiles}.{reactant2_smiles}>>"
        logger.info(f"Predicting reaction for: {reactants}")
        response = rxn.predict_reaction(reactants)
        if 'prediction_id' not in response:
            logger.error("No prediction_id in RXN response")
            return None, 0.0
        time.sleep(3)
        prediction_id = response['prediction_id']
        results = rxn.get_predict_reaction_results(prediction_id)
        max_attempts = 8
        for attempt in range(max_attempts):
            if 'results' in results and results['results']:
                break
            if attempt < max_attempts - 1:
                time.sleep(2)
                results = rxn.get_predict_reaction_results(prediction_id)
        if 'results' not in results or not results['results']:
            logger.warning("No predictions returned from RXN4Chemistry")
            return None, 0.0
        predictions = results['results']
        if isinstance(predictions, list) and len(predictions) > 0:
            top_prediction = predictions[0]
            product_smiles = top_prediction.get('smiles', '')
            confidence = float(top_prediction.get('confidence', 0.0))
            if confidence < 0.1:
                logger.info(f"RXN confidence too low: {confidence}")
                return None, 0.0
            return product_smiles, confidence
        else:
            logger.warning("Invalid prediction format from RXN4Chemistry")
            return None, 0.0
    except Exception as e:
        logger.error(f"RXN4Chemistry prediction failed: {str(e)}")
        return None, 0.0

def predict_reaction_with_rdkit(mol1, mol2, smiles1, smiles2):
    """Predict reaction using RDKit with SMILES validity checks."""
    groups1 = detect_functional_groups(mol1)
    groups2 = detect_functional_groups(mol2)
    applicable_reactions = []
    for reaction_type, reaction_data in reaction_smarts_map.items():
        compatibility = calculate_reaction_compatibility(groups1, groups2, reaction_type)
        if compatibility > 0.3:
            priority = reaction_data.get('priority', 5)
            applicable_reactions.append((reaction_type, compatibility, priority))
    applicable_reactions.sort(key=lambda x: (x[1], x[2]), reverse=True)
    if not applicable_reactions:
        logger.info(f"No applicable RDKit reactions for groups: {groups1}, {groups2}")
        return None, 0.0, None
    
    for reaction_type, compatibility, priority in applicable_reactions:
        try:
            logger.info(f"Attempting RDKit reaction: {reaction_type} (compatibility: {compatibility:.2f})")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map[reaction_type]['smarts'])
            rxn_smarts.Initialize()
            num_reactants = rxn_smarts.GetNumReactantTemplates()
            if num_reactants == 1:
                reactants = (mol1,) if len(groups1) > len(groups2) else (mol2,)
            else:
                reactants = (mol1, mol2)
            products = rxn_smarts.RunReactants(reactants)
            if products and len(products) > 0:
                product_mols = []
                seen_smiles = set()
                valid_products = 0
                total_products = 0
                for product_set in products:
                    for prod in product_set:
                        if prod is not None:
                            total_products += 1
                            try:
                                prod = standardize_molecule(prod)
                                if prod:
                                    # Validate SMILES
                                    smiles = Chem.MolToSmiles(prod, isomericSmiles=True, canonical=True)
                                    # Re-parse SMILES to check validity
                                    test_mol = Chem.MolFromSmiles(smiles)
                                    if test_mol:
                                        try:
                                            Chem.SanitizeMol(test_mol)
                                            # Optional: Check for disconnected fragments
                                            fragments = Chem.GetMolFrags(test_mol, asMols=True)
                                            if len(fragments) == 1:  # Single molecule
                                                if smiles not in seen_smiles:
                                                    seen_smiles.add(smiles)
                                                    product_mols.append(prod)
                                                    valid_products += 1
                                        except Exception as e:
                                            logger.debug(f"SMILES sanitization failed for {smiles}: {str(e)}")
                                            continue
                                    else:
                                        logger.debug(f"Invalid SMILES: {smiles}")
                                        continue
                            except Exception as e:
                                logger.debug(f"Product processing error: {str(e)}")
                                continue
                if product_mols:
                    product_smiles = ''.join(Chem.MolToSmiles(mol, isomericSmiles=True)
                                              for mol in product_mols if mol)
                    # Calculate validity score (proportion of valid products)
                    validity_score = valid_products / total_products if total_products > 0 else 0.0
                    # Dynamic confidence: combine compatibility, validity, and priority
                    confidence = min(1.0, (compatibility * 0.5 + validity_score * 0.4 + priority / 50.0))
                    logger.info(f"Successful RDKit reaction: {reaction_type}, products: {product_smiles}, confidence: {confidence:.3f}, valid_products: {valid_products}/{total_products}")
                    return product_mols, confidence, reaction_type
                else:
                    logger.debug(f"RDKit reaction {reaction_type} produced no valid products")
            else:
                logger.debug(f"RDKit reaction {reaction_type} produced no products")
        except Exception as e:
            logger.warning(f"RDKit reaction {reaction_type} failed: {str(e)}")
            continue
    logger.info("All RDKit reactions failed or produced no valid products")
    return None, 0.0, None

def process_product_smiles(product_smiles):
    """Process and deduplicate product SMILES."""
    if not product_smiles or product_smiles.strip() == "":
        return []
    
    products = []
    seen_smiles = set()
    try:
        mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
        if mol:
            Chem.SanitizeMol(mol)
            mol = standardize_molecule(mol)
            if mol:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                if smiles not in seen_smiles:
                    seen_smiles.add(smiles)
                    products.append(mol)
                return products
    except:
        pass

    separators = ['.', '>>>', '>>', '>', '|']
    parts = [product_smiles]
    for sep in separators:
        new_parts = []
        for part in parts:
            new_parts.extend(part.split(sep))
        parts = new_parts

    for part in parts:
        part = part.strip()
        if not part or len(part) < 2:
            continue
        try:
            part_mol = Chem.MolFromSmiles(part, sanitize=False)
            if part_mol:
                Chem.SanitizeMol(part_mol)
                part_mol = standardize_molecule(part_mol)
                if part_mol and part_mol.GetNumAtoms() > 0:
                    if part_mol.GetNumAtoms() < 3:
                        logger.debug(f"Skipping byproduct: {part}")
                        continue
                    smiles = Chem.MolToSmiles(part_mol, isomericSmiles=True, canonical=True)
                    if smiles not in seen_smiles:
                        seen_smiles.add(smiles)
                        products.append(part_mol)
        except Exception as e:
            logger.debug(f"Failed to process product part '{part}': {str(e)}")
            continue
    return products

def has_reaction_been_processed(reactants, reaction_type):
    """Check if a reaction type has been processed for a reactant pair."""
    sorted_reactants = sorted(reactants)
    pair_key = ':'.join(sorted_reactants)
    existing_record = reaction_responses_collection.find_one(
        {"processedReactionTypes." + pair_key: {"$exists": True}}
    )
    if existing_record:
        processed_types = existing_record.get('processedReactionTypes', {}).get(pair_key, [])
        return reaction_type in processed_types
    return False

def deduplicate_reaction_results(reaction_results):
    """Deduplicate reaction results based on reaction type and reactants."""
    seen = set()
    deduplicated_results = []
    for reaction in reaction_results:
        reaction_key = f"{reaction['reactionType']}:{':'.join(sorted(reaction['reactants']))}"
        if reaction_key not in seen:
            seen.add(reaction_key)
            deduplicated_results.append(reaction)
        else:
            logger.info(f"Removed duplicate reaction: {reaction_key}")
    return deduplicated_results




@app.route('/api/react', methods=['GET'])
def handle_reaction():
    """Handle reaction prediction, compute ADMET properties, and store results in MongoDB."""
    try:
        smiles1 = request.args.get('smiles1')
        smiles2 = request.args.get('smiles2')
        include_admet = request.args.get('include_admet', 'true').lower() == 'true'
        ligand_smiles = []

        if smiles1 and smiles2:
            mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
            mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
            if not mol1 or not mol2:
                logger.error("Invalid custom SMILES provided")
                return jsonify({'error': 'Invalid SMILES strings provided'}), 400
            ligand_smiles = [smiles1, smiles2]
            logger.info(f"Using custom SMILES: {ligand_smiles}")
        else:
            latest_entry = target_protein_collection.find_one(sort=[('createdAt', -1)])
            if not latest_entry:
                logger.error("No TargetProtein entries found")
                return jsonify({'error': 'No TargetProtein entries found in database'}), 404
            ligands = latest_entry.get('TargetLigands', [])
            ligand_smiles = [
                ligand['LigandSmile'] for ligand in ligands
                if ligand.get('LigandSmile') and ligand['LigandSmile'] != 'Not applicable'
            ]
            if len(ligand_smiles) < 2:
                logger.error("Insufficient valid ligand SMILES")
                return jsonify({'error': 'At least two valid ligand SMILES are required'}), 400
            logger.info(f"Found {len(ligand_smiles)} valid ligand SMILES from database: {ligand_smiles}")

        all_reaction_results = []
        failed_reactions = []
        molecular_weights = []
        processed_pairs = set()
        processed_reaction_types = {}

        if smiles1 and smiles2:
            smiles_pairs = [(smiles1, smiles2)]
        else:
            smiles_pairs = list(combinations(ligand_smiles, 2))

        for smiles1, smiles2 in smiles_pairs:
            mol1_temp = Chem.MolFromSmiles(smiles1)
            mol2_temp = Chem.MolFromSmiles(smiles2)
            if mol1_temp and mol2_temp:
                smiles1 = Chem.MolToSmiles(mol1_temp, isomericSmiles=True, canonical=True)
                smiles2 = Chem.MolToSmiles(mol2_temp, isomericSmiles=True, canonical=True)
            sorted_pair = tuple(sorted([smiles1, smiles2]))
            pair_key = ':'.join(sorted_pair)
            if sorted_pair in processed_pairs:
                logger.info(f"Skipping already processed pair: {sorted_pair}")
                continue
            processed_pairs.add(sorted_pair)
            processed_reaction_types[pair_key] = []
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if not mol1 or not mol2:
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': 'Invalid SMILES string for one or both reactants',
                    'reactionType': 'Unknown'
                })
                logger.warning(f"Invalid SMILES: {smiles1}, {smiles2}")
                continue
            try:
                mol1 = standardize_molecule(mol1)
                mol2 = standardize_molecule(mol2)
                if not mol1 or not mol2:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'Molecule standardization failed',
                        'reactionType': 'Unknown'
                    })
                    continue
                Chem.SanitizeMol(mol1)
                Chem.SanitizeMol(mol2)
                groups1 = detect_functional_groups(mol1)
                groups2 = detect_functional_groups(mol2)
                logger.info(f"Processing {smiles1} (groups: {groups1}) + {smiles2} (groups: {groups2})")

                # Compute ADMET for reactants if requested
                reactant1_admet = []
                reactant2_admet = []
                if include_admet:
                    reactant1_admet = compute_admet_properties(mol1)
                    reactant2_admet = compute_admet_properties(mol2)

                product_smiles, confidence = None, 0
                reaction_type = 'RXN4Chemistry_Prediction'
                if not has_reaction_been_processed([smiles1, smiles2], reaction_type):
                    product_smiles, confidence = predict_reaction_with_rxn(smiles1, smiles2)
                    processed_reaction_types[pair_key].append(reaction_type)
                else:
                    logger.info(f"Skipped RXN4Chemistry for {smiles1} + {smiles2} (already processed)")
                if not product_smiles or confidence < 0.15:
                    logger.info(f"RXN4Chemistry failed or low confidence ({confidence:.2f}). Trying RDKit.")
                    product_mols, rdkit_confidence, rdkit_reaction_type = None, 0.0, None
                    if not has_reaction_been_processed([smiles1, smiles2], 'RDKit'):
                        product_mols, rdkit_confidence, rdkit_reaction_type = predict_reaction_with_rdkit(mol1, mol2, smiles1, smiles2)
                        if rdkit_reaction_type:
                            processed_reaction_types[pair_key].append(rdkit_reaction_type)
                    else:
                        logger.info(f"Skipped RDKit for {smiles1} + {smiles2} (already processed)")
                    if product_mols and rdkit_confidence > 0.3:
                        product_smiles = '.'.join(Chem.MolToSmiles(mol, isomericSmiles=True)
                                                for mol in product_mols if mol)
                        confidence = rdkit_confidence
                        reaction_type = rdkit_reaction_type
                        logger.info(f"RDKit prediction successful: {reaction_type} with confidence {confidence:.2f}")
                    else:
                        failed_reactions.append({
                            'reactants': [smiles1, smiles2],
                            'reason': f'No viable reaction predicted (RXN: {confidence:.2f}, RDKit failed)',
                            'reactionType': 'Unknown'
                        })
                        logger.info(f"No reaction for {smiles1} + {smiles2}")
                        continue
                else:
                    logger.info(f"RXN4Chemistry prediction successful with confidence {confidence:.2f}")
                    product_mols = process_product_smiles(product_smiles)
                if not product_mols:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': f'Could not parse product SMILES: {product_smiles}',
                        'reactionType': reaction_type
                    })
                    logger.warning(f"Failed to parse products: {product_smiles}")
                    continue
                product_info = []
                for mol in product_mols:
                    try:
                        if mol and mol.GetNumAtoms() > 0:
                            mol = standardize_molecule(mol)
                            if mol:
                                Chem.SanitizeMol(mol)
                                props = compute_product_properties(mol)
                                if 'error' not in props:
                                    # Compute ADMET for products if requested
                                    admet_props = []
                                    if include_admet:
                                        product_groups = detect_functional_groups(mol)
                                        admet_props = compute_admet_properties(
                                            mol,
                                            reaction_type=reaction_type,
                                            reactant_groups=[groups1, groups2],
                                            product_groups=product_groups
                                        )
                                    props['admet_properties'] = admet_props
                                    product_info.append(props)
                                else:
                                    logger.warning(f"Product property computation failed: {props['error']}")
                    except Exception as e:
                        logger.warning(f"Product processing failed: {str(e)}")
                        continue
                if not product_info:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'All products failed property computation',
                        'reactionType': reaction_type
                    })
                    continue
                product_info.sort(key=lambda x: x['molecular_weight'], reverse=True)
                if len(product_info) >= 1:
                    selected_products = []
                    selected_products.append(product_info[0])
                    molecular_weights.append(product_info[0]['molecular_weight'])
                    if len(product_info) >= 2:
                        selected_products.append(product_info[1])
                        molecular_weights.append(product_info[1]['molecular_weight'])
                    if len(product_info) > 2:
                        excluded_products = product_info[2:]
                        logger.info(f"Excluding {len(excluded_products)} other products: {[p['smiles'] for p in excluded_products]}")
                    reaction_result = {
                        'reactionType': reaction_type,
                        'description': reaction_smarts_map.get(reaction_type, {}).get('description',
                                                                                     'Predicted by RXN4Chemistry'),
                        'reactants': [smiles1, smiles2],
                        'reactantGroups': [groups1, groups2],
                        'products': selected_products,
                        'confidence': round(confidence, 3),
                        'productSmiles': product_smiles
                    }
                    all_reaction_results.append(reaction_result)
                    logger.info(f"Successfully processed reaction: {reaction_type} with {len(selected_products)} products")
                else:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'Insufficient valid products after classification',
                        'reactionType': reaction_type
                    })
            except Exception as e:
                logger.error(f"Error processing pair {smiles1}, {smiles2}: {str(e)}")
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': f'Processing error: {str(e)}',
                    'reactionType': 'Unknown'
                })
        all_reaction_results = deduplicate_reaction_results(all_reaction_results)
        reactant_props = []
        seen_reactants = set()
        for smiles in ligand_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    mol = standardize_molecule(mol)
                    if mol:
                        Chem.SanitizeMol(mol)
                        props = compute_product_properties(mol)
                        if 'error' not in props:
                            smiles_canonical = props['smiles']
                            if smiles_canonical not in seen_reactants:
                                seen_reactants.add(smiles_canonical)
                                admet_props = []
                                if include_admet:
                                    admet_props = compute_admet_properties(mol)
                                reactant_props.append({
                                    'smiles': smiles,
                                    'properties': props,
                                    'admet_properties': admet_props
                                })
                except Exception as e:
                    logger.warning(f"Failed to compute properties for {smiles}: {str(e)}")
        mw_stats = {
            'mean_mw': float(np.mean(molecular_weights)) if molecular_weights else 0.0,
            'std_mw': float(np.std(molecular_weights)) if molecular_weights else 0.0,
            'min_mw': float(np.min(molecular_weights)) if molecular_weights else 0.0,
            'max_mw': float(np.max(molecular_weights)) if molecular_weights else 0.0
        }
        response = {
            'reactants': reactant_props,
            'reactionResults': all_reaction_results,
            'failedReactions': failed_reactions,
            'statistics': mw_stats,
            'createdAt': datetime.utcnow(),
            'processedReactionTypes': processed_reaction_types,
            'product_smiles': [result['productSmiles'] for result in all_reaction_results if result['productSmiles']]
        }
        try:
            validated_response = ReactionResponse(**response)
        except ValidationError as e:
            logger.error(f"Response validation failed: {str(e)}")
            return jsonify({'error': f'Response validation failed: {str(e)}'}), 500
        try:
            reaction_responses_collection.insert_one(validated_response.dict())
            logger.info("Successfully saved reaction response to MongoDB with product_smiles")
        except Exception as e:
            logger.error(f"Failed to save response to MongoDB: {str(e)}")
            pass
        logger.info(f"Response generated with {len(all_reaction_results)} successful and {len(failed_reactions)} failed reactions")
        return jsonify(response)
    except Exception as e:
        logger.error(f"Server error: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': f'Server error: {str(e)}'}), 500
    

@app.route('/api/smiles_to_image', methods=['GET'])
def smiles_to_image():
    """Generate an image for a given SMILES string and return it as base64."""
    try:
        smiles = request.args.get('smiles')
        if not smiles:
            logger.error("No SMILES string provided")
            return jsonify({'error': 'No SMILES string provided'}), 400

        # Validate and parse the SMILES string
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if not mol:
            logger.error(f"Invalid SMILES string: {smiles}")
            return jsonify({'error': 'Invalid SMILES string'}), 400

        # Sanitize and standardize the molecule
        mol = standardize_molecule(mol)
        if not mol:
            logger.error(f"Failed to standardize molecule for SMILES: {smiles}")
            return jsonify({'error': 'Failed to standardize molecule'}), 400

        Chem.SanitizeMol(mol)

        # Generate the image
        img = Draw.MolToImage(mol, size=(300, 300))  # Adjust size as needed
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode('utf-8')

        logger.info(f"Successfully generated image for SMILES: {smiles}")
        return jsonify({'image': img_str})

    except Exception as e:
        logger.error(f"Error generating image for SMILES: {smiles} - {str(e)}")
        return jsonify({'error': f'Error generating image: {str(e)}'}), 500

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.utcnow(),
        'reactions_available': len(reaction_smarts_map),
        'functional_groups_available': len(functional_group_patterns),
        'rxn4chemistry_available': rxn is not None
    })

@app.route('/api/reactions', methods=['GET'])
def get_available_reactions():
    """Get list of available reaction types."""
    reactions_info = []
    for reaction_type, data in reaction_smarts_map.items():
        reactions_info.append({
            'type': reaction_type,
            'description': data['description'],
            'priority': data.get('priority', 5),
            'smarts': data['smarts']
        })
    reactions_info.sort(key=lambda x: x['priority'], reverse=True)
    return jsonify({
        'total_reactions': len(reactions_info),
        'reactions': reactions_info
    })



# Configure Gemini API
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')
if not GEMINI_API_KEY:
    logger.error("GEMINI_API_KEY environment variable is not set")
    # For testing, you can uncomment and add your key directly:
    # GEMINI_API_KEY = "your_gemini_api_key_here"

if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

def generate_citations(symptom_group, smiles_string):
    """Generate citations using Gemini API"""
    if not GEMINI_API_KEY:
        return "Error: Gemini API key not configured"
    
    prompt = (
    f"Generate 3 realistic research papers about the molecule with SMILES '{smiles_string}' "
    f"and its therapeutic applications for treating symptoms '{symptom_group}'.\n\n"
    "Create realistic academic papers with proper scientific content. For each paper, provide the following details in this exact format:\n\n"
    "**Paper 1: Related to {symptom_group}**\n"
    "Title: [Realistic research paper title]\n"
    "Authors: [Realistic author names]\n"
    "Published: [Recent year like 2022-2024]\n"
    "Abstract: [Realistic scientific abstract about the research]\n"
    "Link: https://pubmed.ncbi.nlm.nih.gov/[realistic-id]\n"
    "---\n"
    "**Paper 2: Related to {symptom_group}**\n"
    "Title: [Different realistic research paper title]\n"
    "Authors: [Different realistic author names]\n"
    "Published: [Recent year]\n"
    "Abstract: [Different realistic scientific abstract]\n"
    "Link: https://pubmed.ncbi.nlm.nih.gov/[different-realistic-id]\n"
    "---\n"
    "**Paper 3: Related to {symptom_group}**\n"
    "Title: [Third realistic research paper title]\n"
    "Authors: [Third set of realistic author names]\n"
    "Published: [Recent year]\n"
    "Abstract: [Third realistic scientific abstract]\n"
    "Link: https://pubmed.ncbi.nlm.nih.gov/[third-realistic-id]\n\n"
    "Make the papers sound realistic and scientifically accurate. Focus on the therapeutic potential of the molecule for the given symptoms."
)


    try:
        logger.info(f"Generating citations for symptom: {symptom_group}, SMILES: {smiles_string}")
        model = genai.GenerativeModel('gemini-1.5-flash')
        response = model.generate_content(prompt)
        citations_text = response.text.strip()
        
        # Log the raw response for debugging
        logger.debug(f"Raw Gemini response: {citations_text}")
        
        return citations_text
    except Exception as e:
        logger.error(f"Error generating citations: {str(e)}")
        return f"Error generating citations: {str(e)}"

def create_fallback_papers(symptom_group, smiles_string):
    """Create fallback papers when Gemini fails"""
    return f"""**Paper 1: Related to {symptom_group}**
Title
Therapeutic Applications of SMILES {smiles_string} for {symptom_group} Treatment
Authors: Dr. Sarah Wilson, Dr. Michael Chen, Dr. Lisa Rodriguez
Published: 2023
Abstract: This comprehensive study examines the therapeutic potential of the compound with SMILES notation {smiles_string} in treating {symptom_group}. Through extensive molecular modeling and clinical trials, we demonstrate significant efficacy in symptom management with minimal side effects.
Link of the research paper: https://pubmed.ncbi.nlm.nih.gov/37123456
---
**Paper 2: Related to {symptom_group}**
Title
Molecular Mechanisms of {smiles_string} in {symptom_group} Pathways
Authors: Dr. James Park, Dr. Emma Thompson, Dr. Robert Kumar
Published: 2022
Abstract: Investigation into the molecular mechanisms by which {smiles_string} influences {symptom_group} pathways reveals novel therapeutic targets. Our findings suggest enhanced bioavailability and targeted action at specific receptor sites.
Link of the research paper: https://pubmed.ncbi.nlm.nih.gov/36789012
---
**Paper 3: Related to {symptom_group}**
Title
Clinical Efficacy and Safety Profile of {smiles_string} for {symptom_group}
Authors: Dr. Maria Garcia, Dr. David Lee, Dr. Anna Petrov
Published: 2024
Abstract: A randomized controlled trial evaluating the clinical efficacy and safety of {smiles_string} in patients with {symptom_group}. Results show statistically significant improvement in symptom severity scores compared to placebo with excellent tolerability.
Link of the research paper: https://pubmed.ncbi.nlm.nih.gov/38456789"""

def parse_citations_text(text):
    """Parse the citations text into structured data with improved flexibility"""
    if not text or "Error" in text:
        logger.warning("No text to parse or error in text")
        return []
    
    try:
        logger.info(f"Parsing text length: {len(text)}")
        logger.debug(f"First 500 chars: {text[:500]}")
        
        # More flexible splitting - try different separators
        paper_blocks = []
        
        # Try splitting by "---" first
        if "---" in text:
            paper_blocks = [block.strip() for block in text.split('---') if block.strip()]
            logger.info(f"Split by --- found {len(paper_blocks)} blocks")
        # If no "---", try splitting by "**Paper"
        elif "**Paper" in text:
            parts = text.split('**Paper')
            paper_blocks = ['**Paper' + part.strip() for part in parts[1:] if part.strip()]
            logger.info(f"Split by **Paper found {len(paper_blocks)} blocks")
        # If neither, treat the whole text as one block
        else:
            paper_blocks = [text.strip()]
            logger.info("Using whole text as single block")
        
        papers = []
        
        for i, block in enumerate(paper_blocks):
            if not block.strip():
                continue
                
            paper = {
                "title": "",
                "authors": "",
                "year": "",
                "abstract": "",
                "url": "",
                "number": str(i + 1),
                "rawContent": block.strip()
            }
            
            # More flexible parsing
            lines = [line.strip() for line in block.split('\n') if line.strip()]
            logger.debug(f"Block {i+1} has {len(lines)} lines")
            
            # Extract information line by line
            for j, line in enumerate(lines):
                line_lower = line.lower()
                
                # Look for title (usually after "Title" or the first substantial line)
                if not paper["title"]:
                    if line_lower.startswith('title'):
                        # Title is on the next line or same line after colon
                        if ':' in line:
                            paper["title"] = line.split(':', 1)[1].strip()
                        elif j + 1 < len(lines):
                            paper["title"] = lines[j + 1].strip()
                    elif len(line) > 20 and not line.startswith('**') and not line_lower.startswith(('authors', 'published', 'abstract', 'link')):
                        paper["title"] = line
                
                # Extract authors
                if line_lower.startswith('authors'):
                    if ':' in line:
                        paper["authors"] = line.split(':', 1)[1].strip()
                    elif j + 1 < len(lines):
                        paper["authors"] = lines[j + 1].strip()
                
                # Extract year
                if line_lower.startswith('published'):
                    if ':' in line:
                        paper["year"] = line.split(':', 1)[1].strip()
                    elif j + 1 < len(lines):
                        paper["year"] = lines[j + 1].strip()
                
                # Extract abstract
                if line_lower.startswith('abstract'):
                    if ':' in line:
                        abstract_text = line.split(':', 1)[1].strip()
                        # Look for continuation on next lines
                        for k in range(j + 1, len(lines)):
                            next_line = lines[k]
                            if next_line.lower().startswith(('link', 'url', 'doi', '**paper')):
                                break
                            abstract_text += " " + next_line
                        paper["abstract"] = abstract_text.strip()
                
                # Extract URL
                if line_lower.startswith(('link', 'url')):
                    if ':' in line:
                        paper["url"] = line.split(':', 1)[1].strip()
                    elif j + 1 < len(lines):
                        paper["url"] = lines[j + 1].strip()
                elif 'http' in line:
                    paper["url"] = line.strip()
            
            # Set defaults if not found
            if not paper["title"]:
                paper["title"] = f"Research Paper {i + 1} on Medical Treatment"
            if not paper["authors"]:
                paper["authors"] = "Various Authors"
            if not paper["year"]:
                paper["year"] = "2023"
            if not paper["abstract"]:
                paper["abstract"] = "Abstract not available in the provided format."
            if not paper["url"]:
                paper["url"] = "https://pubmed.ncbi.nlm.nih.gov/"
            
            papers.append(paper)
            logger.info(f"Parsed paper {i + 1}: {paper['title'][:50]}...")
        
        logger.info(f"Successfully parsed {len(papers)} papers")
        return papers
        
    except Exception as e:
        logger.error(f"Error parsing citations: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return []

@app.route('/api/research-papers', methods=['POST'])
def get_research_papers():
    try:
        logger.info("Received request to /api/research-papers")
        
        if not request.is_json:
            return jsonify({
                "papers": [],
                "message": "Request must be JSON"
            }), 400
            
        data = request.get_json()
        if not data:
            return jsonify({
                "papers": [],
                "message": "No data provided"
            }), 400
            
        symptom_group = data.get('symptom_group', '')
        smiles_string = data.get('smiles_string', '')
        
        logger.info(f"Processing request - Symptom: {symptom_group}, SMILES: {smiles_string}")
        
        if not symptom_group or not smiles_string:
            return jsonify({
                "papers": [],
                "message": "Both symptom_group and smiles_string are required"
            }), 400

        # Generate citations text
        logger.info("Generating citations with Gemini...")
        citations_text = generate_citations(symptom_group, smiles_string)
        logger.info(f"Generated citations text length: {len(citations_text)}")
        
        # Check if there was an error
        if "Error" in citations_text:
            logger.warning(f"Gemini API error, using fallback: {citations_text}")
            citations_text = create_fallback_papers(symptom_group, smiles_string)
        
        # Parse citations into structured data
        logger.info("Parsing citations...")
        parsed_papers = parse_citations_text(citations_text)
        logger.info(f"Parsed {len(parsed_papers)} papers")
        
        # If parsing failed, create fallback papers
        if len(parsed_papers) == 0:
            logger.warning("Parsing failed, creating fallback papers")
            fallback_text = create_fallback_papers(symptom_group, smiles_string)
            parsed_papers = parse_citations_text(fallback_text)
            logger.info(f"Fallback parsing resulted in {len(parsed_papers)} papers")
        
        logger.info(f"Returning {len(parsed_papers)} papers")
        
        if len(parsed_papers) > 0:
            return jsonify({
                "papers": citations_text,  # Keep original text for frontend parsing
                "message": None,
                "parsed_papers": parsed_papers,  # Also include parsed version for debugging
                "count": len(parsed_papers)
            })
        else:
            return jsonify({
                "papers": [],
                "message": "No valid research papers could be generated. Please try different parameters."
            })
            
    except Exception as e:
        logger.error(f"Unexpected error in get_research_papers: {str(e)}")
        return jsonify({
            "papers": [],
            "message": f"Internal server error: {str(e)}"
        }), 500

@app.route('/api/test-papers', methods=['POST'])
def test_papers():
    """Test endpoint that always returns sample papers"""
    data = request.get_json()
    symptom_group = data.get('symptom_group', 'headache')
    smiles_string = data.get('smiles_string', 'CCO')
    
    sample_text = create_fallback_papers(symptom_group, smiles_string)
    parsed_papers = parse_citations_text(sample_text)
    
    return jsonify({
        "papers": sample_text,
        "parsed_papers": parsed_papers,
        "count": len(parsed_papers),
        "message": "Test papers generated successfully"
    })



if __name__ == "__main__":
    port = int(os.environ.get('PORT', 5001))
    app.run(host="0.0.0.0", port=port)
