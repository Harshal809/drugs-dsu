
from rdkit import Chem
from rdkit.Chem import rdChemReactions, MolStandardize
import logging
from modules.fun_grp_pattern import functional_group_patterns

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def detect_functional_groups(mol):
    """Detect functional groups in a molecule using substructure matching."""
    if not mol:
        logger.warning("Invalid molecule provided for functional group detection")
        return []
    try:
        Chem.SanitizeMol(mol)
        groups = []
        for group_name, pattern in functional_group_patterns.items():
            if pattern and mol.HasSubstructMatch(pattern):
                groups.append(group_name)
        logger.debug(f"Detected functional groups: {groups}")
        return sorted(groups)  # Sort for consistent output
    except Exception as e:
        logger.error(f"Functional group detection failed: {str(e)}")
        return []

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
