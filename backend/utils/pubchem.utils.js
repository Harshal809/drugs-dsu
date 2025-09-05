import axios from 'axios';
import pubchem from 'pubchem';

// Cache for storing already fetched data
const pubchemCache = new Map();

/**
 * Fetches PubChem data for a given SMILES string
 * @param {string} smiles - The SMILES string to query
 * @param {boolean} [useCache=true] - Whether to use cached results
 * @returns {Promise<Object>} Molecular data from PubChem
 */
export const getPubChemData = async (smiles, useCache = true) => {
    if (!smiles || typeof smiles !== 'string') {
        throw new Error('Invalid SMILES string');
    }

    // Check cache first
    if (useCache && pubchemCache.has(smiles)) {
        return pubchemCache.get(smiles);
    }

    try {
        // Step 1: Convert SMILES to CID (PubChem Compound ID)
        const compound = await pubchem.lookup(smiles, 'smiles');
        const cid = compound.cid;
        if (!cid) {
            throw new Error('SMILES not found in PubChem');
        }

        // Step 2: Fetch compound properties
        const [properties, synonyms, descriptors] = await Promise.all([
            axios.get(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/MolecularWeight,CanonicalSMILES,InChIKey,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,ExactMass,IsomericSMILES/JSON`),
            axios.get(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/synonyms/JSON`),
            axios.get(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/description/JSON`)
        ]);

        // Extract relevant data
        const result = {
            cid,
            molecularWeight: properties.data.PropertyTable.Properties[0]?.MolecularWeight || 0,
            canonicalSMILES: properties.data.PropertyTable.Properties[0]?.CanonicalSMILES || smiles,
            inchiKey: properties.data.PropertyTable.Properties[0]?.InChIKey || '',
            hBondDonors: properties.data.PropertyTable.Properties[0]?.HBondDonorCount || 0,
            hBondAcceptors: properties.data.PropertyTable.Properties[0]?.HBondAcceptorCount || 0,
            rotatableBonds: properties.data.PropertyTable.Properties[0]?.RotatableBondCount || 0,
            exactMass: properties.data.PropertyTable.Properties[0]?.ExactMass || 0,
            iupacName: descriptors.data.InformationList?.Information[0]?.Description || '',
            synonyms: synonyms.data.InformationList?.Information[0]?.Synonym || [],
            isomericSMILES: properties.data.PropertyTable.Properties[0]?.IsomericSMILES || smiles,
            atomCount: countAtoms(properties.data.PropertyTable.Properties[0]?.CanonicalSMILES || smiles)
        };

        // Cache the result
        pubchemCache.set(smiles, result);
        return result;

    } catch (error) {
        console.error(`PubChem lookup failed for ${smiles}:`, error.message);

        // Fallback: Basic molecular weight calculation from SMILES
        return {
            cid: null,
            molecularWeight: estimateMolecularWeight(smiles),
            canonicalSMILES: smiles,
            inchiKey: '',
            hBondDonors: 0,
            hBondAcceptors: 0,
            rotatableBonds: 0,
            exactMass: 0,
            iupacName: '',
            synonyms: [],
            isomericSMILES: smiles,
            atomCount: countAtoms(smiles)
        };
    }
};

// Helper function to estimate molecular weight from SMILES
const estimateMolecularWeight = (smiles) => {
    if (!smiles) return 0;

    // Very rough estimation (12g/mol per non-hydrogen atom)
    const heavyAtoms = smiles.replace(/[\[\]0-9\+\-Hh\s]/g, '').length;
    return Math.max(50, heavyAtoms * 12);
};

// Count non-hydrogen atoms in SMILES
const countAtoms = (smiles) => {
    if (!smiles) return 0;

    // Remove all non-atom characters
    const atomString = smiles.replace(/[\[\]0-9\+\-=\#\(\)\\\/%@Hh\s]/g, '');

    // Count atoms (handling two-letter elements like Cl, Br)
    let count = 0;
    let i = 0;
    while (i < atomString.length) {
        // Check for two-letter elements
        if (i < atomString.length - 1 &&
            /[A-Z][a-z]/.test(atomString.substring(i, i + 2))) {
            i += 2;
        } else {
            i += 1;
        }
        count += 1;
    }

    return count;
};

// Additional utility functions
export const PubChemService = {
    getPubChemData,
    clearCache: () => pubchemCache.clear(),
    getCacheSize: () => pubchemCache.size
};