import { getRDKit } from '../rdkit.js';

export const handleReaction = async (req, res) => {
  const { smiles1, smiles2, reactionTypes } = req.body;
  const reactionSMARTSMap = {
    Condensation: '[C:1][OH:2].[C:3][OH:4]>>[C:1][O:2][C:3].[O:4]', // Two alcohols form ether + water
    Esterification: '[C:1][OH:2].[C:3](=O)[OH:4]>>[C:1][O:2][C:3](=O).[O:4]', // Alcohol + carboxylic acid form ester + water
    AmideFormation: '[C:1][NH2:2].[C:3](=O)[OH:4]>>[C:1][NH:2][C:3](=O).[O:4]' // Amine + carboxylic acid form amide + water
  };

  let mol1 = null;
  let mol2 = null;
  let reaction = null;
  let molList = null;
  let products = null;
  let currentReactants = [];
  let intermediateProducts = [];
  const deletedObjects = new Set();

  try {
    const rdkit = getRDKit();
    console.log('RDKit instance loaded:', !!rdkit, 'Methods:', Object.keys(rdkit || {}));

    // Validate and initialize molecules
    mol1 = rdkit.get_mol(smiles1);
    mol2 = rdkit.get_mol(smiles2);
    console.log('mol1:', mol1, 'valid:', mol1?.is_valid(), 'SMILES:', smiles1);
    console.log('mol2:', mol2, 'valid:', mol2?.is_valid(), 'SMILES:', smiles2);

    if (!mol1 || !mol1.is_valid() || !mol2 || !mol2.is_valid()) {
      return res.status(400).json({ error: 'Invalid SMILES strings' });
    }

    currentReactants = [mol1, mol2];
    let allProductSMILES = [];
    let finalProducts = [];

    // Handle multiple reactions sequentially
    for (const reactionType of reactionTypes) {
      if (!reactionSMARTSMap[reactionType]) {
        return res.status(400).json({ error: `Invalid reaction type: ${reactionType}` });
      }

      const reactionSMARTS = reactionSMARTSMap[reactionType];
      reaction = rdkit.get_rxn(reactionSMARTS);
      console.log(`Reaction created for ${reactionType}:`, !!reaction, 'Reaction:', reaction);

      if (!reaction) {
        return res.status(400).json({ error: `Invalid reaction SMARTS for ${reactionType}` });
      }

      // Prepare reactants for this step
      molList = new rdkit.MolList();
      currentReactants.forEach(mol => {
        if (mol && !deletedObjects.has(mol) && typeof mol.is_valid === 'function' && mol.is_valid()) {
          molList.append(mol);
        }
      });

      // Run the reaction
      products = reaction.run_reactants(molList);
      console.log(`Products for ${reactionType}:`, products, 'Size:', products?.size ? products.size() : 'N/A');
      console.log(`Products methods:`, Object.keys(products || {}));
      console.log(`Products full structure:`, JSON.stringify(products, null, 2));
      if (!products || !products.size || products.size() === 0) {
        console.warn(`No products formed for ${reactionType} with reactants: ${smiles1}, ${smiles2}`);
        return res.status(400).json({ error: `No products formed for ${reactionType}` });
      }

      // Process products with extreme fallbacks
      let stepProductSMILES = [];
      const size = products.size ? products.size() : 0;
      console.log(`Products size for ${reactionType}: ${size}`);
      for (let i = 0; i < size; i++) {
        let productSet = null;
        try {
          if (typeof products.get === 'function') {
            productSet = products.get(i);
          } else {
            console.warn(`products.get is not a function for ${reactionType}`);
            // Fallback 1: Treat products as a molecule
            if (products && typeof products.is_valid === 'function' && products.is_valid()) {
              const smiles = products.get_smiles();
              console.log(`Fallback 1 - Direct SMILES from products for ${reactionType}: ${smiles}`);
              stepProductSMILES.push(smiles);
              const newMol = rdkit.get_mol(smiles);
              if (newMol && newMol.is_valid()) intermediateProducts.push(newMol);
              continue;
            }
            continue;
          }

          if (!productSet) {
            console.warn(`No product set at index ${i} for ${reactionType}`);
            continue;
          }

          console.log(`Product set for ${reactionType} at index ${i}:`, productSet);
          console.log(`Product set methods:`, Object.keys(productSet || {}));
          console.log(`Product set full structure:`, JSON.stringify(productSet, null, 2));

          // Attempt to access molecules
          let molecules = [];
          if (typeof productSet.size === 'function' && productSet.size() > 0) {
            const setSize = productSet.size();
            console.log(`Product set size for ${reactionType} at index ${i}: ${setSize}`);
            for (let j = 0; j < setSize; j++) {
              let mol = null;
              try {
                if (typeof productSet.get === 'function') {
                  mol = productSet.get(j);
                } else {
                  console.warn(`productSet.get is not a function for ${reactionType} at index ${j}`);
                  // Fallback 2: Treat productSet as a molecule
                  if (productSet && typeof productSet.is_valid === 'function' && productSet.is_valid()) {
                    mol = productSet;
                    j = setSize; // Exit loop
                  } else {
                    // Fallback 3: Check for any SMILES-like property
                    if (productSet && typeof productSet.toString === 'function') {
                      const possibleSmiles = productSet.toString();
                      console.log(`Fallback 3 - toString attempt for ${reactionType} at index ${j}: ${possibleSmiles}`);
                      const testMol = rdkit.get_mol(possibleSmiles);
                      if (testMol && testMol.is_valid()) {
                        mol = testMol;
                        molecules.push({ mol, smiles: possibleSmiles });
                      }
                    }
                    continue;
                  }
                }
                if (mol && typeof mol.is_valid === 'function' && mol.is_valid()) {
                  const smiles = mol.get_smiles();
                  console.log(`Generated SMILES for ${reactionType} at index ${i},${j}: ${smiles}`);
                  molecules.push({ mol, smiles });
                } else {
                  console.warn(`Invalid or null molecule at index ${j} for ${reactionType}`);
                }
              } catch (innerErr) {
                console.warn(`Error processing molecule at index ${j} for ${reactionType}:`, innerErr.message);
              }
            }
          } else {
            // Fallback 2: Treat productSet as a single molecule
            console.log(`Attempting to treat productSet as a molecule for ${reactionType}`);
            if (productSet && typeof productSet.is_valid === 'function' && productSet.is_valid()) {
              const smiles = productSet.get_smiles();
              console.log(`Fallback 2 - Generated SMILES for ${reactionType} at index ${i} (direct): ${smiles}`);
              molecules.push({ mol: productSet, smiles });
            } else {
              // Fallback 3: Check for any SMILES-like property
              console.log(`Fallback 3 - Checking productSet for SMILES-like data for ${reactionType}`);
              if (productSet && typeof productSet.toString === 'function') {
                const possibleSmiles = productSet.toString();
                console.log(`Fallback 3 - toString result for ${reactionType} at index ${i}: ${possibleSmiles}`);
                const testMol = rdkit.get_mol(possibleSmiles);
                if (testMol && testMol.is_valid()) {
                  molecules.push({ mol: testMol, smiles: possibleSmiles });
                } else {
                  console.warn(`productSet is not a valid molecule or SMILES for ${reactionType} at index ${i}`);
                }
              } else {
                console.warn(`No valid fallback for productSet for ${reactionType} at index ${i}`);
              }
            }
          }

          // Process valid molecules
          molecules.forEach(({ mol, smiles }) => {
            stepProductSMILES.push(smiles);
            const newMol = rdkit.get_mol(smiles);
            if (newMol && newMol.is_valid()) {
              intermediateProducts.push(newMol);
            } else {
              console.warn(`Invalid molecule generated for ${reactionType}: ${smiles}`);
            }
            // Delete the original molecule if it’s not the productSet itself
            if (mol !== productSet && mol && typeof mol.delete === 'function' && !deletedObjects.has(mol)) {
              mol.delete();
              deletedObjects.add(mol);
            }
          });
        } catch (setErr) {
          console.warn(`Error processing product set at index ${i} for ${reactionType}:`, setErr.message);
          continue;
        } finally {
          if (productSet && typeof productSet.delete === 'function' && !deletedObjects.has(productSet)) {
            productSet.delete();
            deletedObjects.add(productSet);
          }
        }
      }

      allProductSMILES.push({ reactionType, products: stepProductSMILES });

      // Clean up for this step
      if (molList && typeof molList.delete === 'function' && !deletedObjects.has(molList)) {
        molList.delete();
        deletedObjects.add(molList);
      }
      if (reaction && typeof reaction.delete === 'function' && !deletedObjects.has(reaction)) {
        reaction.delete();
        deletedObjects.add(reaction);
      }
      if (products && typeof products.delete === 'function' && !deletedObjects.has(products)) {
        products.delete();
        deletedObjects.add(products);
      }

      // Prepare for the next reaction: use products as reactants
      currentReactants.forEach(mol => {
        if (mol && typeof mol.delete === 'function' && !deletedObjects.has(mol)) {
          mol.delete();
          deletedObjects.add(mol);
        }
      });
      currentReactants = intermediateProducts;
      intermediateProducts = [];

      if (stepProductSMILES.length === 0) {
        console.warn(`No valid products formed for ${reactionType}`);
        return res.status(400).json({ error: `No valid products formed for ${reactionType}` });
      }
    }

    // Collect final products
    finalProducts = currentReactants
      .map(mol => {
        if (mol && typeof mol.is_valid === 'function' && mol.is_valid() && typeof mol.get_smiles === 'function' && !deletedObjects.has(mol)) {
          return mol.get_smiles();
        }
        return null;
      })
      .filter(smiles => smiles !== null);

    if (finalProducts.length === 0) {
      console.warn('No valid final products formed after all reactions');
      return res.status(400).json({ error: 'No valid final products formed' });
    }

    return res.json({
      reactants: [smiles1, smiles2],
      reactionSteps: allProductSMILES,
      finalProducts,
    });
  } catch (error) {
    console.error('Reaction error:', error);
    return res.status(500).json({ error: `Error processing reaction: ${error.message}` });
  } finally {
    // Clean up all RDKit objects safely
    if (mol1 && typeof mol1.delete === 'function' && !deletedObjects.has(mol1)) {
      mol1.delete();
      deletedObjects.add(mol1);
    }
    if (mol2 && typeof mol2.delete === 'function' && !deletedObjects.has(mol2)) {
      mol2.delete();
      deletedObjects.add(mol2);
    }
    if (reaction && typeof reaction.delete === 'function' && !deletedObjects.has(reaction)) {
      reaction.delete();
      deletedObjects.add(reaction);
    }
    if (molList && typeof molList.delete === 'function' && !deletedObjects.has(molList)) {
      molList.delete();
      deletedObjects.add(molList);
    }
    if (products && typeof products.delete === 'function' && !deletedObjects.has(products)) {
      products.delete();
      deletedObjects.add(products);
    }
    if (currentReactants && Array.isArray(currentReactants)) {
      currentReactants.forEach(mol => {
        if (mol && typeof mol.delete === 'function' && !deletedObjects.has(mol)) {
          mol.delete();
          deletedObjects.add(mol);
        }
      });
    }
    if (intermediateProducts && Array.isArray(intermediateProducts)) {
      intermediateProducts.forEach(mol => {
        if (mol && typeof mol.delete === 'function' && !deletedObjects.has(mol)) {
          mol.delete();
          deletedObjects.add(mol);
        }
      });
    }
  }
};

// Test cases for successful reactions
const testCases = [
  {
    description: 'Condensation: Two alcohols form an ether and water',
    smiles1: 'CCO', // Ethanol
    smiles2: 'CCO', // Ethanol
    reactionTypes: ['Condensation'],
    expectedProducts: ['CCOCC', 'O'] // Diethyl ether + water
  },
  {
    description: 'Esterification: Alcohol and carboxylic acid form an ester and water',
    smiles1: 'CCO', // Ethanol
    smiles2: 'CC(=O)O', // Acetic acid
    reactionTypes: ['Esterification'],
    expectedProducts: ['CCOC(=O)C', 'O'] // Ethyl acetate + water
  },
  {
    description: 'AmideFormation: Amine and carboxylic acid form an amide and water',
    smiles1: 'CCN', // Ethylamine
    smiles2: 'CC(=O)O', // Acetic acid
    reactionTypes: ['AmideFormation'],
    expectedProducts: ['CCNC(=O)C', 'O'] // N-Ethylethanamide + water
  }
];

// Log test cases for reference
console.log('Test Cases for Successful Reactions:');
testCases.forEach(test => {
  console.log(`- ${test.description}`);
  console.log(`  SMILES 1: ${test.smiles1}`);
  console.log(`  SMILES 2: ${test.smiles2}`);
  console.log(`  Reaction Types: ${test.reactionTypes.join(', ')}`);
  console.log(`  Expected Products: ${test.expectedProducts.join(', ')}`);
});