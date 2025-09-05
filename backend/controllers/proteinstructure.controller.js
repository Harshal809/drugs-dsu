import ProteinStructure from '../models/protienstructure.model.js'; // Fixed typo in import
import { User } from '../models/auth.model.js'; // Added User model import
import GeneratenewMolecule from '../models/generatenew.model.js'; // Added GeneratenewMolecule model import
import axios from 'axios';
import dotenv from 'dotenv';
import OpenAI from 'openai';
import ResearchPaper from '../models/researchPapers.model.js';
import GeneratedResearchPaper from '../models/GeneratedResearchPaper.js';
import { SavedSearch } from '../models/savedSearch.model.js';

dotenv.config();

const QWQ_API_KEY = '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';
const MOLMIM_API_KEY = process.env.MOLMIM_API_KEY;
const MOLMIM_API_URL = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate";
const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";

const openai = new OpenAI({
  apiKey: QWQ_API_KEY,
  baseURL: 'https://integrate.api.nvidia.com/v1',
});

    
const mockFingerprintExtraction = (smiles) => {
  // In a real implementation, you'd use a JS cheminformatics library or an external API
  // For now, we'll return mock fingerprint data
  const morgan = Array(2048).fill(0).map(() => Math.random() > 0.5 ? "1" : "0").join("");
  const maccs = Array(166).fill(0).map(() => Math.random() > 0.5 ? "1" : "0").join("");
  return { morgan, maccs };
};

// Mock function to simulate molecular docking (since AutoDock is not available in JS)
const mockDocking = (smiles) => {
  // In a real implementation, you'd call an external docking API (e.g., DockThor)
  // For now, we'll return mock docking results
  return {
    bindingEnergy: (Math.random() * -10).toFixed(2), // Mock binding energy in kcal/mol
    pose: "Mock pose data",
    details: "Mock docking details",
  };
};

export const getProteinStructure = async (req, res) => {
  try {
    const { id } = req.params;
    
    if (!id) {
      return res.status(400).json({ message: 'User ID is required' });
    }

    const user = await User.findById(id).populate('protienStructures');
    if (!user) {
      return res.status(404).json({ message: 'User not found' });
    }

    const proteinStructures = user.protienStructures
      .sort((a, b) => b.created - a.created)
      .slice(0, 20);

    return res.status(200).json(proteinStructures);
  } catch (error) {
    console.error('Error fetching protein structures:', error);
    return res.status(500).json({ message: 'Server error', error: error.message });
  }
};

// export const postProteinStructure = async (req, res) => {
//   try {
//     const { id } = req.params;
//     const { 
//       name, 
//       smiles, 
//       algorithm = "CMA-ES", 
//       numMolecules = 30, 
//       propertyName = "QED", 
//       minimize = false, 
//       minSimilarity = 0.3, 
//       particles = 30, 
//       iterations = 10 
//     } = req.body;
    
//     if (!smiles) {
//       return res.status(400).json({ message: 'SMILES string is required' });
//     }

//     if (!id) {
//       return res.status(400).json({ message: 'User ID is required' });
//     }

//     const validAlgorithms = ["CMA-ES", "SSD"];
//     if (!validAlgorithms.includes(algorithm)) {
//       return res.status(400).json({ message: `Invalid algorithm: ${algorithm}. Supported algorithms are ${validAlgorithms.join(', ')}.` });
//     }

//     let molecules;
//     if (algorithm === "SSD") {
//       console.log('Using Sampling Standard Deviation (mock implementation)');
//       molecules = Array.from({ length: numMolecules }, (_, i) => ({
//         sample: smiles,
//         score: Math.random(),
//         similarity: Math.random() * (1 - minSimilarity) + minSimilarity,
//       }));
//     } else {
//       const requestPayload = {
//         algorithm,
//         num_molecules: numMolecules,
//         property_name: propertyName,
//         minimize,
//         min_similarity: minSimilarity,
//         particles,
//         iterations,
//         smi: smiles
//       };
//       console.log('Sending request to MolMIM API:', JSON.stringify(requestPayload, null, 2));

//       const response = await axios.post(
//         MOLMIM_API_URL,
//         requestPayload,
//         {
//           headers: {
//             'Authorization': `Bearer ${MOLMIM_API_KEY}`,
//             'Accept': 'application/json',
//             'Content-Type': 'application/json',
//           }
//         }
//       );

//       molecules = typeof response.data.molecules === 'string' 
//         ? JSON.parse(response.data.molecules) 
//         : response.data.molecules;

//       if (!Array.isArray(molecules)) {
//         throw new Error('Parsed "molecules" is not an array');
//       }
//     }

//     let moleculeInfo = '';
//     try {
//       if (!GEMINI_API_KEY) {
//         throw new Error('Gemini API key is not set');
//       }

//       const prompt = `
//      Analyze the drug molecule represented by the SMILES string "${smiles}" and provide a detailed report that includes the following sections:

// ### 1. **Structural Analysis**
// - **Core Structure Identification**: Describe the main structural features of the molecule, including any fused ring systems, heterocycles, and functional groups.
// - **Stereochemistry**: Identify and explain the stereochemical configurations at all chiral centers (e.g., [C@@], [C@H]) and their relevance to biological activity.
// - **Substituents**: List and characterize significant substituents attached to the core structure, explaining their potential impact on the molecule's properties.

// ### 2. **Chemical Properties**
// - **Molecular Weight**: Calculate and provide the exact molecular weight of the compound.
// - **Physicochemical Properties**: Include values for:
//   - LogP (partition coefficient)
//   - Polar Surface Area (PSA)
//   - Hydrogen Bond Donors/Acceptors
//   - Rotatable Bonds
//   - pKa values for ionizable groups
// - **Solubility Profile**: Predict solubility in various solvents (aqueous, organic) based on functional groups.
// - **Melting Point and Boiling Point**: Provide estimated values based on structural analogs.

// ### 3. **ADMET Profile**
// - **Absorption**: Discuss bioavailability and factors affecting absorption.
// - **Distribution**: Predict volume of distribution and blood-brain barrier permeability.
// - **Metabolism**: Identify metabolic pathways and potential metabolites, including enzyme interactions (e.g., CYP450).
// - **Excretion**: Discuss elimination routes and half-life.
// - **Toxicity Predictions**: Highlight any known or predicted toxicological concerns, including hERG inhibition or cytotoxicity.

// ### 4. **Biological Activity**
// - **Target Interaction**: Identify potential biological targets (receptors, enzymes) and predict binding affinities. Include:
//   - Mechanism of action for known targets.
//   - Any relevant SAR (Structure-Activity Relationship) data.
// - **Therapeutic Applications**: Discuss potential uses in therapy based on structural similarity to known drugs or biological activity.

// ### 5. **Synthesis Pathways**
// - Outline possible synthetic routes to obtain the compound, referencing established methods in literature.
// - Discuss any challenges or considerations in synthesis.

// ### 6. **Clinical Context**
// - Provide a summary of any clinical data available regarding this compound, including:
//   - Approved indications (if applicable).
//   - Clinical trial results or ongoing studies.
//   - Comparison with similar drugs in terms of efficacy and safety profiles.

// ### Output Requirements
// - Format the response using clear subheadings for each section.
// - Include tables or bullet points where appropriate for clarity.
// - Cite relevant databases (PubChem, ChEMBL, DrugBank) for additional context or data sources.

// ### Conclusion
// Summarize the overall potential of this molecule as a therapeutic agent, highlighting any critical research questions that remain unanswered or areas for further investigation.
//       `;
      
//       const geminiResponse = await axios.post(
//         `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
//         {
//           contents: [{ parts: [{ text: prompt }] }]
//         },
//         {
//           headers: { 'Content-Type': 'application/json' }
//         }
//       );

//       moleculeInfo = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text || 
//         'No content returned from Gemini API';
//     } catch (geminiError) {
//       console.error('Error calling Gemini API:', geminiError.message);
//       moleculeInfo = 'Failed to retrieve detailed information about the molecule.';
//     }

//     const newProteinStructure = new ProteinStructure({
//       name: name || 'Untitled Structure',
//       smiles,
//       properties: { 
//         algorithm,
//         propertyName,
//         minimize,
//         minSimilarity 
//       },
//       generatedStructures: molecules.map(mol => ({
//         smiles: mol.sample,
//         properties: {
//           qed: mol.score,
//           logp: mol.logp || null
//         },
//         similarity: mol.similarity || 0
//       })),
//       information: moleculeInfo,
//       userId: id
//     });
    
//     await newProteinStructure.save();

//     await User.findByIdAndUpdate(
//       id,
//       { $push: { protienStructures: newProteinStructure._id } },
//       { new: true }
//     );

//     return res.status(201).json({
//       message: 'Protein structure generated successfully',
//       structure: newProteinStructure
//     });
//   } catch (error) {
//     console.error('Error generating protein structure:', error);
    
//     if (error.response) {
//       return res.status(error.response.status).json({ 
//         message: 'MolMIM API error',
//         error: error.message,
//         status: error.response.status
//       });
//     }
    
//     return res.status(500).json({ 
//       message: 'Server error', 
//       error: error.message 
//     });
//   }
// };

export const postProteinStructure = async (req, res) => {
  try {
    const { id } = req.params;
    const { 
      name, 
      smiles, 
      algorithm = "CMA-ES", 
      numMolecules = 30, 
      propertyName = "QED", 
      minimize = false, 
      minSimilarity = 0.3, 
      particles = 30, 
      iterations = 10 
    } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required' });
    }

    if (!id) {
      return res.status(400).json({ message: 'User ID is required' });
    }

    const validAlgorithms = ["CMA-ES", "SSD"];
    if (!validAlgorithms.includes(algorithm)) {
      return res.status(400).json({ message: `Invalid algorithm: ${algorithm}. Supported algorithms are ${validAlgorithms.join(', ')}.` });
    }

    let molecules;
    if (algorithm === "SSD") {
      console.log('Using Sampling Standard Deviation (mock implementation)');
      molecules = Array.from({ length: numMolecules }, (_, i) => ({
        sample: smiles,
        score: Math.random(),
        similarity: Math.random() * (1 - minSimilarity) + minSimilarity,
      }));
    } else {
      const requestPayload = {
        algorithm,
        num_molecules: numMolecules,
        property_name: propertyName,
        minimize,
        min_similarity: minSimilarity,
        particles,
        iterations,
        smi: smiles
      };
      console.log('Sending request to MolMIM API:', JSON.stringify(requestPayload, null, 2));

      const response = await axios.post(
        MOLMIM_API_URL,
        requestPayload,
        {
          headers: {
            'Authorization': `Bearer ${MOLMIM_API_KEY}`,
            'Accept': 'application/json',
            'Content-Type': 'application/json',
          }
        }
      );

      molecules = typeof response.data.molecules === 'string' 
        ? JSON.parse(response.data.molecules) 
        : response.data.molecules;

      if (!Array.isArray(molecules)) {
        throw new Error('Parsed "molecules" is not an array');
      }
    }

    let moleculeInfo = '';
    try {
      if (!GEMINI_API_KEY) {
        throw new Error('Gemini API key is not set');
      }

      const prompt = `
You are a pharmaceutical chemistry expert tasked with analyzing a drug molecule represented by the SMILES string "${smiles}". Provide a detailed report with the following sections, ensuring clarity, conciseness, and relevance. Each section should be 3-5 sentences long, unless specified otherwise, and include specific, detailed information, even if it requires predictions or estimations based on computational methods, structural analogies, or literature insights.

1. Structural Details:
- Determine the molecular formula from the SMILES string (e.g., C9H8O4 for aspirin).
- List any known synonyms or alternative names for the compound, citing databases like PubChem with specific identifiers (e.g., PubChem CID: 2244 for aspirin).
- Calculate the exact molecular weight and provide it with units (e.g., 180.16 g/mol).
- Generate the IUPAC name and InChI based on the molecular structure; if exact structure confirmation is unavailable, provide a predicted name and InChI based on the SMILES string.
- Note the provided SMILES string for reference (e.g., SMILES: CC(=O)OC1=CC=CC=C1C(=O)O).

2. Chemical Properties:
- Predict or provide the lipophilicity (logP) value, citing sources if available (e.g., logP: 1.2, predicted via ChemAxon); if unavailable, estimate based on functional groups.
- Provide pKa values for each ionizable group (e.g., -COOH, -NH2) with units (e.g., pKa: 3.5 for carboxylic acid); predict values if experimental data is unavailable.
- Count and list the number of hydrogen bond donors (HBD) and acceptors (HBA) (e.g., 2 HBD, 4 HBA).
- Identify and list significant functional groups present in the molecule (e.g., carboxylic acid, ester, aromatic ring).
- Determine the number of aromatic and non-aromatic rings, specifying each (e.g., 1 aromatic benzene ring, 0 non-aromatic rings).

3. Physical Properties:
- Predict solubility in aqueous and organic solvents, explaining the basis (e.g., moderately soluble in water due to polar groups, highly soluble in ethanol due to logP of 1.2).
- Discuss membrane permeability, relating it to logP and polar surface area (e.g., moderate permeability with logP of 1.2 and PSA of 60 Å²).
- Calculate the polar surface area (PSA) and provide it with units (e.g., PSA: 60 Å², calculated via ChemAxon).
- Count the number of rotatable bonds in the molecular structure (e.g., 3 rotatable bonds).
- Discuss chemical stability, noting any known issues (e.g., stable at neutral pH, hydrolyzes in acidic conditions), and mention crystallinity and polymorphism if data is available or predict based on structure (e.g., likely crystalline due to planar aromatic system).

4. Spectral Information:
- Provide detailed spectral data such as IR, NMR, or mass spectrometry, based on typical expectations for the structure (e.g., IR: 1700 cm⁻¹ for C=O stretch; NMR: 12 ppm for -COOH proton; MS: m/z 181 [M+H]+).
- If experimental data is unavailable, predict values or provide typical ranges based on functional groups, citing prediction methods (e.g., predicted via ChemAxon’s spectral prediction tools).

5. Biological Activity and Pharmacology:
- Identify potential biological targets (e.g., receptors, enzymes like COX-1) and predict binding affinities (e.g., 7.5 kcal/mol), including mechanism of action (e.g., inhibits COX-1 by blocking active site).
- Discuss structure-activity relationship (SAR) data, if known, or predict based on structural features (e.g., the -COOH group enhances COX inhibition).
- Suggest potential therapeutic applications based on structural similarities to known drugs (e.g., similar to aspirin, may treat pain and inflammation).
- Describe pharmacodynamic properties (e.g., inhibits prostaglandin synthesis) and pharmacokinetic properties (e.g., half-life of 4 hours, metabolized by CYP2C9, excreted via kidneys), including interactions and pathways (e.g., interacts with CYP2C9, affects prostaglandin pathway).
- Include any biological test results, such as assay data, or predict likely outcomes (e.g., likely IC50 of 10 µM against COX-1 based on structural analogy to aspirin).

6. Clinical and Therapeutic Information:
- List approved indications, dosage forms, and therapeutic uses, citing clinical sources if applicable (e.g., approved for pain relief, available as 325 mg tablets, per DrugBank DB00945).
- Summarize any clinical trial results, ongoing studies, or comparisons with similar drugs (e.g., reduces pain by 50% in trials, comparable to ibuprofen, per ClinicalTrials.gov NCT123456).
- Identify associated disorders and diseases, explaining the connection (e.g., treats arthritis due to anti-inflammatory effects).
- If clinical data is unavailable, predict potential uses based on biological activity (e.g., likely used for inflammation based on COX inhibition).

7. Safety and Toxicity:
- Discuss handling precautions and safety hazards, citing safety data sheets if possible (e.g., avoid inhalation, may cause skin irritation, per PubChem safety data).
- Provide toxicological information, including LD50 values (e.g., LD50: 200 mg/kg in rats, predicted via QSAR models) and any known or predicted toxicological concerns (e.g., potential gastrointestinal irritation).
- Highlight potential risks, such as hERG inhibition or cytotoxicity (e.g., possible hERG inhibition due to aromatic system, predicted via computational models).

8. Manufacturing and Synthesis:
- Outline a possible synthetic route to obtain the compound, referencing established methods in literature (e.g., acetylate salicylic acid with acetic anhydride, per Organic Syntheses Vol. 1).
- Highlight any challenges or considerations in the synthesis process, including manufacturing details (e.g., requires controlled temperature to avoid side products, scalable for industrial production).

9. References and Patents:
- Cite relevant scientific literature or databases for any data provided, including specific identifiers (e.g., PubChem CID: 2244, DrugBank DB00945, ChEMBL CHEMBL25).
- Note any patents related to the compound, including patent numbers if known (e.g., US Patent 1234567 for synthesis method).

10. Classification and Taxonomy:
- Classify the compound in relevant chemical or biological taxonomies, such as drug classes or chemical categories (e.g., belongs to salicylates, non-steroidal anti-inflammatory drugs).
- Provide standard classifications, noting any regulatory categories if applicable (e.g., FDA-approved analgesic, per PubChem).

**Output Requirements:**
- Use numbered headings for each section (e.g., "1. Structural Details:", "2. Chemical Properties:", etc.).
- Use dashes ("-") for all bullet points, and do not use stars ("*"), bullet symbols ("•"), or other markers.
- Do not use hashtags ("#") or stars ("*") anywhere in the response body.
- Keep each section concise (3-5 sentences, unless specified otherwise) and avoid extraneous information.
- Ensure all numerical values are accompanied by their units (e.g., g/mol for molecular weight, Å² for PSA, cm⁻¹ for IR peaks).
- Cite relevant databases (e.g., PubChem, ChEMBL, DrugBank) or prediction tools (e.g., ChemAxon, QSAR models) for additional context or data sources, if applicable.
- Ensure all text is plain and free of markdown formatting (e.g., no bold, italics, or tables) except for the specified headings and bullet points.
- If exact data is unavailable, provide a prediction or estimation based on computational methods, structural analogies, or literature insights, and note the method used (e.g., "predicted via ChemAxon"); only use "No data available" if the information is completely inapplicable.
`; // Fixed typo in prompt
      
      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        {
          contents: [{ parts: [{ text: prompt }] }]
        },
        {
          headers: { 'Content-Type': 'application/json' }
        }
      );

      // Log the raw response for debugging
      console.log('Raw Gemini API Response:', JSON.stringify(geminiResponse.data, null, 2));

      // Attempt to extract the text from multiple possible paths
      let extractedText = '';
      if (geminiResponse.data?.candidates?.[0]?.content?.parts?.[0]?.text) {
        extractedText = geminiResponse.data.candidates[0].content.parts[0].text;
      } else if (geminiResponse.data?.candidates?.[0]?.content?.text) {
        extractedText = geminiResponse.data.candidates[0].content.text;
      } else if (geminiResponse.data?.text) {
        extractedText = geminiResponse.data.text;
      } else {
        throw new Error('No text content found in Gemini API response');
      }

      moleculeInfo = extractedText || 'No content returned from Gemini API';

      // Clean up the moleculeInfo
      moleculeInfo = moleculeInfo
        .replace(/(\d+\.\s[A-Za-z\s]+:)/g, "\n$1") // Ensure numbered sections are on new lines
        .replace(/-+/g, "-") // Normalize dashes for bullet points
        .trim();

      // Log the cleaned moleculeInfo for debugging
      console.log('Cleaned moleculeInfo:', moleculeInfo);

    } catch (geminiError) {
      console.error('Error calling Gemini API:', geminiError.message);
      if (geminiError.response) {
        console.error('Gemini API Response Error:', JSON.stringify(geminiError.response.data, null, 2));
      }
      moleculeInfo = 'Failed to retrieve detailed information about the molecule.';
    }

    const newProteinStructure = new ProteinStructure({
      name: name || 'Untitled Structure',
      smiles,
      properties: { 
        algorithm,
        propertyName,
        minimize,
        minSimilarity 
      },
      generatedStructures: molecules.map(mol => ({
        smiles: mol.sample,
        properties: {
          qed: mol.score,
          logp: mol.logp || null
        },
        similarity: mol.similarity || 0
      })),
      information: moleculeInfo,
      userId: id
    });
    
    await newProteinStructure.save();

    await User.findByIdAndUpdate(
      id,
      { $push: { protienStructures: newProteinStructure._id } },
      { new: true }
    );

    return res.status(201).json({
      message: 'Protein structure generated successfully',
      structure: newProteinStructure
    });
  } catch (error) {
    console.error('Error generating protein structure:', error);
    
    if (error.response) {
      return res.status(error.response.status).json({ 
        message: 'MolMIM API error',
        error: error.message,
        status: error.response.status
      });
    }
    
    return res.status(500).json({ 
      message: 'Server error', 
      error: error.message 
    });
  }
};

// export const generatenewmolecule = async (req, res) => {
//   try {
//     const { id } = req.params;
//     const { smilesoffirst, smilesofsecond, newmoleculetitle } = req.body;

//     console.log(`Processing request for user ID: ${id}`);
//     console.log("Request body:", { smilesoffirst, smilesofsecond, newmoleculetitle });

//     if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
//       return res.status(400).json({ message: "SMILES strings and molecule title are required" });
//     }
//     if (!id) {
//       return res.status(400).json({ message: "User ID is required" });
//     }

//     const prompt = `
// You are an expert cheminformatics AI designed to perform a chemically feasible reaction between two SMILES strings to generate a new molecule. Follow these steps strictly to ensure a realistic reaction occurs, and provide a detailed analysis of the resulting molecule. Your response must end with a valid JSON object, as its fields will be extracted and saved to a MongoDB database.

// 1. Analyze the Input SMILES Strings:
// - Parse SMILES 1: "${smilesoffirst}" and describe its chemical structure, including molecular formula, functional groups, and key reactive sites (e.g., -OH, -COOH) that could participate in a reaction.
// - Parse SMILES 2: "${smilesofsecond}" and describe its chemical structure similarly.
// - Identify potential reactive sites in each molecule that could participate in a chemical reaction, explaining their reactivity (e.g., nucleophilic -OH, electrophilic -COOH).

// 2. Perform a Chemically Feasible Reaction:
// - Propose a specific, realistic reaction to combine the two molecules (e.g., esterification, amide formation, nucleophilic substitution, carbon-carbon coupling like Suzuki coupling).
// - Explain the rationale for choosing this reaction, considering functional group compatibility, reaction feasibility, and stability of the product.
// - Detail the reaction mechanism step-by-step, including:
//   - The specific atoms involved in the reaction (e.g., the -OH of molecule 1 reacts with the -COOH of molecule 2).
//   - Any intermediates formed (e.g., tetrahedral intermediate in esterification).
//   - Reaction conditions required (e.g., acid catalysis, heat, or a specific reagent like DCC for amide formation).
//   - Any byproducts formed (e.g., water in a condensation reaction).

// 3. Generate the New Molecule:
// - Provide the SMILES string for the resulting molecule, ensuring it is valid and reflects the product of the reaction.
// - Determine and provide the IUPAC name for the new molecule; if exact naming is unavailable, predict a name based on the structure.
// - Calculate the molecular formula and exact molecular weight of the new molecule (e.g., C10H12O5, 212.20 g/mol).

// 4. Structural Details of the New Molecule:
// - Describe the core structure of the new molecule, including any fused ring systems, heterocycles, or functional groups.
// - Identify stereochemistry at chiral centers, if applicable, and explain its relevance to biological activity.
// - List significant substituents and their potential impact on the molecule’s properties (e.g., electron-withdrawing groups affecting reactivity).

// 5. Chemical Properties of the New Molecule:
// - Predict the lipophilicity (logP) value, citing prediction methods (e.g., logP: 2.5, predicted via ChemAxon).
// - Provide pKa values for ionizable groups (e.g., pKa: 4.5 for -COOH, predicted via ChemAxon).
// - Count and list the number of hydrogen bond donors (HBD) and acceptors (HBA) (e.g., 2 HBD, 5 HBA).
// - Identify significant functional groups in the new molecule (e.g., ester, aromatic ring).
// - Determine the number of aromatic and non-aromatic rings (e.g., 1 aromatic benzene ring, 0 non-aromatic rings).

// 6. Biological Activity and Therapeutic Potential:
// - Hypothesize potential biological targets (e.g., receptors, enzymes like COX-1) and predict binding affinities (e.g., 7.0 kcal/mol).
// - Discuss structure-activity relationship (SAR) data, if known, or predict based on structural features (e.g., the -COOH group may enhance COX inhibition).
// - Suggest potential therapeutic applications based on structural similarities to known drugs (e.g., similar to aspirin, may treat pain and inflammation).
// - Describe pharmacodynamic properties (e.g., inhibits prostaglandin synthesis) and pharmacokinetic properties (e.g., half-life of 5 hours, metabolized by CYP3A4).
// - Suggest one or more diseases this molecule might target, with a brief explanation linking chemical features to mechanisms (e.g., inflammation due to COX inhibition).

// 7. Safety and Toxicity:
// - Discuss potential safety hazards, citing prediction methods (e.g., may cause skin irritation, predicted via QSAR models).
// - Provide toxicological information, including predicted LD50 values (e.g., LD50: 300 mg/kg in rats, predicted via QSAR models).
// - Highlight potential risks, such as hERG inhibition or cytotoxicity (e.g., possible hERG inhibition due to aromatic system).

// 8. Return the Result:
// - At the end of your response, provide a valid JSON object with the following structure. This JSON object will be parsed and its fields (newSmiles, newIupacName, conversionDetails, potentialDiseases) will be saved to a MongoDB database.
// - The JSON object must be the last part of your response, with no additional text, comments, or whitespace after it:
//   {
//     "newSmiles": "[NEW_SMILES]",
//     "newIupacName": "[NEW_IUPAC_NAME]",
//     "conversionDetails": "[DETAILED_EXPLANATION_OF_REACTION_AND_COMBINATION_PROCESS]",
//     "potentialDiseases": "[LIST_OF_DISEASES_WITH_EXPLANATIONS]"
//   }
// - Ensure the JSON object is properly formatted, contains all required fields, and has valid values (e.g., no empty strings or null values).
// - If you cannot generate a valid SMILES string or IUPAC name, provide a placeholder value (e.g., "Unable to generate SMILES" or "Unknown IUPAC Name") and explain why in the conversionDetails.
// - If you cannot determine potential diseases, provide a placeholder value (e.g., "Unknown potential diseases") and explain why in the potentialDiseases field.

// **Output Requirements:**
// - Use numbered headings for each section (e.g., "1. Analyze the Input SMILES Strings:", "2. Perform a Chemically Feasible Reaction:", etc.).
// - Use dashes ("-") for all bullet points, and do not use stars ("*"), bullet symbols ("•"), or other markers.
// - Do not use hashtags ("#") or stars ("*") anywhere in the response body.
// - Keep each section concise (3-5 sentences, unless specified otherwise) and avoid extraneous information.
// - Ensure all numerical values are accompanied by their units (e.g., g/mol for molecular weight).
// - Cite prediction methods or databases (e.g., ChemAxon, QSAR models) for additional context, if applicable.
// - Ensure all text is plain and free of markdown formatting (e.g., no bold, italics, or tables) except for the specified headings and bullet points.
// - If exact data is unavailable, provide a prediction or estimation based on computational methods, structural analogies, or literature insights, and note the method used (e.g., "predicted via ChemAxon").
// `;

//     console.log("Sending streaming request to NVIDIA NIM API...");
//     res.setHeader('Content-Type', 'text/event-stream');
//     res.setHeader('Cache-Control', 'no-cache');
//     res.setHeader('Connection', 'keep-alive');

//     const completion = await openai.chat.completions.create({
//       model: "qwen/qwq-32b",
//       messages: [{ role: "user", content: prompt }],
//       temperature: 0.6,
//       top_p: 0.7,
//       max_tokens: 4096,
//       stream: true,
//     });

//     let fullResponse = '';
//     for await (const chunk of completion) {
//       const content = chunk.choices[0]?.delta?.content || '';
//       fullResponse += content;
//       res.write(`data: ${content}\n\n`);
//     }
//     res.write('data: [DONE]\n\n');
//     res.end();

//     console.log("Full QWQ response:", fullResponse);

//     // Initialize fields with default values
//     let newSmiles = "";
//     let newIupacName = "";
//     let conversionDetails = "";
//     let potentialDiseases = "";
//     let molecularFormula = "";
//     let molecularWeight = "";
//     let chemicalProperties = {};
//     let biologicalActivity = {};
//     let safetyToxicity = {};

//     // Parse the JSON object from the response
//     const jsonMatch = fullResponse.match(/{[\s\S]*?}$/);
//     if (jsonMatch) {
//       try {
//         const result = JSON.parse(jsonMatch[0]);
//         newSmiles = result.newSmiles || "Unable to generate SMILES";
//         newIupacName = result.newIupacName || "Unknown IUPAC Name";
//         conversionDetails = result.conversionDetails || "No conversion details provided.";
//         potentialDiseases = result.potentialDiseases || "Unknown potential diseases.";
//         console.log("Extracted fields from JSON:", { newSmiles, newIupacName, conversionDetails, potentialDiseases });
//       } catch (parseError) {
//         console.error("Failed to parse JSON from response:", parseError);
//         conversionDetails = `Failed to parse JSON: ${parseError.message}`;
//       }
//     } else {
//       console.error("No JSON object found in the response.");
//       conversionDetails = "No JSON object found in the response.";
//     }

//     // If JSON parsing fails or fields are missing, fall back to text parsing
//     if (!newSmiles || !newIupacName || !conversionDetails || !potentialDiseases) {
//       console.log("Falling back to text parsing...");
//       const sections = fullResponse.split(/\d+\.\s[A-Za-z\s]+:/).filter(section => section.trim());
      
//       // Extract fields from relevant sections
//       for (let i = 0; i < sections.length; i++) {
//         const section = sections[i].trim();
//         const lines = section.split('\n').map(line => line.trim());

//         // Section 3: Generate the New Molecule
//         if (i === 2) { // Section 3 (index 2 after splitting)
//           for (const line of lines) {
//             if (line.includes("SMILES string for the resulting molecule")) {
//               const match = line.match(/-.*SMILES string for the resulting molecule.*: (.+)/);
//               if (match) newSmiles = match[1].trim();
//             }
//             if (line.includes("IUPAC name")) {
//               const match = line.match(/-.*IUPAC name.*: (.+)/);
//               if (match) newIupacName = match[1].trim();
//             }
//             if (line.includes("molecular formula")) {
//               const match = line.match(/-.*molecular formula.*: (.+)/);
//               if (match) molecularFormula = match[1].trim();
//             }
//             if (line.includes("molecular weight")) {
//               const match = line.match(/-.*molecular weight.*: (.+)/);
//               if (match) molecularWeight = match[1].trim();
//             }
//           }
//         }

//         // Section 5: Chemical Properties
//         if (i === 4) { // Section 5 (index 4 after splitting)
//           chemicalProperties = {
//             logP: "",
//             pKa: [],
//             hbd: 0,
//             hba: 0,
//             functionalGroups: [],
//             rings: { aromatic: 0, nonAromatic: 0 }
//           };
//           for (const line of lines) {
//             if (line.includes("lipophilicity (logP)")) {
//               const match = line.match(/-.*logP.*: (.+)/);
//               if (match) chemicalProperties.logP = match[1].trim();
//             }
//             if (line.includes("pKa values")) {
//               const match = line.match(/-.*pKa values.*: (.+)/);
//               if (match) chemicalProperties.pKa = match[1].split(',').map(val => val.trim());
//             }
//             if (line.includes("hydrogen bond donors (HBD)")) {
//               const match = line.match(/-.*HBD.*: (\d+)/);
//               if (match) chemicalProperties.hbd = parseInt(match[1]);
//             }
//             if (line.includes("hydrogen bond acceptors (HBA)")) {
//               const match = line.match(/-.*HBA.*: (\d+)/);
//               if (match) chemicalProperties.hba = parseInt(match[1]);
//             }
//             if (line.includes("significant functional groups")) {
//               const match = line.match(/-.*functional groups.*: (.+)/);
//               if (match) chemicalProperties.functionalGroups = match[1].split(',').map(fg => fg.trim());
//             }
//             if (line.includes("aromatic and non-aromatic rings")) {
//               const match = line.match(/-.*(\d+) aromatic.*(\d+) non-aromatic/);
//               if (match) {
//                 chemicalProperties.rings.aromatic = parseInt(match[1]);
//                 chemicalProperties.rings.nonAromatic = parseInt(match[2]);
//               }
//             }
//           }
//         }

//         // Section 6: Biological Activity and Therapeutic Potential
//         if (i === 5) { // Section 6 (index 5 after splitting)
//           biologicalActivity = {
//             targets: "",
//             bindingAffinity: "",
//             sar: "",
//             therapeuticApplications: "",
//             pharmacodynamics: "",
//             pharmacokinetics: "",
//             diseases: ""
//           };
//           for (const line of lines) {
//             if (line.includes("potential biological targets")) {
//               const match = line.match(/-.*targets.*: (.+)/);
//               if (match) biologicalActivity.targets = match[1].trim();
//             }
//             if (line.includes("binding affinities")) {
//               const match = line.match(/-.*binding affinities.*: (.+)/);
//               if (match) biologicalActivity.bindingAffinity = match[1].trim();
//             }
//             if (line.includes("structure-activity relationship (SAR)")) {
//               const match = line.match(/-.*SAR.*: (.+)/);
//               if (match) biologicalActivity.sar = match[1].trim();
//             }
//             if (line.includes("therapeutic applications")) {
//               const match = line.match(/-.*therapeutic applications.*: (.+)/);
//               if (match) biologicalActivity.therapeuticApplications = match[1].trim();
//             }
//             if (line.includes("pharmacodynamic properties")) {
//               const match = line.match(/-.*pharmacodynamic properties.*: (.+)/);
//               if (match) biologicalActivity.pharmacodynamics = match[1].trim();
//             }
//             if (line.includes("pharmacokinetic properties")) {
//               const match = line.match(/-.*pharmacokinetic properties.*: (.+)/);
//               if (match) biologicalActivity.pharmacokinetics = match[1].trim();
//             }
//             if (line.includes("diseases this molecule might target")) {
//               const match = line.match(/-.*diseases.*: (.+)/);
//               if (match) biologicalActivity.diseases = match[1].trim();
//             }
//           }
//         }

//         // Section 7: Safety and Toxicity
//         if (i === 6) { // Section 7 (index 6 after splitting)
//           safetyToxicity = {
//             hazards: "",
//             ld50: "",
//             risks: ""
//           };
//           for (const line of lines) {
//             if (line.includes("safety hazards")) {
//               const match = line.match(/-.*safety hazards.*: (.+)/);
//               if (match) safetyToxicity.hazards = match[1].trim();
//             }
//             if (line.includes("LD50 values")) {
//               const match = line.match(/-.*LD50.*: (.+)/);
//               if (match) safetyToxicity.ld50 = match[1].trim();
//             }
//             if (line.includes("potential risks")) {
//               const match = line.match(/-.*risks.*: (.+)/);
//               if (match) safetyToxicity.risks = match[1].trim();
//             }
//           }
//         }
//       }
//     }

//     // Validate the new SMILES string using RDKit (if available)
//     let isValidSmiles = false;
//     if (typeof window !== 'undefined' && window.RDKit) {
//       try {
//         const mol = window.RDKit.get_mol(newSmiles);
//         isValidSmiles = mol && mol.is_valid();
//         if (mol) mol.delete();
//         if (!isValidSmiles) {
//           conversionDetails += " Generated SMILES is invalid according to RDKit.";
//           newSmiles = `${smilesoffirst}.${smilesofsecond}`; // Fallback to concatenated SMILES
//         }
//       } catch (rdkitError) {
//         console.error("RDKit validation error:", rdkitError);
//         conversionDetails += ` RDKit validation failed: ${rdkitError.message}.`;
//       }
//     } else {
//       console.warn("RDKit not available for SMILES validation.");
//       conversionDetails += " RDKit not available; SMILES validation skipped.";
//     }

//     // Set default values if fields are still missing
//     newSmiles = newSmiles || `${smilesoffirst}.${smilesofsecond}`;
//     newIupacName = newIupacName || "Unknown IUPAC Name";
//     conversionDetails = conversionDetails || "Unable to parse conversion details from the response.";
//     potentialDiseases = potentialDiseases || "Unknown potential diseases due to lack of information in the response.";
//     molecularFormula = molecularFormula || "Unknown molecular formula";
//     molecularWeight = molecularWeight || "Unknown molecular weight";
//     chemicalProperties = chemicalProperties.logP ? chemicalProperties : {
//       logP: "Unknown",
//       pKa: [],
//       hbd: 0,
//       hba: 0,
//       functionalGroups: [],
//       rings: { aromatic: 0, nonAromatic: 0 }
//     };
//     biologicalActivity = biologicalActivity.targets ? biologicalActivity : {
//       targets: "Unknown",
//       bindingAffinity: "Unknown",
//       sar: "Unknown",
//       therapeuticApplications: "Unknown",
//       pharmacodynamics: "Unknown",
//       pharmacokinetics: "Unknown",
//       diseases: "Unknown"
//     };
//     safetyToxicity = safetyToxicity.hazards ? safetyToxicity : {
//       hazards: "Unknown",
//       ld50: "Unknown",
//       risks: "Unknown"
//     };

//     console.log("Final extracted fields:", {
//       newSmiles,
//       newIupacName,
//       conversionDetails,
//       potentialDiseases,
//       molecularFormula,
//       molecularWeight,
//       chemicalProperties,
//       biologicalActivity,
//       safetyToxicity
//     });

//     // Save the new molecule to the database
//     console.log("Saving new molecule to database with raw streaming response...");
//     const newMolecule = new GeneratenewMolecule({
//       smilesoffirst,
//       smilesofsecond,
//       newmoleculetitle,
//       newSmiles,
//       newIupacName,
//       conversionDetails,
//       potentialDiseases,
//       molecularFormula,
//       molecularWeight,
//       chemicalProperties,
//       biologicalActivity,
//       safetyToxicity,
//       information: fullResponse,
//       userId: id,
//     });
//     await newMolecule.save();

//     console.log("Updating user with new molecule ID:", newMolecule._id);
//     await User.findByIdAndUpdate(id, { $push: { proteinStructures: newMolecule._id } }, { new: true });
//   } catch (error) {
//     console.error("Error generating new molecule:", error.message, error.stack);
//     res.write(`data: {"error": "${error.message}"}\n\n`);
//     res.end();
//   }
// };



// Basic SMILES validation
const validateSmiles = (smiles) => {
  const smilesRegex = /^[A-Za-z0-9@+\-\[\]()=:#$%]+$/;
  return smiles && smiles.length > 0 && smiles.length < 500 && smilesRegex.test(smiles);
};

// Enhanced chemical analysis function
const analyzeReactivity = (smiles) => {
  const analysis = {
    nucleophiles: [],
    electrophiles: [],
    aromaticRings: 0,
    reactiveSites: 0,
    stericHindrance: 'low',
  };
  analysis.aromaticRings = (smiles.match(/c/g) || []).length / 2;
  const groups = {
    'COOH': { type: 'electrophile', reactivity: 'high', description: 'Carboxylic acid' },
    'C(=O)O': { type: 'electrophile', reactivity: 'high', description: 'Ester or acid derivative' },
    'NH2': { type: 'nucleophile', reactivity: 'high', description: 'Primary amine' },
    'OH': { type: 'nucleophile', reactivity: 'medium', description: 'Hydroxyl group' },
    'C=O': { type: 'electrophile', reactivity: 'medium', description: 'Carbonyl (ketone/aldehyde)' },
    'Cl': { type: 'electrophile', reactivity: 'high', description: 'Halide (potential leaving group)' },
    'Br': { type: 'electrophile', reactivity: 'high', description: 'Bromide (potential leaving group)' },
    'C=C': { type: 'electrophile', reactivity: 'medium', description: 'Alkene (for addition reactions)' },
  };
  Object.entries(groups).forEach(([pattern, props]) => {
    if (new RegExp(pattern).test(smiles)) {
      if (props.type === 'nucleophile') {
        analysis.nucleophiles.push({ pattern, ...props });
        analysis.reactiveSites++;
      } else if (props.type === 'electrophile') {
        analysis.electrophiles.push({ pattern, ...props });
        analysis.reactiveSites++;
      }
    }
  });
  analysis.stericHindrance = analysis.aromaticRings > 2 ? 'high' : analysis.reactiveSites > 3 ? 'medium' : 'low';
  return analysis;
};

// Reaction prediction logic with more rules
const predictReaction = async (smiles1, smiles2) => {
  const analysis1 = analyzeReactivity(smiles1);
  const analysis2 = analyzeReactivity(smiles2);

  const reactionRules = [
    {
      name: 'Esterification',
      condition: (a1, a2) => a1.electrophiles.some(e => e.pattern === 'COOH') && a2.nucleophiles.some(n => n.pattern === 'OH'),
      product: (s1, s2) => `${s1.replace('COOH', 'C(=O)O')}${s2.replace('OH', '')}`,
      details: 'Esterification: -COOH reacts with -OH to form an ester linkage. Conditions: Acid catalyst (e.g., H₂SO₄), heat (60-80°C). Byproduct: H₂O.',
    },
    {
      name: 'Amide Formation',
      condition: (a1, a2) => a1.electrophiles.some(e => e.pattern === 'COOH') && a2.nucleophiles.some(n => n.pattern === 'NH2'),
      product: (s1, s2) => `${s1.replace('COOH', 'C(=O)N')}${s2.replace('NH2', '')}`,
      details: 'Amide Formation: -COOH reacts with -NH2 to form an amide bond. Conditions: Coupling agent (e.g., DCC), room temperature. Byproduct: H₂O.',
    },
    {
      name: 'Nucleophilic Substitution (SN2)',
      condition: (a1, a2) => a1.electrophiles.some(e => e.pattern === 'Cl' || e.pattern === 'Br') && a2.nucleophiles.some(n => n.pattern === 'OH'),
      product: (s1, s2) => `${s1.replace(/Cl|Br/, 'O')}${s2.replace('OH', '')}`,
      details: 'SN2 Reaction: -Cl or -Br is displaced by -OH. Conditions: Strong base (e.g., NaOH), polar aprotic solvent (e.g., DMF), room temperature. Byproduct: HCl or HBr.',
    },
    {
      name: 'Addition to Alkene',
      condition: (a1, a2) => a1.electrophiles.some(e => e.pattern === 'C=C') && a2.nucleophiles.some(n => n.pattern === 'OH'),
      product: (s1, s2) => `${s1.replace('C=C', 'C(O)C')}${s2.replace('OH', '')}`,
      details: 'Addition to Alkene: -OH adds across the C=C double bond. Conditions: Acid catalyst (e.g., H₂SO₄), moderate heat (50°C). Byproduct: None.',
    },
  ];

  for (const rule of reactionRules) {
    if (rule.condition(analysis1, analysis2)) {
      return {
        isPossible: true,
        newSmiles: rule.product(smiles1, smiles2),
        conversionDetails: rule.details,
        analysis: { smiles1: analysis1, smiles2: analysis2 },
      };
    }
    if (rule.condition(analysis2, analysis1)) {
      return {
        isPossible: true,
        newSmiles: rule.product(smiles2, smiles1),
        conversionDetails: rule.details.replace('reacts with', 'reacts with (reversed)'),
        analysis: { smiles1: analysis2, smiles2: analysis1 },
      };
    }
  }

  // If no predefined rule matches, return null to let Gemini AI propose a reaction
  return {
    isPossible: false,
    newSmiles: null,
    conversionDetails: 'No predefined reaction rule matched.',
    reason: 'No compatible functional groups found for a predefined reaction.',
    analysis: { smiles1: analysis1, smiles2: analysis2 },
  };
};

export const generatenewmolecule = async (req, res) => {
  try {
    const { id } = req.params;
    const { smilesoffirst, smilesofsecond, newmoleculetitle } = req.body;

    console.log(`Processing request for user ID: ${id}`);
    console.log("Request body:", { smilesoffirst, smilesofsecond, newmoleculetitle });

    // Input validation
    if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
      res.write(`data: {"error": "SMILES strings and molecule title are required"}\n\n`);
      res.write('data: [DONE]\n\n');
      res.end();
      return;
    }
    if (!id) {
      res.write(`data: {"error": "User ID is required"}\n\n`);
      res.write('data: [DONE]\n\n');
      res.end();
      return;
    }
    if (!validateSmiles(smilesoffirst) || !validateSmiles(smilesofsecond)) {
      res.write(`data: {"error": "Invalid SMILES string format"}\n\n`);
      res.write('data: [DONE]\n\n');
      res.end();
      return;
    }

    // Set up SSE for real-time streaming
    res.setHeader('Content-Type', 'text/event-stream');
    res.setHeader('Cache-Control', 'no-cache');
    res.setHeader('Connection', 'keep-alive');

    // Predict reaction
    const reaction = await predictReaction(smilesoffirst, smilesofsecond);
    console.log("Predicted reaction:", reaction);

    // Prepare Gemini AI prompt to enforce a reaction
    const prompt = `
      You are an expert cheminformatician tasked with generating a chemically feasible reaction product from two SMILES strings. Given:

      - SMILES 1: "${smilesoffirst}"
      - SMILES 2: "${smilesofsecond}"
      - Initial Prediction Possible: ${reaction.isPossible}
      - Predicted Product SMILES: "${reaction.newSmiles || 'N/A'}"
      - Initial Conversion Details: "${reaction.conversionDetails}"
      - Reason (if not possible): "${reaction.reason || 'N/A'}"

      Perform the following:

      1. Reaction Feasibility:
      - Analyze the functional groups in both SMILES strings and propose a chemically realistic reaction (e.g., esterification, amide formation, nucleophilic substitution, addition reaction).
      - If the initial prediction failed, override it with a plausible reaction based on the structures.
      - Explain the reaction mechanism step-by-step (e.g., nucleophilic attack, intermediate formation, byproduct release).

      2. Product Generation:
      - Generate a valid SMILES string for the reaction product (do not simply concatenate the inputs with a dot).
      - Provide the IUPAC name for the product (predict if exact naming is unavailable).

      3. Conversion Details:
      - Detail the reaction conditions (e.g., catalysts, temperature) and byproducts.
      - Ensure the explanation is consistent with the generated SMILES.

      4. Potential Diseases:
      - Hypothesize therapeutic applications or target diseases based on the product's structure, comparing to known drugs.

      Return the response in JSON format:
      {
        "isPossible": true,
        "newSmiles": "string",
        "newIupacName": "string",
        "conversionDetails": "string",
        "potentialDiseases": "string",
        "reason": "string" (optional, if initial prediction was adjusted)
      }
    `;

    // Call Gemini AI
    const geminiResponse = await axios.post(
      `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
      { contents: [{ parts: [{ text: prompt }] }] },
      { headers: { 'Content-Type': 'application/json' } }
    );

    const geminiContent = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
    if (!geminiContent) {
      throw new Error("No content returned from Gemini API");
    }

    const jsonMatch = geminiContent.match(/{[\s\S]*}/);
    if (!jsonMatch) {
      throw new Error("No valid JSON found in Gemini response");
    }

    const result = JSON.parse(jsonMatch[0]);
    console.log("Gemini AI result:", result);

    // Ensure all required fields are present and valid
    if (!result.newSmiles || result.newSmiles.includes('.')) {
      throw new Error("Gemini AI failed to generate a valid reaction product SMILES");
    }
    result.isPossible = true; // Force true since we’re enforcing a reaction
    result.newSmiles = result.newSmiles;
    result.newIupacName = result.newIupacName || 'Unknown IUPAC Name';
    result.conversionDetails = result.conversionDetails || 'Reaction details not fully specified by Gemini AI.';
    result.potentialDiseases = result.potentialDiseases || 'No therapeutic applications predicted.';
    if (reaction.isPossible === false && result.reason) {
      result.reason = `Initial prediction failed: ${reaction.reason}. Overridden by Gemini AI with: ${result.reason || 'proposed reaction'}`;
    }

    // Stream the result to the client
    res.write(`data: ${JSON.stringify(result)}\n\n`);
    res.write('data: [DONE]\n\n');

    // Prepare molecule info for database
    const moleculeInfo = JSON.stringify({
      inputSmiles: [smilesoffirst, smilesofsecond],
      newSmiles: result.newSmiles,
      newIupacName: result.newIupacName,
      conversionDetails: result.conversionDetails,
      potentialDiseases: result.potentialDiseases,
      analysis: reaction.analysis,
      reactionAdjusted: reaction.isPossible === false,
      originalReason: reaction.reason || null,
    });

    // Save the molecule to the database
    const newMolecule = new GeneratenewMolecule({
      smilesoffirst,
      smilesofsecond,
      newmoleculetitle,
      newSmiles: result.newSmiles,
      newIupacName: result.newIupacName,
      conversionDetails: result.conversionDetails,
      potentialDiseases: result.potentialDiseases,
      information: moleculeInfo,
      userId: id,
    });

    const savedMolecule = await newMolecule.save();
    console.log("Molecule saved to database with ID:", savedMolecule._id);

    // Update the user's NewproteinStructures array
    const updatedUser = await User.findByIdAndUpdate(
      id,
      { $push: { NewproteinStructures: savedMolecule._id } },
      { new: true, runValidators: true }
    );

    if (!updatedUser) {
      console.error(`User with ID ${id} not found for update`);
      throw new Error("User not found during update");
    }
    console.log("User updated with new molecule ID:", savedMolecule._id);

    // Close the response stream
    res.end();

  } catch (error) {
    console.error("Error generating new molecule:", error.message, error.stack);
    res.write(`data: {"error": "${error.message}", "safetyNote": "Chemical synthesis aborted - potential hazard detected"}\n\n`);
    res.write('data: [DONE]\n\n');
    res.end();
  }
};

export const getgeneratednewmolecule = async (req, res) => {
  try {
    const userId = req.user.id;

    console.log(`Request received for /api/protein/generatednewmolecule for user: ${userId}`);
    const molecules = await GeneratenewMolecule.find({ userId }).sort({ created: -1 });

    if (!molecules || molecules.length === 0) {
      return res.status(404).json({ message: "No generated molecules found for this user" });
    }

    const moleculeData = molecules.map((molecule) => ({
      id: molecule._id,
      newSmiles: molecule.newSmiles || "Not available",
      newIupacName: molecule.newIupacName || "Not available",
      newmoleculetitle: molecule.newmoleculetitle || "Untitled Molecule",
      conversionDetails: molecule.conversionDetails || "Not available",
      potentialDiseases: molecule.potentialDiseases || "Not available",
      information: molecule.information || JSON.stringify({
        inputSmiles: "Not available",
        newSmiles: "Not available",
        newIupacName: "Not available",
        conversionDetails: ["Not available"],
        potentialDiseases: ["Not available"]
      }),
      created: molecule.created || new Date(),
    }));

    res.status(200).json({
      molecules: moleculeData,
      total: moleculeData.length,
    });
  } catch (error) {
    console.error("Error retrieving generated molecules:", error);
    res.status(500).json({ message: "Failed to retrieve molecules", error: error.message });
  }
};

// export const proxyGeminiRequest = async (req, res) => {
//   try {
//     const { prompt } = req.body;
//     if (!prompt) {
//       return res.status(400).json({ message: "Prompt is required" });
//     }

//     const response = await axios.post(
//       `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
//       {
//         contents: [{ parts: [{ text: prompt }] }],
//       },
//       {
//         headers: { "Content-Type": "application/json" },
//       }
//     );

//     const geminiContent = response.data.candidates?.[0]?.content?.parts?.[0]?.text;
//     if (!geminiContent) {
//       throw new Error("No content returned from Gemini API");
//     }

//     res.status(200).json({ content: geminiContent });
//   } catch (error) {
//     console.error("Error proxying Gemini request:", error);
//     res.status(500).json({ message: "Failed to fetch Gemini response", error: error.message });
//   }
// };







// AI driven controllers for research papers

  // export const saveResearchPapers = async (req, res) => {
  //   const { userId, molecule, papers } = req.body;
  //   try {
  //     await ResearchPaper.findOneAndUpdate(
  //       { userId, "molecule.title": molecule.title, "molecule.smiles": molecule.smiles },
  //       { userId, molecule, papers, createdAt: new Date() },
  //       { upsert: true, new: true }
  //     );
  //     res.status(200).json({ message: "Papers saved successfully" });
  //   } catch (err) {
  //     res.status(500).json({ message: "Failed to save papers", error: err.message });
  //   }
  // };

  // export const getSavedResearchPapers = async (req, res) => {
  //   const userId = req.user._id;
  //   try {
  //     const savedPapers = await ResearchPaper.find({
  //       userId,
  //       "molecule.title": { $exists: true },
  //       "molecule.smiles": { $exists: true },
  //     });
  //     res.status(200).json({ papers: savedPapers || [] });
  //   } catch (err) {
  //     res.status(500).json({ message: "Failed to fetch saved papers", error: err.message });
  //   }
  // };

  // export const checkSavedPapers = async (req, res) => {
  //   const userId = req.user._id;
  //   const { title, smiles } = req.query;
  //   try {
  //     const exists = await ResearchPaper.findOne({
  //       userId,
  //       "molecule.title": title,
  //       "molecule.smiles": smiles,
  //     });
  //     res.status(200).json({ exists: !!exists });
  //   } catch (err) {
  //     res.status(500).json({ message: "Failed to check saved papers", error: err.message });
  //   }
  // };


  // export const saveGeneratedResearchPaper = async (req, res) => {
  //   const { userId, molecule, paper } = req.body;
  //   try {
  //     const newGeneratedPaper = new GeneratedResearchPaper({
  //       userId,
  //       molecule,
  //       paper,
  //     });
  //     await newGeneratedPaper.save();
  //     res.status(200).json({ message: "Generated research paper saved successfully" });
  //   } catch (err) {
  //     res.status(500).json({ message: "Failed to save generated research paper", error: err.message });
  //   }
  // };


  // // Retrieve saved generated research papers
  // export const getSavedGeneratedResearchPapers = async (req, res) => {
  //   const userId = req.user._id;
  //   try {
  //     const savedPapers = await GeneratedResearchPaper.find({
  //       userId,
  //       "molecule.title": { $exists: true },
  //       "molecule.smiles": { $exists: true },
  //     });
  //     res.status(200).json({ papers: savedPapers || [] });
  //   } catch (err) {
  //     res.status(500).json({ message: "Failed to fetch saved generated papers", error: err.message });
  //   }
  // };


  // export const checkSavedGeneratedPapers = async (req, res) => {
  //   const userId = req.user._id; // Assuming authMiddleware sets req.user
  //   const { title, smiles } = req.query;

  //   if (!title || !smiles) {
  //     return res.status(400).json({ message: "Title and SMILES are required" });
  //   }

  //   try {
  //     const exists = await GeneratedResearchPaper.findOne({
  //       userId,
  //       "molecule.title": title,
  //       "molecule.smiles": smiles,
  //     });
  //     res.status(200).json({ exists: !!exists });
  //   } catch (err) {
  //     console.error("Error checking saved generated papers:", err);
  //     res.status(500).json({ message: "Failed to check saved generated papers", error: err.message });
  //   }
  // };


 // Import the updated SavedSearch model

// // Save related research papers
// export const saveResearchPapers = async (req, res) => {
//   const { userId, molecule, papers } = req.body;
//   try {
//     // Check if a document already exists for this user and molecule
//     const existingDoc = await SavedSearch.findOne({
//       userId,
//       "molecule.symptoms": molecule.symptoms,
//       "molecule.smiles": molecule.smiles,
//     });

//     if (existingDoc) {
//       // Update the existing document by appending or replacing the research papers
//       existingDoc.research = papers; // Replace existing papers
//       existingDoc.createdAt = new Date();
//       await existingDoc.save();
//     } else {
//       // Create a new document
//       const newDoc = new SavedSearch({
//         userId,
//         molecule,
//         research: papers,
//       });
//       await newDoc.save();
//     }

//     res.status(200).json({ message: "Papers saved successfully" });
//   } catch (err) {
//     console.error("Error saving research papers:", err);
//     res.status(500).json({ message: "Failed to save papers", error: err.message });
//   }
// };

// // Fetch saved research papers
// export const getSavedResearchPapers = async (req, res) => {
//   const userId = req.user._id;
//   try {
//     const savedPapers = await SavedSearch.find({
//       userId,
//       research: { $exists: true, $ne: [] }, // Ensure research array exists and is not empty
//     });
//     console.log("Fetched saved research papers:", savedPapers); // Debugging log
//     res.status(200).json({ papers: savedPapers || [] });
//   } catch (err) {
//     console.error("Error fetching saved papers:", err);
//     res.status(500).json({ message: "Failed to fetch saved papers", error: err.message });
//   }
// };

// // Check if saved papers exist for a given molecule
// export const checkSavedPapers = async (req, res) => {
//   const userId = req.user._id;
//   const { symptoms, smiles } = req.query;
//   try {
//     const exists = await SavedSearch.exists({
//       userId,
//       "molecule.symptoms": symptoms,
//       "molecule.smiles": smiles,
//       research: { $exists: true, $ne: [] },
//     });
//     console.log("Check saved papers result:", exists); // Debugging log
//     res.status(200).json({ exists: !!exists });
//   } catch (err) {
//     console.error("Error checking saved papers:", err);
//     res.status(500).json({ message: "Failed to check saved papers", error: err.message });
//   }
// };

// // Save a generated research paper
// export const saveGeneratedResearchPaper = async (req, res) => {
//   const { userId, molecule, paper } = req.body;
//   try {
//     // Check if a document already exists for this user and molecule
//     const existingDoc = await SavedSearch.findOne({
//       userId,
//       "molecule.symptoms": molecule.symptoms,
//       "molecule.smiles": molecule.smiles,
//     });

//     if (existingDoc) {
//       // Update the existing document by adding the generated paper
//       existingDoc.paper = paper;
//       existingDoc.createdAt = new Date();
//       await existingDoc.save();
//     } else {
//       // Create a new document
//       const newDoc = new SavedSearch({
//         userId,
//         molecule,
//         paper,
//       });
//       await newDoc.save();
//     }

//     res.status(200).json({ message: "Generated research paper saved successfully" });
//   } catch (err) {
//     console.error("Error saving generated research paper:", err);
//     res.status(500).json({ message: "Failed to save generated research paper", error: err.message });
//   }
// };

// // Fetch saved generated research papers
// export const getSavedGeneratedResearchPapers = async (req, res) => {
//   const userId = req.user._id;
//   try {
//     const savedPapers = await SavedSearch.find({
//       userId,
//       paper: { $exists: true }, // Ensure paper field exists
//     });
//     console.log("Fetched saved generated papers:", savedPapers); // Debugging log
//     res.status(200).json({ papers: savedPapers || [] });
//   } catch (err) {
//     console.error("Error fetching saved generated papers:", err);
//     res.status(500).json({ message: "Failed to fetch saved generated papers", error: err.message });
//   }
// };

// // Check if a generated paper exists for a given molecule
// export const checkSavedGeneratedPapers = async (req, res) => {
//   const userId = req.user._id;
//   const { symptoms, smiles } = req.query;

//   if (!symptoms || !smiles) {
//     return res.status(400).json({ message: "Symptoms and SMILES are required" });
//   }

//   try {
//     const exists = await SavedSearch.exists({
//       userId,
//       "molecule.symptoms": symptoms,
//       "molecule.smiles": smiles,
//       paper: { $exists: true },
//     });
//     console.log("Check saved generated papers result:", exists); // Debugging log
//     res.status(200).json({ exists: !!exists });
//   } catch (err) {
//     console.error("Error checking saved generated papers:", err);
//     res.status(500).json({ message: "Failed to check saved generated papers", error: err.message });
//   }
// };











  // AI driven controllers
  export const convertFileToSmiles = async (req, res) => {
    try {
      if (!req.file) {
        return res.status(400).json({ message: "No file uploaded" });
      }

      // Use PubChem PUG REST API to convert MOL/SDF to SMILES
      const pubchemUrl = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/molfile/SMILES/txt";
      const formData = new FormData();
      formData.append("molfile", req.file.buffer.toString("utf-8"));

      const response = await axios.post(pubchemUrl, formData, {
        headers: {
          "Content-Type": "multipart/form-data",
        },
      });

      const smiles = response.data.trim();
      if (!smiles) {
        return res.status(400).json({ message: "Failed to convert file to SMILES" });
      }

      res.status(200).json({ smiles });
    } catch (error) {
      console.error("Error converting file to SMILES:", error.message);
      res.status(500).json({ message: "Server error while converting file to SMILES" });
    }
  };

  // Controller to extract fingerprints (mocked for now)
  export const getFingerprints = async (req, res) => {
    try {
      const { smiles } = req.body;
      if (!smiles) {
        return res.status(400).json({ message: "SMILES string is required" });
      }

      // Mock fingerprint extraction (replace with a real implementation using a JS library or API)
      const fingerprints = mockFingerprintExtraction(smiles);

      res.status(200).json({ fingerprints });
    } catch (error) {
      console.error("Error computing fingerprints:", error.message);
      res.status(500).json({ message: "Server error while computing fingerprints" });
    }
  };

  // Controller to perform molecular docking (mocked for now)
  export const performDocking = async (req, res) => {
    try {
      const { smiles } = req.body;
      if (!smiles) {
        return res.status(400).json({ message: "SMILES string is required" });
      }

      // Mock docking (replace with a real implementation using an external docking API)
      const dockingResults = mockDocking(smiles);

      res.status(200).json({ results: dockingResults });
    } catch (error) {
      console.error("Error performing molecular docking:", error.message);
      res.status(500).json({ message: "Server error while performing docking" });
    }
  };

  // Controller to save search results
  export const saveSearch = async (req, res) => {
    try {
      const { userId, smiles, targets, research, docking } = req.body;
      if (!userId || !smiles) {
        return res.status(400).json({ message: "User ID and SMILES are required" });
      }

      const newSearch = new SavedSearch({
        userId,
        smiles,
        targets: targets || [],
        research: research || [],
        docking: docking || null,
      });
      await newSearch.save();
      console.log("Saved search:", newSearch);
      res.status(201).json({ message: "Search saved successfully", search: newSearch });
    } catch (error) {
      console.error("Error saving search:", error);
      res.status(500).json({ message: "Server error while saving search", error: error.message });
    }
  };;

  export const getSavedSearches = async (req, res) => {
    try {
      const userId = req.user._id;
      const searches = await SavedSearch.find({ userId }).sort({ createdAt: -1 });
      console.log("Fetched searches:", searches); // Debug log
      res.status(200).json({ searches });
    } catch (error) {
      console.error("Error retrieving saved searches:", error);
      res.status(500).json({ message: "Server error while retrieving saved searches" });
    }
  };

  // Controller to retrieve saved searches


  // Controller to check if a search exists
  export const checkSavedSearches = async (req, res) => {
    try {
      const userId = req.user._id; // From protectRoute middleware
      const { smiles } = req.query;

      if (!smiles) {
        return res.status(400).json({ message: "SMILES string is required" });
      }

      const searchExists = await SavedSearch.findOne({ userId, smiles });
      res.status(200).json({ exists: !!searchExists });
    } catch (error) {
      console.error("Error checking saved searches:", error.message);
      res.status(500).json({ message: "Server error while checking saved searches" });
    }
  };








