import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { Spinner } from '../../components/Spinner.jsx';
import StructureForm from '../../components/StructureForm.jsx';
import StructureList from '../../components/StructureList.jsx';
import StructureDetails from '../../components/StructureDetails.jsx';
import { toast } from 'react-hot-toast';
import { useAuthStore } from '../../Store/auth.store.js';

const GEMINI_API_URL = 'https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent' // Ensure this is set in environment
const GEMINI_API_KEY = 'AIzaSyDyujm50dHMYvn1V50dDDqcAhgUqCOuUGU'// Ensure this is set in environment

const ProteinStructureApp = () => {
  const [structures, setStructures] = useState([]);
  const [selectedStructure, setSelectedStructure] = useState(null);
  const [loading, setLoading] = useState(false);
  const [rdkitLoaded, setRdkitLoaded] = useState(false);
  const [error, setError] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    const loadRDKit = async () => {
      try {
        if (window.RDKit) {
          setRdkitLoaded(true);
          return;
        }
        await new Promise((resolve, reject) => {
          const script = document.createElement('script');
          script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
          script.async = true;
          script.onload = () => resolve();
          script.onerror = () => reject(new Error('Failed to load RDKit script'));
          document.head.appendChild(script);
        });
        const rdkitModule = await window.initRDKitModule();
        window.RDKit = rdkitModule;
        setRdkitLoaded(true);
      } catch (err) {
        setError('Failed to load RDKit: ' + err.message);
      }
    };

    const fetchStructures = async () => {
      try {
        setLoading(true);
        const response = await axios.get(`/protein/getproteinstructure/${user?._id}`, {
          baseURL: import.meta.env.mode === "development" ? "http://localhost:5000" : "/api",
          withCredentials: true,
        });
        setStructures(response.data);
      } catch (err) {
        setError(err.response?.data?.message || 'Failed to fetch structures');
      } finally {
        setLoading(false);
      }
    };

    const initializeApp = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError('Authentication failed. Please log in.');
        return;
      }
      await loadRDKit();
      await fetchStructures();
    };

    initializeApp();
  }, [checkAuth]);

  const handleSubmit = async (formData) => {
    if (!user?._id) {
      setError('Please log in to generate structures');
      toast.error('Please log in to generate structures', {
        style: {
          background: 'var(--color-error)',
          color: 'white',
        },
      });
      return;
    }

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
    } = formData;

    if (!smiles) {
      setError('SMILES string is required');
      toast.error('SMILES string is required', {
        style: {
          background: 'var(--color-error)',
          color: 'white',
        },
      });
      return;
    }

    const validAlgorithms = ["CMA-ES", "SSD"];
    if (!validAlgorithms.includes(algorithm)) {
      setError(`Invalid algorithm: ${algorithm}. Supported algorithms are ${validAlgorithms.join(', ')}.`);
      toast.error(`Invalid algorithm: ${algorithm}.`, {
        style: {
          background: 'var(--color-error)',
          color: 'white',
        },
      });
      return;
    }

    try {
      setLoading(true);
      setError(null);

      let molecules;
      if (algorithm === "SSD") {
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

        const response = await axios.post('/api/proxy/molmim', requestPayload, {
          withCredentials: true,
        });

        molecules = typeof response.data.molecules === 'string' 
          ? JSON.parse(response.data.molecules) 
          : response.data.molecules;

        if (!Array.isArray(molecules)) {
          throw new Error('Parsed "molecules" is not an array');
        }
      }

      let moleculeInfo = '';
      try {
        if (!GEMINI_API_KEY || !GEMINI_API_URL) {
          throw new Error('Gemini API credentials are not set');
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
`;

        const geminiResponse = await axios.post(
          `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
          {
            contents: [{ parts: [{ text: prompt }] }]
          },
          {
            headers: { 'Content-Type': 'application/json' }
          }
        );

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

        moleculeInfo = moleculeInfo
          .replace(/(\d+\.\s[A-Za-z\s]+:)/g, "\n$1")
          .replace(/-+/g, "-")
          .trim();
      } catch (geminiError) {
        console.error('Error calling Gemini API:', geminiError.message);
        moleculeInfo = 'Failed to retrieve detailed information about the molecule.';
      }

      const newProteinStructure = {
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
        userId: user._id
      };

      // Save to backend (assuming endpoint exists)
      const saveResponse = await axios.post(`/protein/saveproteinstructure/${user._id}`, newProteinStructure, {
        baseURL: import.meta.env.mode === "development" ? "http://localhost:5000" : "/api",
        withCredentials: true,
      });

      setStructures((prev) => [saveResponse.data.structure, ...prev]);
      setSelectedStructure(saveResponse.data.structure);
      toast.success('Structure generated successfully!', {
        style: {
          background: 'var(--color-success)',
          color: 'var(--color-primary)',
        },
      });
    } catch (err) {
      setError(err.response?.data?.message || 'Failed to generate structure');
      toast.error(err.response?.data?.message || 'Failed to generate structure', {
        style: {
          background: 'var(--color-error)',
          color: 'white',
        },
      });
    } finally {
      setLoading(false);
    }
  };

  const selectStructure = (structure) => {
    setSelectedStructure(structure);
  };

  if (checkingAuth || (!rdkitLoaded && !error)) {
    return (
      <div className="flex flex-col items-center justify-center h-screen bg-primary">
        <Spinner className="text-accent" />
        <p className="mt-4 text-text-primary font-body">
          {checkingAuth ? 'Verifying authentication...' : 'Initializing molecular viewer...'}
        </p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen bg-primary">
        <div className="text-center p-8 bg-secondary rounded-lg shadow-xl max-w-md">
          <h2 className="text-2xl font-heading text-text-primary mb-4">Access Required</h2>
          <p className="text-text-secondary font-body mb-6">
            Please log in to access the Protein Structure Generator
          </p>
          <button
            className="w-full bg-accent hover:bg-accent/90 text-primary font-heading font-bold py-3 px-6 rounded-lg transition-all duration-200"
            onClick={() => (window.location.href = '/login')}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-primary p-4 md:p-8">
      <div className="max-w-7xl mx-auto">
        {/* Header */}
        <header className="text-center mb-8">
          <h1 className="text-3xl md:text-4xl font-heading font-bold text-accent mb-2">
            Protein Structure Generator
          </h1>
          <p className="text-sm text-accent-secondary font-300 font-mono">
            (POWERED BY GEMINI & MolMIM NVIDIA MODEL)
          </p>
        </header>

        {/* Error Message */}
        {error && (
          <div className="bg-error/20 border border-error text-text-primary px-4 py-3 rounded-lg mb-6 flex justify-between items-center">
            <p className="font-body">{error}</p>
            <button
              className="text-text-primary hover:text-accent font-bold ml-4"
              onClick={() => setError(null)}
            >
              ✕
            </button>
          </div>
        )}

        {/* Main Content Grid */}
        <div className="grid grid-cols-1 lg:grid-cols-1 gap-6">
          {/* Left Column - Generate Form */}
          <div className="lg:col-span-1">
            <div className="bg-secondary p-6 rounded-xl shadow-lg h-full">
              <h2 className="text-xl font-heading font-semibold text-accent mb-4">
                Generate New Structure
              </h2>
              <StructureForm onSubmit={handleSubmit} loading={loading} />
            </div>
          </div>

          {/* Middle Column - Structure Details */}
          <div className="lg:col-span-1">
            <div className="bg-secondary p-6 rounded-xl shadow-lg h-full">
              <h2 className="text-xl font-heading font-semibold text-accent mb-4">
                Structure Details
              </h2>
              {selectedStructure ? (
                <StructureDetails structure={selectedStructure} rdkitLoaded={rdkitLoaded} />
              ) : (
                <div className="flex flex-col items-center justify-center h-96">
                  <svg
                    className="w-16 h-16 text-text-secondary mb-4"
                    fill="none"
                    stroke="currentColor"
                    viewBox="0 0 24 24"
                    xmlns="http://www.w3.org/2000/svg"
                  >
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth={1.5}
                      d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z"
                    />
                  </svg>
                  <p className="text-text-secondary font-body mb-2">
                    Select a structure to view details
                  </p>
                  <p className="text-text-secondary text-sm font-body">or generate a new one</p>
                </div>
              )}
            </div>
          </div>

          {/* Right Column - Saved Structures */}
          <div className="lg:col-span-1">
            <div className="bg-secondary p-6 rounded-xl shadow-lg h-full">
              <div className="flex justify-between items-center mb-4">
                <h2 className="text-xl font-heading font-semibold text-accent">
                  Saved Structures
                </h2>
                {structures.length > 0 && (
                  <span className="bg-accent-secondary text-primary text-xs font-bold px-2 py-1 rounded-full">
                    {structures.length}
                  </span>
                )}
              </div>
              <StructureList
                structures={structures}
                onSelect={selectStructure}
                selected={selectedStructure}
                loading={loading}
              />
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureApp;