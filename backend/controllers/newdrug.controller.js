import { PredictDisease, TargetProtein } from "../models/newdrug.model.js";
import fs from "fs";
import path from "path";
import mongoose from "mongoose";
import {MongoClient} from 'mongodb';
import { fileURLToPath } from "url";
import { GoogleGenAI } from "@google/genai";

// Get __dirname in ES modules
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Initialize Google Gemini AI
const ai = new GoogleGenAI({
  apiKey: process.env.GEMINI_API_KEY,
});
const config = {
  responseMimeType: "application/json", // Changed to JSON to match the new prompt output
};
const model = "gemini-2.5-flash-preview-05-20";

// Helper function to call Gemini API
const callGemini = async (prompt) => {
  const contents = [
    {
      role: "user",
      parts: [
        {
          text: prompt,
        },
      ],
    },
  ];

  const response = await ai.models.generateContentStream({
    model,
    config,
    contents,
  });

  let fullText = "";
  for await (const chunk of response) {
    fullText += chunk.text || "";
  }
  return fullText;
};

// Helper function to parse CSV data
const parseCSV = (filePath) => {
  const fileContent = fs.readFileSync(filePath, "utf-8");
  const lines = fileContent.split("\n").filter((line) => line.trim() !== "");
  const headers = lines[0].split(",").map((header) => header.trim());
  const rows = lines.slice(1).map((line) => {
    const values = line.split(",").map((value) => value.trim());
    const row = {};
    headers.forEach((header, index) => {
      row[header] = values[index] || "";
    });
    return row;
  });
  return rows;
};

// Helper function to parse structured Gemini JSON response
const parseStructuredResponse = (response, category) => {
  try {
    const jsonResponse = JSON.parse(response);
    const items = jsonResponse[category] || [];

    if (category === "knownOrPotentialTargetProteins") {
      return items.map((item) => ({
        protein: item.protein || "Unknown Protein",
        source: item.source || "N/A",
        summary: item.summary || "No summary available.",
        url: item.url || "",
      }));
    }

    return items.map((item) => ({
      source: item.source || "N/A",
      summary: item.summary || "No summary available.",
      url: item.url || "",
    }));
  } catch (error) {
    console.error(`Error parsing ${category}:`, error);
    return [
      {
        source: "N/A",
        summary: "No data available due to parsing error.",
        url: "",
      },
    ];
  }
};

// export const predictDisease = async (req, res) => {
//   try {
//     const { symptoms } = req.body;

//     if (!symptoms || !Array.isArray(symptoms) || symptoms.length === 0) {
//       return res.status(400).json({ error: "Symptoms array is required" });
//     }

//     // Load and parse DiseaseAndSymptoms.csv
//     const diseaseSymptomPath = path.join(
//       __dirname,
//       "../Datasets/DiseaseAndSymptoms.csv"
//     );
//     if (!fs.existsSync(diseaseSymptomPath)) {
//       return res
//         .status(500)
//         .json({ error: "DiseaseAndSymptoms.csv file not found" });
//     }
//     const diseaseSymptomData = parseCSV(diseaseSymptomPath);

//     // Group symptom sets by disease (handle multiple rows for the same disease)
//     const diseaseSymptomSets = {};
//     diseaseSymptomData.forEach((row) => {
//       const disease = row.Disease;
//       if (!disease) return;

//       const diseaseSymptoms = Object.values(row)
//         .slice(1)
//         .filter((symptom) => symptom && symptom !== "");

//       if (diseaseSymptoms.length === 0) return;

//       if (!diseaseSymptomSets[disease]) {
//         diseaseSymptomSets[disease] = [];
//       }
//       diseaseSymptomSets[disease].push(diseaseSymptoms);
//     });

//     // Convert input symptoms to a set for comparison
//     const inputSymptomSet = new Set(symptoms);

//     // Find exact matches and calculate match percentages
//     const diseaseMatches = [];
//     for (const [disease, symptomSets] of Object.entries(diseaseSymptomSets)) {
//       let bestMatch = { percentage: 0, matchedSymptoms: [], allSymptoms: [] };

//       // Check each symptom set for the disease
//       for (const diseaseSymptoms of symptomSets) {
//         const diseaseSymptomSet = new Set(diseaseSymptoms);

//         // Check for exact match (same symptoms, same count, ignoring order)
//         if (
//           inputSymptomSet.size === diseaseSymptomSet.size &&
//           [...inputSymptomSet].every((symptom) => diseaseSymptomSet.has(symptom))
//         ) {
//           bestMatch = {
//             percentage: 100,
//             matchedSymptoms: diseaseSymptoms,
//             allSymptoms: diseaseSymptoms,
//           };
//           break; // Exact match found, no need to check other sets for this disease
//         }

//         // Calculate partial match
//         const matchedSymptoms = diseaseSymptoms.filter((symptom) =>
//           symptoms.includes(symptom)
//         );
//         // Match percentage: proportion of input symptoms matched and disease symptoms matched
//         const inputMatchRatio = matchedSymptoms.length / symptoms.length;
//         const diseaseMatchRatio = matchedSymptoms.length / diseaseSymptoms.length;
//         // Use the minimum of the two ratios to ensure balanced matching
//         const matchPercentage = Math.min(inputMatchRatio, diseaseMatchRatio) * 100;

//         if (matchPercentage > bestMatch.percentage) {
//           bestMatch = {
//             percentage: matchPercentage,
//             matchedSymptoms,
//             allSymptoms: diseaseSymptoms,
//           };
//         }
//       }

//       diseaseMatches.push({
//         disease,
//         matchPercentage: bestMatch.percentage.toFixed(2) + "%",
//         matchedSymptoms: bestMatch.matchedSymptoms,
//         allSymptoms: bestMatch.allSymptoms,
//       });
//     }

//     // Sort diseases by match percentage (descending)
//     const sortedDiseases = diseaseMatches.sort(
//       (a, b) => parseFloat(b.matchPercentage) - parseFloat(a.matchPercentage)
//     );

//     const topDisease = sortedDiseases[0];
//     if (!topDisease) {
//       return res.status(404).json({ error: "No matching disease found" });
//     }

//     // Fetch precautions if 100% match
//     let precautions = [];
//     if (parseFloat(topDisease.matchPercentage) === 100) {
//       const precautionPath = path.join(
//         __dirname,
//         "../Datasets/DiseasePrecaution.csv"
//       );
//       if (!fs.existsSync(precautionPath)) {
//         console.warn("DiseasePrecaution.csv file not found");
//       } else {
//         const precautionData = parseCSV(precautionPath);
//         const precautionRow = precautionData.find(
//           (row) => row.Disease === topDisease.disease
//         );
//         if (precautionRow) {
//           precautions = [
//             precautionRow.Precaution_1,
//             precautionRow.Precaution_2,
//             precautionRow.Precaution_3,
//             precautionRow.Precaution_4,
//           ].filter((precaution) => precaution && precaution !== "");
//         }
//       }
//     }

//     // If no precautions found, set a default to satisfy the schema
//     if (precautions.length === 0) {
//       precautions = ["No specific precautions available"];
//     }

//     // Gemini prompt and response (unchanged for brevity)
//     const geminiPrompt = `Prompt Title: AI-Assisted Biomedical Disease Mapping and Target Discovery Using Web-Scraped Evidence

// Objective:
// You are an advanced biomedical AI assistant integrated within an AI drug discovery platform. Your task is to analyze a set of symptoms and a predicted disease and produce a **comprehensive biological and molecular analysis**. Your analysis should be **deeply researched and backed by web scraping or web search**, using authoritative sources like PubMed, UniProt, KEGG, Reactome, DisGeNET, DrugBank, Ensembl, GO Ontology, Human Cell Atlas, and FDA biomarker databases.

// Input:
// - Comma-separated list of patient symptoms: "${symptoms.join(", ")}"
// - Top-matched disease based on local ML model: "${topDisease.disease}"

// Your output will be used to:
// - Drive downstream protein-ligand mapping
// - Identify potential targets for AI-based SMILES generation
// - Guide molecular docking and drug design simulations

// ### MANDATORY TASKS (All Based on Live Web Scraping or Smart Web Search):

// For **each of the following categories**, you must:

// - **Search and scrape real biomedical sources**
// - **Write long, detailed, medically sound summaries of the scraped data as you can **
// - **Always include the source name and exact working URL**
// - If no credible data is found after scraping, say: "No reliable data found after attempted web search." — but only after trying.

// ### Output Format (JSON, No Markdown or Extra Text):

// {
//   "likelyDiseasesOrDisorders": [
//     {
//       "source": "PubMed / WHO / CDC / Mayo Clinic / MedlinePlus",
//       "summary": "Describe diseases or disorders linked to the given symptoms. Include prevalence, etiology, and disease mechanism. Write 2–4 medical-grade sentences summarizing current research.",
//       "url": "Direct source link"
//     }
//   ],
//   "associatedBiologicalPathways": [
//     {
//       "source": "KEGG / Reactome / WikiPathways",
//       "summary": "Describe cellular or metabolic pathways affected by the disease. Include pathway name, function, and how it is disrupted. Minimum 3 lines.",
//       "url": "Exact working link"
//     }
//   ],
//   "affectedBiologicalProcesses": [
//     {
//       "source": "GO Ontology / Reactome / PubMed",
//       "summary": "Explain key biological processes affected—like inflammation, apoptosis, synaptic signaling, etc. Describe how these are impaired or dysregulated in the disease context. Use proper terminology.",
//       "url": "Relevant scientific link"
//     }
//   ],
//   "involvedOrgansAndTissues": [
//     {
//       "source": "Human Protein Atlas / Biomedical Literature",
//       "summary": "List and describe affected organs or tissues. Explain how symptoms map to these body systems (e.g., vestibular apparatus in vertigo). Be anatomical and precise.",
//       "url": "Source link to anatomical/biomedical data"
//     }
//   ],
//   "relevantCellTypes": [
//     {
//       "source": "Human Cell Atlas / Research Articles",
//       "summary": "Name immune, neural, or structural cell types affected. Describe their normal role and what goes wrong in the disease. Include 2+ sentences.",
//       "url": "Working link to source"
//     }
//   ],
//   "molecularCellularMechanismDisruptions": [
//     {
//       "source": "PubMed / Medline / Molecular Biology Databases",
//       "summary": "Describe in detail the molecular mechanisms disrupted (e.g., ion channel malfunction, oxidative stress, neurotransmitter imbalance). Be mechanistic, not vague.",
//       "url": "Precise research link"
//     }
//   ],
//   "associatedGenes": [
//     {
//       "source": "DisGeNET / Ensembl / NCBI Gene / OMIM",
//       "summary": "List genes implicated in the disease. For each, explain gene function and how mutations or dysregulation relate to symptoms or pathogenesis. 2–3 sentences minimum.",
//       "url": "Gene-specific reference link"
//     }
//   ],
//   "knownOrPotentialTargetProteins": [
//     {
//       "source": "UniProt / ChEMBL / DrugBank / STRING DB",
//       "summary": "List known or predicted druggable proteins. For each, explain its function, role in the disease, and why it is therapeutically relevant. Link to official database entry.",
//       "url": "Exact UniProt/DrugBank link"
//     }
//   ],
//   "biomarkers": [
//     {
//       "source": "FDA / PubMed / ClinicalTrials / Biomarker Databases",
//       "summary": "List validated or investigational biomarkers. Describe what is being measured (gene/protein/metabolite), its clinical significance, and current research status.",
//       "url": "Direct link to biomarker data"
//     }
//   ],
//   "relevantTherapeuticClassesOrDrugExamples": [
//     {
//       "source": "DrugBank / RxNorm / Medscape / PubMed",
//       "summary": "Provide examples of drugs or drug classes used to treat the disease. For each, explain its mechanism of action and clinical use. Include generics and brand names if relevant.",
//       "url": "Working medical or pharmacological source"
//     }
//   ]
// }

// ### FINAL INSTRUCTIONS:
// - All responses MUST be accurate and verifiable with scientific citations.
// - Do NOT use placeholder text like "N/A" unless web search **explicitly fails**.
// - DO NOT hallucinate knowledge—back it up with web content.
// - JSON output only. No extra text, markdown, or explanation outside JSON.
// - This will feed into a biomedical target-ligand pipeline. Scientific precision is critical.
// `;

//     const geminiResponse = (await callGemini(geminiPrompt)) || "{}";

//     // Parse the Gemini response to match the DiseaseAnalysis schema
//     const diseaseAnalysis = {
//       likelyDiseasesOrDisorders: parseStructuredResponse(
//         geminiResponse,
//         "likelyDiseasesOrDisorders"
//       ),
//       associatedBiologicalPathways: parseStructuredResponse(
//         geminiResponse,
//         "associatedBiologicalPathways"
//       ),
//       affectedBiologicalProcesses: parseStructuredResponse(
//         geminiResponse,
//         "affectedBiologicalProcesses"
//       ),
//       involvedOrgansAndTissues: parseStructuredResponse(
//         geminiResponse,
//         "involvedOrgansAndTissues"
//       ),
//       relevantCellTypes: parseStructuredResponse(
//         geminiResponse,
//         "relevantCellTypes"
//       ),
//       molecularCellularMechanismDisruptions: parseStructuredResponse(
//         geminiResponse,
//         "molecularCellularMechanismDisruptions"
//       ),
//       associatedGenes: parseStructuredResponse(
//         geminiResponse,
//         "associatedGenes"
//       ),
//       knownOrPotentialTargetProteins: parseStructuredResponse(
//         geminiResponse,
//         "knownOrPotentialTargetProteins"
//       ),
//       biomarkers: parseStructuredResponse(geminiResponse, "biomarkers"),
//       relevantTherapeuticClassesOrDrugExamples: parseStructuredResponse(
//         geminiResponse,
//         "relevantTherapeuticClassesOrDrugExamples"
//       ),
//     };

//     // Ensure all DiseaseAnalysis fields have at least one entry
//     Object.keys(diseaseAnalysis).forEach((key) => {
//       if (diseaseAnalysis[key].length === 0) {
//         diseaseAnalysis[key] = [
//           {
//             source: "N/A",
//             summary: "No data available from web sources.",
//             url: "",
//           },
//         ];
//       }
//     });

//     // Save to MongoDB
//     const prediction = new PredictDisease({
//       symptoms,
//       predictedDiseases: [
//         {
//           diseaseName: topDisease.disease,
//           DiseaseMatchness: topDisease.matchPercentage,
//           diseaseCautions: precautions,
//         },
//       ],
//       DiseaseAnalysis: diseaseAnalysis,
//     });

//     await prediction.save();

//     // Respond with the result
//     res.status(200).json({
//       symptoms,
//       predictedDiseases: [
//         {
//           diseaseName: topDisease.disease,
//           DiseaseMatchness: topDisease.matchPercentage,
//           diseaseCautions: precautions,
//         },
//       ],
//       DiseaseAnalysis: diseaseAnalysis,
//     });
//   } catch (error) {
//     console.error("Error in predictDisease:", error);
//     res.status(500).json({ error: "Internal server error" });
//   }
// };

// Controller for predictDisease route
// export const predictDisease = async (req, res) => {
//   try {
//     const { symptoms } = req.body;

//     if (!symptoms || !Array.isArray(symptoms) || symptoms.length === 0) {
//       return res.status(400).json({ error: "Symptoms array is required" });
//     }

//     // Load and parse DiseaseAndSymptoms.csv
//     const diseaseSymptomPath = path.join(
//       __dirname,
//       "../Datasets/DiseaseAndSymptoms.csv"
//     );
//     if (!fs.existsSync(diseaseSymptomPath)) {
//       return res
//         .status(500)
//         .json({ error: "DiseaseAndSymptoms.csv file not found" });
//     }
//     const diseaseSymptomData = parseCSV(diseaseSymptomPath);

//     // Group symptom sets by disease (handle multiple rows for the same disease)
//     const diseaseSymptomSets = {};
//     diseaseSymptomData.forEach((row) => {
//       const disease = row.Disease;
//       if (!disease) return;

//       const diseaseSymptoms = Object.values(row)
//         .slice(1)
//         .filter((symptom) => symptom && symptom !== "");

//       if (diseaseSymptoms.length === 0) return;

//       if (!diseaseSymptomSets[disease]) {
//         diseaseSymptomSets[disease] = [];
//       }
//       diseaseSymptomSets[disease].push(diseaseSymptoms);
//     });

//     // Convert input symptoms to a set for comparison
//     const inputSymptomSet = new Set(symptoms);

//     // Find exact matches and calculate match percentages
//     const diseaseMatches = [];
//     for (const [disease, symptomSets] of Object.entries(diseaseSymptomSets)) {
//       let bestMatch = { percentage: 0, matchedSymptoms: [], allSymptoms: [] };

//       // Check each symptom set for the disease
//       for (const diseaseSymptoms of symptomSets) {
//         const diseaseSymptomSet = new Set(diseaseSymptoms);

//         // Check for exact match (same symptoms, same count, ignoring order)
//         if (
//           inputSymptomSet.size === diseaseSymptomSet.size &&
//           [...inputSymptomSet].every((symptom) => diseaseSymptomSet.has(symptom))
//         ) {
//           bestMatch = {
//             percentage: 100,
//             matchedSymptoms: diseaseSymptoms,
//             allSymptoms: diseaseSymptoms,
//           };
//           break; // Exact match found, no need to check other sets for this disease
//         }

//         // Calculate partial match
//         const matchedSymptoms = diseaseSymptoms.filter((symptom) =>
//           symptoms.includes(symptom)
//         );
//         // Match percentage: proportion of input symptoms matched and disease symptoms matched
//         const inputMatchRatio = matchedSymptoms.length / symptoms.length;
//         const diseaseMatchRatio = matchedSymptoms.length / diseaseSymptoms.length;
//         // Use the minimum of the two ratios to ensure balanced matching
//         const matchPercentage = Math.min(inputMatchRatio, diseaseMatchRatio) * 100;

//         if (matchPercentage > bestMatch.percentage) {
//           bestMatch = {
//             percentage: matchPercentage,
//             matchedSymptoms,
//             allSymptoms: diseaseSymptoms,
//           };
//         }
//       }

//       diseaseMatches.push({
//         disease,
//         matchPercentage: bestMatch.percentage.toFixed(2) + "%",
//         matchedSymptoms: bestMatch.matchedSymptoms,
//         allSymptoms: bestMatch.allSymptoms,
//       });
//     }

//     // Sort diseases by match percentage (descending)
//     const sortedDiseases = diseaseMatches.sort(
//       (a, b) => parseFloat(b.matchPercentage) - parseFloat(a.matchPercentage)
//     );

//     const topDisease = sortedDiseases[0];
//     if (!topDisease) {
//       return res.status(404).json({ error: "No matching disease found" });
//     }

//     // Check if the top match is 0.00%
//     let predictedDisease;
//     let diseaseAnalysis = {};
//     let precautions = [];

//     if (parseFloat(topDisease.matchPercentage) === 0) {
//       // No match found in the dataset
//       predictedDisease = {
//         diseaseName: "NO disease found on datasets its a new disease",
//         DiseaseMatchness: "0.00%",
//         diseaseCautions: ["No precautions available as this may be a new disease"],
//       };

//       // Populate DiseaseAnalysis with a placeholder
//       diseaseAnalysis = {
//         likelyDiseasesOrDisorders: [
//           {
//             source: "N/A",
//             summary: "No matching disease found in the dataset. The symptoms may indicate a new or uncharacterized disease.",
//             url: "",
//           },
//         ],
//         associatedBiologicalPathways: [
//           {
//             source: "N/A",
//             summary: "No pathways identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         affectedBiologicalProcesses: [
//           {
//             source: "N/A",
//             summary: "No biological processes identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         involvedOrgansAndTissues: [
//           {
//             source: "N/A",
//             summary: "No organs or tissues identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         relevantCellTypes: [
//           {
//             source: "N/A",
//             summary: "No cell types identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         molecularCellularMechanismDisruptions: [
//           {
//             source: "N/A",
//             summary: "No molecular mechanisms identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         associatedGenes: [
//           {
//             source: "N/A",
//             summary: "No genes identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         knownOrPotentialTargetProteins: [
//           {
//             source: "N/A",
//             summary: "No target proteins identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         biomarkers: [
//           {
//             source: "N/A",
//             summary: "No biomarkers identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         relevantTherapeuticClassesOrDrugExamples: [
//           {
//             source: "N/A",
//             summary: "No therapeutic classes or drugs identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//       };
//     } else {
//       // Proceed with normal matching
//       // Fetch precautions if 100% match
//       if (parseFloat(topDisease.matchPercentage) === 100) {
//         const precautionPath = path.join(
//           __dirname,
//           "../Datasets/DiseasePrecaution.csv"
//         );
//         if (!fs.existsSync(precautionPath)) {
//           console.warn("DiseasePrecaution.csv file not found");
//         } else {
//           const precautionData = parseCSV(precautionPath);
//           const precautionRow = precautionData.find(
//             (row) => row.Disease === topDisease.disease
//           );
//           if (precautionRow) {
//             precautions = [
//               precautionRow.Precaution_1,
//               precautionRow.Precaution_2,
//               precautionRow.Precaution_3,
//               precautionRow.Precaution_4,
//             ].filter((precaution) => precaution && precaution !== "");
//           }
//         }
//       }

//       // If no precautions found, set a default to satisfy the schema
//       if (precautions.length === 0) {
//         precautions = ["No specific precautions available"];
//       }

//       predictedDisease = {
//         diseaseName: topDisease.disease,
//         DiseaseMatchness: topDisease.matchPercentage,
//         diseaseCautions: precautions,
//       };

//       // Gemini prompt for DiseaseAnalysis
//       const geminiPrompt = `Prompt Title: AI-Assisted Biomedical Disease Mapping and Target Discovery Using Web-Scraped Evidence

// Objective:
// You are an advanced biomedical AI assistant integrated within an AI drug discovery platform. Your task is to analyze a set of symptoms and a predicted disease and produce a **comprehensive biological and molecular analysis**. Your analysis should be **deeply researched and backed by web scraping or web search**, using authoritative sources like PubMed, UniProt, KEGG, Reactome, DisGeNET, DrugBank, Ensembl, GO Ontology, Human Cell Atlas, and FDA biomarker databases.

// Input:
// - Comma-separated list of patient symptoms: "${symptoms.join(", ")}"
// - Top-matched disease based on local ML model: "${topDisease.disease}"

// Your output will be used to:
// - Drive downstream protein-ligand mapping
// - Identify potential targets for AI-based SMILES generation
// - Guide molecular docking and drug design simulations

// ### MANDATORY TASKS (All Based on Live Web Scraping or Smart Web Search):

// For **each of the following categories**, you must:

// - **Search and scrape real biomedical sources**
// - **Write long, detailed, medically sound summaries of the scraped data as you can **
// - **Always include the source name and exact working URL**
// - If no credible data is found after scraping, say: "No reliable data found after attempted web search." — but only after trying.

// ### Output Format (JSON, No Markdown or Extra Text):

// {
//   "likelyDiseasesOrDisorders": [
//     {
//       "source": "PubMed / WHO / CDC / Mayo Clinic / MedlinePlus",
//       "summary": "Describe diseases or disorders linked to the given symptoms. Include prevalence, etiology, and disease mechanism. Write 2–4 medical-grade sentences summarizing current research.",
//       "url": "Direct source link"
//     }
//   ],
//   "associatedBiologicalPathways": [
//     {
//       "source": "KEGG / Reactome / WikiPathways",
//       "summary": "Describe cellular or metabolic pathways affected by the disease. Include pathway name, function, and how it is disrupted. Minimum 3 lines.",
//       "url": "Exact working link"
//     }
//   ],
//   "affectedBiologicalProcesses": [
//     {
//       "source": "GO Ontology / Reactome / PubMed",
//       "summary": "Explain key biological processes affected—like inflammation, apoptosis, synaptic signaling, etc. Describe how these are impaired or dysregulated in the disease context. Use proper terminology.",
//       "url": "Relevant scientific link"
//     }
//   ],
//   "involvedOrgansAndTissues": [
//     {
//       "source": "Human Protein Atlas / Biomedical Literature",
//       "summary": "List and describe affected organs or tissues. Explain how symptoms map to these body systems (e.g., vestibular apparatus in vertigo). Be anatomical and precise.",
//       "url": "Source link to anatomical/biomedical data"
//     }
//   ],
//   "relevantCellTypes": [
//     {
//       "source": "Human Cell Atlas / Research Articles",
//       "summary": "Name immune, neural, or structural cell types affected. Describe their normal role and what goes wrong in the disease. Include 2+ sentences.",
//       "url": "Working link to source"
//     }
//   ],
//   "molecularCellularMechanismDisruptions": [
//     {
//       "source": "PubMed / Medline / Molecular Biology Databases",
//       "summary": "Describe in detail the molecular mechanisms disrupted (e.g., ion channel malfunction, oxidative stress, neurotransmitter imbalance). Be mechanistic, not vague.",
//       "url": "Precise research link"
//     }
//   ],
//   "associatedGenes": [
//     {
//       "source": "DisGeNET / Ensembl / NCBI Gene / OMIM",
//       "summary": "List genes implicated in the disease. For each, explain gene function and how mutations or dysregulation relate to symptoms or pathogenesis. 2–3 sentences minimum.",
//       "url": "Gene-specific reference link"
//     }
//   ],
//   "knownOrPotentialTargetProteins": [
//     {
//       "source": "UniProt / ChEMBL / DrugBank / STRING DB",
//       "summary": "List known or predicted druggable proteins. For each, explain its function, role in the disease, and why it is therapeutically relevant. Link to official database entry.",
//       "url": "Exact UniProt/DrugBank link"
//     }
//   ],
//   "biomarkers": [
//     {
//       "source": "FDA / PubMed / ClinicalTrials / Biomarker Databases",
//       "summary": "List validated or investigational biomarkers. Describe what is being measured (gene/protein/metabolite), its clinical significance, and current research status.",
//       "url": "Direct link to biomarker data"
//     }
//   ],
//   "relevantTherapeuticClassesOrDrugExamples": [
//     {
//       "source": "DrugBank / RxNorm / Medscape / PubMed",
//       "summary": "Provide examples of drugs or drug classes used to treat the disease. For each, explain its mechanism of action and clinical use. Include generics and brand names if relevant.",
//       "url": "Working medical or pharmacological source"
//     }
//   ]
// }

// ### FINAL INSTRUCTIONS:
// - All responses MUST be accurate and verifiable with scientific citations.
// - Do NOT use placeholder text like "N/A" unless web search **explicitly fails**.
// - DO NOT hallucinate knowledge—back it up with web content.
// - JSON output only. No extra text, markdown, or explanation outside JSON.
// - This will feed into a biomedical target-ligand pipeline. Scientific precision is critical.
// `;

//       const geminiResponse = (await callGemini(geminiPrompt)) || "{}";

//       // Parse the Gemini response to match the DiseaseAnalysis schema
//       diseaseAnalysis = {
//         likelyDiseasesOrDisorders: parseStructuredResponse(
//           geminiResponse,
//           "likelyDiseasesOrDisorders"
//         ),
//         associatedBiologicalPathways: parseStructuredResponse(
//           geminiResponse,
//           "associatedBiologicalPathways"
//         ),
//         affectedBiologicalProcesses: parseStructuredResponse(
//           geminiResponse,
//           "affectedBiologicalProcesses"
//         ),
//         involvedOrgansAndTissues: parseStructuredResponse(
//           geminiResponse,
//           "involvedOrgansAndTissues"
//         ),
//         relevantCellTypes: parseStructuredResponse(
//           geminiResponse,
//           "relevantCellTypes"
//         ),
//         molecularCellularMechanismDisruptions: parseStructuredResponse(
//           geminiResponse,
//           "molecularCellularMechanismDisruptions"
//         ),
//         associatedGenes: parseStructuredResponse(
//           geminiResponse,
//           "associatedGenes"
//         ),
//         knownOrPotentialTargetProteins: parseStructuredResponse(
//           geminiResponse,
//           "knownOrPotentialTargetProteins"
//         ),
//         biomarkers: parseStructuredResponse(geminiResponse, "biomarkers"),
//         relevantTherapeuticClassesOrDrugExamples: parseStructuredResponse(
//           geminiResponse,
//           "relevantTherapeuticClassesOrDrugExamples"
//         ),
//       };

//       // Ensure all DiseaseAnalysis fields have at least one entry
//       Object.keys(diseaseAnalysis).forEach((key) => {
//         if (diseaseAnalysis[key].length === 0) {
//           diseaseAnalysis[key] = [
//             {
//               source: "N/A",
//               summary: "No data available from web sources.",
//               url: "",
//             },
//           ];
//         }
//       });
//     }

//     // Save to MongoDB
//     const prediction = new PredictDisease({
//       symptoms,
//       predictedDiseases: [predictedDisease],
//       DiseaseAnalysis: diseaseAnalysis,
//     });

//     await prediction.save();

//     // Respond with the result
//     res.status(200).json({
//       symptoms,
//       predictedDiseases: [predictedDisease],
//       DiseaseAnalysis: diseaseAnalysis,
//     });
//   } catch (error) {
//     console.error("Error in predictDisease:", error);
//     res.status(500).json({ error: "Internal server error" });
//   }
// };

// Controller for predictDisease route

// export const predictDisease = async (req, res) => {
//   try {
//     const { symptoms } = req.body;

//     if (!symptoms || !Array.isArray(symptoms) || symptoms.length === 0) {
//       return res.status(400).json({ error: "Symptoms array is required" });
//     }

//     // Load and parse DiseaseAndSymptoms.csv
//     const diseaseSymptomPath = path.join(
//       __dirname,
//       "../Datasets/DiseaseAndSymptoms.csv"
//     );
//     if (!fs.existsSync(diseaseSymptomPath)) {
//       return res
//         .status(500)
//         .json({ error: "DiseaseAndSymptoms.csv file not found" });
//     }
//     const diseaseSymptomData = parseCSV(diseaseSymptomPath);

//     // Group symptom sets by disease (handle multiple rows for the same disease)
//     const diseaseSymptomSets = {};
//     diseaseSymptomData.forEach((row) => {
//       const disease = row.Disease;
//       if (!disease) return;

//       const diseaseSymptoms = Object.values(row)
//         .slice(1)
//         .filter((symptom) => symptom && symptom !== "");

//       if (diseaseSymptoms.length === 0) return;

//       if (!diseaseSymptomSets[disease]) {
//         diseaseSymptomSets[disease] = [];
//       }
//       diseaseSymptomSets[disease].push(diseaseSymptoms);
//     });

//     // Convert input symptoms to a set for comparison
//     const inputSymptomSet = new Set(symptoms);

//     // Find exact matches and calculate match percentages
//     const diseaseMatches = [];
//     for (const [disease, symptomSets] of Object.entries(diseaseSymptomSets)) {
//       let bestMatch = { percentage: 0, matchedSymptoms: [], allSymptoms: [] };

//       // Check each symptom set for the disease
//       for (const diseaseSymptoms of symptomSets) {
//         const diseaseSymptomSet = new Set(diseaseSymptoms);

//         // Check for exact match (same symptoms, same count, ignoring order)
//         if (
//           inputSymptomSet.size === diseaseSymptomSet.size &&
//           [...inputSymptomSet].every((symptom) => diseaseSymptomSet.has(symptom))
//         ) {
//           bestMatch = {
//             percentage: 100,
//             matchedSymptoms: diseaseSymptoms,
//             allSymptoms: diseaseSymptoms,
//           };
//           break; // Exact match found, no need to check other sets for this disease
//         }

//         // Calculate partial match
//         const matchedSymptoms = diseaseSymptoms.filter((symptom) =>
//           symptoms.includes(symptom)
//         );
//         // Match percentage: proportion of input symptoms matched and disease symptoms matched
//         const inputMatchRatio = matchedSymptoms.length / symptoms.length;
//         const diseaseMatchRatio = matchedSymptoms.length / diseaseSymptoms.length;
//         // Use the minimum of the two ratios to ensure balanced matching
//         const matchPercentage = Math.min(inputMatchRatio, diseaseMatchRatio) * 100;

//         if (matchPercentage > bestMatch.percentage) {
//           bestMatch = {
//             percentage: matchPercentage,
//             matchedSymptoms,
//             allSymptoms: diseaseSymptoms,
//           };
//         }
//       }

//       diseaseMatches.push({
//         disease,
//         matchPercentage: bestMatch.percentage.toFixed(2) + "%",
//         matchedSymptoms: bestMatch.matchedSymptoms,
//         allSymptoms: bestMatch.allSymptoms,
//       });
//     }

//     // Sort diseases by match percentage (descending)
//     const sortedDiseases = diseaseMatches.sort(
//       (a, b) => parseFloat(b.matchPercentage) - parseFloat(a.matchPercentage)
//     );

//     const topDisease = sortedDiseases[0];
//     if (!topDisease) {
//       return res.status(404).json({ error: "No matching disease found" });
//     }

//     // Check if the top match is 0.00%, and adjust match percentage if > 0.00%
//     let predictedDisease;
//     let diseaseAnalysis = {};
//     let precautions = [];

//     const topMatchPercentage = parseFloat(topDisease.matchPercentage);
//     let adjustedMatchPercentage = topMatchPercentage;

//     // Adjust match percentage: if > 0.00%, set to 100.00%
//     if (topMatchPercentage > 0) {
//       adjustedMatchPercentage = 100.00;
//       topDisease.matchPercentage = "100.00%"; // Update for use in Gemini prompt
//     }

//     if (topMatchPercentage === 0) {
//       // No match found in the dataset
//       predictedDisease = {
//         diseaseName: "NO disease found on datasets its a new disease",
//         DiseaseMatchness: "0.00%",
//         diseaseCautions: ["No precautions available as this may be a new disease"],
//       };

//       // Populate DiseaseAnalysis with a placeholder
//       diseaseAnalysis = {
//         likelyDiseasesOrDisorders: [
//           {
//             source: "N/A",
//             summary: "No matching disease found in the dataset. The symptoms may indicate a new or uncharacterized disease.",
//             url: "",
//           },
//         ],
//         associatedBiologicalPathways: [
//           {
//             source: "N/A",
//             summary: "No pathways identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         affectedBiologicalProcesses: [
//           {
//             source: "N/A",
//             summary: "No biological processes identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         involvedOrgansAndTissues: [
//           {
//             source: "N/A",
//             summary: "No organs or tissues identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         relevantCellTypes: [
//           {
//             source: "N/A",
//             summary: "No cell types identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         molecularCellularMechanismDisruptions: [
//           {
//             source: "N/A",
//             summary: "No molecular mechanisms identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         associatedGenes: [
//           {
//             source: "N/A",
//             summary: "No genes identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         knownOrPotentialTargetProteins: [
//           {
//             source: "N/A",
//             summary: "No target proteins identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         biomarkers: [
//           {
//             source: "N/A",
//             summary: "No biomarkers identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//         relevantTherapeuticClassesOrDrugExamples: [
//           {
//             source: "N/A",
//             summary: "No therapeutic classes or drugs identified due to lack of matched disease in the dataset.",
//             url: "",
//           },
//         ],
//       };
//     } else {
//       // Proceed with normal matching, with adjusted match percentage
//       // Fetch precautions (now always fetched since match is treated as 100%)
//       const precautionPath = path.join(
//         __dirname,
//         "../Datasets/DiseasePrecaution.csv"
//       );
//       if (!fs.existsSync(precautionPath)) {
//         console.warn("DiseasePrecaution.csv file not found");
//       } else {
//         const precautionData = parseCSV(precautionPath);
//         const precautionRow = precautionData.find(
//           (row) => row.Disease === topDisease.disease
//         );
//         if (precautionRow) {
//           precautions = [
//             precautionRow.Precaution_1,
//             precautionRow.Precaution_2,
//             precautionRow.Precaution_3,
//             precautionRow.Precaution_4,
//           ].filter((precaution) => precaution && precaution !== "");
//         }
//       }

//       // If no precautions found, set a default to satisfy the schema
//       if (precautions.length === 0) {
//         precautions = ["No specific precautions available"];
//       }

//       predictedDisease = {
//         diseaseName: topDisease.disease,
//         DiseaseMatchness: "100.00%", // Adjusted to 100%
//         diseaseCautions: precautions,
//       };

//       // Gemini prompt for DiseaseAnalysis
//       const geminiPrompt = `Prompt Title: AI-Assisted Biomedical Disease Mapping and Target Discovery Using Web-Scraped Evidence

// Objective:
// You are an advanced biomedical AI assistant integrated within an AI drug discovery platform. Your task is to analyze a set of symptoms and a predicted disease and produce a **comprehensive biological and molecular analysis**. Your analysis should be **deeply researched and backed by web scraping or web search**, using authoritative sources like PubMed, UniProt, KEGG, Reactome, DisGeNET, DrugBank, Ensembl, GO Ontology, Human Cell Atlas, and FDA biomarker databases.

// Input:
// - Comma-separated list of patient symptoms: "${symptoms.join(", ")}"
// - Top-matched disease based on local ML model: "${topDisease.disease}"

// Your output will be used to:
// - Drive downstream protein-ligand mapping
// - Identify potential targets for AI-based SMILES generation
// - Guide molecular docking and drug design simulations

// ### MANDATORY TASKS (All Based on Live Web Scraping or Smart Web Search):

// For **each of the following categories**, you must:

// - **Search and scrape real biomedical sources**
// - **Write long, detailed, medically sound summaries of the scraped data as you can **
// - **Always include the source name and exact working URL**
// - If no credible data is found after scraping, say: "No reliable data found after attempted web search." — but only after trying.

// ### Output Format (JSON, No Markdown or Extra Text):

// {
//   "likelyDiseasesOrDisorders": [
//     {
//       "source": "PubMed / WHO / CDC / Mayo Clinic / MedlinePlus",
//       "summary": "Describe diseases or disorders linked to the given symptoms. Include prevalence, etiology, and disease mechanism. Write 2–4 medical-grade sentences summarizing current research.",
//       "url": "Direct source link"
//     }
//   ],
//   "associatedBiologicalPathways": [
//     {
//       "source": "KEGG / Reactome / WikiPathways",
//       "summary": "Describe cellular or metabolic pathways affected by the disease. Include pathway name, function, and how it is disrupted. Minimum 3 lines.",
//       "url": "Exact working link"
//     }
//   ],
//   "affectedBiologicalProcesses": [
//     {
//       "source": "GO Ontology / Reactome / PubMed",
//       "summary": "Explain key biological processes affected—like inflammation, apoptosis, synaptic signaling, etc. Describe how these are impaired or dysregulated in the disease context. Use proper terminology.",
//       "url": "Relevant scientific link"
//     }
//   ],
//   "involvedOrgansAndTissues": [
//     {
//       "source": "Human Protein Atlas / Biomedical Literature",
//       "summary": "List and describe affected organs or tissues. Explain how symptoms map to these body systems (e.g., vestibular apparatus in vertigo). Be anatomical and precise.",
//       "url": "Source link to anatomical/biomedical data"
//     }
//   ],
//   "relevantCellTypes": [
//     {
//       "source": "Human Cell Atlas / Research Articles",
//       "summary": "Name immune, neural, or structural cell types affected. Describe their normal role and what goes wrong in the disease. Include 2+ sentences.",
//       "url": "Working link to source"
//     }
//   ],
//   "molecularCellularMechanismDisruptions": [
//     {
//       "source": "PubMed / Medline / Molecular Biology Databases",
//       "summary": "Describe in detail the molecular mechanisms disrupted (e.g., ion channel malfunction, oxidative stress, neurotransmitter imbalance). Be mechanistic, not vague.",
//       "url": "Precise research link"
//     }
//   ],
//   "associatedGenes": [
//     {
//       "source": "DisGeNET / Ensembl / NCBI Gene / OMIM",
//       "summary": "List genes implicated in the disease. For each, explain gene function and how mutations or dysregulation relate to symptoms or pathogenesis. 2–3 sentences minimum.",
//       "url": "Gene-specific reference link"
//     }
//   ],
//   "knownOrPotentialTargetProteins": [
//     {
//       "source": "UniProt / ChEMBL / DrugBank / STRING DB",
//       "summary": "List known or predicted druggable proteins. For each, explain its function, role in the disease, and why it is therapeutically relevant. Link to official database entry.",
//       "url": "Exact UniProt/DrugBank link"
//     }
//   ],
//   "biomarkers": [
//     {
//       "source": "FDA / PubMed / ClinicalTrials / Biomarker Databases",
//       "summary": "List validated or investigational biomarkers. Describe what is being measured (gene/protein/metabolite), its clinical significance, and current research status.",
//       "url": "Direct link to biomarker data"
//     }
//   ],
//   "relevantTherapeuticClassesOrDrugExamples": [
//     {
//       "source": "DrugBank / RxNorm / Medscape / PubMed",
//       "summary": "Provide examples of drugs or drug classes used to treat the disease. For each, explain its mechanism of action and clinical use. Include generics and brand names if relevant.",
//       "url": "Working medical or pharmacological source"
//     }
//   ]
// }

// ### FINAL INSTRUCTIONS:
// - All responses MUST be accurate and verifiable with scientific citations.
// - Do NOT use placeholder text like "N/A" unless web search **explicitly fails**.
// - DO NOT hallucinate knowledge—back it up with web content.
// - JSON output only. No extra text, markdown, or explanation outside JSON.
// - This will feed into a biomedical target-ligand pipeline. Scientific precision is critical.
// `;

//       const geminiResponse = (await callGemini(geminiPrompt)) || "{}";

//       // Parse the Gemini response to match the DiseaseAnalysis schema
//       diseaseAnalysis = {
//         likelyDiseasesOrDisorders: parseStructuredResponse(
//           geminiResponse,
//           "likelyDiseasesOrDisorders"
//         ),
//         associatedBiologicalPathways: parseStructuredResponse(
//           geminiResponse,
//           "associatedBiologicalPathways"
//         ),
//         affectedBiologicalProcesses: parseStructuredResponse(
//           geminiResponse,
//           "affectedBiologicalProcesses"
//         ),
//         involvedOrgansAndTissues: parseStructuredResponse(
//           geminiResponse,
//           "involvedOrgansAndTissues"
//         ),
//         relevantCellTypes: parseStructuredResponse(
//           geminiResponse,
//           "relevantCellTypes"
//         ),
//         molecularCellularMechanismDisruptions: parseStructuredResponse(
//           geminiResponse,
//           "molecularCellularMechanismDisruptions"
//         ),
//         associatedGenes: parseStructuredResponse(
//           geminiResponse,
//           "associatedGenes"
//         ),
//         knownOrPotentialTargetProteins: parseStructuredResponse(
//           geminiResponse,
//           "knownOrPotentialTargetProteins"
//         ),
//         biomarkers: parseStructuredResponse(geminiResponse, "biomarkers"),
//         relevantTherapeuticClassesOrDrugExamples: parseStructuredResponse(
//           geminiResponse,
//           "relevantTherapeuticClassesOrDrugExamples"
//         ),
//       };

//       // Ensure all DiseaseAnalysis fields have at least one entry
//       Object.keys(diseaseAnalysis).forEach((key) => {
//         if (diseaseAnalysis[key].length === 0) {
//           diseaseAnalysis[key] = [
//             {
//               source: "N/A",
//               summary: "No data available from web sources.",
//               url: "",
//             },
//           ];
//         }
//       });
//     }

//     // Save to MongoDB
//     const prediction = new PredictDisease({
//       symptoms,
//       predictedDiseases: [predictedDisease],
//       DiseaseAnalysis: diseaseAnalysis,
//     });

//     await prediction.save();

//     // Respond with the result
//     res.status(200).json({
//       symptoms,
//       predictedDiseases: [predictedDisease],
//       DiseaseAnalysis: diseaseAnalysis,
//     });
//   } catch (error) {
//     console.error("Error in predictDisease:", error);
//     res.status(500).json({ error: "Internal server error" });
//   }
// };

export const predictDisease = async (req, res) => {
  try {
    const { symptoms } = req.body;
    const userId = req.params.id; // Use req.params.id as defined in the route // Extract userId (adjust based on your auth middleware or route)

    // Validate userId
    if (!userId || !mongoose.Types.ObjectId.isValid(userId)) {
      return res.status(400).json({ error: 'Invalid or missing user ID' });
    }

    // Validate symptoms
    if (!symptoms || !Array.isArray(symptoms) || symptoms.length === 0) {
      return res.status(400).json({ error: 'Symptoms array is required' });
    }

    // Load and parse DiseaseAndSymptoms.csv
    const diseaseSymptomPath = path.join(
      __dirname,
      '../Datasets/DiseaseAndSymptoms.csv'
    );
    if (!fs.existsSync(diseaseSymptomPath)) {
      return res.status(500).json({ error: 'DiseaseAndSymptoms.csv file not found' });
    }
    const diseaseSymptomData = parseCSV(diseaseSymptomPath);

    // Group symptom sets by disease
    const diseaseSymptomSets = {};
    diseaseSymptomData.forEach((row) => {
      const disease = row.Disease;
      if (!disease) return;

      const diseaseSymptoms = Object.values(row)
        .slice(1)
        .filter((symptom) => symptom && symptom !== '');

      if (diseaseSymptoms.length === 0) return;

      if (!diseaseSymptomSets[disease]) {
        diseaseSymptomSets[disease] = [];
      }
      diseaseSymptomSets[disease].push(diseaseSymptoms);
    });

    // Convert input symptoms to a set for comparison
    const inputSymptomSet = new Set(symptoms);

    // Find exact matches and calculate match percentages
    const diseaseMatches = [];
    for (const [disease, symptomSets] of Object.entries(diseaseSymptomSets)) {
      let bestMatch = { percentage: 0, matchedSymptoms: [], allSymptoms: [] };

      for (const diseaseSymptoms of symptomSets) {
        const diseaseSymptomSet = new Set(diseaseSymptoms);

        // Check for exact match
        if (
          inputSymptomSet.size === diseaseSymptomSet.size &&
          [...inputSymptomSet].every((symptom) => diseaseSymptomSet.has(symptom))
        ) {
          bestMatch = {
            percentage: 100,
            matchedSymptoms: diseaseSymptoms,
            allSymptoms: diseaseSymptoms,
          };
          break;
        }

        // Calculate partial match
        const matchedSymptoms = diseaseSymptoms.filter((symptom) =>
          symptoms.includes(symptom)
        );
        const inputMatchRatio = matchedSymptoms.length / symptoms.length;
        const diseaseMatchRatio = matchedSymptoms.length / diseaseSymptoms.length;
        const matchPercentage = Math.min(inputMatchRatio, diseaseMatchRatio) * 100;

        if (matchPercentage > bestMatch.percentage) {
          bestMatch = {
            percentage: matchPercentage,
            matchedSymptoms,
            allSymptoms: diseaseSymptoms,
          };
        }
      }

      diseaseMatches.push({
        disease,
        matchPercentage: bestMatch.percentage.toFixed(2) + '%',
        matchedSymptoms: bestMatch.matchedSymptoms,
        allSymptoms: bestMatch.allSymptoms,
      });
    }

    // Sort diseases by match percentage
    const sortedDiseases = diseaseMatches.sort(
      (a, b) => parseFloat(b.matchPercentage) - parseFloat(a.matchPercentage)
    );

    const topDisease = sortedDiseases[0];
    if (!topDisease) {
      return res.status(404).json({ error: 'No matching disease found' });
    }

    let predictedDisease;
    let diseaseAnalysis = {};
    let precautions = [];

    const topMatchPercentage = parseFloat(topDisease.matchPercentage);
    let adjustedMatchPercentage = topMatchPercentage;

    if (topMatchPercentage > 0) {
      adjustedMatchPercentage = 100.00;
      topDisease.matchPercentage = '100.00%';
    }

    if (topMatchPercentage === 0) {
      predictedDisease = {
        diseaseName: 'NO disease found on datasets its a new disease',
        DiseaseMatchness: '0.00%',
        diseaseCautions: ['No precautions available as this may be a new disease'],
      };

      diseaseAnalysis = {
        likelyDiseasesOrDisorders: [
          {
            source: 'N/A',
            summary: 'No matching disease found in the dataset. The symptoms may indicate a new or uncharacterized disease.',
            url: '',
          },
        ],
        associatedBiologicalPathways: [
          {
            source: 'N/A',
            summary: 'No pathways identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        affectedBiologicalProcesses: [
          {
            source: 'N/A',
            summary: 'No biological processes identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        involvedOrgansAndTissues: [
          {
            source: 'N/A',
            summary: 'No organs or tissues identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        relevantCellTypes: [
          {
            source: 'N/A',
            summary: 'No cell types identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        molecularCellularMechanismDisruptions: [
          {
            source: 'N/A',
            summary: 'No molecular mechanisms identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        associatedGenes: [
          {
            source: 'N/A',
            summary: 'No genes identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        knownOrPotentialTargetProteins: [
          {
            source: 'N/A',
            summary: 'No target proteins identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        biomarkers: [
          {
            source: 'N/A',
            summary: 'No biomarkers identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
        relevantTherapeuticClassesOrDrugExamples: [
          {
            source: 'N/A',
            summary: 'No therapeutic classes or drugs identified due to lack of matched disease in the dataset.',
            url: '',
          },
        ],
      };
    } else {
      // Fetch precautions
      const precautionPath = path.join(
        __dirname,
        '../Datasets/DiseasePrecaution.csv'
      );
      if (!fs.existsSync(precautionPath)) {
        console.warn('DiseasePrecaution.csv file not found');
      } else {
        const precautionData = parseCSV(precautionPath);
        const precautionRow = precautionData.find(
          (row) => row.Disease === topDisease.disease
        );
        if (precautionRow) {
          precautions = [
            precautionRow.Precaution_1,
            precautionRow.Precaution_2,
            precautionRow.Precaution_3,
            precautionRow.Precaution_4,
          ].filter((precaution) => precaution && precaution !== '');
        }
      }

      if (precautions.length === 0) {
        precautions = ['No specific precautions available'];
      }

      predictedDisease = {
        diseaseName: topDisease.disease,
        DiseaseMatchness: '100.00%',
        diseaseCautions: precautions,
      };

      // Gemini prompt for DiseaseAnalysis (unchanged)
      const geminiPrompt = `Prompt Title: AI-Assisted Biomedical Disease Mapping and Target Discovery Using Web-Scraped Evidence

Objective:
You are an advanced biomedical AI assistant integrated within an AI drug discovery platform. Your task is to analyze a set of symptoms and a predicted disease and produce a *comprehensive biological and molecular analysis. Your analysis should be **deeply researched and backed by web scraping or web search*, using authoritative sources like PubMed, UniProt, KEGG, Reactome, DisGeNET, DrugBank, Ensembl, GO Ontology, Human Cell Atlas, and FDA biomarker databases.

Input:
- Comma-separated list of patient symptoms: "${symptoms.join(", ")}"
- Top-matched disease based on local ML model: "${topDisease.disease}"

Your output will be used to:
- Drive downstream protein-ligand mapping
- Identify potential targets for AI-based SMILES generation
- Guide molecular docking and drug design simulations

### MANDATORY TASKS (All Based on Live Web Scraping or Smart Web Search):

For *each of the following categories*, you must:

- *Search and scrape real biomedical sources*
- *Write long, detailed, medically sound summaries (minimum 2–4 sentences each)*
- *Always include the source name and exact working URL*
- If no credible data is found after scraping, say: "No reliable data found after attempted web search." — but only after trying.

### Output Format (JSON, No Markdown or Extra Text):

{
  "likelyDiseasesOrDisorders": [
    {
      "source": "PubMed / WHO / CDC / Mayo Clinic / MedlinePlus",
      "summary": "Describe diseases or disorders linked to the given symptoms. Include prevalence, etiology, and disease mechanism. Write 2–4 medical-grade sentences summarizing current research.",
      "url": "Direct source link"
    }
  ],
  "associatedBiologicalPathways": [
    {
      "source": "KEGG / Reactome / WikiPathways",
      "summary": "Describe cellular or metabolic pathways affected by the disease. Include pathway name, function, and how it is disrupted. Minimum 3 lines.",
      "url": "Exact working link"
    }
  ],
  "affectedBiologicalProcesses": [
    {
      "source": "GO Ontology / Reactome / PubMed",
      "summary": "Explain key biological processes affected—like inflammation, apoptosis, synaptic signaling, etc. Describe how these are impaired or dysregulated in the disease context. Use proper terminology.",
      "url": "Relevant scientific link"
    }
  ],
  "involvedOrgansAndTissues": [
    {
      "source": "Human Protein Atlas / Biomedical Literature",
      "summary": "List and describe affected organs or tissues. Explain how symptoms map to these body systems (e.g., vestibular apparatus in vertigo). Be anatomical and precise.",
      "url": "Source link to anatomical/biomedical data"
    }
  ],
  "relevantCellTypes": [
    {
      "source": "Human Cell Atlas / Research Articles",
      "summary": "Name immune, neural, or structural cell types affected. Describe their normal role and what goes wrong in the disease. Include 2+ sentences.",
      "url": "Working link to source"
    }
  ],
  "molecularCellularMechanismDisruptions": [
    {
      "source": "PubMed / Medline / Molecular Biology Databases",
      "summary": "Describe in detail the molecular mechanisms disrupted (e.g., ion channel malfunction, oxidative stress, neurotransmitter imbalance). Be mechanistic, not vague.",
      "url": "Precise research link"
    }
  ],
  "associatedGenes": [
    {
      "source": "DisGeNET / Ensembl / NCBI Gene / OMIM",
      "summary": "List genes implicated in the disease. For each, explain gene function and how mutations or dysregulation relate to symptoms or pathogenesis. 2–3 sentences minimum.",
      "url": "Gene-specific reference link"
    }
  ],
  "knownOrPotentialTargetProteins": [
    {
      "source": "UniProt / ChEMBL / DrugBank / STRING DB",
      "summary": "List known or predicted druggable proteins. For each, explain its function, role in the disease, and why it is therapeutically relevant. Link to official database entry.",
      "url": "Exact UniProt/DrugBank link"
    }
  ],
  "biomarkers": [
    {
      "source": "FDA / PubMed / ClinicalTrials / Biomarker Databases",
      "summary": "List validated or investigational biomarkers. Describe what is being measured (gene/protein/metabolite), its clinical significance, and current research status.",
      "url": "Direct link to biomarker data"
    }
  ],
  "relevantTherapeuticClassesOrDrugExamples": [
    {
      "source": "DrugBank / RxNorm / Medscape / PubMed",
      "summary": "Provide examples of drugs or drug classes used to treat the disease. For each, explain its mechanism of action and clinical use. Include generics and brand names if relevant.",
      "url": "Working medical or pharmacological source"
    }
  ]
}

### FINAL INSTRUCTIONS:
- All responses MUST be accurate and verifiable with scientific citations.
- Do NOT use placeholder text like "N/A" unless web search *explicitly fails*.
- DO NOT hallucinate knowledge—back it up with web content.
- JSON output only. No extra text, markdown, or explanation outside JSON.
- This will feed into a biomedical target-ligand pipeline. Scientific precision is critical.
`;
      const geminiResponse = (await callGemini(geminiPrompt)) || '{}';

      diseaseAnalysis = {
        likelyDiseasesOrDisorders: parseStructuredResponse(
          geminiResponse,
          'likelyDiseasesOrDisorders'
        ),
        associatedBiologicalPathways: parseStructuredResponse(
          geminiResponse,
          'associatedBiologicalPathways'
        ),
        affectedBiologicalProcesses: parseStructuredResponse(
          geminiResponse,
          'affectedBiologicalProcesses'
        ),
        involvedOrgansAndTissues: parseStructuredResponse(
          geminiResponse,
          'involvedOrgansAndTissues'
        ),
        relevantCellTypes: parseStructuredResponse(
          geminiResponse,
          'relevantCellTypes'
        ),
        molecularCellularMechanismDisruptions: parseStructuredResponse(
          geminiResponse,
          'molecularCellularMechanismDisruptions'
        ),
        associatedGenes: parseStructuredResponse(
          geminiResponse,
          'associatedGenes'
        ),
        knownOrPotentialTargetProteins: parseStructuredResponse(
          geminiResponse,
          'knownOrPotentialTargetProteins'
        ),
        biomarkers: parseStructuredResponse(geminiResponse, 'biomarkers'),
        relevantTherapeuticClassesOrDrugExamples: parseStructuredResponse(
          geminiResponse,
          'relevantTherapeuticClassesOrDrugExamples'
        ),
      };

      Object.keys(diseaseAnalysis).forEach((key) => {
        if (diseaseAnalysis[key].length === 0) {
          diseaseAnalysis[key] = [
            {
              source: 'N/A',
              summary: 'No data available from web sources.',
              url: '',
            },
          ];
        }
      });
    }

    // Save to MongoDB with userId
    const prediction = new PredictDisease({
      userId: new mongoose.Types.ObjectId(userId), // Ensure userId is stored
      symptoms,
      predictedDiseases: [predictedDisease],
      DiseaseAnalysis: diseaseAnalysis,
    });

    await prediction.save();

    // Respond with the result
    res.status(200).json({
      symptoms,
      predictedDiseases: [predictedDisease],
      DiseaseAnalysis: diseaseAnalysis,
    });
  } catch (error) {
    console.error('Error in predictDisease:', error);
    res.status(500).json({ error: 'Internal server error' });
  }
};

export const predictTargetProtein = async (req, res) => {
  try {
    // Fetch the latest PredictDisease entry from MongoDB
    const latestPrediction = await PredictDisease.findOne()
      .sort({ createdAt: -1 })
      .exec();

    if (!latestPrediction) {
      return res
        .status(404)
        .json({
          error: "No previous disease prediction found in the database",
        });
    }

    const { symptoms, predictedDiseases } = latestPrediction;
    if (!predictedDiseases || predictedDiseases.length === 0) {
      return res
        .status(400)
        .json({ error: "No predicted diseases found in the database entry" });
    }

    const disease = predictedDiseases[0].diseaseName;
    if (!symptoms || !disease) {
      return res
        .status(400)
        .json({
          error:
            "Required fields (symptoms, diseaseName) are missing in the database entry",
        });
    }

    // Use Gemini to identify target proteins and ligands, matching the TargetProteinSchema
   const geminiPrompt = `Act as a biomedical, bioinformatics, and drug discovery expert.

Objective:
For the associated symptoms (${symptoms.join(
  ", "
)}), identify and present **comprehensive**, **evidence-based**, and **multi-entity** information about all relevant target proteins and their associated therapeutic ligands.
Instructions:

1. Perform an exhaustive search using authoritative biomedical databases and public repositories, including but not limited to:
   - UniProt
   - Protein Data Bank (PDB)
   - PubMed
   - DrugBank
   - ChEMBL
   - DisGeNET
   - GeneCards
   - STRING-DB
   - KEGG
   - Open Targets Platform

2. Search and extract accurate, verifiable, and biologically relevant information from these databases using live queries, scraping, or API access. Use only reliable, peer-reviewed, or curated biomedical sources.

3. For each **target protein** involved in the pathophysiology or biological mechanism of the symptom(s)/disease(s), return greater than 4 unique and scientifically validated proteins**. For **each**, provide:
   - Full name of the protein
   - Detailed biological function and its relevance to the entered symptoms/disease mechanism
   - List of associated diseases or pathological conditions
   - Gene symbol or identifier
   - UniProt ID (if available)
   - PDB ID (if available, otherwise state "Not available")
   - Brief protein description (1–2 sentences) as "ProtienDiscription"
   - Detailed functional and structural description (3–5 sentences), including its cellular localization, pathway involvement, molecular interactions, and relevance to the symptom or condition under investigation (field: "proteinDetailedDiscription")
   - Data source links for traceability
   - Note: Proteins **do not** have SMILES; include "Not applicable" under "ProteinSmile"

4. For each identified target protein, return a minimum of **one therapeutic ligand or drug compound**, with a greater than  4 or more distinct ligands**. For each ligand, provide:
   - Ligand name and type (e.g., antagonist, agonist, inhibitor, modulator)
   - Full pharmacological or biological function in relation to the target protein and symptoms
   - Mechanism of action
   - Associated diseases or symptoms the ligand is known to treat
   - LigandDiscription: A concise (1–2 sentence) overview
   - LigandDetailedDiscription: An expanded (3–5 sentence) description covering mechanism of action, therapeutic relevance, pharmacodynamics/pharmacokinetics, and any clinical insights
   - Valid SMILES string (from DrugBank, PubChem, or ChEMBL); if unavailable, state "Not applicable"
   - DrugBank ID (or state "Not available")
   - Source links used

5. Ensure factual scientific clarity, use domain-specific nomenclature, and provide traceable source links for every field. Avoid summarization — prioritize depth, citation, and comprehensiveness.

 Output JSON Format:
{
  "TargetProteins": [
    {
      "proteinName": "...",
      "proteinFunction": "...",
      "associatedDiseases": ["..."],
      "associatedGene": "...",
      "proteinUniport": "...",
      "pdbID": "...",
      "ProtienDiscription": "...",
      "proteinDetailedDiscription": "...",
      "ProteinSmile": "Not applicable",
      "dataSources": ["..."]
    }
    // Minimum of 4 target proteins required
  ],
  "TargetLigands": [
    {
      "ligandName": "...",
      "ligandFunction": "...",
      "ligandType": "...",
      "associatedDiseases": ["..."],
      "LigandDiscription": "...",
      "LigandDetailedDiscription": "...",
      "LigandSmile": "...",
      "ligandDrugBankID": "...",
      "dataSources": ["..."]
    }
    // Minimum of 4 therapeutic ligands required
  ]
}

 Additional Requirements:
- Do not leave fields blank — if data is missing, use "Not available"
- Do not provide SMILES for proteins
- Include complete citations/URLs for verification
- Do not simplify or generalize any descriptions — be thorough
- All entries should be relevant to the symptom(s) provided
- Provide **at least 4 well-validated entries each** for both proteins and ligands to ensure rich scientific value
`;


    // Call Gemini API with the updated prompt
    const geminiResponse = (await callGemini(geminiPrompt)) || "{}";

    // Parse the Gemini JSON response
    let targetProteins = [];
    let targetLigands = [];

    try {
      const jsonResponse = JSON.parse(geminiResponse);

      // Parse TargetProteins
      targetProteins =
        jsonResponse.TargetProteins?.map((item) => ({
          proteinName: item.proteinName || "Unknown Protein",
          proteinFunction: item.proteinFunction || "Not available",
          associatedDiseases: item.associatedDiseases?.length
            ? item.associatedDiseases
            : [disease],
          ProtienDiscription: item.ProtienDiscription || "Not available",
          proteinDetailedDiscription:
            item.proteinDetailedDiscription || "Not available",
          ProteinSmile: item.ProteinSmile || "Not applicable",
          proteinUniport: item.proteinUniport || "Not available",
        })) || [];

      // Parse TargetLigands
      targetLigands =
        jsonResponse.TargetLigands?.map((item) => ({
          ligandName: item.ligandName || "Unknown Ligand",
          ligandFunction: item.ligandFunction || "Not available",
          associatedDiseases: item.associatedDiseases?.length
            ? item.associatedDiseases
            : [disease],
          LigandDiscription: item.LigandDiscription || "Not available",
          LigandDetailedDiscription:
            item.LigandDetailedDiscription || "Not available",
          LigandSmile: item.LigandSmile || "Not applicable",
          ligandDrugBankID: item.ligandDrugBankID || "Not available",
        })) || [];
    } catch (error) {
      console.error("Error parsing Gemini response:", error);
    }

    // Ensure at least one entry for each section to satisfy the schema
    if (targetProteins.length === 0) {
      targetProteins.push({
        proteinName: "Unknown Protein",
        proteinFunction: "Not available",
        associatedDiseases: [disease],
        ProtienDiscription: "No description available from web sources.",
        proteinDetailedDiscription:
          "No detailed description available from web sources.",
        ProteinSmile: "Not applicable",
        proteinUniport: "Not available",
      });
    }

    if (targetLigands.length === 0) {
      targetLigands.push({
        ligandName: "Unknown Ligand",
        ligandFunction: "Not available",
        associatedDiseases: [disease],
        LigandDiscription: "No description available from web sources.",
        LigandDetailedDiscription:
          "No detailed description available from web sources.",
        LigandSmile: "Not applicable",
        ligandDrugBankID: "Not available",
      });
    }

    // Save to MongoDB
    const targetProteinEntry = new TargetProtein({
      TargetProteins: targetProteins,
      TargetLigands: targetLigands,
      userId: req.user ? req.user._id : null, // Assuming user authentication middleware sets req.user
    });

    await targetProteinEntry.save();

    // Respond with the result
    res.status(200).json({
      disease,
      symptoms,
      TargetProteins: targetProteins,
      TargetLigands: targetLigands,
    });
  } catch (error) {
    console.error("Error in predictTargetProtein:", error);
    res.status(500).json({ error: "Internal server error" });
  }
};


// Controller to fetch product_smiles from reaction_responses collection
export const getnewdrug = async (req, res) => {
    try {
        // MongoDB connection (update URI as per your setup)
        const client = new MongoClient(process.env.MONGODB_URI);
        await client.connect();
        
        const db = client.db('test'); // Replace with your database name
        const collection = db.collection('reaction_responses');
        
        // Query to fetch product_smiles field
        const products = await collection
            .find({}, { projection: { product_smiles: 1, _id: 0 } })
            .toArray();
        
        // Extract product_smiles arrays and flatten them
        const productSmiles = products
            .flatMap(doc => doc.product_smiles)
            .filter(smiles => smiles); // Remove any null/undefined entries
        
        // Close the MongoDB connection
        await client.close();
        
        // Send response
        res.status(200).json({
            success: true,
            data: productSmiles
        });
        console.log('Product smiles fetched successfully:', productSmiles);
    } catch (error) {
        console.error('Error fetching product smiles:', error);
        res.status(500).json({
            success: false,
            message: 'Error fetching product smiles',
            error: error.message
        });
    }
};


export const getSymptoms = async (req, res) => {
  try {
    const userId = req.params.id;

    // Validate ObjectId
    if (!mongoose.Types.ObjectId.isValid(userId)) {
      return res.status(400).json({ message: 'Invalid user ID format' });
    }

    // Find documents by userId and select symptoms field
    const predictions = await PredictDisease.find({ userId }, 'symptoms');

    if (!predictions || predictions.length === 0) {
      return res.status(404).json({ message: 'No symptoms found for this user' });
    }

    // Map predictions to an array of symptom arrays and deduplicate
    const seen = new Set();
    const symptoms = predictions
      .map(prediction => prediction.symptoms)
      .filter(symptomArray => {
        // Convert symptom array to a sorted string for comparison
        const key = symptomArray.sort().join(',');
        if (!seen.has(key)) {
          seen.add(key);
          return true;
        }
        return false;
      });

    res.status(200).json({ symptoms });
  } catch (error) {
    console.error('Error fetching symptoms:', error);
    res.status(500).json({ message: 'Server error while fetching symptoms' });
  }
};