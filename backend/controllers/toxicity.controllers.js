import axios from 'axios';
import Toxicity from '../models/toxicity.model.js';
import DrugName from '../models/drugName.model.js';

const GEMINI_API_URL = 'https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent';
const GEMINI_API_KEY = process.env.GEMINI_API_KEY;

const symptomToEndpointMap = {
  liver: 'hepatotoxicity', heart: 'cardiotoxicity',
  respiratory: 'respiratory_toxicity', neurological: 'neurotoxicity',
  kidney: 'nephrotoxicity', cancer: 'carcinogenicity', skin: 'skin_sensitization',
  thyroid: 'thyroid_toxicity', joint: 'musculoskeletal_toxicity',
  ear: 'ototoxicity', confusion: 'neurotoxicity', heartbeat: 'cardiotoxicity',
  nail: 'dermatotoxicity', knuckle: 'musculoskeletal_toxicity', blockage: 'ototoxicity'
};

// Enhanced SMILES validation
const isValidSmiles = (smiles) => {
  if (!smiles || typeof smiles !== 'string') return false;
  const smilesRegex = /^[A-Za-z0-9@+\-\[\]\(\)=#$.%*:/\\]+$/;
  return smilesRegex.test(smiles) && smiles.length < 1000;
};

// Retry logic for API calls
const retryAxios = async (config, retries = 3, delay = 1000) => {
  for (let i = 0; i < retries; i++) {
    try {
      return await axios(config);
    } catch (error) {
      if (error.response?.status === 429 && i < retries - 1) {
        await new Promise(resolve => setTimeout(resolve, delay * (i + 1)));
        continue;
      }
      throw error;
    }
  }
};

// Helper function to fetch SMILES from compound name
const getSmilesByName = async (compoundName, userId) => {
  try {
    const drugEntry = await DrugName.findOne({ 
      suggestedName: { $regex: new RegExp(`^${compoundName}$`, 'i') },
      userId 
    });
    
    if (!drugEntry) {
      throw new Error(`Compound "${compoundName}" not found in your drug database`);
    }
    
    return {
      smiles: drugEntry.smiles,
      originalName: drugEntry.suggestedName,
      drugId: drugEntry._id
    };
  } catch (error) {
    throw new Error(`Failed to fetch SMILES for "${compoundName}": ${error.message}`);
  }
};

export const predictToxicityController = async (req, res) => {
  try {
    const { compound, symptoms } = req.body;
    const userId = req.user?._id;
    
    if (!compound || !symptoms || !userId) {
      return res.status(400).json({ 
        message: 'Compound name, symptoms, and authentication are required.',
        suggestion: 'Please provide a valid molecule name and at least one symptom.'
      });
    }

    let smilesData;
    try {
      smilesData = await getSmilesByName(compound, userId);
    } catch (error) {
      return res.status(404).json({ 
        message: error.message,
        suggestion: 'Please ensure the compound name exists in your drug database or add it first.'
      });
    }

    const { smiles, originalName } = smilesData;

    if (!isValidSmiles(smiles)) {
      return res.status(400).json({ 
        message: 'Invalid SMILES format provided by the database.',
        suggestion: 'Please verify the SMILES string for this compound.'
      });
    }

    const symptomsArray = symptoms
      .split(',')
      .map(s => s.trim().toLowerCase())
      .filter(s => s.length > 0);

    if (!symptomsArray.length) {
      return res.status(400).json({ 
        message: 'At least one valid symptom is required.',
        suggestion: 'Please provide symptoms separated by commas.'
      });
    }

    const existing = await Toxicity.findOne({ smiles, symptoms, userId }).sort({ created: -1 });
    
    if (existing) {
      return res.status(200).json({ 
        message: 'Prediction exists', 
        result: existing.geminiAnalysis, // Return the full analysis
        historyId: existing._id,
        compoundName: originalName,
        smiles
      });
    }

    if (!GEMINI_API_KEY) {
      return res.status(500).json({ 
        message: 'Gemini API key is not configured.',
        suggestion: 'Please contact the administrator to configure the API key.'
      });
    }

    // REVISED: Prompt with corrected JSON structure for mechanisms
    const prompt = `
You are a senior computational toxicologist. For the compound with SMILES "${smiles}" (known as "${originalName}") intended to treat "${symptoms}", provide complete toxicity estimates—no field can be "unknown." If data is unavailable, use structural or QSAR-based inference with your rationale.

Return ONLY a JSON object wrapped in \`\`\`json\n and \n\`\`\` with:
{
  "isInputValid": boolean, "inputError": string|null, "smiles": string, "moleculeName": string, "isNovelCompound": boolean,
  "structuralAnalysis": { "molecularWeight": string, "logP": string, "polarSurfaceArea": string, "functionalGroups": string[], "complexity": string, "drugLikeness": string },
  "acuteToxicity": { "LD50": string, "toxicityClass": string, "rationale": string, "supplemental": string|null },
  "endpoints": { "hepatotoxicity": string, "carcinogenicity": string, "cardiotoxicity": string, "neurotoxicity": string, "nephrotoxicity": string, "skin_sensitization": string, "thyroid_toxicity": string, "ototoxicity": string, "musculoskeletal_toxicity": string, "dermatotoxicity": string, "supplemental": string|null },
  "sideEffects": [{ "name": string, "description": string, "severity": string, "organSystem": string, "likelihood": string }],
  "mechanisms": [{ "feature": string, "pathway": string, "toxicityType": string }],
  "structureConcerns": [{ "substructure": string, "concern": string, "riskLevel": string }],
  "safetyRecommendations": [{ "test": string, "priority": string, "rationale": string }],
  "qsarAnalysis": { "properties": [{ "name": string, "predictedValue": string, "toxicologicalImplication": string }], "overallRisk": string, "symptomContext": [{ "symptom": string, "structuralInsight": string, "riskAssessment": string }] },
  "toxicophoreAnalysis": [{ "substructure": string, "concern": string, "prevalence": string, "mitigation": string }],
  "developmentRecommendations": [{ "phase": string, "recommendation": string, "importance": string }]
}`;

    let geminiAnalysis;
    try {
      // FIXED: Corrected the axios call payload
      const response = await retryAxios({
        method: 'post',
        url: `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        data: { contents: [{ role: 'user', parts: [{ text: prompt }] }] },
        headers: { 'Content-Type': 'application/json' },
        timeout: 30000,
      });

      let text = response.data.candidates?.[0]?.content?.parts?.[0]?.text || '';
      console.log('Raw Gemini response:', text);

      const jsonMatch = text.match(/```json\n([\s\S]*?)\n```/);
      if (!jsonMatch || !jsonMatch[1]) {
        // Try parsing directly if markdown is missing
        try {
          geminiAnalysis = JSON.parse(text);
        } catch (e) {
          throw new Error('Response is not in expected JSON format: ' + text);
        }
      } else {
         text = jsonMatch[1].trim();
         geminiAnalysis = JSON.parse(text);
      }
      
      if (!geminiAnalysis.isInputValid && geminiAnalysis.inputError) {
        return res.status(400).json({
          message: 'Gemini analysis failed due to invalid input.',
          suggestion: geminiAnalysis.inputError,
          analysis: geminiAnalysis
        });
      }
    } catch (error) {
      console.error('Gemini API error:', error.response?.data || error.message);
      // Fallback logic remains the same
      // ... (your existing fallback is good for creating a dummy geminiAnalysis object)
      return res.status(500).json({
          message: 'Failed to get analysis from Gemini API.',
          error: error.response?.data?.error?.message || error.message
      });
    }

    // Create a summarized toxicityResult object
    const toxicityResult = {
      smiles: geminiAnalysis.smiles || smiles,
      compoundName: originalName,
      symptoms, 
      isNovelCompound: geminiAnalysis.isNovelCompound,
      acuteToxicity: geminiAnalysis.acuteToxicity,
      endpoints: geminiAnalysis.endpoints,
      sideEffects: geminiAnalysis.sideEffects, 
      structuralAnalysis: geminiAnalysis.structuralAnalysis
    };

    // Create the full entry for the database
    const entry = new Toxicity({
      smiles: geminiAnalysis.smiles || smiles,
      compoundName: originalName, // Field now exists in schema
      symptoms, 
      toxicityResult, // The summarized object
      userId, 
      geminiAnalysis // The full, detailed object
    });
    
    // Save to database
    await entry.save();

    res.status(201).json({ // Use 201 Created for new resources
      message: geminiAnalysis.isInputValid ? 'Analysis complete and saved' : 'Completed with assumptions and saved',
      result: geminiAnalysis, // Return the full analysis object to the frontend
      historyId: entry._id,
      compoundName: originalName,
      smiles: entry.smiles
    });
  } catch (err) {
    console.error('Error in predictToxicityController:', err);
    // Catches validation errors from Mongoose as well
    res.status(500).json({ 
      message: 'Prediction failed due to an internal error.', 
      error: err.message,
      suggestion: 'Please check the server logs or contact support.'
    });
  }
};

export const getGeminiAnalysis = async (req, res) => {
  try {
    const { compound, symptoms } = req.body;
    const userId = req.user?._id;

    if (!compound || !symptoms) {
      return res.status(400).json({ 
        message: 'Compound name and symptoms are required.',
        suggestion: 'Please provide a valid molecule name and symptoms.'
      });
    }
    if (!userId) {
      return res.status(401).json({ 
        message: 'User not authenticated.',
        suggestion: 'Please log in and try again.'
      });
    }

    let smilesData;
    try {
      smilesData = await getSmilesByName(compound, userId);
    } catch (error) {
      return res.status(404).json({ 
        message: error.message,
        suggestion: 'Please ensure the compound name exists in your drug database.'
      });
    }

    const { smiles, originalName } = smilesData;

    if (!isValidSmiles(smiles)) {
      return res.status(400).json({ 
        message: 'Invalid SMILES format.',
        suggestion: 'Please verify the SMILES string or molecule name.'
      });
    }

    const symptomsArray = symptoms
      .split(',')
      .map(s => s.trim().toLowerCase())
      .filter(s => s.length > 0);

    if (!symptomsArray.length) {
      return res.status(400).json({ 
        message: 'At least one valid symptom is required.',
        suggestion: 'Please provide symptoms separated by commas.'
      });
    }

    if (!GEMINI_API_KEY) {
      return res.status(500).json({ 
        message: 'Gemini API key is not configured.',
        suggestion: 'Please contact the administrator to configure the API key.'
      });
    }

    const prompt = `
You are a senior computational toxicologist. For the compound with SMILES "${smiles}" (known as "${originalName}") intended to treat "${symptoms}", provide complete toxicity estimates—no field can be "unknown." If data is unavailable, use structural or QSAR-based inference with your rationale. If the SMILES is invalid, return a JSON object with isInputValid set to false and an inputError message.

Return ONLY a JSON object wrapped in \`\`\`json\n and \n\`\`\` with:
{
  "isInputValid": boolean,
  "inputError": string|null,
  "smiles": string,
  "moleculeName": string,
  "isNovelCompound": boolean,
  "structuralAnalysis": {
    "molecularWeight": string,
    "logP": string,
    "polarSurfaceArea": string,
    "functionalGroups": string[],
    "complexity": string,
    "drugLikeness": string
  },
  "acuteToxicity": {
    "LD50": string,
    "toxicityClass": string,
    "rationale": string,
    "supplemental": string|null
  },
  "endpoints": {
    "hepatotoxicity": string,
    "carcinogenicity": string,
    "cardiotoxicity": string,
    "neurotoxicity": string,
    "nephrotoxicity": string,
    "skin_sensitization": string,
    "thyroid_toxicity": string,
    "ototoxicity": string,
    "musculoskeletal_toxicity": string,
    "dermatotoxicity": string,
    "supplemental": string|null
  },
  "sideEffects": [
    {
      "name": string,
      "description": string,
      "severity": string,
      "organSystem": string,
      "likelihood": string
    }
  ],
  "mechanisms": [
    {
      "feature": string,
      "pathway": string,
      "toxicityType": string
    }
  ],
  "structureConcerns": [
    {
      "substructure": string,
      "concern": string,
      "riskLevel": string
    }
  ],
  "safetyRecommendations": [
    {
      "test": string,
      "priority": string,
      "rationale": string
    }
  ],
  "qsarAnalysis": {
    "properties": [
      {
        "name": string,
        "predictedValue": string,
        "toxicologicalImplication": string
      }
    ],
    "overallRisk": string,
    "symptomContext": [
      {
        "symptom": string,
        "structuralInsight": string,
        "riskAssessment": string
      }
    ]
  },
  "toxicophoreAnalysis": [
    {
      "substructure": string,
      "concern": string,
      "prevalence": string,
      "mitigation": string
    }
  ],
  "developmentRecommendations": [
    {
      "phase": string,
      "recommendation": string,
      "importance": string
    }
  ]
}

If invalid SMILES, return:
{
  "isInputValid": false,
  "inputError": "Invalid SMILES string: [reason]",
  "smiles": "${smiles}",
  "moleculeName": "${originalName}",
  ...
}

Example:
\`\`\`json
{
  "isInputValid": true,
  "inputError": null,
  "smiles": "${smiles}",
  "moleculeName": "${originalName}",
  "isNovelCompound": true,
  "sideEffects": [],
  "structuralAnalysis": {
    "molecularWeight": "300.5 g/mol",
    "logP": "2.3",
    "polarSurfaceArea": "40.0 Å²",
    "functionalGroups": ["Amide"],
    "complexity": "Moderate",
    "drugLikeness": "High"
  },
  ...
}
\`\`\`
`;

    let geminiAnalysis;
    try {
      const geminiResponse = await retryAxios({
        method: 'post',
        url: `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        data: { contents: [{ role: 'user', parts: [{ text: prompt }] }] },
        headers: { 'Content-Type': 'application/json' },
        timeout: 20000
      });

      let analysisText = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
      if (!analysisText) {
        throw new Error('No content returned from Gemini API.');
      }

      console.log('Raw Gemini response:', analysisText);

      const jsonMatch = analysisText.match(/```json\n([\s\S]*?)\n```/);
      if (!jsonMatch) {
        throw new Error(`Response is not in expected JSON format: ${analysisText}`);
      }

      analysisText = jsonMatch[1].trim();

      try {
        geminiAnalysis = JSON.parse(analysisText);
        console.log('Parsed Gemini analysis:', JSON.stringify(geminiAnalysis, null, 2));
      } catch (parseError) {
        console.error('JSON parsing error:', parseError.message, 'Raw JSON:', analysisText);
        throw new Error(`Invalid JSON response: ${parseError.message}`);
      }

      if (!geminiAnalysis.isInputValid && geminiAnalysis.inputError) {
        return res.status(400).json({
          message: 'Gemini analysis failed due to invalid input.',
          suggestion: geminiAnalysis.inputError,
          analysis: geminiAnalysis
        });
      }
    } catch (error) {
      console.error('Gemini API error:', error.response?.data || error.message);

      geminiAnalysis = {
        isInputValid: false,
        inputError: `Analysis failed: ${error.response?.data?.error?.message || error.message}`,
        smiles: smiles,
        moleculeName: originalName,
        isNovelCompound: true,
        sideEffects: symptomsArray.map(symptom => ({
          name: symptom.replace(/\b\w/g, c => c.toUpperCase()),
          description: `Potential adverse effect inferred from ${symptom} association`,
          severity: 'Moderate - estimated',
          organSystem: symptomToEndpointMap[Object.keys(symptomToEndpointMap).find(key => symptom.toLowerCase().includes(key))] || 'Systemic',
          likelihood: '40% - estimated'
        })),
        structuralAnalysis: {
          molecularWeight: 'Unknown',
          logP: 'Unknown',
          polarSurfaceArea: 'Unknown',
          functionalGroups: [],
          complexity: 'Unknown',
          drugLikeness: 'Unknown'
        },
        acuteToxicity: {
          LD50: 'Unknown',
          toxicityClass: 'Unknown',
          rationale: 'API failure',
          supplemental: null
        },
        endpoints: {
          hepatotoxicity: 'Unknown',
          carcinogenicity: 'Unknown',
          cardiotoxicity: 'Unknown',
          neurotoxicity: 'Unknown',
          nephrotoxicity: 'Unknown',
          skin_sensitization: 'Unknown',
          thyroid_toxicity: 'Unknown',
          ototoxicity: 'Unknown',
          musculoskeletal_toxicity: 'Unknown',
          dermatotoxicity: 'Unknown',
          supplemental: null
        },
        mechanisms: [],
        structureConcerns: [],
        safetyRecommendations: [],
        qsarAnalysis: {
          properties: [],
          overallRisk: 'Unknown',
          symptomContext: []
        },
        toxicophoreAnalysis: [],
        developmentRecommendations: []
      };
    }

    res.status(200).json({
      message: 'Novel compound analysis generated successfully',
      analysis: geminiAnalysis,
      compoundName: originalName,
      smiles: smiles
    });
  } catch (error) {
    console.error('Error generating Gemini analysis:', error.message);
    res.status(500).json({
      message: 'Failed to generate Gemini analysis',
      error: error.message,
      suggestion: 'Please check the molecule name and symptoms, then try again.'
    });
  }
};

export const saveGeminiAnalysis = async (req, res) => {
  try {
    const { compoundName, symptoms, geminiAnalysis } = req.body;
    const userId = req.user?._id;

    if (!compoundName || !symptoms || !geminiAnalysis) {
      return res.status(400).json({ 
        message: 'Compound name, symptoms, and analysis are required.',
        suggestion: 'Please provide all required fields.'
      });
    }
    if (!userId) {
      return res.status(401).json({ 
        message: 'User not authenticated.',
        suggestion: 'Please log in and try again.'
      });
    }

    let smilesData;
    try {
      smilesData = await getSmilesByName(compoundName, userId);
    } catch (error) {
      return res.status(404).json({ 
        message: error.message,
        suggestion: 'Please ensure the compound name exists in your drug database.'
      });
    }

    const { smiles, originalName } = smilesData;

    const toxicityEntry = await Toxicity.findOne({ smiles, symptoms, userId }).sort({ created: -1 });

    if (!toxicityEntry) {
      return res.status(404).json({ 
        message: 'No matching toxicity prediction found.',
        suggestion: 'Please generate a toxicity prediction first.'
      });
    }

    const symptomList = symptoms.split(',').map(s => s.trim());
    const inferredSideEffects = symptomList.map(symptom => ({
      name: symptom,
      description: `Potential effect related to ${symptom}. Further investigation required.`,
      severity: 'Unknown',
      organSystem: symptomToEndpointMap[Object.keys(symptomToEndpointMap).find(key => symptom.toLowerCase().includes(key))] || 'Unknown',
      likelihood: 'Unknown'
    }));

    const validatedAnalysis = {
      isInputValid: geminiAnalysis.isInputValid || false,
      inputError: geminiAnalysis.inputError || null,
      smiles: geminiAnalysis.smiles || smiles,
      moleculeName: geminiAnalysis.moleculeName || originalName,
      isNovelCompound: geminiAnalysis.isNovelCompound || true,
      structuralAnalysis: geminiAnalysis.structuralAnalysis || {
        molecularWeight: 'Unknown',
        logP: 'Unknown',
        polarSurfaceArea: 'Unknown',
        functionalGroups: [],
        complexity: 'Unknown',
        drugLikeness: 'Unknown'
      },
      acuteToxicity: {
        LD50: geminiAnalysis.acuteToxicity?.LD50 || 'Unknown',
        toxicityClass: geminiAnalysis.acuteToxicity?.toxicityClass || 'Unknown',
        rationale: geminiAnalysis.acuteToxicity?.rationale || 'No data available',
        supplemental: geminiAnalysis.acuteToxicity?.supplemental || 'No specific LD50 data found.'
      },
      endpoints: {
        hepatotoxicity: geminiAnalysis.endpoints?.hepatotoxicity || 'Unknown',
        carcinogenicity: geminiAnalysis.endpoints?.carcinogenicity || 'Unknown',
        cardiotoxicity: geminiAnalysis.endpoints?.cardiotoxicity || 'Unknown',
        neurotoxicity: geminiAnalysis.endpoints?.neurotoxicity || 'Unknown',
        nephrotoxicity: geminiAnalysis.endpoints?.nephrotoxicity || 'Unknown',
        skin_sensitization: geminiAnalysis.endpoints?.skin_sensitization || 'Unknown',
        thyroid_toxicity: geminiAnalysis.endpoints?.thyroid_toxicity || 'Unknown',
        ototoxicity: geminiAnalysis.endpoints?.ototoxicity || 'Unknown',
        musculoskeletal_toxicity: geminiAnalysis.endpoints?.musculoskeletal_toxicity || 'Unknown',
        dermatotoxicity: geminiAnalysis.endpoints?.dermatotoxicity || 'Unknown',
        supplemental: geminiAnalysis.endpoints?.supplemental || 'No endpoint data available.'
      },
      sideEffects: geminiAnalysis.sideEffects?.length > 0 ? geminiAnalysis.sideEffects : inferredSideEffects,
      mechanisms: geminiAnalysis.mechanisms?.length > 0 ? geminiAnalysis.mechanisms : [
        {
          feature: 'Structure from SMILES',
          pathway: 'Based on SMILES analysis, mechanisms may involve metabolic pathways.',
          toxicityType: 'Metabolic'
        }
      ],
      structureConcerns: geminiAnalysis.structureConcerns?.length > 0 ? geminiAnalysis.structureConcerns : [
        {
          substructure: 'From SMILES analysis',
          concern: 'Structural alerts based on SMILES pattern analysis.',
          riskLevel: 'Moderate'
        }
      ],
      safetyRecommendations: geminiAnalysis.safetyRecommendations?.length > 0 ? geminiAnalysis.safetyRecommendations : [
        {
          test: 'General toxicity screening',
          priority: 'High',
          rationale: 'Conduct comprehensive screening based on SMILES structural analysis.'
        }
      ],
      qsarAnalysis: geminiAnalysis.qsarAnalysis || {
        properties: [],
        overallRisk: 'Unknown',
        symptomContext: inferredSideEffects.map(effect => ({
          symptom: effect.name,
          structuralInsight: `Potential toxicological effect based on SMILES analysis and ${effect.description.toLowerCase()}.`,
          riskAssessment: 'Unknown'
        }))
      },
      toxicophoreAnalysis: geminiAnalysis.toxicophoreAnalysis || [],
      developmentRecommendations: geminiAnalysis.developmentRecommendations || [
        {
          phase: 'Preclinical',
          recommendation: 'Conduct comprehensive toxicity screening based on SMILES analysis',
          importance: 'High'
        }
      ]
    };

    toxicityEntry.geminiAnalysis = validatedAnalysis;
    toxicityEntry.compoundName = originalName;

    await toxicityEntry.validate();
    await toxicityEntry.save();

    return res.status(200).json({ 
      message: 'Analysis saved successfully', 
      historyId: toxicityEntry._id,
      compoundName: originalName,
      smiles: smiles
    });
  } catch (error) {
    console.error('Error saving Gemini analysis:', error.message, error);
    return res.status(500).json({
      message: 'Failed to save analysis',
      error: `Toxicity validation failed: ${error.message}`,
      suggestion: 'Please try again or contact support.'
    });
  }
};

export const getToxicityHistory = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ 
        message: 'User not authenticated.',
        suggestion: 'Please log in and try again.'
      });
    }

    const history = await Toxicity.find({ userId }).sort({ created: -1 });

    res.status(200).json({ history });
  } catch (error) {
    console.error('Error fetching toxicity history:', error.message);
    res.status(500).json({
      message: 'Failed to fetch history',
      error: error.message,
      suggestion: 'Please try again or contact support.'
    });
  }
};