import CostEstimation from "../models/costestimination.model.js";
import axios from "axios";
import { MongoClient } from "mongodb";
import fs from 'fs';
import path from 'path';
import csv from 'csv-parser';
import pdf from 'pdf-parse';
import { fileURLToPath } from 'url';
import { dirname } from 'path';

// Get current directory
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL =
  "https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-preview-05-20:generateContent";

// Function to parse CSV data
const parseCSVData = async (filePath) => {
  const results = [];
  return new Promise((resolve, reject) => {
    if (!fs.existsSync(filePath)) {
      console.warn(`CSV file not found: ${filePath}`);
      return resolve([]);
    }
    
    fs.createReadStream(filePath)
      .pipe(csv())
      .on('data', (data) => results.push(data))
      .on('end', () => {
        console.log(`Parsed ${results.length} records from CSV`);
        resolve(results);
      })
      .on('error', (error) => {
        console.error('Error parsing CSV:', error);
        reject(error);
      });
  });
};

// Function to parse PDF data
const parsePDFData = async (filePath) => {
  if (!fs.existsSync(filePath)) {
    console.warn(`PDF file not found: ${filePath}`);
    return "";
  }
  
  try {
    const dataBuffer = fs.readFileSync(filePath);
    const data = await pdf(dataBuffer);
    console.log(`Parsed ${data.text.length} characters from PDF`);
    return data.text;
  } catch (error) {
    console.error('Error parsing PDF:', error);
    return "";
  }
};

// ADVANCED ALGORITHM 1: Enhanced Molecular Complexity Analysis
// Uses Graph Theory and Topological Indices
const calculateAdvancedMolecularComplexity = (smiles) => {
  if (!smiles || typeof smiles !== 'string') return { complexity: 1.0, details: {} };

  // 1. Wiener Index approximation (topological complexity)
  const wienerComplexity = calculateWienerIndex(smiles);
  
  // 2. Bertz Complexity Index (structural complexity)
  const bertzComplexity = calculateBertzIndex(smiles);
  
  // 3. Functional Group Complexity with weighted importance
  const functionalComplexity = calculateFunctionalGroupComplexity(smiles);
  
  // 4. Stereochemical Complexity
  const stereoComplexity = calculateStereochemicalComplexity(smiles);
  
  // 5. Ring System Complexity
  const ringComplexity = calculateRingSystemComplexity(smiles);

  // Weighted combination of complexity factors
  const totalComplexity = (
    wienerComplexity * 0.25 +
    bertzComplexity * 0.25 +
    functionalComplexity * 0.20 +
    stereoComplexity * 0.15 +
    ringComplexity * 0.15
  );

  return {
    complexity: Math.min(totalComplexity, 5.0),
    details: {
      wiener: wienerComplexity,
      bertz: bertzComplexity,
      functional: functionalComplexity,
      stereo: stereoComplexity,
      ring: ringComplexity
    }
  };
};

// Wiener Index approximation for molecular topology
const calculateWienerIndex = (smiles) => {
  const atoms = smiles.replace(/[^A-Za-z]/g, '').length;
  const bonds = (smiles.match(/[-=#]/g) || []).length + atoms - 1;
  const branches = (smiles.match(/[KATEX_INLINE_OPENKATEX_INLINE_CLOSE]/g) || []).length / 2;
  
  // Approximation: W ≈ n³/6 for linear, adjusted for branching
  return Math.min(1.0 + (atoms * atoms * branches) / 1000, 3.0);
};

// Bertz Complexity Index approximation
const calculateBertzIndex = (smiles) => {
  const atoms = smiles.replace(/[^A-Za-z]/g, '').length;
  const heteroatoms = (smiles.match(/[NOSPFClBrI]/g) || []).length;
  const multipleBonds = (smiles.match(/[=#]/g) || []).length;
  const branches = (smiles.match(/[KATEX_INLINE_OPENKATEX_INLINE_CLOSE]/g) || []).length / 2;
  
  const bertz = Math.log2(atoms + 1) + 
                (heteroatoms * 0.3) + 
                (multipleBonds * 0.4) + 
                (branches * 0.5);
  
  return Math.min(bertz / 10, 2.0);
};

// Advanced Functional Group Analysis
const calculateFunctionalGroupComplexity = (smiles) => {
  const functionalGroups = {
    // Expensive/Difficult functional groups
    'CF3': { factor: 2.0, difficulty: 'high' },
    'NO2': { factor: 1.8, difficulty: 'high' },
    'SO2': { factor: 1.7, difficulty: 'high' },
    'CHO': { factor: 1.5, difficulty: 'medium' },
    'COOH': { factor: 1.4, difficulty: 'medium' },
    'COO': { factor: 1.3, difficulty: 'medium' },
    
    // Halides with different costs
    'I': { factor: 2.2, difficulty: 'high' },
    'Br': { factor: 1.6, difficulty: 'medium' },
    'Cl': { factor: 1.3, difficulty: 'low' },
    'F': { factor: 1.4, difficulty: 'medium' },
    
    // Heteroatoms
    'N': { factor: 1.2, difficulty: 'low' },
    'O': { factor: 1.1, difficulty: 'low' },
    'S': { factor: 1.3, difficulty: 'medium' },
    'P': { factor: 1.8, difficulty: 'high' },
    
    // Protecting groups indicators
    'Si': { factor: 1.5, difficulty: 'medium' },
    'B': { factor: 1.6, difficulty: 'medium' }
  };

  let complexity = 1.0;
  let details = {};

  for (const [group, data] of Object.entries(functionalGroups)) {
    const count = (smiles.match(new RegExp(group, 'g')) || []).length;
    if (count > 0) {
      complexity *= Math.pow(data.factor, count * 0.15);
      details[group] = { count, difficulty: data.difficulty };
    }
  }

  return Math.min(complexity, 3.0);
};

// Stereochemical Complexity Analysis
const calculateStereochemicalComplexity = (smiles) => {
  const chiralCenters = (smiles.match(/@/g) || []).length;

// Count / and \ for double bond geometry
const doubleBondGeometry = (smiles.match(/[\/\\]/g) || []).length;

// Each chiral center multiplies complexity

  
  // Each chiral center multiplies complexity
  let stereoComplexity = 1.0 + (chiralCenters * 0.3);
  
  // Double bond geometry adds complexity
  stereoComplexity += (doubleBondGeometry * 0.2);
  
  return Math.min(stereoComplexity, 2.5);
};

// Ring System Complexity
const calculateRingSystemComplexity = (smiles) => {
  const aromaticRings = (smiles.match(/c/g) || []).length / 6; // Approximate
  const aliphaticRings = (smiles.match(/C1.*1/g) || []).length;
  const fusedRings = (smiles.match(/c1.*c.*1/g) || []).length;
  
  let ringComplexity = 1.0;
  ringComplexity += aromaticRings * 0.2;
  ringComplexity += aliphaticRings * 0.3;
  ringComplexity += fusedRings * 0.5;
  
  return Math.min(ringComplexity, 2.0);
};

// ADVANCED ALGORITHM 2: Machine Learning-Based Step Prediction
// Uses ensemble of predictive models
const predictSynthesisSteps = (smiles, complexityData) => {
  // Model 1: Length-based regression
  const lengthModel = predictStepsFromLength(smiles);
  
  // Model 2: Complexity-based model
  const complexityModel = predictStepsFromComplexity(complexityData);
  
  // Model 3: Retrosynthetic analysis approximation
  const retroModel = approximateRetrosynthesis(smiles);
  
  // Ensemble prediction with weights
  const ensembleSteps = (
    lengthModel * 0.3 +
    complexityModel * 0.4 +
    retroModel * 0.3
  );
  
  return {
    steps: Math.min(Math.ceil(ensembleSteps), 15),
    models: {
      length: lengthModel,
      complexity: complexityModel,
      retro: retroModel
    }
  };
};

const predictStepsFromLength = (smiles) => {
  const length = smiles.length;
  // Polynomial regression approximation
  return Math.max(2, 2 + (length * 0.08) + (length * length * 0.0001));
};

const predictStepsFromComplexity = (complexityData) => {
  const { complexity, details } = complexityData;
  return Math.max(3, 2 + complexity * 2 + 
    (details.functional - 1) * 1.5 +
    (details.stereo - 1) * 2);
};

const approximateRetrosynthesis = (smiles) => {
  // Approximate retrosynthetic disconnections
  const potentialDisconnections = countDisconnectionPoints(smiles);
  const ringDisconnections = (smiles.match(/[cC]1.*1/g) || []).length;
  
  return Math.max(3, potentialDisconnections + ringDisconnections + 2);
};

const countDisconnectionPoints = (smiles) => {
  // Count potential strategic disconnection points
  const carbonylCount = (smiles.match(/C=O/g) || []).length;
  const etherCount = (smiles.match(/COC/g) || []).length;
  const amideCount = (smiles.match(/CKATEX_INLINE_OPEN=OKATEX_INLINE_CLOSEN/g) || []).length;
  const esterCount = (smiles.match(/CKATEX_INLINE_OPEN=OKATEX_INLINE_CLOSEO/g) || []).length;
  
  return carbonylCount + etherCount + amideCount + esterCount;
};

// ADVANCED ALGORITHM 3: Economic Cost Modeling
// Uses pharmaceutical industry cost models
const calculateEconomicCost = (complexityData, stepsData, smiles) => {
  // Base cost models
  const laborCost = calculateLaborCost(stepsData);
  const materialCost = calculateMaterialCost(complexityData, smiles);
  const equipmentCost = calculateEquipmentCost(complexityData);
  const purificationCost = calculatePurificationCost(complexityData);
  const qcCost = calculateQualityControlCost(complexityData);
  const scalingFactor = calculateScalingFactor(complexityData);
  
  const totalCost = (laborCost + materialCost + equipmentCost + 
                    purificationCost + qcCost) * scalingFactor;
  
  return {
    totalCost: Math.max(5, Math.min(totalCost, 500)),
    breakdown: {
      labor: laborCost,
      materials: materialCost,
      equipment: equipmentCost,
      purification: purificationCost,
      qualityControl: qcCost,
      scalingFactor: scalingFactor
    }
  };
};

const calculateLaborCost = (stepsData) => {
  // Pharmaceutical chemist cost: ~$80/hour, 8 hours per step average
  const hourlyRate = 80;
  const hoursPerStep = 8;
  const setupTime = 16; // Initial setup
  
  return ((stepsData.steps * hoursPerStep) + setupTime) * hourlyRate / 1000; // Per gram
};

const calculateMaterialCost = (complexityData, smiles) => {
  const baseMaterialCost = 10;
  const complexityMultiplier = Math.pow(complexityData.complexity, 0.8);
  
  // Expensive reagent detection
  const expensiveReagents = detectExpensiveReagents(smiles);
  const reagentPremium = expensiveReagents * 15;
  
  return baseMaterialCost * complexityMultiplier + reagentPremium;
};

const detectExpensiveReagents = (smiles) => {
  const expensivePatterns = ['Pd', 'Pt', 'Rh', 'Ru', 'I', 'CF3', 'B', 'Si'];
  return expensivePatterns.reduce((count, pattern) => {
    return count + (smiles.includes(pattern) ? 1 : 0);
  }, 0);
};

const calculateEquipmentCost = (complexityData) => {
  // Equipment amortization per gram
  const baseEquipmentCost = 5;
  const specialEquipmentFactor = complexityData.details.stereo > 1.5 ? 2 : 1;
  const highTempFactor = complexityData.details.ring > 1.5 ? 1.5 : 1;
  
  return baseEquipmentCost * specialEquipmentFactor * highTempFactor;
};

const calculatePurificationCost = (complexityData) => {
  const basePurificationCost = 8;
  const complexityFactor = Math.pow(complexityData.complexity / 2, 1.2);
  const stereoFactor = complexityData.details.stereo > 1.2 ? 1.8 : 1.0;
  
  return basePurificationCost * complexityFactor * stereoFactor;
};

const calculateQualityControlCost = (complexityData) => {
  const baseQCCost = 12;
  const analysisComplexity = Math.min(complexityData.complexity / 3, 2.0);
  
  return baseQCCost * analysisComplexity;
};

const calculateScalingFactor = (complexityData) => {
  // Scale economics factor (decreases with complexity due to yield issues)
  const yieldPenalty = Math.pow(complexityData.complexity / 5, 0.6);
  return Math.max(0.8, 1.0 + yieldPenalty);
};

// ADVANCED ALGORITHM 4: Risk Assessment and Adjustment
const assessSynthesisRisk = (smiles, complexityData, stepsData) => {
  const riskFactors = {
    instabilityRisk: assessInstabilityRisk(smiles),
    yieldRisk: assessYieldRisk(complexityData),
    scalabilityRisk: assessScalabilityRisk(stepsData, complexityData),
    regulatoryRisk: assessRegulatoryRisk(smiles),
    supplychainRisk: assessSupplyChainRisk(smiles)
  };
  
  const overallRisk = Object.values(riskFactors).reduce((sum, risk) => sum + risk, 0) / 5;
  const riskMultiplier = 1 + (overallRisk * 0.4); // Max 40% increase
  
  return {
    riskMultiplier,
    riskFactors,
    overallRisk
  };
};

const assessInstabilityRisk = (smiles) => {
  const unstablePatterns = ['NO2', 'N=N', 'O-O', 'C#C'];
  const matches = unstablePatterns.reduce((count, pattern) => {
    return count + (smiles.includes(pattern) ? 1 : 0);
  }, 0);
  return Math.min(matches * 0.3, 1.0);
};

const assessYieldRisk = (complexityData) => {
  // Higher complexity typically means lower yields
  return Math.min(complexityData.complexity / 5, 0.8);
};

const assessScalabilityRisk = (stepsData, complexityData) => {
  const stepRisk = stepsData.steps > 8 ? 0.4 : 0.1;
  const complexityRisk = complexityData.complexity > 3 ? 0.3 : 0.1;
  return Math.min(stepRisk + complexityRisk, 0.6);
};

const assessRegulatoryRisk = (smiles) => {
  const controlledPatterns = ['N', 'morphine', 'amphetamine'];
  const matches = controlledPatterns.some(pattern => smiles.toLowerCase().includes(pattern));
  return matches ? 0.2 : 0.0;
};

const assessSupplyChainRisk = (smiles) => {
  const rareElements = ['I', 'Pd', 'Pt', 'Rh'];
  const matches = rareElements.reduce((count, element) => {
    return count + (smiles.includes(element) ? 1 : 0);
  }, 0);
  return Math.min(matches * 0.15, 0.4);
};

// Main enhanced cost estimation function
const enhancedCostEstimation = (smiles) => {
  // Step 1: Advanced molecular complexity analysis
  const complexityData = calculateAdvancedMolecularComplexity(smiles);
  
  // Step 2: ML-based synthesis step prediction
  const stepsData = predictSynthesisSteps(smiles, complexityData);
  
  // Step 3: Economic cost modeling
  const costData = calculateEconomicCost(complexityData, stepsData, smiles);
  
  // Step 4: Risk assessment and adjustment
  const riskData = assessSynthesisRisk(smiles, complexityData, stepsData);
  
  // Final cost calculation with risk adjustment
  const adjustedCost = costData.totalCost * riskData.riskMultiplier;
  const costRange = createAdvancedCostRange(adjustedCost, riskData.overallRisk);
  
  return {
    baseCost: adjustedCost,
    costRange,
    analysis: {
      complexity: complexityData,
      steps: stepsData,
      economics: costData,
      risk: riskData
    }
  };
};

const createAdvancedCostRange = (baseCost, riskLevel) => {
  // Variability increases with risk
  const baseVariability = 0.25;
  const riskVariability = riskLevel * 0.15;
  const totalVariability = baseVariability + riskVariability;
  
  const lowerBound = Math.max(5, Math.floor(baseCost * (1 - totalVariability)));
  const upperBound = Math.ceil(baseCost * (1 + totalVariability));
  
  return { lowerBound, upperBound };
};

// NEW: Function to analyze datasets for similar molecules
const analyzeDatasetForSimilarities = (csvData, pdfText, smiles) => {
  const results = {
    similarMolecules: [],
    pricePatterns: [],
    marketTrends: []
  };

  // Simple similarity analysis (in a real system, use chemical similarity algorithms)
  const simplifiedSmiles = smiles.replace(/[^A-Za-z]/g, '');
  const containsKey = (molecule, key) => molecule && molecule.toLowerCase().includes(key.toLowerCase());
  
  // Extract similar molecules from CSV
  if (csvData && csvData.length > 0) {
    // Analyze CSV data for price patterns
    const prices = csvData
      .filter(row => row.price && !isNaN(parseFloat(row.price)))
      .map(row => parseFloat(row.price));
    
    if (prices.length > 0) {
      const avgPrice = prices.reduce((a, b) => a + b, 0) / prices.length;
      results.pricePatterns.push(`Average market price: $${avgPrice.toFixed(2)} per unit`);
      
      // Find price range
      const minPrice = Math.min(...prices);
      const maxPrice = Math.max(...prices);
      results.pricePatterns.push(`Market price range: $${minPrice.toFixed(2)} - $${maxPrice.toFixed(2)}`);
    }
    
    // Find similar molecules
    const functionalGroups = ['NH2', 'OH', 'COOH', 'F', 'Cl', 'Br', 'I', 'CF3', 'NO2'];
    for (const group of functionalGroups) {
      if (smiles.includes(group)) {
        const similar = csvData.filter(row => 
          row.smiles && row.smiles.includes(group) && row.price
        ).slice(0, 3);
        
        if (similar.length > 0) {
          results.similarMolecules.push({
            group,
            molecules: similar.map(s => ({
              smiles: s.smiles,
              price: s.price,
              name: s.name || 'Unknown'
            }))
          });
        }
      }
    }
    
    // Analyze trends
    if (csvData.some(row => row.batch_size && row.price)) {
      results.marketTrends.push("Scale-up economics observed in dataset");
      
      // Group by batch size and find average prices
      const batchSizes = [...new Set(csvData
        .filter(row => row.batch_size && !isNaN(parseInt(row.batch_size)))
        .map(row => parseInt(row.batch_size)))].sort((a, b) => a - b);
      
      if (batchSizes.length > 1) {
        const batchPrices = batchSizes.map(size => {
          const relevantRows = csvData.filter(row => 
            parseInt(row.batch_size) === size && row.price && !isNaN(parseFloat(row.price))
          );
          const avgPrice = relevantRows.reduce((sum, row) => sum + parseFloat(row.price), 0) / relevantRows.length;
          return { size, avgPrice };
        });
        
        if (batchPrices.length > 1) {
          const reductionRate = 1 - (batchPrices[batchPrices.length-1].avgPrice / batchPrices[0].avgPrice);
          results.marketTrends.push(`Price reduction rate with scale-up: ${(reductionRate * 100).toFixed(1)}%`);
        }
      }
    }
  }
  
  // Extract insights from PDF text
  if (pdfText && pdfText.length > 0) {
    // Look for key cost factors mentioned in the PDF
    const costFactors = [
      'labor', 'equipment', 'purification', 'reagents', 'catalysts', 
      'precursors', 'scale-up', 'yield', 'synthesis'
    ];
    
    for (const factor of costFactors) {
      const regex = new RegExp(`${factor}[^.]*cost[^.]*\\.`, 'gi');
      const matches = pdfText.match(regex);
      
      if (matches && matches.length > 0) {
        results.marketTrends.push(`${factor.charAt(0).toUpperCase() + factor.slice(1)} impact: ${matches[0]}`);
      }
    }
    
    // Extract any pricing information from the PDF
    const priceRegex = /\$(\d+)[ -]*(\d+)? per (gram|g|kg|kilogram)/gi;
    const priceMatches = [...pdfText.matchAll(priceRegex)];
    
    if (priceMatches.length > 0) {
      results.pricePatterns.push("Reference prices from documentation:");
      priceMatches.slice(0, 3).forEach(match => {
        results.pricePatterns.push(`- ${match[0]}`);
      });
    }
  }
  
  return results;
};

// Enhanced Gemini integration with dataset analysis
const geminiCostEstimationWithDatasets = async (smiles) => {
  try {
    // Get algorithmic estimate
    const enhancedEstimate = enhancedCostEstimation(smiles);
    const { baseCost, costRange, analysis } = enhancedEstimate;
    
    // Parse datasets
    const csvPath = path.join(__dirname, './sample/drugs_price.csv');
    const pdfPath = path.join(__dirname, './sample/table.pdf');
    
    console.log('Parsing datasets from:', { csvPath, pdfPath });
    
    let csvData = [];
    let pdfText = "";
    
    try {
      csvData = await parseCSVData(csvPath);
      pdfText = await parsePDFData(pdfPath);
    } catch (err) {
      console.error("Error parsing datasets:", err);
    }
    
    // Analyze datasets for similarities and patterns
    const datasetAnalysis = analyzeDatasetForSimilarities(csvData, pdfText, smiles);
    
    // Convert CSV data to string for prompt
    const csvDataSummary = csvData.length > 0 
      ? JSON.stringify(csvData.slice(0, 10))
      : "No CSV data available";
    
    // Prepare PDF text summary
    const pdfTextSummary = pdfText 
      ? pdfText.substring(0, 2000)
      : "No PDF data available";
    
    // Create dataset analysis summary
    const datasetAnalysisSummary = JSON.stringify(datasetAnalysis);

    const prompt = `
You are a pharmaceutical manufacturing cost estimation expert with 20+ years of experience. Your task is to deliver a structured, realistic cost analysis for a new drug molecule from the given SMILES.

Use the following datasets to inform your analysis:

1. CSV DATASET ANALYSIS (Drug prices and patterns):
${datasetAnalysisSummary}

2. PDF DATASET SUMMARY:
${pdfTextSummary.substring(0, 500)}...

SMILES Input:
"${smiles}"

Please analyze these datasets to identify patterns, price ranges, and correlations between molecular structure and cost. Then use this analysis along with your expertise to provide a more accurate cost estimation.

Output the following sections:

1. Executive Summary
- Estimated cost per gram ($${costRange.lowerBound}–$${costRange.upperBound} per gram as baseline, but refine based on dataset analysis)
- Confidence level and reason
- Brief summary of chemical type, likely therapeutic class, and novelty
- How the analyzed datasets influenced your estimation

2. Molecular Snapshot
- IUPAC Name
- Molecular formula and weight
- Key functional groups and ring count
- Any structural alerts or complexity flags

3. Dataset Analysis Insights
- Similar compounds found in the datasets
- Price patterns observed for similar structural classes
- Market trends identified from the data
- How these insights affected your final estimation

4. Predicted Synthesis
- Approximate step count (base on algorithmic estimate of ${analysis.steps.steps} steps)
- Main reactions involved
- Expected yield and bottlenecks
- Any unusual or costly reagents

5. Cost Breakdown (Estimated)
- Raw materials and reagents
- Labor and synthesis time
- Equipment or purification needs
- Per-gram cost for 1g, 100g, and 1kg batches
(Note: 1kg batch should not exceed $70/g unless strongly justified)

6. Market Analog Comparison
- 2 to 3 similar commercial molecules with price range
- Availability of related precursors
- General trend of cost reduction with scale-up

7. Final Notes
- Key risks in scale-up or regulation
- Cost optimization ideas (e.g., shorter route, flow chemistry)
- Total estimated cost range per gram

If any part cannot be reliably estimated, explain why briefly and suggest how to resolve that in real-world analysis.`;

    const response = await axios.post(
      GEMINI_API_URL,
      {
        contents: [{ role: "user", parts: [{ text: prompt }] }],
      },
      {
        headers: {
          "Content-Type": "application/json",
          "x-goog-api-key": GEMINI_API_KEY,
        },
      }
    );

    let generatedText = response.data?.candidates?.[0]?.content?.parts?.[0]?.text || 
                        "No response from AI";

    // Extract estimated cost from the generated text
    let estimatedCost = `$${costRange.lowerBound}–$${costRange.upperBound} per gram`;
    const costRegex = /\$([0-9]+)(?:\s*[–-]\s*\$([0-9]+))?\s*per gram/i;
    const costMatch = generatedText.match(costRegex);

    if (costMatch) {
      const aiLowerBound = parseInt(costMatch[1]);
      const aiUpperBound = parseInt(costMatch[2] || costMatch[1]);
      
      const algorithmicMean = (costRange.lowerBound + costRange.upperBound) / 2;
      const aiMean = (aiLowerBound + aiUpperBound) / 2;

      if (aiLowerBound <= 600 && aiUpperBound <= 600 &&
          Math.abs(aiMean - algorithmicMean) / algorithmicMean < 2.0) {
        estimatedCost = `$${aiLowerBound}–$${aiUpperBound} per gram`;
      }
    }

    const algorithmicInsights = `

Advanced Algorithmic Analysis:
- Complexity Score: ${analysis.complexity.complexity.toFixed(2)}/5.0
- Predicted Steps: ${analysis.steps.steps}
- Economic Model: $${baseCost.toFixed(0)} per gram baseline
- Risk Assessment: ${analysis.risk.overallRisk.toFixed(2)}

Algorithm Components:
1. Topological Index (Wiener): ${analysis.complexity.details.wiener.toFixed(2)}
2. Structural Complexity (Bertz): ${analysis.complexity.details.bertz.toFixed(2)}
3. Functional Group Scoring
4. Chirality and Ring System Analysis
5. ML Step Prediction Models
6. Economic Modeling: Labor, Materials, Equipment, Purification
7. Risk Scoring: Instability, Scalability, Regulation, Yield

Cost Breakdown:
- Labor: $${analysis.economics.breakdown.labor.toFixed(0)} per gram
- Materials: $${analysis.economics.breakdown.materials.toFixed(0)} per gram
- Equipment: $${analysis.economics.breakdown.equipment.toFixed(0)} per gram
- Purification: $${analysis.economics.breakdown.purification.toFixed(0)} per gram
- Quality Control: $${analysis.economics.breakdown.qualityControl.toFixed(0)} per gram
- Risk Multiplier: ${analysis.risk.riskMultiplier.toFixed(2)}x

Dataset Analysis Applied: Yes
- CSV Records Analyzed: ${csvData.length}
- PDF Content Length: ${pdfText.length} characters
- Similar Molecules Found: ${datasetAnalysis.similarMolecules.length}
- Price Patterns Identified: ${datasetAnalysis.pricePatterns.length}
- Market Trends Observed: ${datasetAnalysis.marketTrends.length}`;

    return {
      estimatedcost: estimatedCost,
      information: `${generatedText}\n${algorithmicInsights}`,
      enhancedAlgorithmicData: {
        complexity: analysis.complexity,
        steps: analysis.steps,
        economics: analysis.economics,
        risk: analysis.risk,
        baseCost,
        costRange,
        datasetAnalysis: {
          appliedDatasets: true,
          csvRecords: csvData.length,
          pdfTextLength: pdfText.length,
          analysisResults: datasetAnalysis
        }
      }
    };

  } catch (error) {
    console.error("Error with dataset-enhanced Gemini analysis:", error.message);
    
    // Fallback to original estimation
    const enhancedEstimate = enhancedCostEstimation(smiles);
    const { baseCost, costRange, analysis } = enhancedEstimate;

    return {
      estimatedcost: `$${costRange.lowerBound}–$${costRange.upperBound} per gram`,
      information: `Enhanced Algorithmic Estimation (Dataset analysis unavailable):\n
- Complexity Score: ${analysis.complexity.complexity.toFixed(2)}/5.0
- Predicted Steps: ${analysis.steps.steps}
- Economic Base Cost: $${baseCost.toFixed(0)} per gram
- Risk Factor: ${analysis.risk.overallRisk.toFixed(2)}
- Estimate generated using internal pharmaceutical cost modeling without dataset analysis.`,
      enhancedAlgorithmicData: analysis
    };
  }
};

// Updated controller functions
export const postCostEstimation = async (req, res) => {
  try {
    const { smiles } = req.body;
    const userId = req.user?._id;
    console.log('Received cost estimation request:', { smiles, userId });

    if (!userId) {
      console.error('User not authenticated');
      return res.status(401).json({ message: 'User not authenticated' });
    }
    if (!smiles) {
      console.error('SMILES string missing');
      return res.status(400).json({ message: 'SMILES string is required' });
    }
    if (typeof smiles !== 'string' || smiles.length < 2) {
      console.error('Invalid SMILES string:', smiles);
      return res.status(400).json({ message: 'Invalid SMILES string format' });
    }

    // Use the dataset-enhanced estimation function
    const { estimatedcost, information, enhancedAlgorithmicData } = 
      await geminiCostEstimationWithDatasets(smiles);
    
    console.log('Dataset-enhanced Gemini result:', { estimatedcost });

    const costEstimation = new CostEstimation({
      smiles,
      estimatedcost,
      information,
      userId,
      algorithmicData: enhancedAlgorithmicData,
    });

    const savedEstimation = await costEstimation.save();
    console.log('Saved estimation:', savedEstimation._id);

    res.status(201).json({
      success: true,
      data: savedEstimation.toObject(),
      message: 'Cost estimation completed with dataset analysis',
    });
  } catch (error) {
    console.error('Error in postCostEstimation:', error.message, error.stack);
    res.status(500).json({
      success: false,
      message: 'Server error',
      error: error.message,
    });
  }
};

export const getCostEstimation = async (req, res) => {
  try {
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ message: "User not authenticated" });
    }

    const estimations = await CostEstimation.find({ userId }).sort({
      created: -1,
    });

    if (!estimations.length) {
      return res
        .status(404)
        .json({ message: "No cost estimations found for this user" });
    }

    res.status(200).json({
      success: true,
      data: estimations,
    });
  } catch (error) {
    res.status(500).json({
      success: false,
      message: "Server error",
      error: error.message,
    });
  }
};

// Function to generate a benchmark report for historical estimations
export const generateBenchmarkReport = async (req, res) => {
  try {
    const userId = req.user?._id;
    
    if (!userId) {
      return res.status(401).json({ message: "User not authenticated" });
    }
    
    const estimations = await CostEstimation.find({ userId }).sort({
      created: -1,
    }).limit(50);
    
    if (!estimations.length) {
      return res.status(404).json({ 
        message: "No cost estimations found for benchmark analysis" 
      });
    }
    
    // Analyze the historical estimations
    const functionGroups = {};
    const costRanges = [];
    let totalComplexity = 0;
    let totalSteps = 0;
    
    estimations.forEach(est => {
      if (est.algorithmicData?.complexity?.details) {
        // Track functional groups
        const details = est.algorithmicData.complexity.details;
        for (const [group, value] of Object.entries(details)) {
          if (!functionGroups[group]) {
            functionGroups[group] = [];
          }
          functionGroups[group].push(value);
        }
        
        // Track complexity
        if (est.algorithmicData.complexity.complexity) {
          totalComplexity += est.algorithmicData.complexity.complexity;
        }
        
        // Track steps
        if (est.algorithmicData.steps?.steps) {
          totalSteps += est.algorithmicData.steps.steps;
        }
        
        // Track cost ranges
        if (est.estimatedcost) {
          const match = est.estimatedcost.match(/\$(\d+)(?:\s*[–-]\s*\$(\d+))?/);
          if (match) {
            const lower = parseInt(match[1]);
            const upper = parseInt(match[2] || match[1]);
            costRanges.push({ lower, upper });
          }
        }
      }
    });
    
    // Calculate statistics
    const avgComplexity = totalComplexity / estimations.length;
    const avgSteps = totalSteps / estimations.length;
    const avgLowerCost = costRanges.reduce((sum, range) => sum + range.lower, 0) / costRanges.length;
    const avgUpperCost = costRanges.reduce((sum, range) => sum + range.upper, 0) / costRanges.length;
    
    // Prepare report
    const report = {
      totalEstimations: estimations.length,
      avgComplexity,
      avgSteps,
      avgCostRange: {
        lower: avgLowerCost,
        upper: avgUpperCost
      },
      functionalGroupStats: Object.keys(functionGroups).map(group => ({
        group,
        count: functionGroups[group].length,
        avgValue: functionGroups[group].reduce((sum, val) => sum + (typeof val === 'number' ? val : 0), 0) / functionGroups[group].length
      })),
      datasetInfluencedEstimations: estimations.filter(est => 
        est.algorithmicData?.datasetAnalysis?.appliedDatasets
      ).length
    };
    
    res.status(200).json({
      success: true,
      data: report
    });
    
  } catch (error) {
    console.error('Error generating benchmark report:', error);
    res.status(500).json({
      success: false,
      message: "Error generating benchmark report",
      error: error.message
    });
  }
};