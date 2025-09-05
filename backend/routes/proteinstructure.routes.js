import express from "express";
import {
  getProteinStructure,
  postProteinStructure,
  generatenewmolecule,
  getgeneratednewmolecule,
  convertFileToSmiles, // New controller for converting MOL/SDF to SMILES
  getFingerprints, // New controller for fingerprint extraction (mocked)
  performDocking, // New controller for molecular docking (mocked)
  saveSearch, // New controller for saving searches
  getSavedSearches, // New controller for retrieving saved searches
  checkSavedSearches,

} from "../controllers/proteinstructure.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";
import  ProteinStructure  from "../models/protienstructure.model.js";
import { User } from "../models/auth.model.js";

const router = express.Router();

// Existing Routes
router.get("/getproteinstructure/:id",protectRoute,  getProteinStructure);
router.post("/postproteinstructure/:id",protectRoute,  postProteinStructure);
router.post('/saveproteinstructure/:id', async (req, res) => {
  try {
    const { id } = req.params;
    const { name, smiles, properties, generatedStructures, information } = req.body;

    // Validate required fields
    if (!smiles) {
      return res.status(400).json({ message: 'SMILES string is required for generation' });
    }

    // Create new protein structure
    const newProteinStructure = new ProteinStructure({
      name: name || 'Untitled Structure',
      smiles,
      properties: properties || {},
      generatedStructures: generatedStructures || [],
      information: information || '',
      userId: id,
    });

    // Save to database
    await newProteinStructure.save();

    // Update user's protein structures
    await User.findByIdAndUpdate(
      id,
      { $push: { protienStructures: newProteinStructure._id } }, // Note: typo 'protien' in your original code
      { new: true }
    );

    res.status(201).json({
      message: 'Protein structure saved successfully',
      structure: newProteinStructure,
    });
  } catch (error) {
    console.error('Error saving protein structure:', error.message);
    res.status(500).json({
      message: 'Failed to save protein structure',
      error: error.message,
    });
  }
});

router.post("/generatenewmolecule/:id",protectRoute, generatenewmolecule);
router.get("/generatednewmolecule", protectRoute, getgeneratednewmolecule);

// // research paper
// router.post("/proxy/gemini", proxyGeminiRequest);
// router.post("/save-research-papers", protectRoute, saveResearchPapers);
// router.get("/saved-research-papers", protectRoute, getSavedResearchPapers);
// router.get("/check-saved-papers", protectRoute, checkSavedPapers);
// router.post("/save-generated-research-paper", protectRoute, saveGeneratedResearchPaper);
// router.get("/saved-generated-research-papers", protectRoute, getSavedGeneratedResearchPapers);
// router.get("/check-saved-generated-papers", protectRoute, checkSavedGeneratedPapers);

// New Routes for AI-Driven Target Prediction
router.post("/convert-file-to-smiles",protectRoute,  convertFileToSmiles); // Convert MOL/SDF to SMILES
router.post("/rdkit-fingerprints",protectRoute,  getFingerprints); // Extract fingerprints (mocked)
router.post("/docking",protectRoute, performDocking); // Perform molecular docking (mocked)
router.post("/save-search",protectRoute,  saveSearch); // Save search results
router.get("/saved-searches",protectRoute, getSavedSearches); // Retrieve saved searches
router.get("/check-saved-searches",protectRoute, checkSavedSearches); // Check if a search exists



export default router;