 // Added User model import
import GeneratenewMolecule from '../models/generatenew.model.js'; // Added GeneratenewMolecule model import
import axios from 'axios';

import drugName from '../models/drugName.model.js';
import { User } from '../models/auth.model.js';
import drugNameModel from '../models/drugName.model.js';


const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-preview-05-20:generateContent";

export const generateDrugName = async (req, res) => {
  try {
    const { id } = req.params; // User ID
    const { smiles, symptoms } = req.body;

    if (!id || !smiles || !symptoms) {
      return res.status(400).json({ message: "User ID, SMILES, and symptoms are required" });
    }

    // Normalize SMILES string (trim whitespace, standardize case if needed)
    const normalizedSmiles = smiles.trim();
    console.log(`Generating drug name for userId: ${id}, SMILES: ${normalizedSmiles}, Symptoms: ${symptoms}`);

    // Check if the name already exists for this SMILES and symptoms
    const existingName = await drugName.findOne({ userId: id, smiles: normalizedSmiles, symptoms });
    if (existingName) {
      return res.status(409).json({
        message: "Drug name already generated for this SMILES and symptoms. Check Saved Names.",
        redirectTo: "savedNames",
      });
    }

    // Fetch the molecule details from GeneratenewMolecule using SMILES
    const molecule = await GeneratenewMolecule.findOne({ userId: id, newSmiles: normalizedSmiles });
    if (!molecule) {
      console.warn(`Molecule not found for SMILES: ${normalizedSmiles}, userId: ${id}`);
      // Proceed with default molecule details to allow name generation
      const defaultMolecule = {
        newSmiles: normalizedSmiles,
        newIupacName: "Unknown IUPAC Name",
        conversionDetails: "No conversion details available",
        potentialDiseases: symptoms,
        information: "Generated based on provided SMILES and symptoms",
      };
      console.log("Using default molecule details for name generation");

      const prompt = `
        Analyze the provided SMILES structure and generate 3 candidate drug names that meet these requirements:

        1. **Structural Accuracy**  
          - Include IUPAC-recognized stems reflecting:  
            * Functional groups (e.g., '-mab' for monoclonal antibodies, '-tinib' for kinase inhibitors)  
            * Molecular topology (e.g., 'cyclo-' for cyclic structures, 'naphtha-' for naphthalene-like)  
            * Pharmacological class indicators (e.g., '-vir' for antivirals, '-zole' for antifungals)  

        2. **Ethical & Regulatory Compliance**  
          - Avoid:  
            * Cultural insensitivities/linguistic offensiveness in 10 major languages (English, Spanish, Mandarin, Hindi, Arabic, French, Russian, German, Japanese, Portuguese)  
            * Phonetic similarities to existing drugs (cross-reference FDA Orange Book)  
            * Therapeutic misrepresentation (e.g., no 'cure-' prefixes)  

        3. **Validation Requirements**  
          - Check name availability via OpenFDA API (simulated access, assume a basic uniqueness check)  
          - Confirm ≤50% phonetic similarity to existing names using Levenshtein distance  
          - Verify no trademark conflicts in USPTO database (simulated check)  

        4. **Output Format**  
          | Rank | Name Candidate | Structural Rationale | Compliance Status |  
          |------|----------------|----------------------|-------------------|  
          | 1    | [Name]         | [Matching features]  | [Pass/Fail flags] |  

        Molecule Details:  
        - SMILES: "${defaultMolecule.newSmiles}"  
        - Symptoms: "${symptoms}"  
        - IUPAC Name: "${defaultMolecule.newIupacName}"  
        - Conversion Details: "${defaultMolecule.conversionDetails}"  
        - Potential Diseases: "${defaultMolecule.potentialDiseases}"  
        - Additional Information: "${defaultMolecule.information}"  

        **Fallback Protocol**  
        If no compliant names meet criteria:  
        a) Propose modified stems with structural justification  
        b) Suggest pharmacological class alternatives  
        c) Flag need for human pharmaceutical linguist review  

        Return the response in JSON format:  
        {
          "candidates": [
            {
              "rank": 1,
              "name": "CandidateName1",
              "rationale": "Explanation of structural features and naming",
              "compliance": "Pass/Fail with details"
            },
            {
              "rank": 2,
              "name": "CandidateName2",
              "rationale": "Explanation",
              "compliance": "Pass/Fail with details"
            },
            {
              "rank": 3,
              "name": "CandidateName3",
              "rationale": "Explanation",
              "compliance": "Pass/Fail with details"
            }
          ],
          "fallback": "Optional message if fallback protocol triggered"
        }
      `;

      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: prompt }] }] },
        { headers: { "Content-Type": "application/json" } }
      );

      const geminiContent = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
      if (!geminiContent) {
        throw new Error("No content returned from Gemini API");
      }

      const jsonMatch = geminiContent.match(/{[\s\S]*}/);
      if (!jsonMatch) {
        throw new Error("No valid JSON found in Gemini response");
      }

      const { candidates, fallback } = JSON.parse(jsonMatch[0]);

      if (!candidates || candidates.length === 0) {
        return res.status(500).json({
          message: "No valid drug name candidates generated",
          fallback: fallback || "No fallback provided",
        });
      }

      // Select the top-ranked candidate (Rank 1) to save
      const topCandidate = candidates.find((c) => c.rank === 1);

      const newDrugName = new drugName({
        smiles: normalizedSmiles,
        symptoms,
        suggestedName: topCandidate.name,
        namingDetails: `${topCandidate.rationale} | Compliance: ${topCandidate.compliance}`,
        userId: id,
      });
      await newDrugName.save();

      res.status(201).json({
        message: "Drug name generated and saved successfully",
        drugName: newDrugName,
        allCandidates: candidates,
        fallback: fallback || "Molecule not found in database, used provided SMILES and symptoms",
      });
    } else {
      // Molecule found, proceed as before
      const prompt = `
        Analyze the provided SMILES structure and generate 3 candidate drug names that meet these requirements:

        1. **Structural Accuracy**  
          - Include IUPAC-recognized stems reflecting:  
            * Functional groups (e.g., '-mab' for monoclonal antibodies, '-tinib' for kinase inhibitors)  
            * Molecular topology (e.g., 'cyclo-' for cyclic structures, 'naphtha-' for naphthalene-like)  
            * Pharmacological class indicators (e.g., '-vir' for antivirals, '-zole' for antifungals)  

        2. **Ethical & Regulatory Compliance**  
          - Avoid:  
            * Cultural insensitivities/linguistic offensiveness in 10 major languages (English, Spanish, Mandarin, Hindi, Arabic, French, Russian, German, Japanese, Portuguese)  
            * Phonetic similarities to existing drugs (cross-reference FDA Orange Book)  
            * Therapeutic misrepresentation (e.g., no 'cure-' prefixes)  

        3. **Validation Requirements**  
          - Check name availability via OpenFDA API (simulated access, assume a basic uniqueness check)  
          - Confirm ≤50% phonetic similarity to existing names using Levenshtein distance  
          - Verify no trademark conflicts in USPTO database (simulated check)  

        4. **Output Format**  
          | Rank | Name Candidate | Structural Rationale | Compliance Status |  
          |------|----------------|----------------------|-------------------|  
          | 1    | [Name]         | [Matching features]  | [Pass/Fail flags] |  

        Molecule Details:  
        - SMILES: "${molecule.newSmiles}"  
        - Symptoms: "${symptoms}"  
        - IUPAC Name: "${molecule.newIupacName}"  
        - Conversion Details: "${molecule.conversionDetails}"  
        - Potential Diseases: "${molecule.potentialDiseases}"  
        - Additional Information: "${molecule.information}"  

        **Fallback Protocol**  
        If no compliant names meet criteria:  
        a) Propose modified stems with structural justification  
        b) Suggest pharmacological class alternatives  
        c) Flag need for human pharmaceutical linguist review  

        Return the response in JSON format:  
        {
          "candidates": [
            {
              "rank": 1,
              "name": "CandidateName1",
              "rationale": "Explanation of structural features and naming",
              "compliance": "Pass/Fail with details"
            },
            {
              "rank": 2,
              "name": "CandidateName2",
              "rationale": "Explanation",
              "compliance": "Pass/Fail with details"
            },
            {
              "rank": 3,
              "name": "CandidateName3",
              "rationale": "Explanation",
              "compliance": "Pass/Fail with details"
            }
          ],
          "fallback": "Optional message if fallback protocol triggered"
        }
      `;

      const geminiResponse = await axios.post(
        `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: prompt }] }] },
        { headers: { "Content-Type": "application/json" } }
      );

      const geminiContent = geminiResponse.data.candidates?.[0]?.content?.parts?.[0]?.text;
      if (!geminiContent) {
        throw new Error("No content returned from Gemini API");
      }

      const jsonMatch = geminiContent.match(/{[\s\S]*}/);
      if (!jsonMatch) {
        throw new Error("No valid JSON found in Gemini response");
      }

      const { candidates, fallback } = JSON.parse(jsonMatch[0]);

      if (!candidates || candidates.length === 0) {
        return res.status(500).json({
          message: "No valid drug name candidates generated",
          fallback: fallback || "No fallback provided",
        });
      }

      // Select the top-ranked candidate (Rank 1) to save
      const topCandidate = candidates.find((c) => c.rank === 1);

      const newDrugName = new drugName({
        smiles: normalizedSmiles,
        symptoms,
        suggestedName: topCandidate.name,
        namingDetails: `${topCandidate.rationale} | Compliance: ${topCandidate.compliance}`,
        userId: id,
      });
      await newDrugName.save();

      res.status(201).json({
        message: "Drug name generated and saved successfully",
        drugName: newDrugName,
        allCandidates: candidates,
        fallback: fallback || null,
      });
    }
  } catch (error) {
    console.error("Error generating drug name:", error);
    res.status(500).json({ message: "Server error while generating drug name", error: error.message });
  }
};

export const acceptDrugName = async (req, res) => {
  try {
    if (!req.user) {
      return res.status(401).json({ message: "Unauthorized: No user authenticated" });
    }
    const userId = req.user._id; // Use req.user._id instead of req.params.userId
    const { smiles, symptoms, selectedName, rationale, compliance } = req.body;

    // Verify user exists (optional, since req.user is already validated)
    const user = await User.findById(userId);
    if (!user) {
      return res.status(404).json({ message: "User not found" });
    }

    // Check if an accepted name already exists
    const existingAccepted = await drugNameModel.findOne({ smiles, userId, status: "accepted" });
    if (existingAccepted) {
      return res.status(409).json({ message: "An accepted drug name already exists for this SMILES" });
    }

    // Update or create GeneratenewMolecule document
    const molecule = await GeneratenewMolecule.findOneAndUpdate(
      { userId, $or: [{ originalSmiles: smiles }, { newSmiles: smiles }] },
      {
        $set: {
          newSmiles: selectedName,
          originalSmiles: smiles,
          newmoleculetitle: selectedName,
          potentialDiseases: symptoms || "",
        },
      },
      { upsert: true, new: true }
    );

    // Save to DrugName collection
    const drugName = new drugNameModel({
      userId,
      smiles,
      symptoms: symptoms || "",
      suggestedName: selectedName,
      namingDetails: `Rationale: ${rationale} | Compliance: ${compliance}`,
      status: "accepted",
    });
    await drugName.save();

    // Delete pending drug names for this SMILES
    await drugNameModel.deleteMany({ smiles, userId, status: "pending" });

    return res.status(201).json({
      drugName,
      molecule,
      message: "Drug name accepted successfully",
    });
  } catch (error) {
    console.error("Error in acceptDrugName:", error);
    return res.status(500).json({
      message: "Failed to accept drug name",
      error: error.message,
    });
  }
};

export const savePendingDrugName = async (req, res) => {
  try {
    const { id } = req.params;
    const { smiles, symptoms, candidates } = req.body;

    if (!id || !smiles || !symptoms || !candidates) {
      return res.status(400).json({ message: "User ID, SMILES, symptoms, and candidates are required" });
    }

    const topCandidate = candidates.find((c) => c.rank === 1);
    const existingPendingName = await drugNameModel.findOne({ userId: id, smiles, symptoms, status: "pending" });

    if (!existingPendingName) {
      const newDrugName = new drugNameModel({
        smiles,
        symptoms,
        suggestedName: topCandidate.name,
        namingDetails: `${topCandidate.rationale} | Compliance: ${topCandidate.compliance}`,
        userId: id,
        status: "pending",
      });
      await newDrugName.save();
      res.status(201).json({ message: "Pending drug name saved successfully", drugName: newDrugName });
    } else {
      res.status(200).json({ message: "Pending drug name already exists", drugName: existingPendingName });
    }
  } catch (error) {
    console.error("Error saving pending drug name:", error);
    res.status(500).json({ message: "Server error while saving pending drug name", error: error.message });
  }
};

export const getSavedDrugNames = async (req, res) => {
  try {
    const userId = req.user._id;
    const savedNames = await drugNameModel.find({ userId }).sort({ createdAt: -1 });
    res.status(200).json({ drugNames: savedNames });
  } catch (error) {
    console.error("Error fetching saved drug names:", error);
    res.status(500).json({ message: "Server error while fetching saved drug names", error: error.message });
  }
};

export const checkSavedDrugName = async (req, res) => {
  try {
    const userId = req.user._id;
    const { smiles, symptoms } = req.query;

    if (!smiles || !symptoms) {
      return res.status(400).json({ message: "SMILES and symptoms are required" });
    }

    const normalizedSmiles = smiles.trim();
    const exists = await drugNameModel.findOne({ userId, smiles: normalizedSmiles, symptoms });
    res.status(200).json({ exists: !!exists });
  } catch (error) {
    console.error("Error checking saved drug name:", error);
    res.status(500).json({ message: "Server error while checking saved drug name", error: error.message });
  }
};

export const deleteDrugName = async (req, res) => {
  try {
    const drugNameId = req.params.id;
    // Find the drug name document before deleting
    const drugName = await drugNameModel.findById(drugNameId);
    if (!drugName) {
      return res.status(404).json({ message: "Drug name not found" });
    }

    const { smiles, userId } = drugName;

    // Delete the drug name document
    await drugNameModel.findByIdAndDelete(drugNameId);

    // Update the corresponding GeneratenewMolecule document
    await GeneratenewMolecule.findOneAndUpdate(
      { userId, $or: [{ originalSmiles: smiles }, { newSmiles: smiles }] },
      {
        $set: {
          newSmiles: drugName.smiles, // Reset to original SMILES
          newmoleculetitle: drugName.smiles, // Reset to original SMILES
        },
      }
    );

    res.status(200).json({ message: "Drug name deleted successfully" });
  } catch (error) {
    console.error("Error deleting drug name:", error);
    res.status(500).json({ message: "Error deleting drug name", error: error.message });
  }
};