import dotenv from "dotenv";
import axios from "axios";
import ResearchPaper from "../models/researchPapers.model.js";
import GeneratedResearchPaper from "../models/GeneratedResearchPaper.js";

dotenv.config();

const GEMINI_API_KEY = process.env.GEMINI_API_KEY;
const GEMINI_API_URL = "https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent";

export const proxyGeminiRequest = async (req, res) => {
  try {
    const { prompt } = req.body;
    if (!prompt) {
      return res.status(400).json({ message: "Prompt is required" });
    }

    const response = await axios.post(
      `${GEMINI_API_URL}?key=${GEMINI_API_KEY}`,
      {
        contents: [{ parts: [{ text: prompt }] }],
      },
      {
        headers: { "Content-Type": "application/json" },
      }
    );

    let geminiContent = response.data.candidates?.[0]?.content?.parts?.[0]?.text;
    if (!geminiContent) {
      throw new Error("No content returned from Gemini API");
    }

    // Clean the response to remove Markdown code fences
    geminiContent = geminiContent.replace(/```json|```/g, "").replace(/^`+|`+$/g, "").trim();
    const jsonMatch = geminiContent.match(/(\[.*?\]|\{.*?\})/s);
    geminiContent = jsonMatch ? jsonMatch[0] : geminiContent;

    res.status(200).json({ content: geminiContent });
  } catch (error) {
    console.error("Error proxying Gemini request:", error);
    res.status(500).json({ message: "Failed to fetch Gemini response", error: error.message });
  }
};

export const saveResearchPapers = async (req, res) => {
  const { userId, molecule, papers } = req.body;
  try {
    if (!userId || !molecule?.symptoms || !molecule?.smiles || !papers?.length) {
      return res.status(400).json({ message: "Invalid request data" });
    }

    // Check if a document already exists for this user and molecule
    const existingDoc = await ResearchPaper.findOne({
      userId,
      "molecule.symptoms": molecule.symptoms,
      "molecule.smiles": molecule.smiles,
    });

    if (existingDoc) {
      // Update the existing document by replacing the research papers
      existingDoc.papers = papers;
      existingDoc.createdAt = new Date();
      await existingDoc.save();
    } else {
      // Create a new document
      const newDoc = new ResearchPaper({
        userId,
        molecule,
        papers,
      });
      await newDoc.save();
    }

    res.status(200).json({ message: "Papers saved successfully" });
  } catch (err) {
    console.error("Error saving research papers:", err);
    res.status(500).json({ message: "Failed to save papers", error: err.message });
  }
};

export const getSavedResearchPapers = async (req, res) => {
  const userId = req.user._id;
  try {
    const savedPapers = await ResearchPaper.find({
      userId,
      papers: { $exists: true, $ne: [] },
    });
    console.log("Fetched saved research papers:", savedPapers);
    res.status(200).json({ papers: savedPapers || [] });
  } catch (err) {
    console.error("Error fetching saved papers:", err);
    res.status(500).json({ message: "Failed to fetch saved papers", error: err.message });
  }
};

export const checkSavedPapers = async (req, res) => {
  const userId = req.user._id;
  const { symptoms, smiles } = req.query;
  try {
    if (!symptoms || !smiles) {
      return res.status(400).json({ message: "Symptoms and SMILES are required" });
    }
    const exists = await ResearchPaper.exists({
      userId,
      "molecule.symptoms": symptoms,
      "molecule.smiles": smiles,
      papers: { $exists: true, $ne: [] },
    });
    console.log("Check saved papers result:", exists);
    res.status(200).json({ exists: !!exists });
  } catch (err) {
    console.error("Error checking saved papers:", err);
    res.status(500).json({ message: "Failed to check saved papers", error: err.message });
  }
};

export const saveGeneratedResearchPaper = async (req, res) => {
  const { userId, molecule, paper } = req.body;
  try {
    if (!userId || !molecule?.symptoms || !molecule?.smiles || !paper) {
      return res.status(400).json({ message: "Invalid request data" });
    }

    // Check if a document already exists for this user and molecule
    const existingDoc = await GeneratedResearchPaper.findOne({
      userId,
      "molecule.symptoms": molecule.symptoms,
      "molecule.smiles": molecule.smiles,
    });

    if (existingDoc) {
      // Update the existing document by replacing the generated paper
      existingDoc.paper = paper;
      existingDoc.createdAt = new Date();
      await existingDoc.save();
    } else {
      // Create a new document
      const newDoc = new GeneratedResearchPaper({
        userId,
        molecule,
        paper,
      });
      await newDoc.save();
    }

    res.status(200).json({ message: "Generated research paper saved successfully" });
  } catch (err) {
    console.error("Error saving generated research paper:", err);
    res.status(500).json({ message: "Failed to save generated research paper", error: err.message });
  }
};

export const getSavedGeneratedResearchPapers = async (req, res) => {
  const userId = req.user._id;
  try {
    const savedPapers = await GeneratedResearchPaper.find({
      userId,
      paper: { $exists: true },
    });
    console.log("Fetched saved generated papers:", savedPapers);
    res.status(200).json({ papers: savedPapers || [] });
  } catch (err) {
    console.error("Error fetching saved generated papers:", err);
    res.status(500).json({ message: "Failed to fetch saved generated papers", error: err.message });
  }
};

export const checkSavedGeneratedPapers = async (req, res) => {
  const userId = req.user._id;
  const { symptoms, smiles } = req.query;

  if (!symptoms || !smiles) {
    return res.status(400).json({ message: "Symptoms and SMILES are required" });
  }

  try {
    const exists = await GeneratedResearchPaper.exists({
      userId,
      "molecule.symptoms": symptoms,
      "molecule.smiles": smiles,
      paper: { $exists: true },
    });
    console.log("Check saved generated papers result:", exists);
    res.status(200).json({ exists: !!exists });
  } catch (err) {
    console.error("Error checking saved generated papers:", err);
    res.status(500).json({ message: "Failed to check saved generated papers", error: err.message });
  }
};