import express from "express";

const router = express.Router();
import {
    proxyGeminiRequest,
    saveResearchPapers,
    getSavedResearchPapers,
    checkSavedPapers,
    saveGeneratedResearchPaper,
    getSavedGeneratedResearchPapers,
    checkSavedGeneratedPapers,
} from "../controllers/researchPapers.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";



// research paper
router.post("/proxy/gemini",protectRoute, proxyGeminiRequest);
router.post("/save-research-papers",protectRoute, saveResearchPapers);
router.get("/saved-research-papers",protectRoute, getSavedResearchPapers);
router.get("/check-saved-papers",protectRoute, checkSavedPapers);
router.post("/save-generated-research-paper",protectRoute, saveGeneratedResearchPaper);
router.get("/saved-generated-research-papers",protectRoute, getSavedGeneratedResearchPapers);
router.get("/check-saved-generated-papers",protectRoute,  checkSavedGeneratedPapers);

export default router;