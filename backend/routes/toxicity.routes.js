import express from 'express';
import { protectRoute } from '../middleware/auth.middleware.js';
import {
  getToxicityHistory,
  predictToxicityController,
  saveGeminiAnalysis,
  getGeminiAnalysis,
} from '../controllers/toxicity.controllers.js';

const router = express.Router();

router.post('/predict', protectRoute, predictToxicityController);
router.get('/history', protectRoute, getToxicityHistory);
router.post('/save-analysis', protectRoute, saveGeminiAnalysis);
router.post('/gemini-analysis',protectRoute,  getGeminiAnalysis);

export default router;