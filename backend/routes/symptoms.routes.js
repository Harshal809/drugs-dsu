import express from 'express';
import { getSymptoms, getDiagnosis } from '../controllers/symptoms.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';

const router = express.Router();

router.get('/list',protectRoute, getSymptoms);
router.post('/diagnosis', protectRoute, getDiagnosis);

export default router;