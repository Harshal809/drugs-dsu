import express from 'express';
const router = express.Router();
import {
  submitPrediction,
  getUniprotSummary,
  getUniprotAnnotations,
  getJobStatus,
  getPreviousJobs,
  // getPDBFile
} from '../controllers/alphafold.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';

router.post('/predict',protectRoute, submitPrediction);
router.get('/status/:jobId', protectRoute, getJobStatus);
router.get('/uniprot/summary/:uniprotId',protectRoute,getUniprotSummary);
router.get('/uniprot/annotations/:uniprotId',protectRoute,  getUniprotAnnotations);
router.get('/previous-jobs',protectRoute,  getPreviousJobs);
// router.get('/pdb/:jobId',  getPDBFile);

router.get('/health', (req, res) => {
  res.status(200).json({ status: 'healthy' });
});

export default router;