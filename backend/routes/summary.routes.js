// routes/summary.routes.js
import express from 'express';
import { bulkSaveItems, getMoleculeProgress, getSavedItems, getSummary, saveItem } from '../controllers/summary.controllers.js';
import { protectRoute } from '../middleware/auth.middleware.js';


const router = express.Router();

router.get('/summary',protectRoute, getSummary);
router.post('/save',protectRoute, saveItem);
router.get('/saved',protectRoute,  getSavedItems);
router.post('/bulk-save',protectRoute, bulkSaveItems);
router.get('/molecule-progress',protectRoute,getMoleculeProgress);

export default router;