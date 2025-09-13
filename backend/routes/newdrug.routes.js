import express from 'express';
import { predictDisease,predictTargetProtein, getnewdrug,getSymptoms, compareDrugs } from '../controllers/prediction.controller.js';

const router = express.Router();

router.post("/predictDisease/:id",predictDisease);
router.post("/predictTargetProtein",predictTargetProtein);
router.get('/getnewdrug/:id', getnewdrug);
router.get('/symptoms/:id', getSymptoms);
// Add to your newdrug routes file
router.post('/compareDrugs/:id', compareDrugs);


export default router;