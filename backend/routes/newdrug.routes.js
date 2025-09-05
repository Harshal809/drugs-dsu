import express from 'express';
import { predictDisease,predictTargetProtein, getnewdrug,getSymptoms } from '../controllers/newdrug.controller.js';

const router = express.Router();

router.post("/predictDisease/:id",predictDisease);
router.post("/predictTargetProtein",predictTargetProtein);
router.get('/getnewdrug/:id', getnewdrug);
router.get('/symptoms/:id', getSymptoms);


export default router;