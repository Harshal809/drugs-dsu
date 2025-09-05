import  express from  'express';
import {getTopHeadlines, getEverything, getSources, getSavedNews} from '../controllers/news.controller.js';
import { protectRoute } from '../middleware/auth.middleware.js';
const router = express.Router();


router.get('/top-headlines', protectRoute,getTopHeadlines);
router.get('/everything',protectRoute, getEverything);
router.get('/sources',protectRoute,getSources);
router.get('/saved',protectRoute, getSavedNews);


export default router;