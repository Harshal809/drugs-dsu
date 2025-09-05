import express from 'express';
import { analyzeSentiment, createNote, getNotes, summarizeText, updateNote } from '../controllers/note.controller.js';

const router = express.Router();

router.get('/:userId',getNotes);
router.post('/', createNote);
router.put('/:id', updateNote);
router.post('/analyze/sentiment',analyzeSentiment);
router.post('/analyze/summary',summarizeText);
export default router;