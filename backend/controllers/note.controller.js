import Note from "../models/note.model.js";
import { GoogleGenerativeAI } from "@google/generative-ai";

const genAI = new GoogleGenerativeAI(process.env.GEMINI_API_KEY);

export const getNotes = async (req, res) => {
  try {
    const notes = await Note.find({ userId: req.params.userId });
    res.json(notes);
  } catch (error) {
    console.error('Error fetching notes:', error);
    res.status(500).json({ error: 'Failed to fetch notes: ' + error.message });
  }
};

export const createNote = async (req, res) => {
  try {
    const { userId, content } = req.body;
    if (!userId || !content) {
      return res.status(400).json({ error: 'userId and content are required' });
    }
    const note = new Note({ userId, content });
    await note.save();
    res.status(201).json(note);
  } catch (error) {
    console.error('Error saving note:', error);
    res.status(500).json({ error: 'Failed to save note: ' + error.message });
  }
};

export const updateNote = async (req, res) => {
  try {
    const { content } = req.body;
    if (!content) {
      return res.status(400).json({ error: 'Content is required' });
    }
    const note = await Note.findByIdAndUpdate(req.params.id, { content }, { new: true });
    if (!note) {
      return res.status(404).json({ error: 'Note not found' });
    }
    res.json(note);
  } catch (error) {
    console.error('Error updating note:', error);
    res.status(500).json({ error: 'Failed to update note: ' + error.message });
  }
};

export const analyzeSentiment = async (req, res) => {
  try {
    const { text } = req.body;
    if (!text) {
      return res.status(400).json({ error: 'Text is required' });
    }
    const model = genAI.getGenerativeModel({ model: "gemini-pro" });
    const result = await model.generateContent(`Analyze the sentiment of this text and return 'positive', 'negative', or 'neutral': ${text}`);
    const sentiment = result.response.text().toLowerCase();
    res.json({ sentiment });
  } catch (error) {
    console.error('Error with Gemini sentiment analysis:', error.message);
    const text = req.body.text || '';
    const sentiment = text.includes('happy') || text.includes('good') ? 'positive' : 
                     text.includes('sad') || text.includes('bad') ? 'negative' : 'neutral';
    res.json({ sentiment });
  }
};

export const summarizeText = async (req, res) => {
  try {
    const { text } = req.body;
    if (!text) {
      return res.status(400).json({ error: 'Text is required' });
    }
    const model = genAI.getGenerativeModel({ model: "gemini-pro" });
    const result = await model.generateContent(`Summarize this note in 100 characters or less: ${text}`);
    const summary = result.response.text();
    res.json({ summary });
  } catch (error) {
    console.error('Error with Gemini summarization:', error.message);
    const text = req.body.text || '';
    const summary = text.split(' ').slice(0, 10).join(' ') + '...';
    res.json({ summary });
  }
};