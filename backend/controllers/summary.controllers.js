// controllers/summary.controllers.js
import { v4 as uuidv4 } from 'uuid';
import GeneratenewMolecule from '../models/generatenew.model.js';
import CostEstimation from '../models/costestimination.model.js';
import DrugName from '../models/drugName.model.js';
import ResearchPaper from '../models/researchPapers.model.js';
import { SavedSearch } from '../models/savedSearch.model.js';
import Toxicity from '../models/toxicity.model.js';
import GeneratedResearchPaper from '../models/GeneratedResearchPaper.js';
import SummarySavedItem from '../models/summarySavedItems.js';
import { PredictDisease } from '../models/newdrug.model.js';
import { TargetProtein} from '../models/newdrug.model.js'
import DrugDevelopmentProgress from '../models/DrugDevelopmentProgress.model.js';

const VALID_TYPES = [
  'newDrug',
  'costEstimation',
  'drugName',
  'researchPaper',
  'generatedResearchPaper',
  'toxicityResult',
  'predictDisease',
  'targetProtein',
];

const updateProgress = async (userId, moleculeId, type) => {
  const progress = await DrugDevelopmentProgress.findOne({ userId, moleculeId });
  const allSections = VALID_TYPES;

  if (progress) {
    if (!progress.completedSections.includes(type)) {
      progress.completedSections.push(type);
      progress.progress = Math.round((progress.completedSections.length / allSections.length) * 100);
      await progress.save();
    }
  } else {
    await DrugDevelopmentProgress.create({
      userId,
      moleculeId,
      completedSections: [type],
      progress: Math.round((1 / allSections.length) * 100),
    });
  }
};

export const getSummary = async (req, res) => {
  try {
    const userId = req.user?._id;
    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    const [
      newDrugs,
      costEstimations,
      drugNames,
      researchPapers,
      generatedResearchPapers,
      predictDiseases,
      targetProteins,
    ] = await Promise.all([
      GeneratenewMolecule.find({ userId }).sort({ createdAt: -1 }).lean(),
      CostEstimation.find({ userId }).sort({ createdAt: -1 }).lean(),
      DrugName.find({ userId }).sort({ createdAt: -1 }).lean(),
      ResearchPaper.find({ userId }).sort({ createdAt: -1 }).lean(),
      GeneratedResearchPaper.find({ userId }).sort({ createdAt: -1 }).lean(),
      PredictDisease.find({ userId }).sort({ createdAt: -1 }).lean(),
      TargetProtein.find({ userId }).sort({ createdAt: -1 }).lean(),
    ]);

    res.status(200).json({
      success: true,
      data: {
        newDrugs,
        costEstimations,
        drugNames,
        researchPapers,
        generatedResearchPapers,
        predictDiseases,
        targetProteins,
      },
    });
  } catch (error) {
    console.error('Error fetching summary:', {
      message: error.message,
      stack: error.stack,
      userId: req.user?._id,
    });
    res.status(500).json({
      success: false,
      message: 'Failed to fetch summary data',
      error: error.message,
    });
  }
};

export const saveItem = async (req, res) => {
  try {
    const { item, type, moleculeId } = req.body;
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    if (!item || !type) {
      return res.status(400).json({ success: false, message: 'Item and type are required' });
    }

    if (!VALID_TYPES.includes(type)) {
      return res.status(400).json({ success: false, message: `Invalid type: ${type}` });
    }

    const savedItem = await SummarySavedItem.create({
      userId,
      type,
      moleculeId: moleculeId || uuidv4(),
      data: item,
    });

    await updateProgress(userId, savedItem.moleculeId, type);

    res.status(201).json({ success: true, data: savedItem });
  } catch (error) {
    console.error('Error saving item:', error);
    res.status(500).json({ success: false, message: 'Failed to save item', error: error.message });
  }
};

export const bulkSaveItems = async (req, res) => {
  try {
    const items = req.body;
    const userId = req.user?._id;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    if (!Array.isArray(items) || items.length === 0) {
      return res.status(400).json({ success: false, message: 'Items array is required and cannot be empty' });
    }

    for (const { item, type, moleculeId } of items) {
      if (!item || !type) {
        return res.status(400).json({ success: false, message: 'Each item must have item and type' });
      }
      if (!VALID_TYPES.includes(type)) {
        return res.status(400).json({ success: false, message: `Invalid type: ${type}` });
      }
    }

    const itemsToSave = items.map(({ item, type, moleculeId }) => ({
      userId,
      type,
      moleculeId: moleculeId || uuidv4(),
      data: item,
    }));

    const savedItems = await SummarySavedItem.insertMany(itemsToSave);

    const progressUpdates = itemsToSave.reduce((acc, { type, moleculeId }) => {
      if (!acc[moleculeId]) acc[moleculeId] = new Set();
      acc[moleculeId].add(type);
      return acc;
    }, {});

    for (const [moleculeId, types] of Object.entries(progressUpdates)) {
      for (const type of types) {
        await updateProgress(userId, moleculeId, type);
      }
    }

    res.status(201).json({ success: true, data: savedItems });
  } catch (error) {
    console.error('Error bulk saving items:', error);
    res.status(500).json({ success: false, message: 'Failed to bulk save items', error: error.message });
  }
};

export const getSavedItems = async (req, res) => {
  try {
    const userId = req.user?._id;
    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    const savedItems = await SummarySavedItem.find({ userId }).sort({ createdAt: -1 }).lean();

    const groupedItems = savedItems.reduce((acc, item) => {
      if (!acc[item.moleculeId]) {
        acc[item.moleculeId] = [];
      }
      acc[item.moleculeId].push(item);
      return acc;
    }, {});

    res.status(200).json({
      success: true,
      data: groupedItems,
    });
  } catch (error) {
    console.error('Error fetching saved items:', error);
    res.status(500).json({
      success: false,
      message: 'Failed to fetch saved items',
      error: error.message,
    });
  }
};

export const getMoleculeProgress = async (req, res) => {
  try {
    const userId = req.user?._id;
    const { moleculeId } = req.query;

    if (!userId) {
      return res.status(401).json({ success: false, message: 'User not authenticated' });
    }

    const progress = await DrugDevelopmentProgress.findOne({ userId, moleculeId }).lean();
    const latestMolecule = await GeneratenewMolecule.findOne({ userId }).sort({ createdAt: -1 }).lean();

    res.status(200).json({
      success: true,
      data: {
        progress: progress ? { percentage: progress.progress, completedSections: progress.completedSections } : { percentage: 0, completedSections: [] },
        molecule: latestMolecule || null,
      },
    });
  } catch (error) {
    console.error('Error fetching molecule progress:', error);
    res.status(500).json({
      success: false,
      message: 'Failed to fetch molecule progress',
      error: error.message,
    });
  }
};