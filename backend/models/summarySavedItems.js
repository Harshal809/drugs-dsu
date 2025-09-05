// models/summarySavedItems.js
import mongoose from 'mongoose';

const summarySavedItemSchema = new mongoose.Schema({
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    required: true,
    ref: 'User',
  },
  type: {
    type: String,
    required: true,
    enum: [
      'newDrug',
      'costEstimation',
      'drugName',
      'researchPaper',
      'generatedResearchPaper',
      
      'toxicityResult',
      'predictDisease', // New type
      'targetProtein', // New type
      'reactionResponse', // New type
    ],
  },
  moleculeId: {
    type: String,
    required: true,
  },
  data: {
    type: mongoose.Schema.Types.Mixed,
    required: true,
  },
  createdAt: {
    type: Date,
    default: Date.now,
  },
});

export default mongoose.model('SummarySavedItem', summarySavedItemSchema);