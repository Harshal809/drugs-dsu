// models/DrugDevelopmentProgress.js
import mongoose from 'mongoose';

const drugDevelopmentProgressSchema = new mongoose.Schema({
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: 'User',
    required: true,
  },
  moleculeId: {
    type: String,
    required: true,
  },
  progress: {
    type: Number,
    default: 0,
    min: 0,
    max: 100,
  },
  completedSections: {
    type: [String],
    enum: [
      'newDrug',
      'costEstimation',
      'drugName',
      'researchPaper',
      'generatedResearchPaper',
      'targetPrediction',
      'toxicityResult',
      'predictDisease',
      'targetProtein',
      'reactionResponse',
    ],
    default: [],
  },
  createdAt: {
    type: Date,
    default: Date.now,
  },
  updatedAt: {
    type: Date,
    default: Date.now,
  },
});

drugDevelopmentProgressSchema.pre('save', function (next) {
  this.updatedAt = Date.now();
  next();
});

export default mongoose.model('DrugDevelopmentProgress', drugDevelopmentProgressSchema);