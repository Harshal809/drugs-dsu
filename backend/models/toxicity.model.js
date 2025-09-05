// models/toxicity.model.js

import mongoose from 'mongoose';

// Reusable sub-schemas for consistency
const StructuralAnalysisSchema = new mongoose.Schema({
  molecularWeight: { type: String, default: 'Unknown' },
  logP: { type: String, default: 'Unknown' },
  polarSurfaceArea: { type: String, default: 'Unknown' },
  functionalGroups: [String],
  complexity: { type: String, default: 'Unknown' },
  drugLikeness: { type: String, default: 'Unknown' },
}, { _id: false });

const AcuteToxicitySchema = new mongoose.Schema({
  LD50: { type: String, default: 'Unknown' },
  toxicityClass: { type: String, default: 'Unknown' },
  rationale: { type: String, default: 'No rationale provided.' },
  supplemental: { type: String, default: null },
}, { _id: false });

const EndpointsSchema = new mongoose.Schema({
  hepatotoxicity: { type: String, default: 'Unknown' },
  carcinogenicity: { type: String, default: 'Unknown' },
  cardiotoxicity: { type: String, default: 'Unknown' },
  neurotoxicity: { type: String, default: 'Unknown' },
  nephrotoxicity: { type: String, default: 'Unknown' },
  skin_sensitization: { type: String, default: 'Unknown' },
  thyroid_toxicity: { type: String, default: 'Unknown' },
  ototoxicity: { type: String, default: 'Unknown' },
  musculoskeletal_toxicity: { type: String, default: 'Unknown' },
  dermatotoxicity: { type: String, default: 'Unknown' },
  supplemental: { type: String, default: null },
}, { _id: false });

const SideEffectSchema = new mongoose.Schema({
  name: { type: String, required: true },
  description: { type: String, default: 'No description available.' },
  severity: { type: String, default: 'Unknown' },
  organSystem: { type: String, default: 'Systemic' },
  likelihood: { type: String, default: 'Unknown' },
}, { _id: false });

const MechanismSchema = new mongoose.Schema({
  feature: { type: String, required: true },
  pathway: { type: String, required: true },
  toxicityType: { type: String, required: true },
}, { _id: false });

const StructureConcernSchema = new mongoose.Schema({
  substructure: { type: String, required: true },
  concern: { type: String, required: true },
  riskLevel: { type: String, required: true },
}, { _id: false });

const SafetyRecommendationSchema = new mongoose.Schema({
  test: { type: String, required: true },
  priority: { type: String, required: true },
  rationale: { type: String, required: true },
}, { _id: false });

const QsarPropertySchema = new mongoose.Schema({
  name: { type: String, required: true },
  predictedValue: { type: String, required: true },
  toxicologicalImplication: { type: String, required: true },
}, { _id: false });

const QsarSymptomContextSchema = new mongoose.Schema({
  symptom: { type: String, required: true },
  structuralInsight: { type: String, required: true },
  riskAssessment: { type: String, required: true },
}, { _id: false });

const QsarAnalysisSchema = new mongoose.Schema({
  properties: [QsarPropertySchema],
  overallRisk: { type: String, default: 'Unknown' },
  symptomContext: [QsarSymptomContextSchema],
}, { _id: false });

const ToxicophoreAnalysisSchema = new mongoose.Schema({
  substructure: { type: String, required: true },
  concern: { type: String, required: true },
  prevalence: { type: String, required: true },
  mitigation: { type: String, required: true },
}, { _id: false });

const DevelopmentRecommendationSchema = new mongoose.Schema({
  phase: { type: String, required: true },
  recommendation: { type: String, required: true },
  importance: { type: String, required: true },
}, { _id: false });

// Main Toxicity Schema
const toxicitySchema = new mongoose.Schema({
  smiles: { type: String, required: true, index: true },
  compoundName: { type: String, required: true }, // FIXED: Added compoundName
  symptoms: { type: String, required: true },
  userId: { type: mongoose.Schema.Types.ObjectId, ref: 'User', required: true, index: true },
  created: { type: Date, default: Date.now },

  // This object is a summary saved for quick retrieval
  toxicityResult: {
    smiles: String,
    compoundName: String,
    symptoms: String,
    isNovelCompound: Boolean,
    structuralAnalysis: StructuralAnalysisSchema,
    acuteToxicity: AcuteToxicitySchema,
    endpoints: EndpointsSchema,
    sideEffects: [SideEffectSchema],
  },
  
  // This is the full, detailed analysis from Gemini
  geminiAnalysis: {
    isInputValid: { type: Boolean, default: false },
    inputError: { type: String, default: null },
    smiles: String,
    moleculeName: String,
    isNovelCompound: Boolean,
    structuralAnalysis: StructuralAnalysisSchema,
    acuteToxicity: AcuteToxicitySchema,
    endpoints: EndpointsSchema,
    sideEffects: [SideEffectSchema],
    mechanisms: [MechanismSchema],
    structureConcerns: [StructureConcernSchema],
    safetyRecommendations: [SafetyRecommendationSchema],
    qsarAnalysis: QsarAnalysisSchema,
    toxicophoreAnalysis: [ToxicophoreAnalysisSchema],
    developmentRecommendations: [DevelopmentRecommendationSchema],
  }
});

export default mongoose.model('Toxicity', toxicitySchema);