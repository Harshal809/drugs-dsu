// import mongoose from "mongoose";

// const savedSearchSchema = new mongoose.Schema({
//   userId: { type: mongoose.Schema.Types.ObjectId, ref: "User", required: true },
//   smiles: { type: String, required: true },
//   targets: [{
//     protein: { type: String, required: true },
//     confidence: { type: Number, required: true },
//     moa: { type: String, required: true },
//     pathways: [{ type: String }],
//     diseases: [{ type: String }],
//     knownInteractions: { type: Object, default: null },
//   }],
//   research: [{
//     title: { type: String, required: true },
//     authors: { type: String, required: true },
//     journal: { type: String },
//     year: { type: String },
//     abstract: { type: String },
//     doi: { type: String },
//     url: { type: String },
//   }],
//   docking: {
//     bindingEnergy: { type: String },
//     pose: { type: String },
//     details: { type: String },
//   },
//   createdAt: { type: Date, default: Date.now },
// });

// export const SavedSearch = mongoose.model("SavedSearch", savedSearchSchema);

import mongoose from "mongoose";

const savedSearchSchema = new mongoose.Schema({
  userId: { type: mongoose.Schema.Types.ObjectId, ref: "User", required: true },
  molecule: {
    symptoms: { type: String, required: true }, // Store symptoms as a string (e.g., "fever, cough")
    smiles: { type: String, required: true },   // Store SMILES string
  },
  research: [{
    title: { type: String, required: false },
    authors: { type: String, required: true },
    year: { type: String, required: true },
    abstract: { type: String, default: "Abstract not available" },
    doi: { type: String },
    url: { type: String },
    is_simulated: { type: Boolean, default: false },
  }],
  paper: {
    title: { type: String },
    authors: { type: String },
    abstract: { type: String },
    keywords: [{ type: String }],
    introduction: { type: String },
    methodology: { type: String },
    resultsAndDiscussion: { type: String },
    conclusion: { type: String },
    references: [{ type: String }],
  },
  createdAt: { type: Date, default: Date.now },
});

export const SavedSearch = mongoose.model("SavedSearch", savedSearchSchema);