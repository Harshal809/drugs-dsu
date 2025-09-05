import mongoose from "mongoose";

const researchPaperSchema = new mongoose.Schema({
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: "User",
    required: true,
  },
  molecule: {
    symptoms: {
      type: String,
      required: true,
    },
    smiles: {
      type: String,
      required: true,
    },
    // Remove or make title optional if it exists
    // title: { type: String }, // Commented out or removed
  },
  papers: [
    {
      title: String,
      authors: String,
      year: String,
      abstract: String,
      doi: String,
      url: String,
      is_simulated: Boolean,
    },
  ],
  createdAt: {
    type: Date,
    default: Date.now,
  },
});

// Add index for faster queries
researchPaperSchema.index({ userId: 1, "molecule.symptoms": 1, "molecule.smiles": 1 });

export default mongoose.model("ResearchPaper", researchPaperSchema);