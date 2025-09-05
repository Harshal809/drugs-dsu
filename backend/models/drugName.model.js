import mongoose from "mongoose";

const DrugNameSchema = new mongoose.Schema({
  smiles: {
    type: String,
    required: true,
  },
  symptoms: {
    type: String,
    required: true,
  },
  suggestedName: {
    type: String,
    required: true,
  },
  namingDetails: {
    type: String,
    default: "",
  },
  status: {
    type: String,
    enum: ["pending", "accepted"],
    default: "pending",
  },
  createdAt: {
    type: Date,
    default: Date.now,
  },
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: "User",
    required: true,
  },
});

export default mongoose.model("DrugName", DrugNameSchema);