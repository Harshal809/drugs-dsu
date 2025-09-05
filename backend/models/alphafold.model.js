import mongoose from "mongoose";

const JobSchema = new mongoose.Schema(
  {
    uniprot_id: {
      type: String,
      required: [true, "UniProt ID is required"],
      trim: true,
      validate: {
        validator: (v) => /^[A-Z0-9]{5,10}$/i.test(v),
        message: "Invalid UniProt ID format",
      },
    },
    status: {
      type: String,
      enum: ["pending", "processing", "completed", "failed"],
      default: "pending",
    },
    pythonJobId: {
      type: String,
      index: true,
      sparse: true,
    },
    pdbUrl: {
      type: String,
      default: null,
    },
    error: {
      type: String,
      default: null,
    },
    createdAt: {
      type: Date,
      default: Date.now,
    },
    updatedAt: {
      type: Date,
      default: Date.now,
    },
    completedAt: {
      type: Date,
      default: null,
    },
    userId: {
      type: mongoose.Schema.Types.ObjectId,
      ref: "User",
    },
  },
  {
    timestamps: false,
    versionKey: false,
  }
);

export const Job = mongoose.model("Job", JobSchema);