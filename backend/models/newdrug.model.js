import mongoose from "mongoose";
import { User } from "./auth.model.js";

const predictDiseaseSchema = new mongoose.Schema({
  symptoms: {
    type: [String],
    required: true,
  },
   productSmileId:{
    type: mongoose.Schema.Types.ObjectId,
    ref: "ReactionResponse",
  },
  predictedDiseases: [
    {
      diseaseName: {
        type: String,
        required: true,
      },
      DiseaseMatchness: {
        type: String,
        required: true,
      },
      diseaseCautions: {
        type: [String],
        required: true,
      },
    },
  ],
  DiseaseAnalysis: {
    likelyDiseasesOrDisorders: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    associatedBiologicalPathways: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    affectedBiologicalProcesses: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    involvedOrgansAndTissues: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    relevantCellTypes: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    molecularCellularMechanismDisruptions: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    associatedGenes: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    knownOrPotentialTargetProteins: [
      {
        protein: String,
        source: String,
        summary: String,
        url: String,
      },
    ],
    biomarkers: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
    relevantTherapeuticClassesOrDrugExamples: [
      {
        source: String,
        summary: String,
        url: String,
      },
    ],
  },
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: "User",
  },

  createdAt: { type: Date, default: Date.now },
});

export const PredictDisease = mongoose.model(
  "PredictDisease",
  predictDiseaseSchema
);

// Updated TargetProteinSchema
// Updated TargetProteinSchema
const TargetProteinSchema = new mongoose.Schema({
  TargetProteins: [
    {
      proteinName: {
        type: String,
        required: true,
      },
      proteinFunction: {
        type: String,
        required: true,
      },
      associatedDiseases: {
        type: [String],
        required: true,
      },
      ProtienDiscription: {
        type: String,
        required: true,
      },
      proteinDetailedDiscription: {
        type: String,
        required: true,
      },
      ProteinSmile: {
        type: String,
        required: true,
      },
      proteinUniport: {
        type: String,
        required: true,
      },
    },
  ],
  TargetLigands: [
    {
      ligandName: {
        type: String,
        required: true,
      },
      ligandFunction: {
        type: String,
        required: true,
      },
      associatedDiseases: {
        type: [String],
        required: true,
      },
      LigandDiscription: {
        type: String,
        required: true,
      },
      LigandDetailedDiscription: {
        type: String,
        required: true,
      },
      LigandSmile: {
        type: String,
        required: true,
      },
      ligandDrugBankID: {
        // Renamed from ligandUniport to reflect the correct identifier
        type: String,
        required: true,
      },
    },
  ],
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: "User",
  },
  createdAt: { type: Date, default: Date.now },
});

export const TargetProtein = mongoose.model(
  "TargetProtein",
  TargetProteinSchema
);

// Define the ReactionResponse schema based on your sample data
const ReactionResponseSchema = new mongoose.Schema(
  {
    userId: {
      type: mongoose.Schema.Types.ObjectId,
      ref: "User",
      required: [true, "User ID is required"],
      validate: {
        validator: async function (value) {
          const user = await User.findById(value);
          return !!user;
        },
        message: "User ID must reference an existing User",
      },
    },
    symptomsId: {
      type: mongoose.Schema.Types.ObjectId,
      ref: "PredictDisease",
      validate: {
        validator: async function (value) {
          if (!value) return true; // Allow null/undefined if not required
          const predictDisease = await PredictDisease.findById(value);
          return !!predictDisease;
        },
        message: "Symptoms ID must reference an existing PredictDisease document",
      },
    },
    reactants: [
      {
        smiles: String,
        properties: {
          smiles: String,
          molecular_weight: Number,
          num_atoms: Number,
          logP: Number,
          tpsa: Number,
          num_h_donors: Number,
          num_h_acceptors: Number,
          num_rotatable_bonds: Number,
          has_stereochemistry: Boolean,
          functional_groups: [String],
        },
        admet_properties: [
          {
            section: String,
            properties: [
              {
                name: String,
                prediction: mongoose.Schema.Types.Mixed,
                units: String,
              },
            ],
          },
        ],
      },
    ],
    reactionResults: [
      {
        reactionType: String,
        description: String,
        reactants: [String],
        reactantGroups: [[String]],
        products: [
          {
            smiles: String,
            molecular_weight: Number,
            num_atoms: Number,
            logP: Number,
            tpsa: Number,
            num_h_donors: Number,
            num_h_acceptors: Number,
            num_rotatable_bonds: Number,
            has_stereochemistry: Boolean,
            functional_groups: [String],
            admet_properties: [
              {
                section: String,
                properties: [
                  {
                    name: String,
                    prediction: mongoose.Schema.Types.Mixed,
                    units: String,
                  },
                ],
              },
            ],
          },
        ],
        confidence: Number,
        productSmiles: String,
      },
    ],
    failedReactions: [
      {
        reactants: [String],
        reason: String,
        reactionType: String,
      },
    ],
    statistics: {
      mean_mw: Number,
      std_mw: Number,
      min_mw: Number,
      max_mw: Number,
    },
    createdAt: { type: Date, default: Date.now },
    processedReactionTypes: mongoose.Schema.Types.Mixed,
    product_smiles: [String],
  },
  {
    strict: true, // Enforce schema validation
    indexes: [
      { key: { userId: 1 } }, // Index on userId for faster queries
      { key: { symptomsId: 1 } }, // Index on symptomsId for faster queries
    ],
  }
);

export const ReactionResponse = mongoose.model(
  "ReactionResponse",
  ReactionResponseSchema
);
