import mongoose from 'mongoose';

const ProteinStructureSchema = new mongoose.Schema({
  name: {
    type: String,
    required: true,
    trim: true
  },
  smiles: {
    type: String,
    required: true
  },
  properties: {
    type: Object,
    default: {}
  },
  generatedStructures: [{
    smiles: String,
    properties: Object,
    similarity: Number,
    created: {
      type: Date,
      default: Date.now
    }
  }],
  information: {
    type: String, // Store the detailed information as a string
    default: ''
  },
  created: {
    type: Date,
    default: Date.now
  },
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: 'User',
  }
});

export default mongoose.model('ProteinStructure', ProteinStructureSchema);