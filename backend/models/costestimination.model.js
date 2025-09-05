import mongoose from 'mongoose';

const costEstimationSchema = new mongoose.Schema({
  smiles: {
    type: String,
    required: [true, 'SMILES string is required'],
    trim: true,
  },
  estimatedcost: {
    type: String, // Changed from 'estiminatedcost' to 'estimatedcost'
    required: [true, 'Estimated cost is required'],
  },
  information: {
    type: String,
    required: [true, 'Cost estimation information is required'],
  },
  userId: {
    type: mongoose.Schema.Types.ObjectId,
    ref: 'User',
    required: [true, 'User ID is required'],
  },
  created: {
    type: Date,
    default: Date.now,
  }
});

const CostEstimation = mongoose.model('CostEstimation', costEstimationSchema);

export default CostEstimation;