import mongoose from 'mongoose';

const newsSchema = new mongoose.Schema({
  title: {
    type: String,
    required: true
  },
  description: {
    type: String
  },
  url: {
    type: String,
    required: true
  },
  source: {
    type: String
  },
  publishedAt: {
    type: Date,
    default: Date.now
  },
  category: {
    type: String
  }
}, {
  timestamps: true
});

export default mongoose.model('News', newsSchema);