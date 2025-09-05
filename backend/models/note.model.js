import mongoose from 'mongoose';
const noteSchema = new mongoose.Schema({
  userId: { type: String, required: true },
  content: { type: String, required: true },
  sentiment: { type: String },
  summary: { type: String },
  createdAt: { type: Date, default: Date.now },
});

const Note = mongoose.model('Note', noteSchema);
export default Note;