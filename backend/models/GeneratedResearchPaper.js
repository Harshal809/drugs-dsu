// // GeneratedResearchPaperSchema.js
// import mongoose from "mongoose";

// const GeneratedResearchPaperSchema = new mongoose.Schema({
//   userId: { type: mongoose.Schema.Types.ObjectId, ref: "User", required: true },
//   molecule: {
//     title: { type: String, required: true },
//     smiles: { type: String, required: true },
//   },
//   paper: {
//     title: String,
//     authors: String,
//     abstract: String,
//     keywords: [String],
//     introduction: String,
//     methodology: String,
//     resultsAndDiscussion: String,
//     conclusion: String,
//     references: [String],
//   },
//   createdAt: { type: Date, default: Date.now },
// });

// export default mongoose.model("GeneratedResearchPaper", GeneratedResearchPaperSchema);


import mongoose from "mongoose";

const GeneratedResearchPaperSchema = new mongoose.Schema({
  userId: { type: mongoose.Schema.Types.ObjectId, ref: "User", required: true },
  molecule: {
    symptoms: { type: String, required: true },
    smiles: { type: String, required: true },
  },
  paper: {
    title: String,
    authors: String,
    abstract: String,
    keywords: [String],
    introduction: String,
    methodology: String,
    resultsAndDiscussion: String,
    conclusion: String,
    references: [String],
  },
  createdAt: { type: Date, default: Date.now },
});

export default mongoose.model("GeneratedResearchPaper", GeneratedResearchPaperSchema);