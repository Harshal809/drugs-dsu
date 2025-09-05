import express from "express";
import dotenv from "dotenv";
import cors from "cors";
import cookieParser from "cookie-parser";
import path from "path";
import { connectionDb } from "./lib/db.js";
import multer from "multer";
import axios from "axios"; // Added axios for proxy requests

import proteinRoutes from "./routes/proteinstructure.routes.js";
import authRoutes from "./routes/auth.routes.js";
import costestiminationRoutes from "./routes/costestimination.routes.js";
import newsRoutes from "./routes/news.routes.js";
import alphafoldRoutes from "./routes/alphafold.routes.js";
import toxicityRoutes from "./routes/toxicity.routes.js";
import summaryRoutes from "./routes/summary.routes.js";
import noteRoutes from "./routes/note.routes.js";
import newdrugRoutes from "./routes/newdrug.routes.js";
import getsymptomproductRoutes from "./routes/getsymptomproduct.routes.js";
import drugNameRoutes from "./routes/drugnaming.routes.js";
import researchPaperRoutes from "./routes/researchPapers.routes.js";

dotenv.config();

const app = express();
const PORT = process.env.PORT || 5000;
const __dirname = path.resolve();

const upload = multer({ storage: multer.memoryStorage() });

// Dynamic CORS configuration for Render
const allowedOrigins = [
  process.env.CLIENT_URL,
  'https://drugs-10979.firebaseapp.com',
  'http://localhost:5173',
].filter(Boolean);

const corsOptions = {
  origin: function (origin, callback) {
    if (!origin) return callback(null, true);
    if (
      allowedOrigins.some(
        (allowedOrigin) =>
          origin === allowedOrigin ||
          origin.startsWith(allowedOrigin) ||
          origin.includes('render.com')
      )
    ) {
      callback(null, true);
    } else {
      callback(new Error('Not allowed by CORS'));
    }
  },
  credentials: true,
  methods: ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
  allowedHeaders: ['Content-Type', 'Authorization'],
};

// Apply CORS middleware
app.use(cors(corsOptions));

// Handle preflight requests for all routes
app.options('*', cors(corsOptions));

app.use(express.json());
app.use(express.urlencoded({ extended: true }));
app.use(cookieParser());

// Serve static files from the 'jobs' directory under '/results'
app.use('/results', express.static(path.join(__dirname, 'jobs')));

// Proxy route for NVIDIA MolMIM API
app.post("/api/proxy/molmim", async (req, res) => {
  try {
    const MOLMIM_API_URL = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate";
    const MOLMIM_API_KEY = process.env.MOLMIM_API_KEY;

    if (!MOLMIM_API_KEY) {
      return res.status(500).json({ message: "MolMIM API key is not configured" });
    }

    const response = await axios.post(MOLMIM_API_URL, req.body, {
      headers: {
        'Authorization': `Bearer ${MOLMIM_API_KEY}`,
        'Accept': 'application/json',
        'Content-Type': 'application/json',
      },
    });

    res.status(200).json(response.data);
  } catch (error) {
    console.error('Error proxying MolMIM API request:', error.message);
    res.status(error.response?.status || 500).json({
      message: 'Failed to make proxy MolMIM API request',
      error: error.message,
    });
  }
});

// API Routes
app.use("/api/protein", upload.single("file"), proteinRoutes);
app.use("/api/auth", authRoutes);
app.use("/api/costestimation", costestiminationRoutes);
app.use("/api/news", newsRoutes);
app.use("/api/alphafold", alphafoldRoutes);
app.use("/api/toxicity", toxicityRoutes);
app.use("/api/summary", summaryRoutes);
app.use("/api/notes", noteRoutes);
app.use("/api/newdrug", newdrugRoutes);
app.use("/api/getdata", getsymptomproductRoutes);
app.use("/api/drugname", drugNameRoutes);
app.use("/api/researchPaper", researchPaperRoutes);

app.listen(PORT, () => {
  connectionDb();
  console.log(`Server running on port ${PORT}`);
  console.log(`Allowed origins: ${allowedOrigins.join(', ')}`);
});