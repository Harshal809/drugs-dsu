import express from "express";
import {
  generateDrugName,
  getSavedDrugNames,
  checkSavedDrugName,
  acceptDrugName,
  savePendingDrugName, // New controller for checking if a search exists
  deleteDrugName,
} from "../controllers/drugnaming.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";


const router = express.Router();



// ai-naming
router.post("/generate-drug-name/:id",protectRoute, generateDrugName);
router.post("/accept-drug-name/:id",protectRoute,  acceptDrugName);
router.post("/save-pending-drug-name/:id",protectRoute, savePendingDrugName);
router.get("/saved-drug-names",protectRoute,getSavedDrugNames);
router.get("/check-saved-drug-name",protectRoute, checkSavedDrugName);
router.delete("/delete-drug-name/:id", protectRoute, deleteDrugName);
export default router;