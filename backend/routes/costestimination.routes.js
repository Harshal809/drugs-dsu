// routes/costEstimationRoutes.js
import express from "express";
import { protectRoute } from "../middleware/auth.middleware.js";
import {
  postCostEstimation,
  getCostEstimation,
  // getSymptomsAndAllProducts
} from "../controllers/costestimination.controller.js";

const router = express.Router();

router.post("/cost-estimation", protectRoute, postCostEstimation);
router.get("/getcostestimation/:userId",protectRoute,  getCostEstimation);
// router.get("/getsymptoms-product/:id", protectRoute, getSymptomsAndAllProducts);


export default router;