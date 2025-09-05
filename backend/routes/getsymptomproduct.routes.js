// routes/costEstimationRoutes.js
import express from "express";

import {

  getSymptomsAndAllProducts
} from "../controllers/getsymptomproduct.controller.js";

const router = express.Router();


router.get("/getsymptoms-product/:id",  getSymptomsAndAllProducts);


export default router;