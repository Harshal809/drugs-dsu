import express from "express";
import { sendMessage, getMessages, getUsers } from "../controllers/message.controller.js";
import { protectRoute } from "../middleware/auth.middleware.js";

const router = express.Router();

router.post("/send",  sendMessage);
router.get("/:userId",  getMessages);
router.get("/users/list",  getUsers);

export default router;