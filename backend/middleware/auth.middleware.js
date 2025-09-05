import jwt from "jsonwebtoken";
import { User } from "../models/auth.model.js";

export const protectRoute = async (req, res, next) => {
  let token;

  if (req.cookies?.token) {
    token = req.cookies.token;
  } else if (req.headers.authorization?.startsWith("Bearer")) {
    token = req.headers.authorization.split(" ")[1];
  }

  if (!token) {
    return res.status(401).json({ success: false, message: "Not authorized, please login" });
  }

  try {
    const decoded = jwt.verify(token, process.env.JWT_SECRET);
    const user = await User.findById(decoded.userId).select("-password");

    if (!user) {
      return res.status(401).json({ success: false, message: "User no longer exists" });
    }

    req.user = user;
    next();
  } catch (error) {
    let message = "Not authorized, invalid token";
    if (error.name === "TokenExpiredError") {
      message = "Session expired, please login again";
    }
    return res.status(401).json({ success: false, message });
  }
};
