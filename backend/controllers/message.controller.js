import { Message } from "../models/message.model.js";
import { User } from "../models/auth.model.js";
import { io } from "../server.js";

// Setup socket events
export const setupSocket = (socket) => {
  socket.on("join", (userId) => {
    const room = userId.toString(); // Convert ObjectId to string
    socket.join(room);
    console.log(`Server: User ${userId} joined room ${room}`);
  });

  socket.on("disconnect", () => {
    console.log(`Server: User disconnected: ${socket.id}`);
  });
};

// Send a message
export const sendMessage = async (req, res) => {
  try {
    const { receiverId, content } = req.body;
    const senderId = req.user?._id;

    if (!senderId || !receiverId || !content) {
      return res.status(400).json({ error: "Missing required fields" });
    }

    const message = new Message({
      sender: senderId,
      receiver: receiverId,
      content,
    });

    await message.save();

    const populatedMessage = await Message.findById(message._id)
      .populate("sender", "username _id")
      .populate("receiver", "username _id");

    console.log(`Server: Sending message from ${senderId} to ${receiverId}`);
    console.log("Server: Message data:", populatedMessage);

    // Emit to sender and receiver rooms (convert ObjectId to string)
    io.to(senderId.toString()).emit("newMessage", populatedMessage);
    io.to(receiverId.toString()).emit("newMessage", populatedMessage);

    console.log(`Server: Emitted newMessage to rooms: ${senderId} and ${receiverId}`);

    res.status(201).json(populatedMessage);
  } catch (error) {
    console.error("Server: Error in sendMessage:", error);
    res.status(500).json({ error: error.message });
  }
};

// Get messages between two users
export const getMessages = async (req, res) => {
  try {
    const { userId } = req.params;
    const currentUserId = req.user?._id;

    if (!userId || !currentUserId) {
      return res.status(400).json({ error: "Missing user IDs" });
    }

    const messages = await Message.find({
      $or: [
        { sender: currentUserId, receiver: userId },
        { sender: userId, receiver: currentUserId },
      ],
    })
      .populate("sender", "username _id")
      .populate("receiver", "username _id")
      .sort("timestamp");

    res.status(200).json(messages);
  } catch (error) {
    console.error("Server: Error in getMessages:", error);
    res.status(500).json({ error: error.message });
  }
};

// Get list of users
export const getUsers = async (req, res) => {
  try {
    const currentUserId = req.user?._id;

    if (!currentUserId) {
      return res.status(401).json({ error: "Unauthorized" });
    }

    const users = await User.find({
      _id: { $ne: currentUserId },
      username: { $exists: true, $ne: null },
    }).select("username profilePicture _id");

    res.status(200).json(users);
  } catch (error) {
    console.error("Server: Error in getUsers:", error);
    res.status(500).json({ error: error.message });
  }
};