"use client"

import { useState, useEffect, useRef } from "react";
import { io } from "socket.io-client";
import axios from "axios";
import { useAuthStore } from "../../Store/auth.store.js";
import { toast, Toaster } from "react-hot-toast";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000";

const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

function Message() {
  const [messages, setMessages] = useState([]);
  const [newMessage, setNewMessage] = useState("");
  const [users, setUsers] = useState([]);
  const [selectedUser, setSelectedUser] = useState(null);
  const [loading, setLoading] = useState({
    auth: true,
    users: false,
    messages: false
  });
  const [error, setError] = useState(null);
  const [socketConnected, setSocketConnected] = useState(false);
  const { user, checkAuth } = useAuthStore();
  const socketRef = useRef(null);
  const messageEndRef = useRef(null);

  useEffect(() => {
    const initialize = async () => {
      try {
        await checkAuth();
        const currentUser = useAuthStore.getState().user;
        
        if (!currentUser) {
          setError("Authentication failed. Please log in.");
          setLoading(prev => ({ ...prev, auth: false }));
          return;
        }

        // Initialize socket connection with explicit namespace and path
        socketRef.current = io(`${API_BASE_URL}`, {
          path: '/socket.io', // Match server path
          transports: ['websocket', 'polling'],
          withCredentials: true,
          autoConnect: false,
          reconnectionAttempts: 5,
          reconnectionDelay: 1000,
          query: { userId: currentUser._id }, // Pass user ID for server-side validation
        });

        // Socket event handlers
        socketRef.current.on('connect', () => {
          console.log('Socket connected with ID:', socketRef.current.id);
          setSocketConnected(true);
          socketRef.current.emit("join", currentUser._id);
        });

        socketRef.current.on('disconnect', () => {
          console.log('Socket disconnected');
          setSocketConnected(false);
        });

        socketRef.current.on('connect_error', (err) => {
          console.error('Socket connection error:', err.message);
          setSocketConnected(false);
          toast.error(`Connection to chat server failed: ${err.message}. Trying to reconnect...`);
        });

        socketRef.current.on('newMessage', (message) => {
          handleNewMessage(message, currentUser._id);
        });

        // Connect socket after setting up handlers
        socketRef.current.connect();

        // Load users list
        await loadUsers();

        setLoading(prev => ({ ...prev, auth: false }));
      } catch (err) {
        setError(err.response?.data?.message || "Initialization failed");
        console.error("Initialization error:", err);
        setLoading(prev => ({ ...prev, auth: false }));
      }
    };

    initialize();

    return () => {
      if (socketRef.current) {
        socketRef.current.disconnect();
        socketRef.current.off('connect');
        socketRef.current.off('disconnect');
        socketRef.current.off('connect_error');
        socketRef.current.off('newMessage');
      }
    };
  }, [checkAuth]);

  const handleNewMessage = (message, currentUserId) => {
    if (!selectedUser) return;
    
    const isRelevant = 
      (message.sender._id === currentUserId && message.receiver._id === selectedUser._id) ||
      (message.sender._id === selectedUser._id && message.receiver._id === currentUserId);
    
    if (isRelevant) {
      setMessages(prev => {
        const exists = prev.some(msg => msg._id === message._id);
        return exists ? prev : [...prev, message];
      });
    }
  };

  const loadUsers = async () => {
    try {
      setLoading(prev => ({ ...prev, users: true }));
      const response = await axiosInstance.get("/api/message/users/list"); // Corrected endpoint
      const validUsers = response.data.filter(u => u.username && typeof u.username === "string");
      setUsers(validUsers);
    } catch (err) {
      setError(err.response?.data?.message || "Failed to fetch users");
      console.error("Error fetching users:", err);
    } finally {
      setLoading(prev => ({ ...prev, users: false }));
    }
  };

  useEffect(() => {
    const loadMessages = async () => {
      if (!selectedUser || !user) return;
      try {
        setLoading(prev => ({ ...prev, messages: true }));
        const response = await axiosInstance.get(`/api/message/${selectedUser._id}`);
        setMessages(response.data);
      } catch (err) {
        setError(err.response?.data?.message || "Failed to fetch messages");
        console.error("Error fetching messages:", err);
      } finally {
        setLoading(prev => ({ ...prev, messages: false }));
      }
    };

    loadMessages();
  }, [selectedUser, user]);

  const handleSendMessage = async (e) => {
    e.preventDefault();
    if (!newMessage.trim() || !selectedUser || !user || !socketRef.current) return;

    try {
      const response = await axiosInstance.post("/api/message/send", {
        receiverId: selectedUser._id,
        content: newMessage,
      });

      if (response.status === 201) {
        const message = response.data; // Assuming the server returns the created message
        socketRef.current.emit("sendMessage", message); // Notify other clients
        setNewMessage("");
        toast.success("Message sent successfully!", {
          position: "top-right",
          duration: 2000,
        });
      }
    } catch (error) {
      setError(error.response?.data?.message || "Failed to send message");
      console.error("Error sending message:", error);
      toast.error("Failed to send message", {
        position: "top-right",
        duration: 2000,
      });
    }
  };

  useEffect(() => {
    messageEndRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages]);

  useEffect(() => {
    const storedUser = localStorage.getItem("selectedUser");
    if (storedUser) {
      setSelectedUser(JSON.parse(storedUser));
    }
  }, []);

  useEffect(() => {
    if (selectedUser) {
      localStorage.setItem("selectedUser", JSON.stringify(selectedUser));
      toast(
        "This messaging feature is currently under development.",
        {
          position: "top-right",
          duration: 5000,
          style: {
            background: "#fefcbf",
            color: "#92400e",
            border: "1px solid #f59e0b",
          },
          icon: "⚠️",
        }
      );
    }
  }, [selectedUser]);

  if (loading.auth) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-100">
        <div className="text-center space-y-2">
          <div className="inline-block animate-spin rounded-full h-8 w-8 border-4 border-t-transparent border-blue-500"></div>
          <p className="text-gray-600">Verifying authentication...</p>
        </div>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen overflow-y-hidden bg-gray-100">
        <div className="text-center p-6 bg-white rounded-lg shadow-md max-w-sm space-y-4">
          <h2 className="text-2xl font-semibold text-gray-800">Access Restricted</h2>
          <p className="text-gray-600">Please sign in to use messaging</p>
          <button
            className="w-full px-4 py-2 bg-blue-500 text-white rounded-lg hover:bg-blue-600 transition-colors font-medium"
            onClick={() => (window.location.href = "/login")}
          >
            Continue to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="flex h-[90vh] bg-gray-100">
      <Toaster />
      <div className="w-1/4 min-w-[250px] bg-white border-r border-gray-200 flex flex-col">
        <div className="p-4 border-b border-gray-200 flex items-center justify-between">
          <h2 className="text-lg font-semibold text-gray-800">Messages</h2>
          <div className={`flex items-center ${socketConnected ? 'text-green-500' : 'text-red-500'}`}>
            <div className={`w-2 h-2 rounded-full mr-1 ${socketConnected ? 'bg-green-500' : 'bg-red-500'}`}></div>
            <span className="text-xs">{socketConnected ? 'Online' : 'Offline'}</span>
          </div>
        </div>
        <div className="flex-1 overflow-y-auto">
          {loading.users ? (
            <div className="p-4 space-y-3">
              {[...Array(5)].map((_, i) => (
                <div key={i} className="flex items-center space-x-3">
                  <div className="w-10 h-10 bg-gray-200 rounded-full animate-pulse"></div>
                  <div className="flex-1 space-y-1">
                    <div className="h-4 bg-gray-200 rounded w-3/4 animate-pulse"></div>
                    <div className="h-3 bg-gray-200 rounded w-1/2 animate-pulse"></div>
                  </div>
                </div>
              ))}
            </div>
          ) : users.length === 0 ? (
            <div className="flex flex-col items-center justify-center h-full p-4 text-center">
              <p className="text-gray-500">No users available</p>
              <p className="text-sm text-gray-400 mt-1">Users will appear here once they join</p>
            </div>
          ) : (
            users.map((u) => (
              <div
                key={u._id}
                onClick={() => setSelectedUser(u)}
                className={`flex items-center p-4 cursor-pointer hover:bg-gray-50 transition-colors ${
                  selectedUser?._id === u._id ? "bg-blue-50 border-l-2 border-blue-500" : ""
                }`}
              >
                <div className="w-10 h-10 rounded-full bg-gray-300 flex items-center justify-center text-white font-medium">
                  {u.profilePicture ? (
                    <img
                      src={u.profilePicture}
                      alt={u.username}
                      className="w-full h-full rounded-full object-cover"
                    />
                  ) : (
                    <span>{u.username[0]?.toUpperCase()}</span>
                  )}
                </div>
                <div className="ml-3 flex-1 overflow-hidden">
                  <p className="text-gray-800 font-medium truncate">{u.username}</p>
                  <p className="text-gray-500 text-sm truncate">Click to start chatting</p>
                </div>
                <div className="w-2 h-2 rounded-full bg-green-500 ml-2"></div>
              </div>
            ))
          )}
        </div>
      </div>

      <div className="flex-1 flex flex-col h-full">
        {selectedUser ? (
          <>
            <div className="p-4 bg-white border-b border-gray-200 flex items-center">
              <div className="w-10 h-10 rounded-full bg-gray-300 flex items-center justify-center text-white font-medium">
                {selectedUser.profilePicture ? (
                  <img
                    src={selectedUser.profilePicture}
                    alt={selectedUser.username}
                    className="w-full h-full rounded-full object-cover"
                  />
                ) : (
                  <span>{selectedUser.username[0]?.toUpperCase()}</span>
                )}
              </div>
              <div className="ml-3">
                <h2 className="text-lg font-semibold text-gray-800">{selectedUser.username}</h2>
                <p className="text-xs text-gray-500">Online</p>
              </div>
            </div>

            <div
              id="message-container"
              className="flex-1 p-4 overflow-y-auto bg-gray-50"
            >
              {loading.messages && messages.length === 0 ? (
                <div className="space-y-4">
                  {[...Array(3)].map((_, i) => (
                    <div key={i} className="animate-pulse flex space-x-4">
                      <div className="h-8 w-8 bg-gray-200 rounded-full"></div>
                      <div className="flex-1 space-y-2">
                        <div className="h-4 bg-gray-200 rounded w-1/4"></div>
                        <div className="h-4 bg-gray-200 rounded w-3/4"></div>
                      </div>
                    </div>
                  ))}
                </div>
              ) : error ? (
                <div className="p-4 bg-red-50 rounded-lg flex items-center gap-3 max-w-md mx-auto my-4">
                  <svg className="w-6 h-6 text-red-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth="2"
                      d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"
                    />
                  </svg>
                  <div>
                    <h3 className="font-medium text-red-700">Error</h3>
                    <p className="text-sm text-red-600">{error}</p>
                  </div>
                </div>
              ) : messages.length === 0 ? (
                <div className="flex flex-col items-center justify-center h-full">
                  <p className="text-gray-600 font-medium">
                    Start a conversation with {selectedUser.username}!
                  </p>
                  <p className="text-gray-500 text-sm mt-2">Your messages will appear here</p>
                </div>
              ) : (
                <div className="space-y-3">
                  {messages.map((msg, index) => {
                    const isCurrentUser = msg.sender._id === user._id;
                    const showAvatar = index === 0 || messages[index - 1].sender._id !== msg.sender._id;

                    return (
                      <div key={msg._id} className={`flex ${isCurrentUser ? "justify-end" : "justify-start"}`}>
                        {!isCurrentUser && showAvatar && (
                          <div className="w-8 h-8 rounded-full bg-gray-300 flex items-center justify-center text-white text-xs mr-2">
                            {selectedUser.username[0]?.toUpperCase()}
                          </div>
                        )}

                        <div
                          className={`max-w-xs md:max-w-md p-3 rounded-lg ${
                            isCurrentUser
                              ? "bg-blue-500 text-white rounded-tr-none"
                              : "bg-white text-gray-800 rounded-tl-none border border-gray-200"
                          }`}
                        >
                          <p className="whitespace-pre-wrap break-words">{msg.content}</p>
                          <p className={`text-xs mt-1 text-right ${isCurrentUser ? "text-blue-100" : "text-gray-500"}`}>
                            {new Date(msg.timestamp).toLocaleTimeString([], {
                              hour: "numeric",
                              minute: "2-digit",
                              hour12: true,
                            })}
                          </p>
                        </div>

                        {isCurrentUser && showAvatar && (
                          <div className="w-8 h-8 rounded-full bg-gray-300 flex items-center justify-center text-white text-xs ml-2">
                            {user.username[0]?.toUpperCase()}
                          </div>
                        )}
                      </div>
                    );
                  })}
                  <div ref={messageEndRef} />
                </div>
              )}
            </div>

            <div className="p-4 bg-white border-t border-gray-200">
              <form onSubmit={handleSendMessage} className="flex items-center gap-2">
                <input
                  type="text"
                  value={newMessage}
                  onChange={(e) => setNewMessage(e.target.value)}
                  placeholder="Type a message..."
                  className="flex-1 px-4 py-2 border border-gray-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
                  disabled={!socketConnected}
                />
                <button
                  type="submit"
                  className="p-2 bg-blue-500 text-white rounded-lg hover:bg-blue-600 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
                  disabled={!newMessage.trim() || !socketConnected}
                >
                  <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth="2"
                      d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8"
                    />
                  </svg>
                </button>
              </form>
            </div>
          </>
        ) : (
          <div className="flex flex-col items-center justify-center h-full bg-gray-50 p-6">
            <h3 className="text-lg font-semibold text-gray-800 mb-2">Your Messages</h3>
            <p className="text-gray-500 text-center">Select a user from the sidebar to start chatting</p>
          </div>
        )}
      </div>
    </div>
  );
}

export default Message;