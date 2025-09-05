import React, { useState, useEffect, useRef } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";
import { motion, AnimatePresence } from "framer-motion";
import "./name.css";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
});

const formatSmiles = (smiles) => {
  if (!smiles) return "";
  return smiles.match(/.{1,60}(?=\s|$)/g)?.join('\n') || "";
};

const truncateSmiles = (smiles, maxLength = 30) => {
  if (!smiles) return "";
  if (smiles.length <= maxLength) return smiles;
  return `${smiles.substring(0, maxLength)}...`;
};



const AINamingSuggestion = () => {
  const [activeTab, setActiveTab] = useState("generate");
  const [symptomGroupIndex, setSymptomGroupIndex] = useState("");
  const [selectedSmiles, setSelectedSmiles] = useState("");
  const [symptomGroups, setSymptomGroups] = useState([]);
  const [productSmilesGroups, setProductSmilesGroups] = useState([]);
  const [moleculeDetails, setMoleculeDetails] = useState({});
  const [suggestedNames, setSuggestedNames] = useState([]);
  const [savedNames, setSavedNames] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [copiedText, setCopiedText] = useState(null);
  const [showAcceptModal, setShowAcceptModal] = useState(false);
  const [selectedCandidate, setSelectedCandidate] = useState(null);
  const dropdownRef = useRef(null);
  const { user, checkAuth, checkingAuth, setCheckingAuth } = useAuthStore();
  const [isSmilesDropdownOpen, setIsSmilesDropdownOpen] = useState(false);
  const hasCheckedAuth = useRef(false); // Track if auth has been checked

  const toastOptions = {
    style: {
      background: '#172A45',
      color: '#E0E0E0',
      border: '1px solid #5E81F4',
      borderRadius: '8px',
      padding: '12px',
      fontFamily: 'Roboto, Open Sans, sans-serif',
    },
    success: { style: { borderColor: '#70E000' }, iconTheme: { primary: '#70E000', secondary: '#E0E0E0' } },
    error: { style: { borderColor: '#FF4C4C' } },
  };

  useEffect(() => {
    const initialize = async () => {
      if (hasCheckedAuth.current) return; // Skip if already checked
      hasCheckedAuth.current = true; // Mark as checked
      setCheckingAuth(true);
      try {
        const authTimeout = new Promise((_, reject) =>
          setTimeout(() => reject(new Error("Authentication check timed out")), 5000)
        );
        await Promise.race([checkAuth(), authTimeout]);
        console.log("Auth check completed, user:", user);
        if (!user) {
          setError("Authentication failed. Please log in.");
          return;
        }
        await fetchSymptomsAndProducts();
        await fetchSavedNames();
      } catch (err) {
        console.error("Authentication error:", err.message);
        setError("Authentication failed: " + err.message);
        toast.error("Authentication failed. Redirecting to login...", toastOptions);
        setTimeout(() => (window.location.href = "/login"), 2000);
      } finally {
        setCheckingAuth(false);
      }
    };
    initialize();
  }, [checkAuth, setCheckingAuth]); // Removed user from dependencies

  const fetchSymptomsAndProducts = async (targetOriginalSmiles = null, acceptedName = null) => {
    if (!user?._id) return;

    setLoading(true);
    try {
      const response = await axiosInstance.get(`/getdata/getsymptoms-product/${user._id}`);
      console.log("Symptoms and Products Response:", response.data);
      const { symptoms, productSmiles, moleculeDetails, message } = response.data;

      if (!symptoms || symptoms.length === 0) {
        throw new Error(message || "No symptoms returned from server");
      }

      setSymptomGroups(symptoms);
      setProductSmilesGroups(productSmiles || []);
      setMoleculeDetails(moleculeDetails || {});
      setError(null);

      // If a symptom group is not selected or was reset, default to the first group
      if (symptoms.length > 0 && (symptomGroupIndex === "" || symptomGroupIndex >= symptoms.length)) {
        setSymptomGroupIndex("0");
      }

      // If accepting a name, find and select the updated SMILES (drug name)
      if (targetOriginalSmiles && acceptedName) {
        let foundSmiles = null;
        productSmiles.forEach((group, index) => {
          const match = group.find(
            (smiles) =>
              moleculeDetails[smiles]?.originalSmiles === targetOriginalSmiles ||
              smiles === acceptedName
          );
          if (match) {
            foundSmiles = match;
            setSymptomGroupIndex(index.toString());
          }
        });
        if (foundSmiles) {
          setSelectedSmiles(foundSmiles);
        } else if (productSmiles[0]?.length > 0) {
          setSelectedSmiles(productSmiles[0][0]);
        }
      } else if (symptoms.length > 0 && productSmiles[symptomGroupIndex]?.length > 0) {
        // Preserve existing selection if possible
        setSelectedSmiles(productSmiles[symptomGroupIndex][0]);
      }
    } catch (err) {
      console.error("Error fetching symptoms and products:", err);
      setError(err.response?.data?.message || "Failed to fetch symptoms and products");
      setSymptomGroups([]);
      setProductSmilesGroups([]);
      setMoleculeDetails({});
      toast.error(err.response?.data?.message || "Failed to fetch data", toastOptions);
    } finally {
      setLoading(false);
    }
  };

  const fetchSavedNames = async () => {
    if (!user?._id) return;

    try {
      const response = await axiosInstance.get("/drugname/saved-drug-names");
      const drugNames = response.data.drugNames || [];
      console.log("Saved Drug Names:", drugNames);
      setSavedNames(drugNames.map((drugName) => ({
        ...drugName,
        symptomsGrp: drugName.symptoms || "N/A",
      })));
    } catch (err) {
      console.error("Error fetching saved drug names:", err);
      setError(err.response?.data?.message || "Failed to fetch saved names");
      toast.error("Failed to fetch saved names", toastOptions);
    }
  };

  const checkIfNameExists = async (smiles) => {
    try {
      const response = await axiosInstance.get("/drugname/check-saved-drug-name", {
        params: { smiles },
      });
      return response.data.exists;
    } catch (err) {
      console.error("Error checking saved name:", err);
      return false;
    }
  };

  const handleGenerateName = async () => {
    if (symptomGroupIndex === "" || !selectedSmiles) {
      toast.error("Please select both a symptom group and SMILES string", toastOptions);
      return;
    }

    const symptoms = symptomGroups[symptomGroupIndex]?.join(", ") || "";
    const originalSmiles = moleculeDetails[selectedSmiles]?.originalSmiles || selectedSmiles;
    const nameExists = await checkIfNameExists(originalSmiles);
    if (nameExists) {
      const savedName = savedNames.find((n) => n.smiles === originalSmiles);
      if (savedName?.status === "accepted") {
        toast("An accepted drug name already exists for this SMILES.", {
          ...toastOptions,
          type: "info",
        });
        setActiveTab("saved");
        return;
      }
    }

    setLoading(true);
    setError(null);
    setSuggestedNames([]);

    try {
      const response = await axiosInstance.post(`/drugname/generate-drug-name/${user._id}`, {
        smiles: originalSmiles,
        symptoms,
      });

      if (response.status === 409) {
        toast("An accepted drug name already exists. Redirecting to Saved Names.", {
          ...toastOptions,
          type: "info",
        });
        setActiveTab("saved");
        return;
      }

      setSuggestedNames(response.data.allCandidates);
      toast.success("Drug names generated successfully!", toastOptions);
    } catch (err) {
      console.error("Error generating drug name:", err);
      setError(err.response?.data?.message || "Failed to generate drug name");
      toast.error("Failed to generate drug name", toastOptions);
    } finally {
      setLoading(false);
    }
  };

  const handleAcceptName = async (candidate) => {
    setSelectedCandidate(candidate);
    setShowAcceptModal(true);
  };

  const confirmAcceptName = async () => {
    if (!selectedCandidate) return;

    const symptoms = symptomGroups[symptomGroupIndex]?.join(", ") || "";
    const originalSmiles = moleculeDetails[selectedSmiles]?.originalSmiles || selectedSmiles;
    setLoading(true);
    setShowAcceptModal(false);

    try {
      const response = await axiosInstance.post(`/drugname/accept-drug-name/${user._id}`, {
        smiles: originalSmiles,
        symptoms,
        selectedName: selectedCandidate.name,
        rationale: selectedCandidate.rationale,
        compliance: selectedCandidate.compliance,
      });

      setSuggestedNames([]);
      // Do not reset symptomGroupIndex: setSymptomGroupIndex("");
      // Do not reset selectedSmiles yet: setSelectedSmiles("");
      await fetchSymptomsAndProducts(originalSmiles, selectedCandidate.name); // Pass params to select updated SMILES
      await fetchSavedNames();
      toast.success("Drug name accepted successfully!", toastOptions);
    } catch (err) {
      console.error("Error accepting drug name:", err);
      const errorMessage = err.response?.data?.message || "Failed to accept drug name";
      setError(errorMessage);
      toast.error(errorMessage, toastOptions);
    } finally {
      setLoading(false);
      setSelectedCandidate(null);
    }
  };

  const handleRejectName = () => {
    setSuggestedNames([]);
  };

  const handleTabChange = async (tab) => {
    if (activeTab === "generate" && suggestedNames.length > 0) {
      try {
        const symptoms = symptomGroups[symptomGroupIndex]?.join(", ") || "";
        const originalSmiles = moleculeDetails[selectedSmiles]?.originalSmiles || selectedSmiles;
        await axiosInstance.post(`/drugname/save-pending-drug-name/${user._id}`, {
          smiles: originalSmiles,
          symptoms,
          candidates: suggestedNames,
        });
      } catch (err) {
        console.error("Error saving pending drug name:", err);
        toast.error("Failed to save pending drug name", toastOptions);
      }
    }
    setActiveTab(tab);
    setSuggestedNames([]);
    setError(null);
  };

  const copyToClipboard = (text) => {
    navigator.clipboard.writeText(text);
    setCopiedText(text);
    setTimeout(() => setCopiedText(null), 2000);
  };

  const toggleSmilesDropdown = () => {
    if (!loading && symptomGroupIndex !== "" && productSmilesGroups[symptomGroupIndex]?.length > 0) {
      setIsSmilesDropdownOpen(!isSmilesDropdownOpen);
    }
  };

  const handleSmilesSelect = (smiles) => {
    setSelectedSmiles(smiles);
    setIsSmilesDropdownOpen(false);
  };

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-primary">
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          className="text-center py-10 text-accent"
        >
          <p>Verifying authentication...</p>
          <p className="text-sm text-text-secondary mt-2">Please wait a moment.</p>
        </motion.div>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-primary">
        <motion.div
          initial={{ scale: 0.9, opacity: 0 }}
          animate={{ scale: 1, opacity: 1 }}
          className="text-center p-6 bg-secondary rounded-lg shadow-lg max-w-md w-full"
        >
          <p className="text-text-secondary mb-4">Please log in to access AI Naming Suggestion</p>
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            className="px-4 py-2 bg-accent text-primary rounded-lg hover:bg-accent transition-colors w-full sm:w-auto font-medium"
            onClick={() => (window.location.href = "/login")}
          >
            Go to Login
          </motion.button>
        </motion.div>
      </div>
    );
  }

  return (
    <div className="min-h-screen py-6 px-4 sm:py-8 sm:px-6 lg:py-12 lg:px-8 bg-primary">
      <div className="max-w-6xl mx-auto">
        <motion.div
          initial={{ y: -20, opacity: 0 }}
          animate={{ y: 0, opacity: 1 }}
          transition={{ duration: 0.5 }}
          className="text-center mb-8 sm:mb-12"
        >
          <h1 className="text-3xl sm:text-4xl font-bold text-accent mb-2">
            AI Drug Naming Suggestion
          </h1>
          <p className="text-sm text-accent-secondary font-mono">(POWERED BY GEMINI)</p>
        </motion.div>

        <div className="flex flex-col sm:flex-row justify-center mb-6 sm:mb-10 space-y-2 sm:space-y-0 sm:space-x-4">
          {["generate", "saved"].map((tab) => (
            <motion.button
              key={tab}
              whileHover={{ scale: 1.05 }}
              whileTap={{ scale: 0.95 }}
              className={`px-6 py-3 rounded-lg font-medium transition-all duration-300 ${activeTab === tab
                  ? "bg-accent text-primary shadow-lg"
                  : "bg-secondary text-text-primary hover:bg-opacity-80"
                }`}
              onClick={() => handleTabChange(tab)}
            >
              {tab === "generate" && "Generate Drug Name"}
              {tab === "saved" && "Saved Names"}
            </motion.button>
          ))}
        </div>

        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ delay: 0.2 }}
          className="bg-secondary p-6 sm:p-8 rounded-xl shadow-lg border border-gray-700"
        >
          {activeTab === "generate" && (
            <>
              <h2 className="text-xl sm:text-2xl font-semibold text-accent mb-6">Generate Drug Name</h2>

              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
                <div className="space-y-4">
                  <div>
                    <label className="block text-sm font-medium text-text-secondary mb-2">
                      Select Symptom Group
                    </label>
                    <select
                      value={symptomGroupIndex}
                      onChange={(e) => {
                        const index = e.target.value;
                        setSymptomGroupIndex(index);
                        setSelectedSmiles(productSmilesGroups[index]?.[0] || "");
                        console.log("Selected symptom group:", symptomGroups[index]);
                      }}
                      className="w-full p-3 bg-primary border border-gray-600 rounded-lg focus:outline-none focus:ring-2 focus:ring-accent text-text-primary font-mono"
                      disabled={loading || symptomGroups.length === 0}
                    >
                      {symptomGroups.length === 0 ? (
                        <option value="">No symptom groups available</option>
                      ) : (
                        symptomGroups.map((group, index) => (
                          <option key={index} value={index}>
                            {group.join(", ") || "Unnamed Group"}
                          </option>
                        ))
                      )}
                    </select>
                  </div>

                  <div className="relative" ref={dropdownRef}>
                    <label className="block text-sm font-medium text-text-secondary mb-2">
                      Select Drug Name / SMILES
                    </label>
                    <div
                      className={`w-full p-3 bg-primary border border-gray-600 rounded-lg text-text-primary font-mono flex justify-between items-center cursor-pointer ${loading || symptomGroupIndex === "" || productSmilesGroups[symptomGroupIndex]?.length === 0
                          ? "opacity-50 cursor-not-allowed"
                          : "hover:bg-gray-700"
                        }`}
                      onClick={toggleSmilesDropdown}
                    >
                      <span className="truncate">
                        {selectedSmiles ? truncateSmiles(selectedSmiles) : "Select a drug name / SMILES"}
                      </span>
                      <svg
                        className={`w-4 h-4 transform transition-transform ${isSmilesDropdownOpen ? "rotate-180" : ""}`}
                        fill="none"
                        stroke="currentColor"
                        viewBox="0 0 24 24"
                        xmlns="http://www.w3.org/2000/svg"
                      >
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M19 9l-7 7-7-7" />
                      </svg>
                    </div>
                    {isSmilesDropdownOpen && (
                      <div className="absolute z-10 w-full mt-1 bg-primary border border-gray-600 rounded-lg shadow-lg max-h-60 overflow-y-auto">
                        {productSmilesGroups[symptomGroupIndex]?.length === 0 ? (
                          <div className="p-3 text-text-secondary">No entries available</div>
                        ) : (
                          productSmilesGroups[symptomGroupIndex]?.map((smiles, index) => (
                            <div
                              key={index}
                              className="p-3 text-text-primary font-mono hover:bg-gray-700 cursor-pointer group relative"
                              onClick={() => handleSmilesSelect(smiles)}
                            >
                              <span className="truncate block">{truncateSmiles(smiles)}</span>
                              {moleculeDetails[smiles]?.originalSmiles && (
                                <div className="absolute left-0 top-full mt-1 hidden group-hover:block bg-gray-800 text-white text-xs rounded p-2 z-20 max-w-md whitespace-pre-wrap">
                                  Original SMILES: {moleculeDetails[smiles].originalSmiles}
                                </div>
                              )}
                            </div>
                          ))
                        )}
                      </div>
                    )}
                  </div>

                  <motion.button
                    whileHover={{ scale: 1.02 }}
                    whileTap={{ scale: 0.98 }}
                    onClick={handleGenerateName}
                    disabled={loading || symptomGroupIndex === "" || !selectedSmiles}
                    className={`w-full py-3 px-6 rounded-lg transition-all duration-300 ${loading || symptomGroupIndex === "" || !selectedSmiles
                        ? "bg-gray-600 cursor-not-allowed"
                        : "bg-accent hover:bg-accent-secondary text-primary"
                      } font-medium`}
                  >
                    {loading ? (
                      <span className="flex items-center justify-center">
                        <svg
                          className="animate-spin -ml-1 mr-3 h-5 w-5 text-white"
                          xmlns="http://www.w3.org/2000/svg"
                          fill="none"
                          viewBox="0 0 24 24"
                        >
                          <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                          <path
                            className="opacity-75"
                            fill="currentColor"
                            d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                          ></path>
                        </svg>
                        Generating...
                      </span>
                    ) : (
                      "Start Prediction"
                    )}
                  </motion.button>
                </div>

                {selectedSmiles && (
                  <div className="bg-primary p-4 rounded-lg border border-gray-600">
                    <div className="flex justify-between items-center mb-3">
                      <h3 className="text-lg font-medium text-accent">Selected Drug Name / SMILES</h3>
                      <motion.button
                        whileHover={{ scale: 1.05 }}
                        whileTap={{ scale: 0.95 }}
                        onClick={() => copyToClipboard(moleculeDetails[selectedSmiles]?.originalSmiles || selectedSmiles)}
                        className="text-xs bg-accent-secondary text-primary px-2 py-1 rounded"
                      >
                        {copiedText === (moleculeDetails[selectedSmiles]?.originalSmiles || selectedSmiles) ? "Copied!" : "Copy"}
                      </motion.button>
                    </div>
                    <div className="bg-gray-900 p-3 rounded font-mono text-sm text-green-400 whitespace-pre-wrap overflow-y-auto max-h-60">
                      {moleculeDetails[selectedSmiles]?.originalSmiles ? (
                        <>
                          <div>Drug Name: {selectedSmiles}</div>
                          <div>Original SMILES: {moleculeDetails[selectedSmiles].originalSmiles}</div>
                        </>
                      ) : (
                        selectedSmiles
                      )}
                    </div>
                  </div>
                )}
              </div>

              {error && (
                <motion.div
                  initial={{ opacity: 0, y: -10 }}
                  animate={{ opacity: 1, y: 0 }}
                  className="mt-4 bg-error bg-opacity-20 border border-error rounded-lg p-4 flex justify-between items-center"
                >
                  <p className="text-text-primary">{error}</p>
                  <button
                    className="text-error underline hover:text-opacity-80"
                    onClick={() => setError(null)}
                  >
                    Dismiss
                  </button>
                </motion.div>
              )}

              <AnimatePresence>
                {suggestedNames.length > 0 && (
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -20 }}
                    className="mt-8 bg-primary p-6 rounded-xl border border-gray-600"
                  >
                    <h3 className="text-xl font-semibold text-accent mb-6">Suggested Drug Names</h3>
                    <div className="space-y-6">
                      {suggestedNames.map((candidate) => (
                        <motion.div
                          key={candidate.rank}
                          initial={{ opacity: 0 }}
                          animate={{ opacity: 1 }}
                          transition={{ delay: candidate.rank * 0.1 }}
                          className="border-b border-gray-600 pb-6 last:border-b-0"
                        >
                          <div className="flex items-center mb-3">
                            <span className="text-sm font-medium text-text-secondary mr-3">
                              Rank {candidate.rank}:
                            </span>
                            <h4 className="text-xl font-bold text-accent">{candidate.name}</h4>
                          </div>
                          <div className="mb-3">
                            <p className="text-sm text-text-secondary mb-1">Structural Rationale</p>
                            <p className="text-text-primary bg-gray-800 p-3 rounded font-mono text-sm">{candidate.rationale}</p>
                          </div>
                          <div className="flex space-x-4">
                            <motion.button
                              whileHover={{ scale: 1.05 }}
                              whileTap={{ scale: 0.95 }}
                              onClick={() => handleAcceptName(candidate)}
                              className="px-4 py-2 bg-success text-primary rounded-lg hover:bg-opacity-90 transition-colors font-medium"
                              disabled={loading}
                            >
                              Accept
                            </motion.button>
                            {/* <motion.button
                              whileHover={{ scale: 1.05 }}
                              whileTap={{ scale: 0.95 }}
                              onClick={async () => {
                                try {
                                  await axiosInstance.delete(`/drugname/delete-drug-name/${drugName._id}`);
                                  await Promise.all([fetchSavedNames(), fetchSymptomsAndProducts()]); // Refresh both
                                  toast.success("Drug name deleted successfully!", toastOptions);
                                } catch (err) {
                                  console.error("Error deleting drug name:", err);
                                  toast.error("Failed to delete drug name", toastOptions);
                                }
                              }}
                              className="px-4 py-2 bg-error text-primary rounded-lg hover:bg-opacity-90 transition-colors font-medium"
                              disabled={loading}
                            >
                              Reject
                            </motion.button> */}
                          </div>
                        </motion.div>
                      ))}
                    </div>
                  </motion.div>
                )}
              </AnimatePresence>
            </>
          )}

          {activeTab === "saved" && (
            <>
              <h2 className="text-xl sm:text-2xl font-semibold text-accent mb-6">Saved Drug Names</h2>
              {savedNames.length > 0 ? (
                <div className="space-y-6">
                  <AnimatePresence>
                    {savedNames.map((drugName) => (
                      <motion.div
                        key={drugName._id}
                        initial={{ opacity: 0, y: 20 }}
                        animate={{ opacity: 1, y: 0 }}
                        exit={{ opacity: 0, x: -20 }}
                        className="bg-primary p-6 rounded-xl border border-gray-600 transition-all duration-200 hover:shadow-lg"
                      >
                        <div className="flex justify-between items-start mb-4">
                          <div>
                            <h3 className="text-lg font-semibold text-accent mb-1">
                              {drugName.suggestedName}
                            </h3>
                            <span
                              className={`text-xs px-2 py-1 rounded ${drugName.status === "accepted"
                                  ? "bg-success bg-opacity-20 text-success"
                                  : "bg-accent-secondary bg-opacity-20 text-accent-secondary"
                                }`}
                            >
                              {drugName.status.charAt(0).toUpperCase() + drugName.status.slice(1)}
                            </span>
                          </div>
                          <span className="text-xs text-text-secondary">
                            {new Date(drugName.createdAt).toLocaleDateString()}
                          </span>
                        </div>
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
                          <div>
                            <p className="text-sm text-text-secondary mb-1">Symptom Group</p>
                            <p className="text-text-primary">{drugName.symptomsGrp}</p>
                          </div>
                          <div>
                            <p className="text-sm text-text-secondary mb-1">SMILES</p>
                            <div className="flex items-center">
                              <div className="bg-gray-900 p-2 rounded font-mono text-xs text-green-400 overflow-x-auto max-w-xs">
                                {truncateSmiles(drugName.smiles)}
                              </div>
                              <motion.button
                                whileHover={{ scale: 1.1 }}
                                whileTap={{ scale: 0.9 }}
                                onClick={() => copyToClipboard(drugName.smiles)}
                                className="ml-2 text-xs bg-accent-secondary text-primary px-2 py-1 rounded"
                              >
                                {copiedText === drugName.smiles ? "Copied!" : "Copy"}
                              </motion.button>
                            </div>
                          </div>
                        </div>
                        <div className="mb-4">
                          <p className="text-sm text-text-secondary mb-1">Naming Details</p>
                          <p className="text-text-primary bg-gray-800 p-3 rounded font-mono text-sm whitespace-pre-wrap">
                            {drugName.namingDetails}
                          </p>
                        </div>
                        {drugName.status === "pending" && (
                          <div className="flex space-x-4">
                            <motion.button
                              whileHover={{ scale: 1.05 }}
                              whileTap={{ scale: 0.95 }}
                              onClick={() =>
                                handleAcceptName({
                                  name: drugName.suggestedName,
                                  rationale: drugName.namingDetails.split(" | Compliance: ")[0].replace("Rationale: ", ""),
                                  compliance: drugName.namingDetails.split(" | Compliance: ")[1],
                                })
                              }
                              className="px-4 py-2 bg-success text-primary rounded-lg hover:bg-opacity-90 transition-colors font-medium"
                              disabled={loading}
                            >
                              Accept
                            </motion.button>
                            <motion.button
                              whileHover={{ scale: 1.05 }}
                              whileTap={{ scale: 0.95 }}
                              onClick={async () => {
                                await axiosInstance.delete(`/drugname/delete-drug-name/${drugName._id}`);
                                await fetchSavedNames();
                              }}
                              className="px-4 py-2 bg-error text-primary rounded-lg hover:bg-opacity-90 transition-colors font-medium"
                              disabled={loading}
                            >
                              Reject
                            </motion.button>
                          </div>
                        )}
                      </motion.div>
                    ))}
                  </AnimatePresence>
                </div>
              ) : (
                <motion.div
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  className="text-center py-10"
                >
                  <p className="text-text-secondary">No saved drug names found.</p>
                </motion.div>
              )}
            </>
          )}
        </motion.div>

        <AnimatePresence>
          {showAcceptModal && selectedCandidate && (
            <motion.div
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              exit={{ opacity: 0 }}
              className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50"
            >
              <motion.div
                initial={{ scale: 0.9, y: 20 }}
                animate={{ scale: 1, y: 0 }}
                exit={{ scale: 0.9, y: 20 }}
                className="bg-secondary p-6 rounded-xl shadow-lg max-w-md w-full border border-gray-700"
              >
                <h3 className="text-xl font-semibold text-accent mb-4">Confirm Acceptance</h3>
                <p className="text-text-primary mb-6">
                  By accepting "<span className="font-bold">{selectedCandidate.name}</span>", this name will become the final title for this SMILES across the database. Proceed?
                </p>
                <div className="flex justify-end space-x-4">
                  <motion.button
                    whileHover={{ scale: 1.05 }}
                    whileTap={{ scale: 0.95 }}
                    onClick={() => {
                      setShowAcceptModal(false);
                      setSelectedCandidate(null);
                    }}
                    className="px-4 py-2 bg-error text-primary rounded-lg hover:bg-opacity-90 transition-colors font-medium"
                  >
                    Cancel
                  </motion.button>
                  <motion.button
                    whileHover={{ scale: 1.05 }}
                    whileTap={{ scale: 0.95 }}
                    onClick={confirmAcceptName}
                    className="px-4 py-2 bg-success text-primary rounded-lg hover:bg-opacity-90 transition-colors font-medium"
                  >
                    Confirm
                  </motion.button>
                </div>
              </motion.div>
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </div>
  );
};

export default AINamingSuggestion;