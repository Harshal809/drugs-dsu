"use client";

import { useState, useEffect, useRef } from "react";
import { useAuthStore } from "../../Store/auth.store.js";
import { postCostEstimation, getCostEstimations } from "../../api/costestimination.jsx";
import axios from "axios";
import { AlertCircle, DollarSign, Clock, Database, X, Info, RefreshCw, LogIn, FileText, ChevronDown, ChevronUp, Download, FlaskConical, Pill, TestTube2, Syringe, Atom, Leaf } from "lucide-react";
import { jsPDF } from "jspdf";
import domtoimage from "dom-to-image";
import Loader from "../../components/Loader.jsx";
import { motion, AnimatePresence } from "framer-motion";
import { toast } from "react-hot-toast";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
});


const floatingIcons = [
  { icon: FlaskConical, size: 24, delay: 0, duration: 5, color: "text-accent/20" },
  { icon: Pill, size: 20, delay: 0.5, duration: 6, color: "text-accent-secondary/20" },
  { icon: TestTube2, size: 22, delay: 0.8, duration: 4.5, color: "text-success/20" },
  { icon: Syringe, size: 18, delay: 1.2, duration: 5.5, color: "text-error/20" },
  { icon: Atom, size: 26, delay: 0.3, duration: 6.5, color: "text-accent/20" },
  { icon: Leaf, size: 20, delay: 0.7, duration: 5, color: "text-success/20" },
];

const CostEstimationForm = () => {
  const [symptomGroupIndex, setSymptomGroupIndex] = useState("");
  const [selectedSmiles, setSelectedSmiles] = useState("");
  const [symptomGroups, setSymptomGroups] = useState([]);
  const [productSmilesGroups, setProductSmilesGroups] = useState([]);
  const [result, setResult] = useState(null);
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [historyLoading, setHistoryLoading] = useState(false);
  const [isResultInfoOpen, setIsResultInfoOpen] = useState(false);
  const [openHistoryItems, setOpenHistoryItems] = useState({});

  const { user, checkAuth, checkingAuth } = useAuthStore();
  const infoRef = useRef(null);
  const userId = user?._id;

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

  const fetchSymptomsAndProducts = async () => {
    if (!userId) return;

    try {
      const response = await axiosInstance.get(`/getdata/getsymptoms-product/${userId}`);
      const { symptoms, productSmiles } = response.data;
      setSymptomGroups(symptoms || []);
      setProductSmilesGroups(productSmiles || []);
      if (symptoms?.length > 0) {
        setSymptomGroupIndex("0");
        if (productSmiles?.[0]?.length > 0) {
          setSelectedSmiles(productSmiles[0][0]);
        }
      }
    } catch (err) {
      console.error("Error fetching symptoms and products:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch symptoms and products");
      setSymptomGroups([]);
      setProductSmilesGroups([]);
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!userId) {
      setError("Please log in to estimate costs");
      return;
    }
    if (symptomGroupIndex === "") {
      setError("Please select a symptom group");
      return;
    }
    if (!selectedSmiles) {
      setError("Please select a SMILES string");
      return;
    }
    setLoading(true);
    setError(null);
    try {
      const symptomsGrp = symptomGroups[symptomGroupIndex]?.join(", ") || "N/A";
      const data = await postCostEstimation(selectedSmiles);
      setResult(data.data);

      const estimationId = data.data._id;
      // Debug: Log the values being stored
      console.log("Storing in localStorage - userId:", userId, "estimationId:", estimationId, "symptomsGrp:", symptomsGrp);
      
      // Store the symptom group in localStorage
      const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${userId}`) || "{}");
      storedSymptomGroups[estimationId] = symptomsGrp;
      localStorage.setItem(`symptomGroups_${userId}`, JSON.stringify(storedSymptomGroups));
      
      // Debug: Verify the stored data
      console.log("Stored symptom groups in localStorage:", JSON.parse(localStorage.getItem(`symptomGroups_${userId}`)));

      toast.success("Estimating the cost", toastOptions);
      await fetchHistory();
    } catch (error) {
      setError("Failed to estimate cost. Please try again.");
      console.error("Error estimating cost:", error);
    } finally {
      setLoading(false);
    }
  };

  const fetchHistory = async () => {
    if (!user?._id) {
      setError("Please log in to view history");
      return;
    }
    setHistoryLoading(true);
    try {
      const data = await getCostEstimations(user._id);
      const validHistory = Array.isArray(data.data) ? data.data.filter((item) => item && typeof item === "object") : [];

      // Debug: Log the userId and history data
      console.log("Fetching history for userId:", user._id);
      console.log("Raw history data:", validHistory);

      // Retrieve symptom groups from localStorage
      const storedSymptomGroupsRaw = localStorage.getItem(`symptomGroups_${user._id}`);
      console.log("Raw symptom groups from localStorage:", storedSymptomGroupsRaw);
      
      const storedSymptomGroups = storedSymptomGroupsRaw ? JSON.parse(storedSymptomGroupsRaw) : {};
      console.log("Parsed symptom groups from localStorage:", storedSymptomGroups);

      // Map through history and attach the symptom group
      const updatedHistory = validHistory.map((item) => {
        const symptomsGrp = storedSymptomGroups[item._id] || "N/A";
        console.log(`History item ID: ${item._id}, Symptom Group: ${symptomsGrp}`);
        return {
          ...item,
          symptomsGrp,
        };
      });

      setHistory(updatedHistory);
    } catch (error) {
      setError("No previous estimations found.");
      console.error("Error fetching history:", error);
    } finally {
      setHistoryLoading(false);
    }
  };

  const toggleResultInfo = () => {
    setIsResultInfoOpen(!isResultInfoOpen);
  };

  const toggleHistoryItem = (id) => {
    setOpenHistoryItems((prev) => ({
      ...prev,
      [id]: !prev[id],
    }));
  };

  const renderInformation = (info) => {
    if (!info) return <p className="text-text-secondary">No additional information available</p>;

    const lines = info.split("\n").filter((line) => line.trim() !== "");
    let introParagraph = "";
    const sections = [];
    let currentSection = null;

    lines.forEach((line, index) => {
      if (index === 0 && !line.match(/^\d+\./)) {
        introParagraph = line.trim();
      } else if (line.match(/^\d+\.\s/)) {
        if (currentSection) sections.push(currentSection);
        currentSection = { title: line.trim(), bullets: [] };
      } else if (line.trim().startsWith("-")) {
        if (currentSection) currentSection.bullets.push(line.trim().replace("-", "").trim());
      } else if (currentSection) {
        currentSection.bullets.push(line.trim());
      }
    });
    if (currentSection) sections.push(currentSection);

    return (
      <div className="font-body">
        {introParagraph && <p className="mb-2 text-text-primary">{introParagraph}</p>}
        {sections.length > 0 && (
          <ul className="list-disc pl-5 space-y-2">
            {sections.map((section, index) => (
              <li key={index} className="font-semibold text-accent">
                {section.title}
                {section.bullets.length > 0 && (
                  <ul className="list-disc pl-5 mt-1 space-y-1">
                    {section.bullets.map((bullet, bulletIndex) => (
                      <li key={bulletIndex} className="text-text-primary">{bullet}</li>
                    ))}
                  </ul>
                )}
              </li>
            ))}
          </ul>
        )}
      </div>
    );
  };

  const exportToPDF = async () => {
    const pdf = new jsPDF("p", "mm", "a4");
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const margin = 10;
    let yPosition = margin;

    pdf.setFontSize(16);
    pdf.text("Drug Cost Estimation Report", margin, yPosition);
    yPosition += 10;

    pdf.setFontSize(12);
    pdf.text("Symptom Group:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(symptomGroups[symptomGroupIndex]?.join(", ") || "N/A", margin, yPosition);
    yPosition += 8;

    pdf.setFontSize(12);
    pdf.text("SMILES:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(result?.smiles || selectedSmiles || "N/A", margin, yPosition);
    yPosition += 8;

    pdf.setFontSize(12);
    pdf.text("Estimated Cost:", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.text(result?.estimatedcost || "N/A", margin, yPosition);
    yPosition += 10;

    if (infoRef.current) {
      pdf.setFontSize(12);
      pdf.text("Information:", margin, yPosition);
      yPosition += 8;

      const infoImgData = await domtoimage.toPng(infoRef.current, { quality: 1 });
      const infoImgProps = pdf.getImageProperties(infoImgData);
      const infoImgWidth = pageWidth - 2 * margin;
      let infoImgHeight = (infoImgProps.height * infoImgWidth) / infoImgProps.width;

      let remainingHeight = infoImgHeight;
      let yOffset = 0;

      while (remainingHeight > 0) {
        const spaceLeft = pageHeight - yPosition - margin;
        const heightToRender = Math.min(remainingHeight, spaceLeft);

        const tempCanvas = document.createElement("canvas");
        const tempCtx = tempCanvas.getContext("2d");
        const img = new Image();
        img.src = infoImgData;
        await new Promise((resolve) => {
          img.onload = resolve;
        });
        tempCanvas.width = img.width;
        tempCanvas.height = (heightToRender / infoImgHeight) * img.height;
        tempCtx.drawImage(img, 0, yOffset, img.width, tempCanvas.height, 0, 0, img.width, tempCanvas.height);
        const croppedImgData = tempCanvas.toDataURL("image/png");

        pdf.addImage(croppedImgData, "PNG", margin, yPosition, infoImgWidth, heightToRender);
        yPosition += heightToRender;
        remainingHeight -= heightToRender;
        yOffset += tempCanvas.height;

        if (remainingHeight > 0) {
          pdf.addPage();
          yPosition = margin;
        }
      }
    }

    pdf.save("drug-cost-estimation-report.pdf");
  };

  useEffect(() => {
    const initialize = async () => {
      setLoading(true);
      try {
        const checkedUser = await checkAuth();
        if (!checkedUser || !useAuthStore.getState().user?._id) {
          setError("Authentication failed. Please log in.");
          setLoading(false);
          return;
        }

        await fetchSymptomsAndProducts();
        await fetchHistory();
      } catch (err) {
        console.error("Initialization error:", err);
        setError("Failed to verify authentication. Please try refreshing the page or logging in again.");
      } finally {
        setLoading(false);
      }
    };

    initialize();
  }, [checkAuth]);

  if (!user || !user._id) {
    return (
      <div className="flex flex-col items-center justify-center h-screen bg-primary relative overflow-hidden">
        {/* Floating drug discovery icons */}
        {floatingIcons.map((iconData, index) => (
          <motion.div
            key={index}
            initial={{ y: 0, opacity: 0 }}
            animate={{
              y: [0, -50, 0],
              opacity: [0, 1, 0],
            }}
            transition={{
              duration: iconData.duration,
              delay: iconData.delay,
              repeat: Infinity,
              repeatType: "reverse",
              ease: "easeInOut",
            }}
            className={`absolute ${iconData.color}`}
            style={{
              left: `${Math.random() * 90 + 5}%`,
              top: `${Math.random() * 80 + 10}%`,
            }}
          >
            <iconData.icon size={iconData.size} />
          </motion.div>
        ))}

        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          className="text-center max-w-md p-8 bg-secondary rounded-xl shadow-lg border border-secondary relative z-10"
        >
          <motion.div
            animate={{ rotate: 360 }}
            transition={{ duration: 8, repeat: Infinity, ease: "linear" }}
            className="absolute -top-16 -right-16 opacity-10"
          >
            <Atom size={120} />
          </motion.div>
          <LogIn size={48} className="mx-auto text-accent mb-4" />
          <h2 className="text-2xl font-bold text-text-primary mb-3 font-heading">Authentication Required</h2>
          <p className="text-text-secondary mb-6 font-body">
            Please log in to access the Drug Cost Estimator tool and view your estimation history.
          </p>
          <motion.button
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
            className="w-full px-4 py-3 bg-accent text-primary font-medium rounded-lg hover:bg-accent/90 transition-colors duration-200 flex items-center justify-center font-heading"
            onClick={() => (window.location.href = "/login")}
          >
            <LogIn size={18} className="mr-2" />
            Go to Login
          </motion.button>
        </motion.div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-primary relative overflow-hidden">
      {/* Floating drug discovery icons */}
      {floatingIcons.map((iconData, index) => (
        <motion.div
          key={index}
          initial={{ y: 0, opacity: 0 }}
          animate={{
            y: [0, -50, 0],
            opacity: [0, 1, 0],
          }}
          transition={{
            duration: iconData.duration,
            delay: iconData.delay,
            repeat: Infinity,
            repeatType: "reverse",
            ease: "easeInOut",
          }}
          className={`absolute ${iconData.color}`}
          style={{
            left: `${Math.random() * 90 + 5}%`,
            top: `${Math.random() * 80 + 10}%`,
          }}
        >
          <iconData.icon size={iconData.size} />
        </motion.div>
      ))}

      <div className="container mx-auto px-4 max-w-6xl relative z-10">
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          className="mb-12 text-center"
        >
          <h1 className="text-4xl font-bold text-text-primary mb-2 font-heading tracking-tight">
            Drug Cost Estimation
            <p className="text-sm text-accent-secondary font-300 font-mono">(POWERED BY GEMINI)</p>
          </h1>
          <p className="text-text-secondary max-w-2xl mx-auto font-body text-lg">
            Select a symptom group and a corresponding SMILES string to estimate the cost of drug synthesis and production.
          </p>
        </motion.div>

        <div className="space-y-10">
          {/* Form Section */}
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.5, delay: 0.2 }}
            className="bg-secondary p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg border border-secondary"
          >
            <div className="flex items-center mb-6">
              <motion.div
                animate={{ y: [-5, 5, -5] }}
                transition={{ duration: 4, repeat: Infinity, ease: "easeInOut" }}
              >
                <FlaskConical className="h-7 w-7 text-accent mr-2" />
              </motion.div>
              <h2 className="text-2xl font-semibold text-text-primary font-heading">Estimate New Cost</h2>
            </div>

            <form onSubmit={handleSubmit} className="space-y-5">
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ duration: 0.5, delay: 0.3 }}
              >
                <label htmlFor="symptomGroup" className="block text-sm font-medium text-text-secondary mb-1 font-body">
                  Select Symptom Group
                </label>
                <select
                  id="symptomGroup"
                  value={symptomGroupIndex}
                  onChange={(e) => {
                    setSymptomGroupIndex(e.target.value);
                    setSelectedSmiles(productSmilesGroups[e.target.value]?.[0] || "");
                  }}
                  className="w-full px-4 py-3 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent transition-all duration-200 font-body"
                  disabled={loading || symptomGroups.length === 0}
                >
                  {symptomGroups.length === 0 ? (
                    <option value="" className="bg-primary">No symptom groups available</option>
                  ) : (
                    symptomGroups.map((group, index) => (
                      <option key={index} value={index} className="bg-primary">
                        {group.join(", ")}
                      </option>
                    ))
                  )}
                </select>
                <p className="mt-1 text-xs text-text-secondary font-body">
                  Select a group of symptoms
                </p>
              </motion.div>

              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ duration: 0.5, delay: 0.4 }}
              >
                <label htmlFor="smiles" className="block text-sm font-medium text-text-secondary mb-1 font-body">
                  Select SMILES String
                </label>
                <select
                  id="smiles"
                  value={selectedSmiles}
                  onChange={(e) => setSelectedSmiles(e.target.value)}
                  className="w-full px-4 py-3 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent transition-all duration-200 font-body"
                  disabled={loading || !symptomGroupIndex || productSmilesGroups[symptomGroupIndex]?.length === 0}
                >
                  {productSmilesGroups[symptomGroupIndex]?.length === 0 ? (
                    <option value="" className="bg-primary">No SMILES available</option>
                  ) : (
                    productSmilesGroups[symptomGroupIndex]?.map((smiles, index) => (
                      <option key={index} value={smiles} className="bg-primary">
                        {smiles}
                      </option>
                    ))
                  )}
                </select>
                <p className="mt-1 text-xs text-text-secondary font-body">
                  Select a SMILES string corresponding to the symptom group
                </p>
              </motion.div>

              <motion.button
                type="submit"
                disabled={loading || !selectedSmiles || symptomGroupIndex === ""}
                className={`w-full px-4 py-3 text-primary rounded-lg transition-all duration-200 flex items-center justify-center font-heading ${
                  loading || !selectedSmiles || symptomGroupIndex === ""
                    ? "bg-gray-400 cursor-not-allowed"
                    : "bg-accent hover:bg-accent/90"
                }`}
                whileHover={!loading && selectedSmiles && symptomGroupIndex !== "" ? { scale: 1.02 } : {}}
                whileTap={!loading && selectedSmiles && symptomGroupIndex !== "" ? { scale: 0.98 } : {}}
              >
                <>
                  {loading ? (
                    <motion.div
                      animate={{ rotate: 20, y: [0, -2, 0] }}
                      transition={{
                        duration: 0.5,
                        repeat: Infinity,
                        repeatType: "reverse",
                        ease: "easeInOut",
                      }}
                    >
                      <DollarSign size={20} className="mr-2" />
                    </motion.div>
                  ) : (
                    <DollarSign size={20} className="mr-2" />
                  )}
                  {loading ? "Estimating..." : "Estimate Cost"}
                </>
              </motion.button>
            </form>

            {error && (
              <motion.div
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                className="mt-6 bg-error/10 border border-error text-text-primary px-4 py-3 rounded-lg flex justify-between items-center font-body"
              >
                <p>{error}</p>
                <button
                  className="text-error underline hover:text-error/80"
                  onClick={() => setError(null)}
                >
                  Dismiss
                </button>
              </motion.div>
            )}

            {result && (
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ duration: 0.5 }}
                className="mt-8"
              >
                <div className="flex items-center mb-4">
                  <motion.div
                    animate={{ scale: [1, 1.1, 1] }}
                    transition={{ duration: 2, repeat: Infinity }}
                  >
                    <DollarSign className="h-6 w-6 text-success mr-2" />
                  </motion.div>
                  <h3 className="text-xl font-semibold text-text-primary font-heading">Estimation Result</h3>
                </div>
                <div className="bg-gradient-to-r from-primary/20 to-secondary/20 p-6 rounded-lg border border-accent/20">
                  <div className="space-y-4">
                    <div>
                      <p className="text-sm text-text-secondary mb-1 font-body">Symptom Group</p>
                      <p className="text-sm bg-primary p-2 text-text-secondary rounded border border-secondary font-body">
                        {symptomGroups[symptomGroupIndex].join(", ") || "N/A"}
                      </p>
                    </div>
                    <div>
                      <p className="text-sm text-text-secondary mb-1 font-body">SMILES</p>
                      <p className="font-mono text-sm bg-primary p-2 text-text-secondary rounded border border-secondary overflow-x-auto font-code">
                        {result.smiles || "N/A"}
                      </p>
                    </div>
                    <div>
                      <p className="text-sm text-text-secondary mb-1 font-body">Estimated Cost</p>
                      <motion.p
                        className="text-3xl font-bold text-success font-heading"
                        initial={{ scale: 0.5 }}
                        animate={{ scale: 1 }}
                        transition={{ type: "spring", stiffness: 300 }}
                      >
                        {result.estimatedcost || "N/A"}
                      </motion.p>
                    </div>
                    <div className="space-y-3">
                      <motion.button
                        onClick={toggleResultInfo}
                        whileHover={{ scale: 1.01 }}
                        whileTap={{ scale: 0.99 }}
                        className="flex items-center justify-between w-full p-3 bg-accent-secondary/10 hover:bg-accent-secondary/20 rounded-lg transition-colors duration-200 border border-accent-secondary/30 focus:outline-none focus:ring-2 focus:ring-accent-secondary/50 focus:ring-offset-2 focus:ring-offset-secondary font-body"
                      >
                        <div className="flex items-center space-x-2">
                          <Info className="h-5 w-5 text-accent-secondary" />
                          <span className="text-base font-semibold text-text-primary">
                            Detailed Analysis <span className="text-sm text-accent font-label">(powered by Gemini)</span>
                          </span>
                        </div>
                        {isResultInfoOpen ? (
                          <ChevronUp className="h-5 w-5 text-accent-secondary transform transition-transform duration-300" />
                        ) : (
                          <ChevronDown className="h-5 w-5 text-accent-secondary transform transition-transform duration-300" />
                        )}
                      </motion.button>
                      <AnimatePresence>
                        {isResultInfoOpen && (
                          <motion.div
                            initial={{ opacity: 0, height: 0 }}
                            animate={{ opacity: 1, height: "auto" }}
                            exit={{ opacity: 0, height: 0 }}
                            transition={{ duration: 0.3 }}
                            className="overflow-hidden"
                          >
                            <div className="space-y-4">
                              <div
                                ref={infoRef}
                                className="p-4 bg-secondary rounded-lg shadow-lg border border-secondary"
                              >
                                <div className="prose prose-sm text-text-primary max-w-none font-body">
                                  {renderInformation(result.information)}
                                </div>
                              </div>
                              <div className="flex justify-end">
                                <motion.button
                                  onClick={exportToPDF}
                                  whileHover={{ scale: 1.03 }}
                                  whileTap={{ scale: 0.97 }}
                                  className="flex items-center px-4 py-2.5 bg-gradient-to-br from-success to-success/90 text-primary font-medium rounded-md hover:from-success/90 hover:to-success transition-all duration-200 shadow-sm hover:shadow-md focus:outline-none focus:ring-2 focus:ring-success focus:ring-offset-2 focus:ring-offset-secondary font-heading"
                                >
                                  <Download size={18} className="mr-2" />
                                  Export Report
                                </motion.button>
                              </div>
                            </div>
                          </motion.div>
                        )}
                      </AnimatePresence>
                    </div>
                  </div>
                </div>
              </motion.div>
            )}
          </motion.div>

          {/* History Section */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.5, delay: 0.2 }}
            className="bg-secondary p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg border border-secondary"
          >
            <div className="flex items-center justify-between mb-6">
              <div className="flex items-center">
                <motion.div
                  animate={{ rotate: 360 }}
                  transition={{ duration: 8, repeat: Infinity, ease: "linear" }}
                >
                  <Clock className="h-6 w-6 text-accent mr-2" />
                </motion.div>
                <h2 className="text-2xl font-semibold text-text-primary font-heading">Estimation History</h2>
              </div>
              <motion.button
                onClick={fetchHistory}
                disabled={historyLoading}
                className="text-accent hover:text-accent/80 transition-colors duration-200"
                title="Refresh history"
                whileHover={{ rotate: 90 }}
                whileTap={{ scale: 0.9 }}
              >
                <RefreshCw size={20} className={historyLoading ? "animate-spin" : ""} />
              </motion.button>
            </div>

            {historyLoading ? (
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                className="flex flex-col items-center justify-center py-12"
              >
                <RefreshCw size={24} className="animate-spin mb-4 text-accent" />
                <p className="text-text-secondary font-body">Loading history...</p>
              </motion.div>
            ) : history.length > 0 ? (
              <div className="space-y-4 max-h-[500px] overflow-y-auto pr-2 custom-scrollbar">
                <AnimatePresence>
                  {history.map((item) => (
                    <motion.div
                      key={item._id}
                      initial={{ opacity: 0, y: 20 }}
                      animate={{ opacity: 1, y: 0 }}
                      exit={{ opacity: 0, x: -20 }}
                      transition={{ duration: 0.3 }}
                      className="border border-secondary p-4 rounded-lg hover:border-accent-secondary transition-all duration-200 hover:shadow-sm"
                    >
                      <div className="flex justify-between items-start mb-2">
                        <div className="flex items-center">
                          <DollarSign size={16} className="text-success mr-1" />
                          <span className="font-bold text-lg text text-success font-heading">
                            {item.estimatedcost || "N/A"}
                          </span>
                        </div>
                        <span className="text-xs text-text-secondary bg-primary px-2 py-1 rounded-full font-body">
                          {item.created ? new Date(item.created).toLocaleDateString() : "N/A"}
                        </span>
                      </div>

                      <div className="mt-2">
                        <p className="text-xs text-text-secondary mb-1 font-body">Symptom Group</p>
                        <p className="text-sm bg-primary p-2 rounded border text-text-secondary border-secondary font-body">
                          {item.symptomsGrp || "N/A"}
                        </p>
                      </div>

                      <div className="mt-2">
                        <p className="text-xs text-text-secondary mb-1 font-body">SMILES</p>
                        <p className="font-mono text-sm bg-primary p-2 rounded border text-text-secondary border-secondary overflow-x-auto font-code">
                          {item.smiles || "N/A"}
                        </p>
                      </div>

                      {item.information && (
                        <div className="mt-3 space-y-2">
                          <motion.button
                            onClick={() => toggleHistoryItem(item._id)}
                            whileHover={{ scale: 1.01 }}
                            whileTap={{ scale: 0.99 }}
                            className="flex items-center justify-between w-full p-2.5 bg-accent-secondary/10 hover:bg-accent-secondary/20 rounded-md transition-colors duration-200 border border-accent-secondary/30 focus:outline-none focus:ring-2 focus:ring-accent-secondary/50 focus:ring-offset-1 focus:ring-offset-secondary font-body"
                          >
                            <div className="flex items-center space-x-2">
                              <Info className="h-4 w-4 text-accent-secondary" />
                              <span className="text-xs font-medium text-text-primary uppercase tracking-wide">
                                Details <span className="text-xs text-accent font-label">(powered by Gemini)</span>
                              </span>
                            </div>
                            {openHistoryItems[item._id] ? (
                              <ChevronUp className="h-4 w-4 text-accent-secondary transform transition-transform duration-300" />
                            ) : (
                              <ChevronDown className="h-4 w-4 text-accent-secondary transform transition-transform duration-300" />
                            )}
                          </motion.button>

                          <AnimatePresence>
                            {openHistoryItems[item._id] && (
                              <motion.div
                                initial={{ opacity: 0, height: 0 }}
                                animate={{ opacity: 1, height: "auto" }}
                                exit={{ opacity: 0, height: 0 }}
                                transition={{ duration: 0.3 }}
                                className="overflow-hidden"
                              >
                                <div className="space-y-3">
                                  <div className="p-3 bg-secondary rounded-md shadow-sm border border-secondary">
                                    <div className="prose prose-sm text-text-primary max-w-none font-body">
                                      {renderInformation(item.information)}
                                    </div>
                                  </div>
                                </div>
                              </motion.div>
                            )}
                          </AnimatePresence>
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
                className="flex flex-col items-center justify-center py-12 text-center"
              >
                <Database size={32} className="text-text-secondary mb-4" />
                <p className="text-text-secondary mb-2 font-body">No previous estimations found.</p>
                <p className="text-sm text-text-secondary font-body">
                  Your estimation history will appear here after you submit your first request.
                </p>
              </motion.div>
            )}
          </motion.div>
        </div>
      </div>

      <style jsx>{`
        .animate-fadeIn {
          animation: fadeIn 0.5s ease-in-out;
        }
        
        @keyframes fadeIn {
          from { opacity: 0; transform: translateY(10px); }
          to { opacity: 1; transform: translateY(0); }
        }
        
        .custom-scrollbar::-webkit-scrollbar {
          width: 6px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-track {
          background: #172A45;
          border-radius: 10px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-thumb {
          background: #5E81F4;
          border-radius: 10px;
        }
        
        .custom-scrollbar::-webkit-scrollbar-thumb:hover {
          background: #00F5D4;
        }
      `}</style>
    </div>
  );
};

export default CostEstimationForm;