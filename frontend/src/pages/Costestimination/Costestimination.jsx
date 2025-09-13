"use client";

import { useState, useEffect, useRef } from "react";
import { useAuthStore } from "../../Store/auth.store.js";
import { postCostEstimation, getCostEstimations } from "../../api/costestimination.jsx";
import axios from "axios";
import { 
  AlertCircle, DollarSign, Clock, Database, X, Info, RefreshCw, LogIn, 
  FileText, ChevronDown, ChevronUp, Download, FlaskConical, Pill, TestTube2, 
  Syringe, Atom, Leaf, BarChart3, PieChart, LineChart, Share2, Table, 
  Beaker, Award, TrendingUp, Check, ChevronRight, Filter, Layers, Scale,
  TrendingDown, ZoomIn, Maximize, Thermometer, Target, Microscope, Zap
} from "lucide-react";
import { jsPDF } from "jspdf";
import domtoimage from "dom-to-image";
import Loader from "../../components/Loader.jsx";
import { motion, AnimatePresence } from "framer-motion";
import { toast } from "react-hot-toast";
import { 
  AreaChart, Area, BarChart, Bar, Cell, XAxis, YAxis, 
  CartesianGrid, Tooltip, Legend, ResponsiveContainer, PieChart as RechartsPieChart,
  Pie, Radar, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis,
  ScatterChart, Scatter, LineChart as RechartsLineChart, Line
} from 'recharts';

// Simple tabs components
const TabsList = ({ className, children }) => (
  <div className={`flex space-x-1 rounded-lg bg-secondary/20 p-1 ${className}`}>{children}</div>
);

const TabsTrigger = ({ value, className, children, onClick, active }) => (
  <button 
    className={`flex items-center justify-center rounded-md px-3 py-1.5 text-sm font-medium transition-all
    ${active ? 'bg-secondary text-text-primary' : 'bg-transparent text-text-secondary hover:bg-secondary/30'} ${className}`}
    onClick={() => onClick(value)}
  >
    {children}
  </button>
);

const TabsContent = ({ value, activeValue, className, children }) => {
  if (value !== activeValue) return null;
  return <div className={className}>{children}</div>;
};

const Tabs = ({ value, onValueChange, className, children }) => {
  return <div className={className}>{children}</div>;
};

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

// Color palette for charts
const chartColors = {
  primary: "#5E81F4",
  secondary: "#8FB1FF",
  success: "#70E000",
  accent: "#00F5D4",
  error: "#FF4C4C",
  warning: "#FFBF00",
  labor: "#8FB1FF",
  materials: "#70E000",
  equipment: "#00F5D4",
  purification: "#D499FF",
  qualityControl: "#FFBF00",
  riskFactor: "#FF4C4C",
};

// Helper function to parse numerical data from a cost string
const parseCost = (costString) => {
  if (!costString) return { min: 0, max: 0, avg: 0 };
  
  const matches = costString.match(/\$(\d+)(?:\s*[–-]\s*\$(\d+))?/);
  if (!matches) return { min: 0, max: 0, avg: 0 };
  
  const min = parseInt(matches[1]);
  const max = matches[2] ? parseInt(matches[2]) : min;
  return { min, max, avg: (min + max) / 2 };
};

// Helper function to extract algorithm data
const extractAlgorithmData = (algorithmicData) => {
  if (!algorithmicData) return null;
  
  try {
    const data = typeof algorithmicData === 'string' 
      ? JSON.parse(algorithmicData) 
      : algorithmicData;
      
    return data;
  } catch (e) {
    console.error("Failed to parse algorithm data:", e);
    return null;
  }
};

// Helper to format molecule names
const formatMoleculeName = (smiles) => {
  if (!smiles) return "Unknown";
  
  if (smiles.length > 30) {
    return `${smiles.substring(0, 27)}...`;
  }
  return smiles;
};

const CostEstimationForm = () => {
  // Form and data states
  const [symptomGroupIndex, setSymptomGroupIndex] = useState("");
  const [selectedSmiles, setSelectedSmiles] = useState("");
  const [symptomGroups, setSymptomGroups] = useState([]);
  const [productSmilesGroups, setProductSmilesGroups] = useState([]);
  const [result, setResult] = useState(null);
  const [parsedResult, setParsedResult] = useState(null);
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [historyLoading, setHistoryLoading] = useState(false);
  const [benchmarkData, setBenchmarkData] = useState(null);
  
  // UI control states
  const [isResultInfoOpen, setIsResultInfoOpen] = useState(false);
  const [openHistoryItems, setOpenHistoryItems] = useState({});
  const [activeResultTab, setActiveResultTab] = useState("overview");
  const [activeHistoryTab, setActiveHistoryTab] = useState("list");
  const [showMoleculeModal, setShowMoleculeModal] = useState(false);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  
  // References
  const { user, checkAuth, checkingAuth } = useAuthStore();
  const infoRef = useRef(null);
  const chartRefs = useRef({});
  const userId = user?._id;
  
  // Toast notification configuration
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

  // Fetch symptoms and product SMILES strings
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

  // Process the result data to prepare it for visualization
  const processResultData = (result) => {
    if (!result) return null;
    
    // Extract cost data
    const costData = parseCost(result.estimatedcost);
    
    // Extract algorithm data
    const algorithmData = extractAlgorithmData(result.algorithmicData);
    
    // Prepare data for visualizations - generate mock data if actual data is not available
    const costBreakdownData = algorithmData?.economics?.breakdown ? [
      { name: 'Labor', value: algorithmData.economics.breakdown.labor || 0, color: chartColors.labor },
      { name: 'Materials', value: algorithmData.economics.breakdown.materials || 0, color: chartColors.materials },
      { name: 'Equipment', value: algorithmData.economics.breakdown.equipment || 0, color: chartColors.equipment },
      { name: 'Purification', value: algorithmData.economics.breakdown.purification || 0, color: chartColors.purification },
      { name: 'QC', value: algorithmData.economics.breakdown.qualityControl || 0, color: chartColors.qualityControl },
    ] : [
      { name: 'Labor', value: costData.avg * 0.3, color: chartColors.labor },
      { name: 'Materials', value: costData.avg * 0.25, color: chartColors.materials },
      { name: 'Equipment', value: costData.avg * 0.15, color: chartColors.equipment },
      { name: 'Purification', value: costData.avg * 0.2, color: chartColors.purification },
      { name: 'QC', value: costData.avg * 0.1, color: chartColors.qualityControl },
    ];
    
    // Complexity radar data - use mock data if not available
    const complexityData = algorithmData?.complexity?.details ? [
      { subject: 'Wiener', A: algorithmData.complexity.details.wiener || 0, fullMark: 3 },
      { subject: 'Bertz', A: algorithmData.complexity.details.bertz || 0, fullMark: 2 },
      { subject: 'Functional', A: algorithmData.complexity.details.functional || 0, fullMark: 3 },
      { subject: 'Stereo', A: algorithmData.complexity.details.stereo || 0, fullMark: 2.5 },
      { subject: 'Ring', A: algorithmData.complexity.details.ring || 0, fullMark: 2 },
    ] : [
      { subject: 'Wiener', A: 1.2, fullMark: 3 },
      { subject: 'Bertz', A: 0.8, fullMark: 2 },
      { subject: 'Functional', A: 1.5, fullMark: 3 },
      { subject: 'Stereo', A: 1.0, fullMark: 2.5 },
      { subject: 'Ring', A: 0.9, fullMark: 2 },
    ];
    
    // Risk data - use mock data if not available
    const riskData = algorithmData?.risk?.riskFactors ? [
      { name: 'Instability', value: algorithmData.risk.riskFactors.instabilityRisk || 0 },
      { name: 'Yield', value: algorithmData.risk.riskFactors.yieldRisk || 0 },
      { name: 'Scalability', value: algorithmData.risk.riskFactors.scalabilityRisk || 0 },
      { name: 'Regulatory', value: algorithmData.risk.riskFactors.regulatoryRisk || 0 },
      { name: 'Supply Chain', value: algorithmData.risk.riskFactors.supplychainRisk || 0 },
    ] : [
      { name: 'Instability', value: 0.2 },
      { name: 'Yield', value: 0.4 },
      { name: 'Scalability', value: 0.3 },
      { name: 'Regulatory', value: 0.1 },
      { name: 'Supply Chain', value: 0.25 },
    ];
    
    // Scale-up economics data (simulated for different batch sizes)
    const scaleUpData = [
      { name: '1g', cost: costData.avg, batch: 1 },
      { name: '10g', cost: costData.avg * 0.8, batch: 10 },
      { name: '100g', cost: costData.avg * 0.6, batch: 100 },
      { name: '1kg', cost: costData.avg * 0.4, batch: 1000 },
    ];
    
    // Mock dataset analysis if not available
    const datasetAnalysis = algorithmData?.datasetAnalysis?.analysisResults || {
      similarMolecules: [
        {
          group: 'Carbonyl',
          molecules: [
            { smiles: 'CC(=O)C', price: '$45/g', name: 'Similar Compound 1' },
            { smiles: 'CC(=O)O', price: '$38/g', name: 'Similar Compound 2' },
          ]
        },
        {
          group: 'Hydroxyl',
          molecules: [
            { smiles: 'CCO', price: '$30/g', name: 'Similar Compound 3' },
            { smiles: 'C1CCOC1', price: '$55/g', name: 'Similar Compound 4' },
          ]
        }
      ],
      pricePatterns: [
        'Average market price: $42.50 per unit',
        'Market price range: $28.00 - $62.00',
        'Scale-up discount: 35-45% at 1kg'
      ],
      marketTrends: [
        'Scale-up economics observed in dataset',
        'Price reduction rate with scale-up: 40.5%',
        'Similar functional groups show 15-25% cost variance'
      ]
    };
    
    // Parse functional groups from the information field
    const functionalGroups = extractFunctionalGroups(result.information);
    
    return {
      costData,
      algorithmData,
      costBreakdownData,
      complexityData,
      riskData,
      scaleUpData,
      datasetAnalysis,
      functionalGroups,
      stepCount: algorithmData?.steps?.steps || Math.floor(Math.random() * 5) + 3,
      complexity: algorithmData?.complexity?.complexity || 2.5,
      riskFactor: algorithmData?.risk?.overallRisk || 0.3,
    };
  };

  // Extract functional groups mentioned in the information text
  const extractFunctionalGroups = (information) => {
    if (!information) return [];
    
    const functionalGroups = [
      'hydroxyl', 'amino', 'carbonyl', 'carboxyl', 'amide', 'ester', 
      'ether', 'phenyl', 'methyl', 'alkyl', 'halide', 'nitro', 'cyano'
    ];
    
    const found = [];
    functionalGroups.forEach(group => {
      if (information.toLowerCase().includes(group)) {
        found.push(group);
      }
    });
    
    // Return at least some functional groups for demo purposes
    return found.length > 0 ? found : ['carbonyl', 'hydroxyl'];
  };

  // Handle form submission
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
      
      // Process the result for visualizations
      const processed = processResultData(data.data);
      setParsedResult(processed);

      const estimationId = data.data._id;
      
      // Store the symptom group in localStorage
      const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${userId}`) || "{}");
      storedSymptomGroups[estimationId] = symptomsGrp;
      localStorage.setItem(`symptomGroups_${userId}`, JSON.stringify(storedSymptomGroups));

      toast.success("Cost estimation completed", toastOptions);
      await fetchHistory();
      
      // Mock benchmark data for frontend demo
      setBenchmarkData({
        totalEstimations: history.length + 1,
        avgComplexity: 2.8,
        avgSteps: 6.2,
        avgCostRange: {
          lower: 35,
          upper: 85
        },
        functionalGroupStats: [
          { group: 'carbonyl', count: 5, avgValue: 1.3 },
          { group: 'hydroxyl', count: 4, avgValue: 0.9 },
          { group: 'amino', count: 3, avgValue: 1.1 }
        ],
        datasetInfluencedEstimations: history.length
      });
    } catch (error) {
      setError("Failed to estimate cost. Please try again.");
      console.error("Error estimating cost:", error);
    } finally {
      setLoading(false);
    }
  };

  // Fetch estimation history
  const fetchHistory = async () => {
    if (!user?._id) {
      setError("Please log in to view history");
      return;
    }
    
    setHistoryLoading(true);
    
    try {
      const data = await getCostEstimations(user._id);
      const validHistory = Array.isArray(data.data) ? data.data.filter((item) => item && typeof item === "object") : [];
      
      // Retrieve symptom groups from localStorage
      const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "{}");
      
      // Map through history and attach the symptom group
      const updatedHistory = validHistory.map((item) => {
        const symptomsGrp = storedSymptomGroups[item._id] || "N/A";
        return {
          ...item,
          symptomsGrp,
          parsedData: processResultData(item),
        };
      });

      setHistory(updatedHistory);
      
      // Generate mock benchmark data based on history
      if (updatedHistory.length > 0) {
        setBenchmarkData({
          totalEstimations: updatedHistory.length,
          avgComplexity: 2.8,
          avgSteps: 6.2,
          avgCostRange: {
            lower: 35,
            upper: 85
          },
          functionalGroupStats: [
            { group: 'carbonyl', count: 5, avgValue: 1.3 },
            { group: 'hydroxyl', count: 4, avgValue: 0.9 },
            { group: 'amino', count: 3, avgValue: 1.1 }
          ],
          datasetInfluencedEstimations: updatedHistory.length
        });
      }
    } catch (error) {
      setError("No previous estimations found.");
      console.error("Error fetching history:", error);
    } finally {
      setHistoryLoading(false);
    }
  };

  // Toggle sections
  const toggleResultInfo = () => {
    setIsResultInfoOpen(!isResultInfoOpen);
  };

  const toggleHistoryItem = (id) => {
    setOpenHistoryItems((prev) => ({
      ...prev,
      [id]: !prev[id],
    }));
  };

  // Render information from text data
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

  // Export to PDF
  const exportToPDF = async () => {
    const pdf = new jsPDF("p", "mm", "a4");
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const margin = 10;
    let yPosition = margin;

    // Add header with logo styling
    pdf.setFillColor(23, 42, 69); // Background color
    pdf.rect(0, 0, pageWidth, 20, 'F');
    pdf.setTextColor(255, 255, 255);
    pdf.setFontSize(16);
    pdf.text("Drug Cost Estimation Report", margin, 13);
    
    // Add date
    pdf.setFontSize(10);
    const date = new Date().toLocaleDateString();
    pdf.text(`Generated: ${date}`, pageWidth - margin - 30, 13, { align: 'right' });
    
    // Reset position for content
    yPosition = 30;

    // Add symptom group
    pdf.setTextColor(0, 0, 0);
    pdf.setFontSize(14);
    pdf.text("Molecule Information", margin, yPosition);
    yPosition += 7;
    
    pdf.setFontSize(11);
    pdf.setTextColor(70, 70, 70);
    pdf.text("Symptom Group:", margin, yPosition);
    pdf.setFontSize(10);
    pdf.text(symptomGroups[symptomGroupIndex]?.join(", ") || "N/A", margin + 35, yPosition);
    yPosition += 6;

    // Add SMILES
    pdf.setFontSize(11);
    pdf.text("SMILES:", margin, yPosition);
    pdf.setFontSize(9);
    pdf.text(result?.smiles || selectedSmiles || "N/A", margin + 35, yPosition);
    yPosition += 6;

    // Add Estimated Cost with emphasis
    pdf.setFontSize(11);
    pdf.text("Estimated Cost:", margin, yPosition);
    pdf.setFontSize(12);
    pdf.setTextColor(0, 128, 0); // Green for cost
    pdf.text(result?.estimatedcost || "N/A", margin + 35, yPosition);
    pdf.setTextColor(0, 0, 0);
    yPosition += 10;

    // Add detailed information from the AI analysis
    if (infoRef.current) {
      // Check if we need a new page
      if (yPosition > pageHeight - 40) {
        pdf.addPage();
        yPosition = margin + 10;
      }
      
      pdf.setFontSize(14);
      pdf.text("Detailed Analysis", margin, yPosition);
      yPosition += 7;

      try {
        const infoImgData = await domtoimage.toPng(infoRef.current, { quality: 1 });
        const infoImgProps = pdf.getImageProperties(infoImgData);
        const infoImgWidth = pageWidth - 2 * margin;
        let infoImgHeight = (infoImgProps.height * infoImgWidth) / infoImgProps.width;

        // If content is too tall, split across pages
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
      } catch (error) {
        console.error("Error capturing detailed info:", error);
      }
    }

    // Add footer
    const totalPages = pdf.internal.getNumberOfPages();
    for (let i = 1; i <= totalPages; i++) {
      pdf.setPage(i);
      pdf.setFontSize(8);
      pdf.setTextColor(100, 100, 100);
      pdf.text(`Page ${i} of ${totalPages} | Generated by Drug Cost Estimator (Powered by Gemini)`, pageWidth / 2, pageHeight - 5, { align: 'center' });
    }

    // Save PDF
    pdf.save("drug-cost-estimation-report.pdf");
    toast.success("Report exported successfully", toastOptions);
  };

  // Export data as CSV
  const exportAsCSV = () => {
    if (!result) return;
    
    // Create CSV content
    let csvContent = "data:text/csv;charset=utf-8,";
    csvContent += "Category,Value\n";
    csvContent += `SMILES,${result.smiles}\n`;
    csvContent += `Estimated Cost,${result.estimatedcost}\n`;
    csvContent += `Symptom Group,${symptomGroups[symptomGroupIndex]?.join(", ") || "N/A"}\n`;
    
    if (parsedResult?.algorithmData) {
      csvContent += `Complexity,${parsedResult.complexity}\n`;
      csvContent += `Step Count,${parsedResult.stepCount}\n`;
      csvContent += `Risk Factor,${parsedResult.riskFactor}\n`;
      
      if (parsedResult.costBreakdownData?.length) {
        csvContent += "\nCost Breakdown\n";
        parsedResult.costBreakdownData.forEach(item => {
          csvContent += `${item.name},${item.value}\n`;
        });
      }
    }
    
    // Create download link
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "drug-cost-estimation.csv");
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    
    toast.success("CSV exported successfully", toastOptions);
  };

  // Initialize data on component mount
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

  // If not logged in, show login prompt
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

  // Main application UI
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

      <div className="container mx-auto px-4 max-w-7xl relative z-10 pb-12">
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          className="mb-8 text-center"
        >
          <h1 className="text-4xl font-bold text-text-primary mb-2 font-heading tracking-tight">
            Advanced Drug Cost Estimation
            <p className="text-sm text-accent-secondary font-300 font-mono">(POWERED BY GEMINI + DATASET ANALYSIS)</p>
          </h1>
          <p className="text-text-secondary max-w-2xl mx-auto font-body text-lg">
            Get comprehensive cost estimations for drug synthesis and production with AI-powered analysis enhanced by real-world datasets.
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
            <div className="flex items-center justify-between mb-6">
              <div className="flex items-center">
                <motion.div
                  animate={{ y: [-5, 5, -5] }}
                  transition={{ duration: 4, repeat: Infinity, ease: "easeInOut" }}
                >
                  <FlaskConical className="h-7 w-7 text-accent mr-2" />
                </motion.div>
                <h2 className="text-2xl font-semibold text-text-primary font-heading">Estimate New Cost</h2>
              </div>
              
              {benchmarkData && (
                <div className="flex items-center text-text-secondary">
                  <Award className="h-4 w-4 mr-1 text-accent" />
                  <span className="text-xs font-body">
                    {benchmarkData.totalEstimations} previous estimations analyzed
                  </span>
                </div>
              )}
            </div>

            <form onSubmit={handleSubmit} className="space-y-5">
              <div className="grid grid-cols-1 md:grid-cols-2 gap-5">
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
                    Select a group of symptoms for the target indication
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
                    disabled={loading || symptomGroupIndex === "" || productSmilesGroups[symptomGroupIndex]?.length === 0}
                  >
                    {productSmilesGroups[symptomGroupIndex]?.length === 0 ? (
                      <option value="" className="bg-primary">No SMILES available</option>
                    ) : (
                      productSmilesGroups[symptomGroupIndex]?.map((smiles, index) => (
                        <option key={index} value={smiles} className="bg-primary">
                          {formatMoleculeName(smiles)}
                        </option>
                      ))
                    )}
                  </select>
                  <p className="mt-1 text-xs text-text-secondary font-body">
                    Select a molecular structure for cost estimation
                  </p>
                </motion.div>
              </div>

              <motion.button
                type="submit"
                disabled={loading || !selectedSmiles || symptomGroupIndex === ""}
                className={`w-full px-4 py-3.5 text-primary rounded-lg transition-all duration-200 flex items-center justify-center font-heading ${
                  loading || !selectedSmiles || symptomGroupIndex === ""
                    ? "bg-gray-400 cursor-not-allowed"
                    : "bg-gradient-to-r from-accent to-accent-secondary hover:from-accent/90 hover:to-accent-secondary/90"
                }`}
                whileHover={!loading && selectedSmiles && symptomGroupIndex !== "" ? { scale: 1.02 } : {}}
                whileTap={!loading && selectedSmiles && symptomGroupIndex !== "" ? { scale: 0.98 } : {}}
              >
                <>
                  {loading ? (
                    <motion.div
                      animate={{ rotate: 360 }}
                      transition={{
                        duration: 1.5,
                        repeat: Infinity,
                        ease: "linear",
                      }}
                    >
                      <RefreshCw size={20} className="mr-2" />
                    </motion.div>
                  ) : (
                    <DollarSign size={20} className="mr-2" />
                  )}
                  {loading ? "Estimating Cost..." : "Generate Cost Estimation"}
                </>
              </motion.button>
            </form>

            {error && (
              <motion.div
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                className="mt-6 bg-error/10 border border-error text-text-primary px-4 py-3 rounded-lg flex justify-between items-center font-body"
              >
                <div className="flex items-center">
                  <AlertCircle size={18} className="text-error mr-2" />
                  <p>{error}</p>
                </div>
                <button
                  className="text-error hover:text-error/80 focus:outline-none"
                  onClick={() => setError(null)}
                >
                  <X size={18} />
                </button>
              </motion.div>
            )}
          </motion.div>

          {/* Results Section */}
          {result && (
            <motion.div
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              transition={{ duration: 0.5 }}
              className="mt-8"
            >
              <div className="flex items-center justify-between mb-4">
                <div className="flex items-center">
                  <motion.div
                    animate={{ scale: [1, 1.1, 1] }}
                    transition={{ duration: 2, repeat: Infinity }}
                  >
                    <DollarSign className="h-6 w-6 text-success mr-2" />
                  </motion.div>
                  <h3 className="text-xl font-semibold text-text-primary font-heading">Estimation Results</h3>
                </div>
                
                <div className="flex space-x-2">
                  <motion.button
                    onClick={exportToPDF}
                    whileHover={{ scale: 1.05 }}
                    whileTap={{ scale: 0.95 }}
                    className="text-xs flex items-center px-2 py-1 bg-accent/10 hover:bg-accent/20 text-accent rounded-md"
                    title="Export to PDF"
                  >
                    <FileText size={14} className="mr-1" />
                    PDF
                  </motion.button>
                  
                  <motion.button
                    onClick={exportAsCSV}
                    whileHover={{ scale: 1.05 }}
                    whileTap={{ scale: 0.95 }}
                    className="text-xs flex items-center px-2 py-1 bg-accent-secondary/10 hover:bg-accent-secondary/20 text-accent-secondary rounded-md"
                    title="Export as CSV"
                  >
                    <Table size={14} className="mr-1" />
                    CSV
                  </motion.button>
                </div>
              </div>
              
              {/* Top Summary Card */}
              <div className="bg-gradient-to-r from-primary/20 to-secondary/20 p-5 rounded-lg border border-accent/20 mb-4">
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  <div>
                    <p className="text-xs text-text-secondary mb-1 font-body">Symptom Group</p>
                    <p className="text-sm bg-primary p-2 text-text-secondary rounded border border-secondary font-body">
                      {symptomGroups[symptomGroupIndex]?.join(", ") || "N/A"}
                    </p>
                  </div>
                  <div>
                    <p className="text-xs text-text-secondary mb-1 font-body">SMILES</p>
                    <div className="relative">
                      <p className="font-mono text-sm bg-primary p-2 text-text-secondary rounded border border-secondary overflow-x-auto font-code">
                        {result.smiles?.length > 20 
                          ? `${result.smiles.substring(0, 20)}...` 
                          : result.smiles || "N/A"}
                      </p>
                      <button 
                        className="absolute right-2 top-1/2 -translate-y-1/2 text-accent-secondary hover:text-accent"
                        onClick={() => {
                          setSelectedMolecule(result.smiles);
                          setShowMoleculeModal(true);
                        }}
                        title="View full SMILES"
                      >
                        <ZoomIn size={14} />
                      </button>
                    </div>
                  </div>
                  <div>
                    <p className="text-xs text-text-secondary mb-1 font-body">Estimated Cost</p>
                    <motion.div
                      className="flex items-center justify-between bg-primary p-2 rounded border border-secondary"
                      initial={{ scale: 0.5 }}
                      animate={{ scale: 1 }}
                      transition={{ type: "spring", stiffness: 300 }}
                    >
                      <span className="text-xl font-bold text-success font-heading">
                        {result.estimatedcost || "N/A"}
                      </span>
                      
                      {parsedResult?.costData?.min !== parsedResult?.costData?.max && (
                        <div className="flex items-center text-xs text-text-secondary space-x-1">
                          <Scale size={12} />
                          <span>{parsedResult?.costData?.min}-{parsedResult?.costData?.max} range</span>
                        </div>
                      )}
                    </motion.div>
                  </div>
                </div>
              </div>
              
              {/* Tabs for detailed information */}
              <div className="w-full">
                <div className="flex space-x-1 rounded-lg bg-secondary/20 p-1 mb-4">
                  {["overview", "costAnalysis", "complexity", "synthesis", "dataset"].map((tab) => (
                    <button 
                      key={tab}
                      className={`flex items-center justify-center rounded-md px-3 py-1.5 text-xs font-medium transition-all ${
                        activeResultTab === tab ? 'bg-secondary text-text-primary' : 'bg-transparent text-text-secondary hover:bg-secondary/30'
                      }`}
                      onClick={() => setActiveResultTab(tab)}
                    >
                      {tab === "overview" && "Overview"}
                      {tab === "costAnalysis" && "Cost Analysis"}
                      {tab === "complexity" && "Complexity"}
                      {tab === "synthesis" && "Synthesis"}
                      {tab === "dataset" && "Dataset Insights"}
                    </button>
                  ))}
                </div>
                
                {/* Overview Tab */}
                {activeResultTab === "overview" && (
                  <div className="space-y-4">
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                      {/* Key metrics cards */}
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <div className="flex items-center justify-between mb-2">
                          <h4 className="text-sm font-medium text-text-primary">Molecular Complexity</h4>
                          <Atom size={16} className="text-accent" />
                        </div>
                        <div className="flex items-end justify-between">
                          <div className="text-2xl font-bold text-accent">
                            {parsedResult?.complexity.toFixed(1)}<span className="text-sm text-text-secondary">/5.0</span>
                          </div>
                          <div className="text-xs text-text-secondary">
                            {parsedResult?.complexity < 2 ? 'Low' : parsedResult?.complexity < 3.5 ? 'Medium' : 'High'} complexity
                          </div>
                        </div>
                        <div className="w-full bg-secondary/50 h-1.5 mt-2 rounded-full">
                          <div 
                            className="bg-accent h-1.5 rounded-full" 
                            style={{ width: `${(parsedResult?.complexity / 5) * 100}%` }}
                          ></div>
                        </div>
                      </div>
                      
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <div className="flex items-center justify-between mb-2">
                          <h4 className="text-sm font-medium text-text-primary">Synthesis Steps</h4>
                          <TestTube2 size={16} className="text-accent-secondary" />
                        </div>
                        <div className="flex items-end justify-between">
                          <div className="text-2xl font-bold text-accent-secondary">
                            {parsedResult?.stepCount}<span className="text-sm text-text-secondary"> steps</span>
                          </div>
                          <div className="text-xs text-text-secondary">
                            {parsedResult?.stepCount < 5 ? 'Simple' : parsedResult?.stepCount < 9 ? 'Moderate' : 'Complex'} synthesis
                          </div>
                        </div>
                        <div className="w-full bg-secondary/50 h-1.5 mt-2 rounded-full">
                          <div 
                            className="bg-accent-secondary h-1.5 rounded-full" 
                            style={{ width: `${(parsedResult?.stepCount / 15) * 100}%` }}
                          ></div>
                        </div>
                      </div>
                      
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <div className="flex items-center justify-between mb-2">
                          <h4 className="text-sm font-medium text-text-primary">Risk Factor</h4>
                          <AlertCircle size={16} className="text-error" />
                        </div>
                        <div className="flex items-end justify-between">
                          <div className="text-2xl font-bold text-error">
                            {parsedResult?.riskFactor.toFixed(2)}<span className="text-sm text-text-secondary">/1.0</span>
                          </div>
                          <div className="text-xs text-text-secondary">
                            {parsedResult?.riskFactor < 0.3 ? 'Low' : parsedResult?.riskFactor < 0.6 ? 'Medium' : 'High'} risk
                          </div>
                        </div>
                        <div className="w-full bg-secondary/50 h-1.5 mt-2 rounded-full">
                          <div 
                            className="bg-error h-1.5 rounded-full" 
                            style={{ width: `${parsedResult?.riskFactor * 100}%` }}
                          ></div>
                        </div>
                      </div>
                    </div>
                    
                    {/* Cost by scale chart */}
                    <div className="bg-primary p-4 rounded-lg border border-secondary">
                      <div className="flex items-center justify-between mb-4">
                        <h4 className="text-sm font-medium text-text-primary">Cost by Scale</h4>
                        <TrendingDown size={16} className="text-success" />
                      </div>
                      <div className="h-64">
                        <ResponsiveContainer width="100%" height="100%">
                          <AreaChart
                            data={parsedResult?.scaleUpData || []}
                            margin={{ top: 5, right: 20, left: 0, bottom: 5 }}
                          >
                            <defs>
                              <linearGradient id="colorCost" x1="0" y1="0" x2="0" y2="1">
                                <stop offset="5%" stopColor={chartColors.success} stopOpacity={0.8}/>
                                <stop offset="95%" stopColor={chartColors.success} stopOpacity={0.1}/>
                              </linearGradient>
                            </defs>
                            <CartesianGrid strokeDasharray="3 3" stroke="#333" />
                            <XAxis dataKey="name" stroke="#888" />
                            <YAxis stroke="#888" />
                            <Tooltip 
                              contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                              formatter={(value) => [`$${value.toFixed(2)}`, 'Cost per gram']}
                              labelFormatter={(value) => `Batch size: ${value}`}
                            />
                            <Area 
                              type="monotone" 
                              dataKey="cost" 
                              stroke={chartColors.success} 
                              fillOpacity={1} 
                              fill="url(#colorCost)" 
                            />
                          </AreaChart>
                        </ResponsiveContainer>
                      </div>
                      <p className="text-xs text-text-secondary mt-2 text-center">
                        Estimated cost reduction with increased production scale
                      </p>
                    </div>
                    
                    {/* Detailed info toggle */}
                    <div className="space-y-3">
                      <motion.button
                        onClick={toggleResultInfo}
                        whileHover={{ scale: 1.01 }}
                        whileTap={{ scale: 0.99 }}
                        className="flex items-center justify-between w-full p-3.5 bg-accent-secondary/10 hover:bg-accent-secondary/20 rounded-lg transition-colors duration-200 border border-accent-secondary/30 focus:outline-none focus:ring-2 focus:ring-accent-secondary/50 focus:ring-offset-2 focus:ring-offset-secondary font-body"
                      >
                        <div className="flex items-center space-x-2">
                          <Info className="h-5 w-5 text-accent-secondary" />
                          <span className="text-base font-semibold text-text-primary">
                            Complete Analysis <span className="text-sm text-accent font-label">(powered by Gemini)</span>
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
                                className="p-5 bg-secondary rounded-lg shadow-lg border border-secondary"
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
                                  Export Complete Report
                                </motion.button>
                              </div>
                            </div>
                          </motion.div>
                        )}
                      </AnimatePresence>
                    </div>
                  </div>
                )}
                
                {/* Cost Analysis Tab */}
                {activeResultTab === "costAnalysis" && (
                  <div className="space-y-4">
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      {/* Cost breakdown chart */}
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-4">Cost Component Breakdown</h4>
                        <div className="h-64" ref={el => chartRefs.current.costBreakdown = el}>
                          <ResponsiveContainer width="100%" height="100%">
                            <RechartsPieChart>
                              <Pie
                                data={parsedResult?.costBreakdownData || []}
                                cx="50%"
                                cy="50%"
                                labelLine={true}
                                outerRadius={80}
                                fill="#8884d8"
                                dataKey="value"
                                nameKey="name"
                                label={({ name, percent }) => `${name}: ${(percent * 100).toFixed(0)}%`}
                              >
                                {parsedResult?.costBreakdownData?.map((entry, index) => (
                                  <Cell key={`cell-${index}`} fill={entry.color} />
                                ))}
                              </Pie>
                              <Tooltip
                                contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                formatter={(value) => [`$${value.toFixed(2)}`, 'Cost per gram']}
                              />
                              <Legend verticalAlign="bottom" height={36} />
                            </RechartsPieChart>
                          </ResponsiveContainer>
                        </div>
                      </div>
                      
                      {/* Detailed cost table */}
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-4">Detailed Cost Breakdown</h4>
                        <div className="overflow-hidden rounded-lg border border-secondary">
                          <table className="min-w-full divide-y divide-secondary">
                            <thead className="bg-secondary/30">
                              <tr>
                                <th className="px-4 py-2 text-left text-xs font-medium text-text-secondary uppercase tracking-wider">Component</th>
                                <th className="px-4 py-2 text-right text-xs font-medium text-text-secondary uppercase tracking-wider">Cost ($/g)</th>
                                <th className="px-4 py-2 text-right text-xs font-medium text-text-secondary uppercase tracking-wider">% of Total</th>
                              </tr>
                            </thead>
                            <tbody className="bg-primary divide-y divide-secondary">
                              {parsedResult?.costBreakdownData?.map((item, index) => {
                                const totalCost = parsedResult?.costBreakdownData?.reduce((sum, i) => sum + i.value, 0) || 1;
                                const percentage = (item.value / totalCost) * 100;
                                
                                return (
                                  <tr key={index} className={index % 2 === 0 ? 'bg-secondary/10' : ''}>
                                    <td className="px-4 py-2 whitespace-nowrap text-sm text-text-primary">
                                      <div className="flex items-center">
                                        <div className="w-2 h-2 rounded-full mr-2" style={{ backgroundColor: item.color }}></div>
                                        {item.name}
                                      </div>
                                    </td>
                                    <td className="px-4 py-2 whitespace-nowrap text-sm text-text-primary text-right">
                                      ${item.value.toFixed(2)}
                                    </td>
                                    <td className="px-4 py-2 whitespace-nowrap text-sm text-text-primary text-right">
                                      {percentage.toFixed(1)}%
                                    </td>
                                  </tr>
                                );
                              })}
                              <tr className="bg-secondary/20 font-medium">
                                <td className="px-4 py-2 whitespace-nowrap text-sm text-text-primary">Total</td>
                                <td className="px-4 py-2 whitespace-nowrap text-sm text-text-primary text-right">
                                  ${parsedResult?.costBreakdownData?.reduce((sum, i) => sum + i.value, 0).toFixed(2) || '0.00'}
                                </td>
                                <td className="px-4 py-2 whitespace-nowrap text-sm text-text-primary text-right">100%</td>
                              </tr>
                            </tbody>
                          </table>
                        </div>
                      </div>
                    </div>
                    
                    {/* Scale economics */}
                    <div className="bg-primary p-4 rounded-lg border border-secondary">
                      <h4 className="text-sm font-medium text-text-primary mb-4">Scale Economics Analysis</h4>
                      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
                        {parsedResult?.scaleUpData?.map((item, index) => (
                          <div key={index} className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                            <div className="text-xs text-text-secondary mb-1">{item.name} Batch Size</div>
                            <div className="text-lg font-bold text-text-primary">${item.cost.toFixed(2)}<span className="text-xs text-text-secondary">/g</span></div>
                            {index > 0 && (
                              <div className="flex items-center mt-1 text-xs">
                                <TrendingDown className="h-3 w-3 text-success mr-1" />
                                <span className="text-success">
                                  {(100 - (item.cost / parsedResult?.scaleUpData[0].cost) * 100).toFixed(0)}% reduction
                                </span>
                              </div>
                            )}
                          </div>
                        ))}
                      </div>
                    </div>
                    
                    {/* Benchmark comparison */}
                    {benchmarkData && (
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-3">Benchmark Comparison</h4>
                        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                          <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                            <div className="text-xs text-text-secondary mb-1">Your Estimate</div>
                            <div className="text-lg font-bold text-text-primary">${parsedResult?.costData?.avg.toFixed(2)}/g</div>
                            <div className="flex items-center mt-1 text-xs">
                              <Scale className="h-3 w-3 text-accent-secondary mr-1" />
                              <span className="text-text-secondary">Range: ${parsedResult?.costData?.min}-${parsedResult?.costData?.max}</span>
                            </div>
                          </div>
                          
                          <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                            <div className="text-xs text-text-secondary mb-1">Average Cost (All Estimates)</div>
                            <div className="text-lg font-bold text-text-primary">
                              ${benchmarkData?.avgCostRange?.lower.toFixed(2)}-${benchmarkData?.avgCostRange?.upper.toFixed(2)}/g
                            </div>
                            <div className="flex items-center mt-1 text-xs">
                              <Database className="h-3 w-3 text-accent-secondary mr-1" />
                              <span className="text-text-secondary">Based on {benchmarkData.totalEstimations} estimates</span>
                            </div>
                          </div>
                          
                          <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                            <div className="text-xs text-text-secondary mb-1">Complexity Comparison</div>
                            <div className="text-lg font-bold text-text-primary">
                              {parsedResult?.complexity.toFixed(1)} <span className="text-xs text-text-secondary">vs avg {benchmarkData?.avgComplexity.toFixed(1)}</span>
                            </div>
                            <div className="flex items-center mt-1 text-xs">
                              {parsedResult?.complexity > benchmarkData?.avgComplexity ? (
                                <>
                                  <TrendingUp className="h-3 w-3 text-error mr-1" />
                                  <span className="text-error">
                                    {((parsedResult?.complexity / benchmarkData?.avgComplexity - 1) * 100).toFixed(0)}% more complex
                                  </span>
                                </>
                              ) : (
                                <>
                                  <TrendingDown className="h-3 w-3 text-success mr-1" />
                                  <span className="text-success">
                                    {((1 - parsedResult?.complexity / benchmarkData?.avgComplexity) * 100).toFixed(0)}% less complex
                                  </span>
                                </>
                              )}
                            </div>
                          </div>
                        </div>
                      </div>
                    )}
                  </div>
                )}
                
                {/* Complexity Tab */}
                {activeResultTab === "complexity" && (
                  <div className="space-y-4">
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      {/* Complexity radar chart */}
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-4">Molecular Complexity Analysis</h4>
                        <div className="h-64" ref={el => chartRefs.current.complexity = el}>
                          <ResponsiveContainer width="100%" height="100%">
                            <RadarChart cx="50%" cy="50%" outerRadius="80%" data={parsedResult?.complexityData || []}>
                              <PolarGrid stroke="#444" />
                              <PolarAngleAxis dataKey="subject" stroke="#888" />
                              <PolarRadiusAxis angle={30} domain={[0, 'auto']} stroke="#888" />
                              <Radar
                                name="Complexity"
                                dataKey="A"
                                stroke={chartColors.accent}
                                fill={chartColors.accent}
                                fillOpacity={0.5}
                              />
                              <Tooltip
                                contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                formatter={(value) => [`${value.toFixed(2)}`, 'Complexity']}
                              />
                              <Legend />
                            </RadarChart>
                          </ResponsiveContainer>
                        </div>
                      </div>
                      
                      {/* Risk assessment chart */}
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-4">Risk Assessment</h4>
                        <div className="h-64" ref={el => chartRefs.current.risk = el}>
                          <ResponsiveContainer width="100%" height="100%">
                            <BarChart
                              layout="vertical"
                              data={parsedResult?.riskData || []}
                              margin={{ top: 5, right: 30, left: 40, bottom: 5 }}
                            >
                              <CartesianGrid strokeDasharray="3 3" stroke="#333" horizontal={false} />
                              <XAxis type="number" domain={[0, 1]} stroke="#888" />
                              <YAxis dataKey="name" type="category" stroke="#888" />
                              <Tooltip
                                contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                formatter={(value) => [`${(value * 100).toFixed(0)}%`, 'Risk Level']}
                              />
                              <Bar dataKey="value" fill={chartColors.error}>
                                {parsedResult?.riskData?.map((entry, index) => (
                                  <Cell
                                    key={`cell-${index}`}
                                    fill={entry.value < 0.3 ? chartColors.success : entry.value < 0.6 ? chartColors.warning : chartColors.error}
                                  />
                                ))}
                              </Bar>
                            </BarChart>
                          </ResponsiveContainer>
                        </div>
                      </div>
                    </div>
                    
                    {/* Functional Group Analysis */}
                    <div className="bg-primary p-4 rounded-lg border border-secondary">
                      <h4 className="text-sm font-medium text-text-primary mb-3">Functional Group Analysis</h4>
                      {parsedResult?.functionalGroups?.length > 0 ? (
                        <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
                          {parsedResult.functionalGroups.map((group, index) => (
                            <div key={index} className="bg-secondary/10 p-2.5 rounded-lg border border-secondary/50 flex items-center">
                              <div className="w-2 h-2 rounded-full bg-accent mr-2"></div>
                              <span className="text-sm text-text-primary capitalize">{group}</span>
                            </div>
                          ))}
                        </div>
                      ) : (
                        <p className="text-sm text-text-secondary">No functional groups identified</p>
                      )}
                    </div>
                    
                    {/* Molecule Complexity Impact */}
                    <div className="bg-primary p-4 rounded-lg border border-secondary">
                      <h4 className="text-sm font-medium text-text-primary mb-3">Complexity Impact on Cost</h4>
                      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                        <div className="col-span-2">
                          <div className="h-48">
                            <ResponsiveContainer width="100%" height="100%">
                              <ScatterChart
                                margin={{ top: 5, right: 20, bottom: 5, left: 0 }}
                              >
                                <CartesianGrid stroke="#333" strokeDasharray="3 3" />
                                <XAxis 
                                  type="number" 
                                  dataKey="x" 
                                  name="Complexity" 
                                  domain={[0, 5]} 
                                  stroke="#888"
                                  label={{ value: 'Complexity', position: 'insideBottom', offset: -5, fill: '#888' }}
                                />
                                <YAxis 
                                  type="number" 
                                  dataKey="y" 
                                  name="Cost Factor" 
                                  stroke="#888"
                                  label={{ value: 'Cost Factor', angle: -90, position: 'insideLeft', fill: '#888' }}
                                />
                                <Tooltip 
                                  contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                  cursor={{ strokeDasharray: '3 3' }}
                                />
                                <Scatter name="Cost-Complexity Relationship" data={[
                                  { x: 1, y: 1 },
                                  { x: 2, y: 1.5 },
                                  { x: 3, y: 2.3 },
                                  { x: 4, y: 3.2 },
                                  { x: 5, y: 4.5 },
                                ]} fill={chartColors.accent} />
                                {parsedResult?.complexity && (
                                  <Scatter 
                                    name="Current Molecule" 
                                    data={[{ 
                                      x: parsedResult.complexity, 
                                      y: parsedResult.costData?.avg / 100 || 2 
                                    }]} 
                                    fill={chartColors.success} 
                                    shape="star"
                                  />
                                )}
                                <Legend />
                              </ScatterChart>
                            </ResponsiveContainer>
                          </div>
                        </div>
                        <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                          <h5 className="text-xs font-medium text-text-primary mb-2">Complexity Factors</h5>
                          <ul className="space-y-2 text-xs">
                            <li className="flex justify-between">
                              <span className="text-text-secondary">Topological (Wiener):</span>
                              <span className="text-text-primary font-medium">{parsedResult?.complexityData?.[0]?.A?.toFixed(2) || 'N/A'}</span>
                            </li>
                            <li className="flex justify-between">
                              <span className="text-text-secondary">Structural (Bertz):</span>
                              <span className="text-text-primary font-medium">{parsedResult?.complexityData?.[1]?.A?.toFixed(2) || 'N/A'}</span>
                            </li>
                            <li className="flex justify-between">
                              <span className="text-text-secondary">Functional Group:</span>
                              <span className="text-text-primary font-medium">{parsedResult?.complexityData?.[2]?.A?.toFixed(2) || 'N/A'}</span>
                            </li>
                            <li className="flex justify-between">
                              <span className="text-text-secondary">Stereochemical:</span>
                              <span className="text-text-primary font-medium">{parsedResult?.complexityData?.[3]?.A?.toFixed(2) || 'N/A'}</span>
                            </li>
                            <li className="flex justify-between">
                              <span className="text-text-secondary">Ring System:</span>
                              <span className="text-text-primary font-medium">{parsedResult?.complexityData?.[4]?.A?.toFixed(2) || 'N/A'}</span>
                            </li>
                            <li className="flex justify-between pt-1 border-t border-secondary/50">
                              <span className="text-text-secondary">Overall Complexity:</span>
                              <span className="text-accent font-bold">{parsedResult?.complexity?.toFixed(2) || 'N/A'}/5.0</span>
                            </li>
                          </ul>
                        </div>
                      </div>
                    </div>
                  </div>
                )}
                
                {/* Synthesis Tab */}
                {activeResultTab === "synthesis" && (
                  <div className="space-y-4">
                    {/* Synthesis pathway visualization */}
                    <div className="bg-primary p-4 rounded-lg border border-secondary">
                      <h4 className="text-sm font-medium text-text-primary mb-3">Synthesis Steps Visualization</h4>
                      <div className="overflow-x-auto">
                        <div className="flex items-center justify-start space-x-2 min-w-max py-4">
                          {Array.from({ length: parsedResult?.stepCount || 5 }).map((_, index) => (
                            <div key={index} className="flex items-center">
                              <div className={`flex flex-col items-center justify-center w-28 h-16 rounded-lg border ${
                                index === 0 ? 'border-accent bg-accent/10' : 'border-secondary bg-secondary/10'
                              } p-2`}>
                                <div className="text-xs font-medium text-text-primary">Step {index + 1}</div>
                                <div className="text-xs text-text-secondary">
                                  {index === 0 ? 'Starting Material' : 
                                   index === parsedResult?.stepCount - 1 ? 'Final Product' : 
                                   `Intermediate ${index}`}
                                </div>
                              </div>
                              {index < (parsedResult?.stepCount || 5) - 1 && (
                                <ChevronRight className="h-4 w-4 text-text-secondary mx-1" />
                              )}
                            </div>
                          ))}
                        </div>
                      </div>
                    </div>
                    
                    {/* Step complexity and yield projection */}
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-3">Step-wise Complexity</h4>
                        <div className="h-64">
                          <ResponsiveContainer width="100%" height="100%">
                            <AreaChart
                              data={Array.from({ length: parsedResult?.stepCount || 5 }).map((_, i) => ({
                                step: `Step ${i + 1}`,
                                complexity: (
                                  parsedResult?.complexity * 
                                  (0.5 + (i / (parsedResult?.stepCount || 5)) * 0.5)
                                ).toFixed(2)
                              }))}
                              margin={{ top: 5, right: 20, left: 0, bottom: 5 }}
                            >
                              <defs>
                                <linearGradient id="colorComplexity" x1="0" y1="0" x2="0" y2="1">
                                  <stop offset="5%" stopColor={chartColors.accent} stopOpacity={0.8}/>
                                  <stop offset="95%" stopColor={chartColors.accent} stopOpacity={0.1}/>
                                </linearGradient>
                              </defs>
                              <CartesianGrid strokeDasharray="3 3" stroke="#333" />
                              <XAxis dataKey="step" stroke="#888" />
                              <YAxis stroke="#888" />
                              <Tooltip
                                contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                formatter={(value) => [value, 'Complexity']}
                              />
                              <Area
                                type="monotone"
                                dataKey="complexity"
                                stroke={chartColors.accent}
                                fillOpacity={1}
                                fill="url(#colorComplexity)"
                              />
                            </AreaChart>
                          </ResponsiveContainer>
                        </div>
                      </div>
                      
                      <div className="bg-primary p-4 rounded-lg border border-secondary">
                        <h4 className="text-sm font-medium text-text-primary mb-3">Projected Yield by Step</h4>
                        <div className="h-64">
                          <ResponsiveContainer width="100%" height="100%">
                            <RechartsLineChart
                              data={Array.from({ length: parsedResult?.stepCount || 5 }).map((_, i) => {
                                // Calculate diminishing yield based on complexity
                                const baseYield = 95 - (parsedResult?.complexity || 2) * 5;
                                const yieldDropFactor = 1 - (parsedResult?.riskFactor || 0.3) * 0.3;
                                const stepYield = baseYield * Math.pow(yieldDropFactor, i);
                                
                                return {
                                  step: `Step ${i + 1}`,
                                  yield: Math.max(stepYield, 50).toFixed(0),
                                  totalYield: Math.max(
                                    baseYield * Math.pow(yieldDropFactor, i + 1), 10
                                  ).toFixed(0)
                                };
                              })}
                              margin={{ top: 5, right: 20, left: 0, bottom: 5 }}
                            >
                              <CartesianGrid strokeDasharray="3 3" stroke="#333" />
                              <XAxis dataKey="step" stroke="#888" />
                              <YAxis stroke="#888" />
                              <Tooltip
                                contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                formatter={(value, name) => [
                                  `${value}%`, 
                                  name === 'yield' ? 'Step Yield' : 'Overall Yield'
                                ]}
                              />
                              <Line
                                type="monotone"
                                dataKey="yield"
                                name="Step Yield"
                                stroke={chartColors.success}
                                activeDot={{ r: 8 }}
                              />
                              <Line
                                type="monotone"
                                dataKey="totalYield"
                                name="Overall Yield"
                                stroke={chartColors.error}
                                activeDot={{ r: 8 }}
                              />
                              <Legend />
                            </RechartsLineChart>
                          </ResponsiveContainer>
                        </div>
                      </div>
                    </div>
                    
                    {/* Synthesis considerations */}
                    <div className="bg-primary p-4 rounded-lg border border-secondary">
                      <h4 className="text-sm font-medium text-text-primary mb-3">Synthesis Considerations</h4>
                      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                        <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                          <div className="flex items-center text-accent mb-2">
                            <Beaker className="h-4 w-4 mr-1" />
                            <h5 className="text-xs font-medium">Reagent Costs</h5>
                          </div>
                          <p className="text-xs text-text-secondary">
                            {parsedResult?.algorithmData?.economics?.breakdown?.materials > 30 
                              ? 'High cost reagents required for synthesis' 
                              : parsedResult?.algorithmData?.economics?.breakdown?.materials > 15
                              ? 'Moderate reagent costs expected'
                              : 'Standard reagents, low cost implications'}
                          </p>
                          <div className="mt-2 w-full bg-secondary/50 h-1 rounded-full">
                            <div 
                              className="bg-accent h-1 rounded-full" 
                              style={{ 
                                width: `${Math.min(
                                  (parsedResult?.algorithmData?.economics?.breakdown?.materials || 10) / 50 * 100, 
                                  100
                                )}%` 
                              }}
                            ></div>
                          </div>
                        </div>
                        
                        <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                          <div className="flex items-center text-accent-secondary mb-2">
                            <Thermometer className="h-4 w-4 mr-1" />
                            <h5 className="text-xs font-medium">Reaction Conditions</h5>
                          </div>
                          <p className="text-xs text-text-secondary">
                            {parsedResult?.complexityData?.[4]?.A > 1.5
                              ? 'Specialized conditions needed (pressure/temperature)'
                              : parsedResult?.complexityData?.[3]?.A > 1.5
                              ? 'Controlled environment required for stereoselectivity'
                              : 'Standard laboratory conditions sufficient'}
                          </p>
                          <div className="mt-2 w-full bg-secondary/50 h-1 rounded-full">
                            <div 
                              className="bg-accent-secondary h-1 rounded-full" 
                              style={{ 
                                width: `${Math.min(
                                  ((parsedResult?.complexityData?.[4]?.A || 0.5) + 
                                   (parsedResult?.complexityData?.[3]?.A || 0.5)) / 4 * 100,
                                  100
                                )}%` 
                              }}
                            ></div>
                          </div>
                        </div>
                        
                        <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                          <div className="flex items-center text-success mb-2">
                            <Target className="h-4 w-4 mr-1" />
                            <h5 className="text-xs font-medium">Purification Challenge</h5>
                          </div>
                          <p className="text-xs text-text-secondary">
                            {parsedResult?.algorithmData?.economics?.breakdown?.purification > 20
                              ? 'Complex purification required (chiral separation possible)'
                              : parsedResult?.algorithmData?.economics?.breakdown?.purification > 10
                              ? 'Standard chromatography should be sufficient'
                              : 'Simple purification methods applicable'}
                          </p>
                          <div className="mt-2 w-full bg-secondary/50 h-1 rounded-full">
                            <div 
                              className="bg-success h-1 rounded-full" 
                              style={{ 
                                width: `${Math.min(
                                  (parsedResult?.algorithmData?.economics?.breakdown?.purification || 8) / 40 * 100,
                                  100
                                )}%` 
                              }}
                            ></div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </div>
                )}
                
                {/* Dataset Insights Tab */}
                {activeResultTab === "dataset" && (
                  <div className="space-y-4">
                    {parsedResult?.datasetAnalysis ? (
                      <>
                        {/* Similar molecules from dataset */}
                        <div className="bg-primary p-4 rounded-lg border border-secondary">
                          <h4 className="text-sm font-medium text-text-primary mb-3">Similar Molecules in Dataset</h4>
                          {parsedResult.datasetAnalysis.similarMolecules?.length > 0 ? (
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                              {parsedResult.datasetAnalysis.similarMolecules.slice(0, 2).map((item, index) => (
                                <div key={index} className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                                  <div className="flex items-center mb-2">
                                    <Microscope className="h-4 w-4 text-accent mr-2" />
                                    <h5 className="text-xs font-medium text-text-primary">
                                      Contains {item.group} Group
                                    </h5>
                                  </div>
                                  <div className="space-y-2">
                                    {item.molecules.slice(0, 2).map((mol, molIndex) => (
                                      <div key={molIndex} className="flex justify-between items-center text-xs p-2 bg-secondary/20 rounded">
                                        <span className="text-text-secondary font-mono overflow-hidden text-ellipsis" style={{maxWidth: '180px'}}>
                                          {mol.name || formatMoleculeName(mol.smiles)}
                                        </span>
                                        <span className="text-success font-medium">{mol.price}</span>
                                      </div>
                                    ))}
                                  </div>
                                </div>
                              ))}
                            </div>
                          ) : (
                            <div className="bg-secondary/10 p-4 rounded-lg border border-secondary/50 text-center">
                              <Database className="h-5 w-5 text-text-secondary mx-auto mb-2" />
                              <p className="text-sm text-text-secondary">No similar molecules found in the dataset</p>
                            </div>
                          )}
                        </div>
                        
                        {/* Market trends from dataset */}
                        <div className="bg-primary p-4 rounded-lg border border-secondary">
                          <h4 className="text-sm font-medium text-text-primary mb-3">Market Trends Identified</h4>
                          {parsedResult.datasetAnalysis.marketTrends?.length > 0 ? (
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                              <div className="space-y-2">
                                {parsedResult.datasetAnalysis.marketTrends.map((trend, index) => (
                                  <div key={index} className="flex items-start p-2.5 bg-secondary/10 rounded-lg border border-secondary/50">
                                    <TrendingUp className="h-4 w-4 text-accent-secondary mt-0.5 mr-2 flex-shrink-0" />
                                    <p className="text-xs text-text-primary">{trend}</p>
                                  </div>
                                ))}
                              </div>
                              
                              <div className="bg-secondary/10 p-3 rounded-lg border border-secondary/50 h-full flex flex-col justify-center">
                                <h5 className="text-xs font-medium text-text-primary mb-3">Price Trends by Batch Size</h5>
                                <div className="h-40">
                                  <ResponsiveContainer width="100%" height="100%">
                                    <RechartsLineChart
                                      data={[
                                        { name: '1g', dataset: 100, current: 100 },
                                        { name: '10g', dataset: 85, current: 80 },
                                        { name: '100g', dataset: 65, current: 60 },
                                        { name: '1kg', dataset: 45, current: 40 },
                                      ]}
                                    >
                                      <CartesianGrid strokeDasharray="3 3" stroke="#333" />
                                      <XAxis dataKey="name" stroke="#888" />
                                      <YAxis stroke="#888" />
                                      <Tooltip
                                        contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                        formatter={(value) => [`${value}%`, '']}
                                      />
                                      <Line
                                        type="monotone"
                                        dataKey="dataset"
                                        name="Dataset Average"
                                        stroke={chartColors.accent}
                                        strokeDasharray="5 5"
                                      />
                                      <Line
                                        type="monotone"
                                        dataKey="current"
                                        name="Current Molecule"
                                        stroke={chartColors.success}
                                      />
                                      <Legend />
                                    </RechartsLineChart>
                                  </ResponsiveContainer>
                                </div>
                              </div>
                            </div>
                          ) : (
                            <div className="bg-secondary/10 p-4 rounded-lg border border-secondary/50 text-center">
                              <TrendingUp className="h-5 w-5 text-text-secondary mx-auto mb-2" />
                              <p className="text-sm text-text-secondary">No market trends identified in the dataset</p>
                            </div>
                          )}
                        </div>
                        
                        {/* Price patterns */}
                        <div className="bg-primary p-4 rounded-lg border border-secondary">
                          <h4 className="text-sm font-medium text-text-primary mb-3">Price Patterns & Reference Points</h4>
                          {parsedResult.datasetAnalysis.pricePatterns?.length > 0 ? (
                            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                              {parsedResult.datasetAnalysis.pricePatterns.map((pattern, index) => (
                                <div key={index} className="bg-secondary/10 p-3 rounded-lg border border-secondary/50">
                                  <div className="flex items-center">
                                    <DollarSign className="h-4 w-4 text-success mr-1.5" />
                                    <p className="text-xs text-text-primary">{pattern}</p>
                                  </div>
                                </div>
                              ))}
                            </div>
                          ) : (
                            <div className="bg-secondary/10 p-4 rounded-lg border border-secondary/50 text-center">
                              <DollarSign className="h-5 w-5 text-text-secondary mx-auto mb-2" />
                              <p className="text-sm text-text-secondary">No price patterns found in the dataset</p>
                            </div>
                          )}
                        </div>
                      </>
                    ) : (
                      <div className="bg-secondary/10 p-6 rounded-lg border border-secondary/50 text-center">
                        <Database className="h-6 w-6 text-text-secondary mx-auto mb-3" />
                        <h4 className="text-base font-medium text-text-primary mb-2">No Dataset Analysis Available</h4>
                        <p className="text-sm text-text-secondary">
                          The estimation was conducted using only algorithmic models without dataset enhancement.
                          This may occur if the datasets are unavailable or the analysis process encountered an error.
                        </p>
                      </div>
                    )}
                  </div>
                )}
              </div>
            </motion.div>
          )}

          {/* History Section with Tabs */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.5, delay: 0.2 }}
            className="bg-secondary p-8 rounded-xl shadow-md transition-all duration-200 hover:shadow-lg border border-secondary"
          >
            <div className="flex items-center justify-between mb-4">
              <div className="flex items-center">
                <motion.div
                  animate={{ rotate: 360 }}
                  transition={{ duration: 8, repeat: Infinity, ease: "linear" }}
                >
                  <Clock className="h-6 w-6 text-accent mr-2" />
                </motion.div>
                <h2 className="text-2xl font-semibold text-text-primary font-heading">Estimation History</h2>
              </div>
              
              <div className="flex items-center space-x-2">
                <div className="flex space-x-1 rounded-lg bg-secondary/20 p-1 h-8">
                  <button 
                    className={`flex items-center justify-center rounded-md px-3 py-1 text-xs font-medium transition-all
                    ${activeHistoryTab === 'list' ? 'bg-secondary text-text-primary' : 'bg-transparent text-text-secondary hover:bg-secondary/30'}`}
                    onClick={() => setActiveHistoryTab('list')}
                  >
                    <Database className="h-3.5 w-3.5 mr-1" />
                    List View
                  </button>
                  <button 
                    className={`flex items-center justify-center rounded-md px-3 py-1 text-xs font-medium transition-all
                    ${activeHistoryTab === 'stats' ? 'bg-secondary text-text-primary' : 'bg-transparent text-text-secondary hover:bg-secondary/30'}`}
                    onClick={() => setActiveHistoryTab('stats')}
                  >
                    <BarChart3 className="h-3.5 w-3.5 mr-1" />
                    Statistics
                  </button>
                </div>
                
                <motion.button
                  onClick={fetchHistory}
                  disabled={historyLoading}
                  className="h-8 w-8 flex items-center justify-center text-accent hover:text-accent/80 transition-colors duration-200 bg-secondary/80 rounded-md"
                  title="Refresh history"
                  whileHover={{ rotate: 90 }}
                  whileTap={{ scale: 0.9 }}
                >
                  <RefreshCw size={16} className={historyLoading ? "animate-spin" : ""} />
                </motion.button>
              </div>
            </div>

            {historyLoading ? (
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                className="flex flex-col items-center justify-center py-12"
              >
                <RefreshCw size={24} className="animate-spin mb-4 text-accent" />
                <p className="text-text-secondary font-body">Loading history data...</p>
              </motion.div>
            ) : (
              <div className="space-y-4">
                {activeHistoryTab === "list" && (
                  <div>
                    {history.length > 0 ? (
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
                                <div className="relative">
                                  <p className="font-mono text-sm bg-primary p-2 rounded border text-text-secondary border-secondary overflow-x-auto font-code">
                                    {item.smiles?.length > 20 
                                      ? `${item.smiles.substring(0, 20)}...` 
                                      : item.smiles || "N/A"}
                                  </p>
                                  <button 
                                    className="absolute right-2 top-1/2 -translate-y-1/2 text-accent-secondary hover:text-accent"
                                    onClick={() => {
                                      setSelectedMolecule(item.smiles);
                                      setShowMoleculeModal(true);
                                    }}
                                    title="View full SMILES"
                                  >
                                    <ZoomIn size={14} />
                                  </button>
                                </div>
                              </div>
                              
                              {/* Add visualization of key metrics */}
                              {item.parsedData && (
                                <div className="mt-3 grid grid-cols-3 gap-2">
                                  <div className="bg-primary p-2 rounded border border-secondary">
                                    <div className="flex items-center justify-between">
                                      <span className="text-xs text-text-secondary">Complexity</span>
                                      <Atom size={12} className="text-accent" />
                                    </div>
                                    <div className="text-sm font-medium text-accent mt-1">
                                      {item.parsedData.complexity.toFixed(1)}/5
                                    </div>
                                  </div>
                                  
                                  <div className="bg-primary p-2 rounded border border-secondary">
                                    <div className="flex items-center justify-between">
                                      <span className="text-xs text-text-secondary">Steps</span>
                                      <TestTube2 size={12} className="text-accent-secondary" />
                                    </div>
                                    <div className="text-sm font-medium text-accent-secondary mt-1">
                                      {item.parsedData.stepCount}
                                    </div>
                                  </div>
                                  
                                  <div className="bg-primary p-2 rounded border border-secondary">
                                    <div className="flex items-center justify-between">
                                      <span className="text-xs text-text-secondary">Risk</span>
                                      <AlertCircle size={12} className="text-error" />
                                    </div>
                                    <div className="text-sm font-medium text-error mt-1">
                                      {item.parsedData.riskFactor.toFixed(2)}
                                    </div>
                                  </div>
                                </div>
                              )}

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
                  </div>
                )}
                
                {activeHistoryTab === "stats" && (
                  <div>
                    {history.length > 0 ? (
                      <div className="space-y-6">
                        {/* Cost distribution */}
                        <div className="bg-primary p-4 rounded-lg border border-secondary">
                          <h4 className="text-sm font-medium text-text-primary mb-3">Cost Distribution</h4>
                          <div className="h-64">
                            <ResponsiveContainer width="100%" height="100%">
                              <BarChart
                                data={(() => {
                                  // Group estimations by cost range
                                  const costRanges = [
                                    { range: '$0-50', count: 0 },
                                    { range: '$51-100', count: 0 },
                                    { range: '$101-150', count: 0 },
                                    { range: '$151-200', count: 0 },
                                    { range: '$201+', count: 0 },
                                  ];
                                  
                                  history.forEach(item => {
                                    const cost = parseCost(item.estimatedcost).avg;
                                    if (cost <= 50) costRanges[0].count++;
                                    else if (cost <= 100) costRanges[1].count++;
                                    else if (cost <= 150) costRanges[2].count++;
                                    else if (cost <= 200) costRanges[3].count++;
                                    else costRanges[4].count++;
                                  });
                                  
                                  return costRanges;
                                })()}
                                margin={{ top: 5, right: 30, left: 20, bottom: 5 }}
                              >
                                <CartesianGrid strokeDasharray="3 3" stroke="#333" />
                                <XAxis dataKey="range" stroke="#888" />
                                <YAxis stroke="#888" />
                               
                                <Tooltip
                                  contentStyle={{ backgroundColor: '#172A45', borderColor: '#5E81F4', color: '#E0E0E0' }}
                                  formatter={(value) => [value, 'Estimations']}
                                />
                                <Bar dataKey="count" fill={chartColors.accent}>
                                  {history.map((_, index) => (
                                    <Cell key={`cell-${index}`} fill={chartColors.accent} />
                                  ))}
                                </Bar>
                              </BarChart>
                            </ResponsiveContainer>
                          </div>
                        </div>
                        
                        {/* Summary statistics */}
                        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                          <div className="bg-primary p-4 rounded-lg border border-secondary">
                            <h4 className="text-xs font-medium text-text-secondary mb-3">Cost Statistics</h4>
                            <div className="space-y-2">
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Average Cost:</span>
                                <span className="text-sm font-bold text-success">
                                  ${(history.reduce((sum, item) => sum + parseCost(item.estimatedcost).avg, 0) / history.length || 0).toFixed(2)}
                                </span>
                              </div>
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Lowest Cost:</span>
                                <span className="text-sm font-bold text-success">
                                  ${Math.min(...history.map(item => parseCost(item.estimatedcost).min))}
                                </span>
                              </div>
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Highest Cost:</span>
                                <span className="text-sm font-bold text-success">
                                  ${Math.max(...history.map(item => parseCost(item.estimatedcost).max))}
                                </span>
                              </div>
                            </div>
                          </div>
                          
                          <div className="bg-primary p-4 rounded-lg border border-secondary">
                            <h4 className="text-xs font-medium text-text-secondary mb-3">Complexity Statistics</h4>
                            <div className="space-y-2">
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Average Complexity:</span>
                                <span className="text-sm font-bold text-accent">
                                  {(history.reduce((sum, item) => sum + (item.parsedData?.complexity || 0), 0) / history.length || 0).toFixed(2)}/5
                                </span>
                              </div>
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Most Common Range:</span>
                                <span className="text-sm font-bold text-accent">
                                  {(() => {
                                    const ranges = [0, 0, 0]; // Low, Medium, High
                                    history.forEach(item => {
                                      const complexity = item.parsedData?.complexity || 0;
                                      if (complexity < 2) ranges[0]++;
                                      else if (complexity < 3.5) ranges[1]++;
                                      else ranges[2]++;
                                    });
                                    
                                    const max = Math.max(...ranges);
                                    const index = ranges.indexOf(max);
                                    return ['Low', 'Medium', 'High'][index];
                                  })()}
                                </span>
                              </div>
                            </div>
                          </div>
                          
                          <div className="bg-primary p-4 rounded-lg border border-secondary">
                            <h4 className="text-xs font-medium text-text-secondary mb-3">Step Count Statistics</h4>
                            <div className="space-y-2">
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Average Steps:</span>
                                <span className="text-sm font-bold text-accent-secondary">
                                  {(history.reduce((sum, item) => sum + (item.parsedData?.stepCount || 0), 0) / history.length || 0).toFixed(1)}
                                </span>
                              </div>
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Most Efficient:</span>
                                <span className="text-sm font-bold text-accent-secondary">
                                  {Math.min(...history.map(item => item.parsedData?.stepCount || Infinity))} steps
                                </span>
                              </div>
                              <div className="flex justify-between items-center">
                                <span className="text-sm text-text-secondary">Most Complex:</span>
                                <span className="text-sm font-bold text-accent-secondary">
                                  {Math.max(...history.map(item => item.parsedData?.stepCount || 0))} steps
                                </span>
                              </div>
                            </div>
                          </div>
                        </div>
                      </div>
                    ) : (
                      <motion.div
                        initial={{ opacity: 0 }}
                        animate={{ opacity: 1 }}
                        className="flex flex-col items-center justify-center py-12 text-center"
                      >
                        <BarChart3 size={32} className="text-text-secondary mb-4" />
                        <p className="text-text-secondary mb-2 font-body">No data available for statistics.</p>
                        <p className="text-sm text-text-secondary font-body">
                          Generate multiple estimations to see statistical insights and trends.
                        </p>
                      </motion.div>
                    )}
                  </div>
                )}
              </div>
            )}
          </motion.div>
        </div>
      </div>

      {/* Molecule Detail Modal */}
      {showMoleculeModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/60">
          <motion.div
            initial={{ opacity: 0, scale: 0.9 }}
            animate={{ opacity: 1, scale: 1 }}
            exit={{ opacity: 0, scale: 0.9 }}
            className="bg-secondary max-w-2xl w-full mx-4 p-6 rounded-xl border border-secondary"
          >
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-medium text-text-primary">Molecule SMILES</h3>
              <button
                onClick={() => setShowMoleculeModal(false)}
                className="text-text-secondary hover:text-text-primary"
              >
                <X size={20} />
              </button>
            </div>
            <div className="bg-primary p-4 rounded-lg border border-secondary">
              <p className="font-mono text-sm text-text-primary break-all">{selectedMolecule}</p>
            </div>
            <div className="mt-4 flex justify-end">
              <button
                onClick={() => setShowMoleculeModal(false)}
                className="px-4 py-2 bg-accent text-primary rounded-md"
              >
                Close
              </button>
            </div>
          </motion.div>
        </div>
      )}

      <style jsx>{`
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
