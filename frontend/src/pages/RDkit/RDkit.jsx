

import { useState, useEffect } from "react"
import {
  Container,
  TextField,
  Button,
  Typography,
  Box,
  CircularProgress,
  Alert,
  Paper,
  Grid,
  Collapse,
  IconButton,
  Tooltip,
  Menu,
  MenuItem,
  Chip,
  Divider,
  LinearProgress,
  Card,
  CardContent,
  CardHeader,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableContainer,
  FormControl,
  InputLabel,
  Select,
} from "@mui/material"
import {
  ExpandMore,
  ExpandLess,
  Science,
  CheckCircle,
  LocalHospital,
  Biotech,
  Warning,
  Info,
  Visibility,
  Sync,
  MedicalInformation,
  MedicalServices,
  GetApp,
} from "@mui/icons-material"
import axios from "axios"
import {
  Chart as ChartJS,
  RadialLinearScale,
  PointElement,
  LineElement,
  Filler,
  Tooltip as ChartTooltip,
  Legend,
} from "chart.js"
import { Radar } from "react-chartjs-2"
import toast, { Toaster } from "react-hot-toast"
import { useAuthStore } from "../../Store/auth.store.js"
import "./rdkit.css"
import jsPDF from "jspdf"

ChartJS.register(RadialLinearScale, PointElement, LineElement, Filler, Legend, ChartTooltip)

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
});


const timelineSteps = [
  "Fetching Symptoms",
  "Predicting Diseases",
  "Analyzing Disease Patterns",
  "Identifying Target Proteins",
  "Retrieving Ligands",
  "Fetching Ligands from Database",
  "Detecting Functional Groups",
  "Validating Reactants",
  "Running Chemical Reactions",
  "Computing ADMET Properties",
  "Finalizing Analysis",
]

const loadingMessages = [
  "Processing your symptoms to identify potential diseases...",
  "Analyzing symptom patterns to predict diseases...",
  "Performing detailed disease analysis...",
  "Identifying target proteins for therapeutic intervention...",
  "Retrieving ligand structures from the database...",
  "Fetching ligand SMILES from MongoDB...",
  "Analyzing ligand structures for functional groups...",
  "Validating reactants for chemical compatibility...",
  "Simulating reactions with RDChiral and RDKit...",
  "Computing ADMET properties for reactants and products...",
  "Generating and scoring reaction products...",
]

const toastWarn = (message, options = {}) => {
  toast(message, {
    ...options,
    style: {
      background: "#fff3e0",
      color: "#e65100",
      border: "1px solid #fb8c00",
    },
  })
}

// Storage keys for persistence
// const STORAGE_KEYS = {
//   SYMPTOMS: "rdkit_symptoms",
//   RESULT: "rdkit_result",
//   REACTION_RESULT: "rdkit_reaction_result",
//   SMILES_IMAGES: "rdkit_smiles_images",
//   EXPANDED: "rdkit_expanded",
//   SAVED_SYMPTOMS: "rdkit_saved_symptoms",
// }

// Helper functions for localStorage
// const saveToStorage = (key, data) => {
//   try {
//     localStorage.setItem(key, JSON.stringify(data))
//   } catch (error) {
//     console.warn("Failed to save to localStorage:", error)
//   }
// }

// const loadFromStorage = (key, defaultValue = null) => {
//   try {
//     const stored = localStorage.getItem(key)
//     return stored ? JSON.parse(stored) : defaultValue
//   } catch (error) {
//     console.warn("Failed to load from localStorage:", error)
//     return defaultValue
//   }
// }

// const clearStorage = () => {
//   try {
//     Object.values(STORAGE_KEYS).forEach((key) => {
//       localStorage.removeItem(key)
//     })
//   } catch (error) {
//     console.warn("Failed to clear localStorage:", error)
//   }
// }

const RDkit = () => {
  // Initialize state with persisted data
  // const [symptoms, setSymptoms] = useState(() => loadFromStorage(STORAGE_KEYS.SYMPTOMS, []))
  // const [loading, setLoading] = useState(false)
  // const [error, setError] = useState("")
  // const [result, setResult] = useState(() => loadFromStorage(STORAGE_KEYS.RESULT, null))
  // const [reactionResult, setReactionResult] = useState(() => loadFromStorage(STORAGE_KEYS.REACTION_RESULT, null))
  // // const [expanded, setExpanded] = useState(() => loadFromStorage(STORAGE_KEYS.EXPANDED, {}))
  // const [currentStep, setCurrentStep] = useState(0)
  // const [availableReactions, setAvailableReactions] = useState([])
  // const [reactionsLoading, setReactionsLoading] = useState(false)
  // const [anchorEl, setAnchorEl] = useState(null)
  // const [smilesImages, setSmilesImages] = useState(() => loadFromStorage(STORAGE_KEYS.SMILES_IMAGES, {}))
  // const [retryAttempts, setRetryAttempts] = useState({ disease: 0, protein: 0, reaction: 0, noReactions: 0 })
  // const [maxRetries] = useState(5)
  // const [savedSymptoms, setSavedSymptoms] = useState(() => loadFromStorage(STORAGE_KEYS.SAVED_SYMPTOMS, []))
  // const [symptomsLoading, setSymptomsLoading] = useState(false)
  // const [expanded, setExpanded] = useState(() => loadFromStorage(STORAGE_KEYS.EXPANDED, {
  //   "disease-prediction": false,
  //   "detailed-disease-analysis": false,
  //   "target-proteins": false,
  //   "target-ligands": false,
  //   "chemical-reactions": false,
  //   "failed-reactions": false,
  // }));



  const [hasGeneratedPDF, setHasGeneratedPDF] = useState(false);

  const [symptoms, setSymptoms] = useState([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState("")
  const [result, setResult] = useState(null)
  const [reactionResult, setReactionResult] = useState(null)
  const [currentStep, setCurrentStep] = useState(0)
  const [availableReactions, setAvailableReactions] = useState([])
  const [reactionsLoading, setReactionsLoading] = useState(false)
  const [anchorEl, setAnchorEl] = useState(null)
  const [smilesImages, setSmilesImages] = useState({})
  const [retryAttempts, setRetryAttempts] = useState({ disease: 0, protein: 0, reaction: 0, noReactions: 0 })
  const [maxRetries] = useState(5)
  const [savedSymptoms, setSavedSymptoms] = useState([])
  const [symptomsLoading, setSymptomsLoading] = useState(false)
  const [expanded, setExpanded] = useState({
    "disease-prediction": false,
    "detailed-disease-analysis": false,
    "target-proteins": false,
    "target-ligands": false,
    "chemical-reactions": false,
    "failed-reactions": false,
  })
  const { user } = useAuthStore()
  const userId = user?._id

const [comparisonResult, setComparisonResult] = useState(null);
const [comparisonLoading, setComparisonLoading] = useState(false);
const [showComparison, setShowComparison] = useState(false);

// Add this state to track the raw input value
const [symptomInput, setSymptomInput] = useState("")

const handleCompareWithDataset = async () => {
  if (!reactionResult || !result) {
    toast.error("No drug prediction results available for comparison.", { 
      id: "compare-error", 
      ...toastOptions 
    });
    return;
  }

  setComparisonLoading(true);
  setError("");
  
  try {
    toast("Analyzing generated compounds against reference dataset...", {
      id: "compare-start",
      duration: 3000,
      ...toastOptions,
    });

    const response = await axiosInstance.post(
      `/newdrug/compareDrugs/${userId}`,
      { reactionResult },
      {
        headers: { 
          Authorization: `Bearer ${localStorage.getItem("token")}` 
        },
        timeout: 120000, // 2 minutes timeout for complex analysis
      }
    );

    setComparisonResult(response.data.comparison);
    setShowComparison(true);
    
    toast.success(
      `Comparison completed! Overall accuracy: ${response.data.comparison.overall_accuracy}`, 
      { 
        id: "compare-success", 
        duration: 5000,
        ...toastOptions 
      }
    );

  } catch (err) {
    console.error("Comparison Error:", err.message);
    setError("Failed to compare drugs with dataset: " + (err.response?.data?.error || err.message));
    toast.error("Failed to compare drugs with dataset.", { 
      id: "compare-error", 
      ...toastOptions 
    });
  } finally {
    setComparisonLoading(false);
  }
};

  // Save state to localStorage whenever it changes
  // useEffect(() => {
  //   saveToStorage(STORAGE_KEYS.SYMPTOMS, symptoms)
  // }, [symptoms])

  // useEffect(() => {
  //   saveToStorage(STORAGE_KEYS.RESULT, result)
  // }, [result])

  // useEffect(() => {
  //   saveToStorage(STORAGE_KEYS.REACTION_RESULT, reactionResult)
  // }, [reactionResult])

  // useEffect(() => {
  //   saveToStorage(STORAGE_KEYS.SMILES_IMAGES, smilesImages)
  // }, [smilesImages])

  // useEffect(() => {
  //   saveToStorage(STORAGE_KEYS.EXPANDED, expanded)
  // }, [expanded])

  // useEffect(() => {
  //   saveToStorage(STORAGE_KEYS.SAVED_SYMPTOMS, savedSymptoms)
  // }, [savedSymptoms])

  // Show toast notification on component mount if data was restored
  useEffect(() => {
    const hasPersistedData = result || reactionResult || symptoms.length > 0
    if (hasPersistedData) {
      toast.success("Previous analysis data restored from your last session!", {
        id: "data-restored",
        duration: 4000,
        icon: "🔄", toastOptions
      })
    }
  }, []) // Only run on mount



  // Add this useEffect after your existing useEffect hooks


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

  // Fetch available reactions
  useEffect(() => {
    const fetchAvailableReactions = async () => {
      setReactionsLoading(true)
      try {
        const response = await axios.get("http://127.0.0.1:5001/api/reactions", { timeout: 30000 })
        setAvailableReactions(response.data.reactions || [])
        toast.success("Available reactions fetched successfully!", { id: "fetch-reactions", ...toastOptions })
      } catch (err) {
        console.error("Error fetching available reactions:", err.message)
        setError(err.response?.data?.error || "Error fetching available reactions.")
        toast.error("Failed to fetch available reactions.", {
          id: "fetch-reactions-error",
          ...toastOptions,
        });

      } finally {
        setReactionsLoading(false)
      }
    }

    fetchAvailableReactions()
  }, [])

  // Fetch saved symptoms
  useEffect(() => {
    if (!userId) return

    const fetchSavedSymptoms = async () => {
      setSymptomsLoading(true)
      try {
        const response = await axiosInstance.get(`/newdrug/symptoms/${userId}`, {
          headers: { Authorization: `Bearer ${localStorage.getItem("token")}` },
          timeout: 30000,
        })
        const fetchedSymptoms = response.data.symptoms || []
        // Validate saved symptoms
        const validSymptoms = fetchedSymptoms.filter(
          (group) =>
            Array.isArray(group) && group.length > 0 && group.every((s) => typeof s === "string" && s.trim() !== ""),
        )
        setSavedSymptoms(validSymptoms)
        if (fetchedSymptoms.length !== validSymptoms.length) {
          toastWarn("Some saved symptom groups were invalid and filtered out.", { id: "invalid-symptoms" })
        }
        toast.success("Saved symptoms fetched successfully!", { id: "fetch-symptoms", ...toastOptions })
      } catch (err) {
        console.error("Error fetching saved symptoms:", err.message)
        setError(err.response?.data?.message || "Error fetching saved symptoms.")
        toast.error("Failed to fetch saved symptoms.", { id: "fetch-symptoms-error" , ...toastOptions})
      } finally {
        setSymptomsLoading(false)
      }
    }

    fetchSavedSymptoms()
  }, [userId])

  useEffect(() => {
    if (!loading) return

    const interval = setInterval(() => {
      setCurrentStep((prev) => {
        if (prev < timelineSteps.length - 1) {
          return prev + 1
        }
        return prev
      })
    }, 2000)

    return () => clearInterval(interval)
  }, [loading])

  useEffect(() => {
    if (!loading) {
      if (error) {
        setCurrentStep(0)
      } else if (currentStep > 0) {
        setCurrentStep(timelineSteps.length - 1)
      }
    }
  }, [loading, error, currentStep])

  useEffect(() => {
    const fetchSmilesImages = async () => {
      if (!reactionResult) return

      const allSmiles = new Set()

      reactionResult.reactants?.forEach((reactant) => {
        if (reactant.smiles && reactant.smiles !== "Not applicable" && reactant.smiles.length > 1) {
          allSmiles.add(reactant.smiles)
        }
      })

      reactionResult.reactionResults?.forEach((reaction) => {
        reaction.products?.forEach((product) => {
          if (product.smiles && product.smiles.length > 1) {
            allSmiles.add(product.smiles)
          }
        })
      })

      const newImages = { ...smilesImages }
      for (const smiles of allSmiles) {
        if (!newImages[smiles]) {
          try {
            const response = await axios.get(
              `http://127.0.0.1:5001/api/smiles_to_image?smiles=${encodeURIComponent(smiles)}`,
              { timeout: 15000 },
            )
            newImages[smiles] = response.data.image
          } catch (err) {
            console.warn(`Failed to fetch image for SMILES: ${smiles} - ${err.message}`)
            newImages[smiles] = null
          }
        }
      }
      setSmilesImages(newImages)
    }

    fetchSmilesImages()
  }, [reactionResult, smilesImages])

  useEffect(() => {
    if (
      reactionResult &&
      reactionResult.reactionResults?.length === 0 &&
      retryAttempts.noReactions < maxRetries &&
      !loading
    ) {
      toastWarn(
        `No successful reactions found. Retrying drug prediction (Attempt ${retryAttempts.noReactions + 1}/${maxRetries})...`,
        { id: "no-reactions-retry" },
      )
      setTimeout(() => {
        handleSubmit({ preventDefault: () => { } })
        setRetryAttempts((prev) => ({ ...prev, noReactions: prev.noReactions + 1 }))
      }, 2000)
    } else if (
      reactionResult &&
      reactionResult.reactionResults?.length === 0 &&
      retryAttempts.noReactions >= maxRetries
    ) {
      toast.error(
        "Max retry attempts reached for finding reactions. Please try different symptoms or check the backend.",
        { id: "no-reactions-max-retries" , ...toastOptions}
      )
    }
  }, [reactionResult, retryAttempts.noReactions, loading, maxRetries])



const renderComparisonResults = () => {
  if (!comparisonResult) return null;

  return (
    <Collapse in={showComparison} timeout="auto" unmountOnExit>
      <Box sx={{ mt: 2 }}>
        <Card
          sx={{
            bgcolor: "#0A192F",
            border: "2px solid #70E000",
            color: "#E0E0E0",
            mb: 2,
          }}
        >
          <CardHeader
            title={
              <Box display="flex" alignItems="center" gap={1}>
                <Science sx={{ color: "#70E000" }} />
                <Typography variant="h6" sx={{ color: "#70E000" }}>
                  Drug Similarity Analysis Results
                </Typography>
                <Chip
                  label={`Accuracy: ${comparisonResult.overall_accuracy}`}
                  sx={{
                    bgcolor: "#70E000",
                    color: "#0A192F",
                    fontWeight: 600,
                  }}
                />
              </Box>
            }
            subheader={
              <Typography variant="body2" sx={{ color: "#A0A0A0", mt: 1 }}>
                Comprehensive comparison with reference drug dataset
              </Typography>
            }
          />
          
          <CardContent>
            {/* Overall Summary */}
            <Box sx={{ mb: 3 }}>
              <Typography variant="h6" sx={{ color: "#00F5D4", mb: 1 }}>
                Analysis Summary
              </Typography>
              <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 2 }}>
                {comparisonResult.comparison_summary}
              </Typography>
            </Box>

            {/* Individual Compound Analysis */}
            <Box sx={{ mb: 3 }}>
              <Typography variant="h6" sx={{ color: "#00F5D4", mb: 2 }}>
                Individual Compound Analysis ({comparisonResult.individual_compounds?.length || 0} compounds)
              </Typography>
              
              {comparisonResult.individual_compounds?.map((compound, index) => (
                <Card
                  key={index}
                  sx={{
                    bgcolor: "#172A45",
                    border: "1px solid #5E81F4",
                    mb: 2,
                  }}
                >
                  <CardContent>
                    <Box display="flex" justifyContent="space-between" alignItems="center" mb={1}>
                      <Typography variant="subtitle1" sx={{ color: "#70E000", fontWeight: 600 }}>
                        {compound.compound_id}
                      </Typography>
                      <Chip
                        label={`Similarity: ${compound.similarity_score}`}
                        size="small"
                        sx={{
                          bgcolor: compound.similarity_score >= "45%" ? "#70E000" : "#FFD700",
                          color: "#0A192F",
                          fontWeight: 600,
                        }}
                      />
                    </Box>
                    
                    <Typography variant="body2" sx={{ color: "#A0A0A0", mb: 1, fontFamily: 'monospace' }}>
                      SMILES: {compound.smiles}
                    </Typography>

                    <Grid container spacing={2}>
                      <Grid item xs={12} md={6}>
                        <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                          <strong>Property Match:</strong> {compound.property_match}
                        </Typography>
                        <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                          <strong>Structural Analysis:</strong> {compound.structural_analysis}
                        </Typography>
                      </Grid>
                      <Grid item xs={12} md={6}>
                        <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                          <strong>Drug-Likeness:</strong> {compound.drug_likeness_assessment}
                        </Typography>
                      </Grid>
                    </Grid>

                    {compound.key_similarities?.length > 0 && (
                      <Box sx={{ mt: 2 }}>
                        <Typography variant="body2" sx={{ color: "#00F5D4", mb: 1 }}>
                          <strong>Key Similarities:</strong>
                        </Typography>
                        <Box display="flex" flexWrap="wrap" gap={0.5}>
                          {compound.key_similarities.map((similarity, idx) => (
                            <Chip
                              key={idx}
                              label={similarity}
                              size="small"
                              sx={{
                                bgcolor: "#5E81F4",
                                color: "#E0E0E0",
                                fontSize: "0.75rem",
                              }}
                            />
                          ))}
                        </Box>
                      </Box>
                    )}
                  </CardContent>
                </Card>
              ))}
            </Box>

            {/* Dataset Analysis */}
            {comparisonResult.dataset_analysis && (
              <Box sx={{ mb: 3 }}>
                <Typography variant="h6" sx={{ color: "#00F5D4", mb: 2 }}>
                  Dataset Analysis
                </Typography>
                <Grid container spacing={2}>
                  <Grid item xs={12} md={6}>
                    <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                      <strong>Property Distributions:</strong> {comparisonResult.dataset_analysis.property_distributions}
                    </Typography>
                    <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                      <strong>Chemical Space Coverage:</strong> {comparisonResult.dataset_analysis.chemical_space_coverage}
                    </Typography>
                  </Grid>
                  <Grid item xs={12} md={6}>
                    <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                      <strong>Novelty Assessment:</strong> {comparisonResult.dataset_analysis.novelty_assessment}
                    </Typography>
                    <Typography variant="body2" sx={{ color: "#E0E0E0", mb: 1 }}>
                      <strong>Drug-Likeness Comparison:</strong> {comparisonResult.dataset_analysis.drug_likeness_comparison}
                    </Typography>
                  </Grid>
                </Grid>
              </Box>
            )}

            {/* Recommendations */}
            {comparisonResult.recommendations?.length > 0 && (
              <Box>
                <Typography variant="h6" sx={{ color: "#00F5D4", mb: 2 }}>
                  Optimization Recommendations
                </Typography>
                <Box component="ul" sx={{ color: "#E0E0E0", pl: 2 }}>
                  {comparisonResult.recommendations.map((recommendation, index) => (
                    <Typography component="li" key={index} variant="body2" sx={{ mb: 1 }}>
                      {recommendation}
                    </Typography>
                  ))}
                </Box>
              </Box>
            )}
          </CardContent>
        </Card>
      </Box>
    </Collapse>
  );
};

  const validateSymptoms = (symptomList) => {
    if (!Array.isArray(symptomList) || symptomList.length === 0) {
      return { valid: false, message: "Please enter at least three symptoms." }
    }

    if (symptomList.length < 3) {
      return { valid: false, message: "Please enter at least three symptoms to proceed with drug discovery." }
    }

    if (symptomList.some((s) => typeof s !== "string" || s.trim() === "")) {
      return { valid: false, message: "All symptoms must be non-empty strings." }
    }

    // Check for duplicate symptoms (case-insensitive)
    const normalizedSymptoms = symptomList.map((s) => s.trim().toLowerCase())
    const uniqueSymptoms = new Set(normalizedSymptoms)

    if (uniqueSymptoms.size !== symptomList.length) {
      return { valid: false, message: "Please remove duplicate symptoms. Each symptom should be entered only once." }
    }

    return { valid: true }
  }

  const axiosWithRetry = async (config, type, retriesLeft = maxRetries) => {
    try {
      const response = await axios(config)
      setRetryAttempts((prev) => ({ ...prev, [type]: 0 }))
      return response
    } catch (err) {
      if (retriesLeft <= 1 || !err.message.includes("timeout")) {
        throw err
      }

      setRetryAttempts((prev) => ({ ...prev, [type]: maxRetries - retriesLeft + 1 }))
      toastWarn(`Retrying ${type} request... Attempt ${maxRetries - retriesLeft + 1} of ${maxRetries}`, {
        id: `retry-${type}`,
      })
      await new Promise((resolve) => setTimeout(resolve, 2000))
      return axiosWithRetry(config, type, retriesLeft - 1)
    }
  }

  const handleSubmit = async (e) => {
    e.preventDefault()
    setError("")
    setResult(null)
    setReactionResult(null)
    setLoading(true)
    setCurrentStep(0)
    setRetryAttempts((prev) => ({ ...prev, disease: 0, protein: 0, reaction: 0 }))

    try {
      if (!userId) {
        throw new Error("User ID is missing. Please log in.")
      }

      const validationResult = validateSymptoms(symptoms)
      if (!validationResult.valid) {
        setError(validationResult.message)
        toast.error(validationResult.message, { id: "submit-error", ...toastOptions })
        setLoading(false)
        return
      }

      toast("This process may go for more than 1 Attempt to predict new drug due to complex computations.", {
        id: "process-info", duration: 5000,
        ...toastOptions,
      })

      let diseaseResponse
      try {
        diseaseResponse = await axiosWithRetry(
          {
            method: "post",
            url: `http://localhost:5000/api/newdrug/predictDisease/${userId}`,
            data: { symptoms },
            headers: { Authorization: `Bearer ${localStorage.getItem("token")}` },
            timeout: 60000,
          },
          "disease",
        )
      } catch (err) {
        console.error("Disease Prediction Error:", err.message)
        throw new Error("Failed to predict disease: " + (err.response?.data?.error || err.message))
      }

      const predictedDiseases = diseaseResponse.data?.predictedDiseases || []
      if (predictedDiseases.length > 0 && predictedDiseases[0].DiseaseMatchness === "0%") {
        setResult({
          disease: { predictedDiseases: [{ diseaseName: "No Match", DiseaseMatchness: "0%" }], DiseaseAnalysis: {} },
          proteins: { TargetProteins: [], TargetLigands: [] },
        })
        setLoading(false)
        return
      }

      if (predictedDiseases.length > 0 && predictedDiseases[0].DiseaseMatchness === "100%") {
        setResult({
          disease: diseaseResponse.data || { predictedDiseases: [], DiseaseAnalysis: {} },
          proteins: { TargetProteins: [], TargetLigands: [] },
        })
        setLoading(false)
        return
      }

      let proteinResponse
      try {
        proteinResponse = await axiosWithRetry(
          {
            method: "post",
            url: "http://localhost:5000/api/newdrug/predictTargetProtein",
            data: {},
            headers: { Authorization: `Bearer ${localStorage.getItem("token")}` },
            timeout: 60000,
          },
          "protein",
        )
      } catch (err) {
        console.error("Protein Prediction Error:", err.message)
        if (err.message.includes("timeout")) {
          throw new Error(
            "Protein prediction is taking longer than expected. Please check if the server is running and try again.",
          )
        }
        throw new Error("Failed to predict target proteins: " + (err.response?.data?.error || err.message))
      }

      setResult({
        disease: diseaseResponse.data || { predictedDiseases: [], DiseaseAnalysis: {} },
        proteins: proteinResponse.data || { TargetProteins: [], TargetLigands: [] },
      })

      let reactionResponse
      try {
        reactionResponse = await axiosWithRetry(
          {
            method: "get",
            url: "http://127.0.0.1:5001/api/react?include_admet=true",
            timeout: 30000,
          },
          "reaction",
        )
      } catch (err) {
        console.error("Reaction Error:", err.message)
        throw new Error("Failed to run reactions: " + (err.response?.data?.error || err.message))
      }

      setReactionResult(reactionResponse.data || { reactionResults: [], failedReactions: [], statistics: {} })
      toast.success("Analysis and reactions completed successfully!", { id: "submit-success",...toastOptions, })
    } catch (err) {
      console.error("Submission Error:", err.message)
      setError(err.message || "Error processing analysis. Please try again or check the backend server.")
      toast.error(err.message || "Error processing analysis.", { id: "submit-error",...toastOptions, })
    } finally {
      setLoading(false)
    }
  }

  const handleReset = () => {
    // Clear all state
    setSymptoms([]);
    setError("");
    setResult(null);
    setReactionResult(null);
    setLoading(false);
    setCurrentStep(0);
    setExpanded({});
    setAnchorEl(null);
    setSmilesImages({});
    setRetryAttempts({ disease: 0, protein: 0, reaction: 0, noReactions: 0 });
    setHasGeneratedPDF(false); // Allow new PDF generation after reset

    // Clear localStorage
    // clearStorage();

    toast.success("Form reset successfully! All data cleared from storage.", {
      id: "reset-success",
      icon: "🗑️",
      ...toastOptions,
    });
  };

  const handleToggle = (key) => {
    setExpanded((prev) => ({ ...prev, [key]: !prev[key] }))
  }

  const getFunctionalGroups = (ligandSmiles) => {
    if (ligandSmiles === "Not applicable") {
      return "N/A (Protein)"
    }
    if (!reactionResult || !reactionResult.reactants) {
      return "Run reactions to analyze functional groups"
    }
    const reactant = reactionResult.reactants.find((r) => r.smiles === ligandSmiles)
    if (!reactant || !Array.isArray(reactant.properties?.functional_groups)) {
      console.warn(`No functional groups found for ligand SMILES: ${ligandSmiles}`)
      return "Not identified"
    }
    return reactant.properties.functional_groups.join(", ") || "None"
  }

  const handleMenuOpen = (event) => {
    setAnchorEl(event.currentTarget)
  }

  const handleMenuClose = () => {
    setAnchorEl(null)
  }

  const deduplicateProducts = (products) => {
    const seenSmiles = new Set()
    return products.filter((product) => {
      if (seenSmiles.has(product.smiles)) {
        return false
      }
      seenSmiles.add(product.smiles)
      return true
    })
  }

  const prepareAdmetGraphData = (reactionIndex) => {
    if (!reactionResult || !reactionResult.reactionResults[reactionIndex]) return null

    const reaction = reactionResult.reactionResults[reactionIndex]
    const products = reaction.products

    const propertiesToShow = ["Non-Toxic", "Soluble", "Bioavailable", "hERG Safe"]
    const datasets = []

    products.forEach((product, idx) => {
      const admetSections = product.admet_properties?.find((section) => section.section === "ADMET")
      if (!admetSections) return

      const data = propertiesToShow.map((prop) => {
        const property = admetSections.properties.find((p) => p.name === prop)
        return property ? Number(property.prediction) : 0
      })

      datasets.push({
        label: `Product ${idx + 1} (${product.smiles.slice(0, 10)}...)`,
        data,
        backgroundColor: `rgba(255, 99, 132, 0.2)`,
        borderColor: `rgba(255, 99, 132, 1)`,
        borderWidth: 2,
        pointBackgroundColor: `rgba(255, 99, 132, 1)`,
        pointBorderColor: "#fff",
        pointHoverBackgroundColor: "#fff",
        pointHoverBorderColor: `rgba(255, 99, 132, 1)`,
      })
    })

    return {
      labels: propertiesToShow,
      datasets,
    }
  }

  const handleSymptomGroupSelect = (selectedSymptoms) => {
    if (Array.isArray(selectedSymptoms) && selectedSymptoms.length > 0) {
      setSymptoms(selectedSymptoms)
      toast.success(`Loaded symptom group: ${selectedSymptoms.join(", ")}`, {
        id: "symptom-group-loaded",...toastOptions,
      })
    } else {
      setSymptoms([])
      toast.error("Selected symptom group is invalid or empty.", { id: "symptom-group-error",...toastOptions, })
    }
  }

  const handleExample = (exampleSymptoms) => {
    setSymptoms(exampleSymptoms)
  }

  const generatePDF = async () => {
    const doc = new jsPDF({ orientation: "portrait", unit: "mm", format: "a4" });
    let yPosition = 20;
    const pageHeight = 280;
    const margin = 20;

    // Helper functions
    const checkPageBreak = (requiredSpace = 20) => {
      if (yPosition + requiredSpace > pageHeight) {
        doc.addPage();
        yPosition = 20;
      }
    };

    const addTitle = (title, fontSize = 16, color = [0, 100, 255]) => {
      checkPageBreak(15);
      doc.setFontSize(fontSize);
      doc.setTextColor(...color);
      doc.text(title, margin, yPosition);
      yPosition += fontSize === 20 ? 15 : 10;
    };

    const addText = (text, fontSize = 12, color = [40, 40, 40], indent = 0) => {
      checkPageBreak(10);
      doc.setFontSize(fontSize);
      doc.setTextColor(...color);
      const splitText = doc.splitTextToSize(text, 170 - indent);
      doc.text(splitText, margin + indent, yPosition);
      yPosition += splitText.length * 6;
    };

    const addTable = (headers, rows) => {
      checkPageBreak(30);
      const startY = yPosition;
      const cellHeight = 8;
      const cellWidth = (170 - margin) / headers.length;

      // Draw headers
      doc.setFontSize(10);
      doc.setTextColor(255, 255, 255);
      doc.setFillColor(94, 129, 244);
      doc.rect(margin, startY, 170 - margin, cellHeight, "F");

      headers.forEach((header, i) => {
        doc.text(header, margin + i * cellWidth + 2, startY + 6);
      });

      // Draw rows
      doc.setTextColor(40, 40, 40);
      rows.forEach((row, rowIndex) => {
        const rowY = startY + (rowIndex + 1) * cellHeight;
        if (rowIndex % 2 === 0) {
          doc.setFillColor(245, 245, 245);
          doc.rect(margin, rowY, 170 - margin, cellHeight, "F");
        }

        row.forEach((cell, cellIndex) => {
          doc.text(String(cell), margin + cellIndex * cellWidth + 2, rowY + 6);
        });
      });

      yPosition = startY + cellHeight * (rows.length + 1) + 10;
    };

    const addImage = async (base64Image, caption, width = 80, height = 60) => {
      if (!base64Image) {
        addText(`[Image: ${caption}]`, 10, [150, 150, 150]);
        return;
      }
      try {
        checkPageBreak(height + 15);
        await loadImage(base64Image); // Verify image is loadable
        doc.addImage(
          `data:image/png;base64,${base64Image}`,
          "PNG",
          margin,
          yPosition,
          width,
          height,
          undefined,
          "FAST"
        );
        yPosition += height + 5;
        if (caption) {
          addText(caption, 10, [100, 100, 100]);
        }
        yPosition += 5;
      } catch (error) {
        console.warn("Failed to add image to PDF:", error);
        addText(`[Failed to load image: ${caption}]`, 10, [150, 150, 150]);
      }
    };

    const addAdmetChart = async (chartData, title) => {
      try {
        checkPageBreak(120);
        const canvas = document.createElement("canvas");
        canvas.width = 500;
        canvas.height = 400;
        const ctx = canvas.getContext("2d");

        ctx.fillStyle = "#2d3748";
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        const chartInstance = new ChartJS(ctx, {
          type: "radar",
          data: chartData,
          options: {
            responsive: false,
            animation: false,
            plugins: {
              legend: {
                position: "top",
                labels: {
                  color: "#ffffff",
                  font: { size: 14 },
                },
              },
              title: {
                display: true,
                text: title,
                color: "#ffffff",
                font: { size: 16 },
              },
            },
            scales: {
              r: {
                beginAtZero: true,
                min: 0,
                max: 100,
                ticks: {
                  color: "#ffffff",
                  backdropColor: "transparent",
                },
                pointLabels: {
                  color: "#ffffff",
                  font: { size: 12 },
                },
                grid: { color: "#4a5568" },
                angleLines: { color: "#4a5568" },
              },
            },
            elements: {
              line: { borderWidth: 3 },
              point: { radius: 6, borderWidth: 2 },
            },
          },
        });

        await new Promise((resolve) => setTimeout(resolve, 100));
        const chartImage = canvas.toDataURL("image/png").split(",")[1];
        await addImage(chartImage, title, 140, 110);

        chartInstance.destroy();
        canvas.remove();
      } catch (error) {
        console.warn("Failed to add ADMET chart to PDF:", error);
        addText(`[ADMET Chart: ${title}]`, 10, [150, 150, 150]);
      }
    };

    // Title Page
    doc.setFontSize(24);
    doc.setTextColor(40, 40, 40);
    doc.text("Drug Discovery Analysis Report", margin, yPosition);
    yPosition += 20;

    doc.setFontSize(14);
    doc.text(`Generated on: ${new Date().toLocaleString()}`, margin, yPosition);
    yPosition += 10;

    doc.setFontSize(12);
    doc.text(`User: ${user?.name || "Unknown"}`, margin, yPosition);
    yPosition += 20;

    // INPUT SYMPTOMS SECTION
    if (symptoms.length > 0) {
      addTitle("Input Symptoms", 20, [0, 150, 0]);
      addText(`Total Symptoms Entered: ${symptoms.length}`, 12, [40, 40, 40], 5);
      yPosition += 5;

      addTitle("Symptoms List:", 14, [0, 100, 0]);
      symptoms.forEach((symptom, index) => {
        addText(`${index + 1}. ${symptom}`, 12, [40, 40, 40], 10);
      });

      const validation = validateSymptoms(symptoms);
      if (validation.valid) {
        addText("✓ Symptoms validation: PASSED (3+ unique symptoms)", 11, [0, 150, 0], 5);
      } else {
        addText(`✗ Symptoms validation: ${validation.message}`, 11, [200, 0, 0], 5);
      }

      yPosition += 15;
    }

    // DISEASE PREDICTION RESULTS
    if (result?.disease?.predictedDiseases) {
      addTitle("Disease Prediction Results", 20, [200, 0, 0]);

      result.disease.predictedDiseases.forEach((disease, index) => {
        addTitle(`${index + 1}. ${disease.diseaseName || "Unknown Disease"}`, 16, [0, 150, 0]);
        addText(`Match Confidence: ${disease.DiseaseMatchness || "N/A"}`, 12, [40, 40, 40], 10);

        if (disease.diseaseCautions?.length > 0) {
          addText(`Cautions: ${disease.diseaseCautions.join(", ")}`, 11, [200, 0, 0], 10);
        }

        if (disease.diseaseSymptoms?.length > 0) {
          addText(`Symptoms: ${disease.diseaseSymptoms.join(", ")}`, 11, [40, 40, 40], 10);
        }

        if (disease.diseaseTreatments?.length > 0) {
          addText(`Treatments: ${disease.diseaseTreatments.join(", ")}`, 11, [0, 100, 0], 10);
        }

        if (disease.diseaseDescription) {
          addText(`Description: ${disease.diseaseDescription}`, 10, [60, 60, 60], 10);
        }

        yPosition += 10;
      });

      // Detailed Disease Analysis
      if (
        result.disease.predictedDiseases[0]?.DiseaseMatchness !== "0%" &&
        Object.entries(result.disease?.DiseaseAnalysis || {}).length > 0
      ) {
        addTitle("Detailed Disease Analysis", 18, [0, 100, 150]);

        Object.entries(result.disease.DiseaseAnalysis).forEach(([key, items]) => {
          addTitle(key.replace(/([A-Z])/g, " $1").trim(), 14, [0, 100, 150]);

          items.forEach((item, idx) => {
            addText(`${idx + 1}. ${item.summary || "No summary available"}`, 11, [40, 40, 40], 15);
            if (item.source) {
              addText(`Source: ${item.source}`, 10, [100, 100, 100], 20);
            }
            if (item.url) {
              addText(`URL: ${item.url}`, 10, [0, 0, 200], 20);
            }
            yPosition += 5;
          });
        });
      }
    }

    // TARGET PROTEINS
    if (result?.proteins?.TargetProteins?.length > 0) {
      checkPageBreak(30);
      addTitle("Target Proteins", 20, [150, 0, 150]);

      result.proteins.TargetProteins.forEach((protein, index) => {
        addTitle(`${index + 1}. ${protein.proteinName || "Unknown Protein"}`, 16, [150, 0, 150]);
        addText(`Function: ${protein.proteinFunction || "N/A"}`, 12, [40, 40, 40], 10);
        addText(`UniProt ID: ${protein.proteinUniport || "N/A"}`, 12, [40, 40, 40], 10);
        addText(`Description: ${protein.ProtienDiscription || "N/A"}`, 12, [40, 40, 40], 10);

        if (protein.proteinDetailedDiscription) {
          addText(`Detailed Information: ${protein.proteinDetailedDiscription}`, 11, [60, 60, 60], 10);
        }

        if (protein.proteinStructure) {
          addText(`Structure: ${protein.proteinStructure}`, 11, [40, 40, 40], 10);
        }

        if (protein.proteinPathways?.length > 0) {
          addText(`Pathways: ${protein.proteinPathways.join(", ")}`, 11, [40, 40, 40], 10);
        }

        if (protein.proteinInteractions?.length > 0) {
          addText(`Interactions: ${protein.proteinInteractions.join(", ")}`, 11, [40, 40, 40], 10);
        }

        yPosition += 15;
      });
    }

    // TARGET LIGANDS
    if (result?.proteins?.TargetLigands?.length > 0) {
      checkPageBreak(30);
      addTitle("Target Ligands", 20, [200, 100, 0]);

      for (const [index, ligand] of result.proteins.TargetLigands.entries()) {
        addTitle(`${index + 1}. ${ligand.ligandName || "Unknown Ligand"}`, 16, [200, 100, 0]);
        addText(`Function: ${ligand.ligandFunction || "N/A"}`, 12, [40, 40, 40], 10);
        addText(`DrugBank ID: ${ligand.ligandDrugBankID || "N/A"}`, 12, [40, 40, 40], 10);

        if (ligand.LigandSmile && ligand.LigandSmile !== "Not applicable") {
          addText(`SMILES: ${ligand.LigandSmile}`, 11, [40, 40, 40], 10);

          // Add molecular structure image
          if (smilesImages[ligand.LigandSmile]) {
            await addImage(
              smilesImages[ligand.LigandSmile],
              `Structure of ${ligand.ligandName || "Ligand"}`,
              100,
              75
            );
          } else {
            addText(`[No image available for ${ligand.ligandName}]`, 10, [150, 150, 150]);
          }

          const functionalGroups = getFunctionalGroups(ligand.LigandSmile);
          addText(`Functional Groups: ${functionalGroups}`, 11, [40, 40, 40], 10);

          // Add ADMET Properties for ligands
          if (reactionResult?.reactants) {
            const reactant = reactionResult.reactants.find((r) => r.smiles === ligand.LigandSmile);
            if (reactant?.admet_properties?.length > 0) {
              addTitle("ADMET Properties", 14, [0, 100, 100]);

              reactant.admet_properties.forEach((section) => {
                addText(section.section, 12, [0, 80, 80], 15);

                const headers = ["Property", "Prediction", "Units"];
                const rows = section.properties.map((prop) => [
                  prop.name,
                  prop.prediction,
                  prop.units || "-",
                ]);

                addTable(headers, rows);
              });
            }
          }
        } else {
          addText("SMILES: Not applicable (Protein)", 11, [150, 150, 150], 10);
        }

        if (ligand.LigandDiscription) {
          addText(`Description: ${ligand.LigandDiscription}`, 11, [60, 60, 60], 10);
        }

        if (ligand.ligandMechanism) {
          addText(`Mechanism: ${ligand.ligandMechanism}`, 11, [40, 40, 40], 10);
        }

        if (ligand.ligandSideEffects?.length > 0) {
          addText(`Side Effects: ${ligand.ligandSideEffects.join(", ")}`, 11, [200, 0, 0], 10);
        }

        yPosition += 15;
      }
    }

    // CHEMICAL REACTION RESULTS
    if (reactionResult?.reactionResults?.length > 0) {
      checkPageBreak(30);
      addTitle("Chemical Reaction Results", 20, [200, 0, 0]);

      for (const [index, reaction] of reactionResult.reactionResults.entries()) {
        addTitle(`Reaction ${index + 1}: ${reaction.reactionType || "Unknown Reaction"}`, 16, [200, 0, 0]);
        addText(`Confidence: ${(reaction.confidence * 100).toFixed(1)}%`, 12, [40, 40, 40], 10);
        addText(`Description: ${reaction.description || "No description available"}`, 12, [40, 40, 40], 10);
        addText(`Reactants: ${reaction.reactants?.join(", ") || "N/A"}`, 12, [40, 40, 40], 10);

        // Products
        if (reaction.products?.length > 0) {
          const mainProducts = deduplicateProducts(reaction.products.slice(0, 5));

          for (const [prodIndex, product] of mainProducts.entries()) {
            addTitle(`Product ${prodIndex + 1}`, 14, [0, 150, 0]);

            if (product.smiles) {
              addText(`SMILES: ${product.smiles}`, 11, [40, 40, 40], 15);

              // Add product structure image
              if (smilesImages[product.smiles]) {
                await addImage(
                  smilesImages[product.smiles],
                  `Product ${prodIndex + 1} Structure`,
                  100,
                  75
                );
              } else {
                addText(`[No image available for Product ${prodIndex + 1}]`, 10, [150, 150, 150]);
              }
            }

            const properties = [
              ["Molecular Weight", `${product.molecular_weight?.toFixed(2) || "N/A"} g/mol`],
              ["LogP", product.logP?.toFixed(2) || "N/A"],
              ["TPSA", `${product.tpsa?.toFixed(2) || "N/A"} Å²`],
              ["H-Donors", product.num_h_donors || 0],
              ["H-Acceptors", product.num_h_acceptors || 0],
              ["Rotatable Bonds", product.num_rotatable_bonds || 0],
              ["Aromatic Rings", product.num_aromatic_rings || 0],
            ];

            addTitle("Molecular Properties", 12, [0, 100, 100]);
            properties.forEach(([prop, value]) => {
              addText(`${prop}: ${value}`, 10, [40, 40, 40], 20);
            });

            if (product.functional_groups?.length > 0) {
              addText(`Functional Groups: ${product.functional_groups.join(", ")}`, 11, [40, 40, 40], 15);
            }

            if (product.drug_likeness) {
              addTitle("Drug-Likeness Properties", 12, [0, 100, 100]);
              Object.entries(product.drug_likeness).forEach(([key, value]) => {
                addText(`${key}: ${value}`, 10, [40, 40, 40], 20);
              });
            }

            if (product.admet_properties?.length > 0) {
              addTitle("ADMET Properties", 12, [0, 100, 100]);

              product.admet_properties.forEach((section) => {
                addText(section.section, 11, [0, 80, 80], 20);

                const headers = ["Property", "Prediction", "Units"];
                const rows = section.properties.map((prop) => [
                  prop.name,
                  prop.prediction,
                  prop.units || "-",
                ]);

                addTable(headers, rows);
              });
            }
            yPosition += 10;
          }

          const chartData = prepareAdmetGraphData(index);
          if (chartData) {
            await addAdmetChart(chartData, `ADMET Properties Comparison - Reaction ${index + 1}`);
          }
        }

        if (reaction.statistics) {
          addTitle("Reaction Statistics", 14, [100, 0, 100]);
          Object.entries(reaction.statistics).forEach(([key, value]) => {
            addText(`${key}: ${value}`, 11, [40, 40, 40], 15);
          });
        }

        yPosition += 15;
      }

      // Failed Reactions Summary
      if (reactionResult.failedReactions?.length > 0) {
        checkPageBreak(30);
        addTitle("Failed Reactions Summary", 16, [200, 0, 0]);

        reactionResult.failedReactions.forEach((fail, index) => {
          addText(`${index + 1}. ${fail.reactionType || "Unknown Reaction"}`, 12, [200, 0, 0], 5);
          addText(`Reactants: ${fail.reactants?.join(", ") || "N/A"}`, 11, [40, 40, 40], 15);
          addText(`Reason: ${fail.reason || "Unknown error"}`, 11, [150, 0, 0], 15);
          yPosition += 8;
        });
      }

      // Overall Statistics
      if (reactionResult.statistics) {
        checkPageBreak(40);
        addTitle("Overall Reaction Statistics", 16, [100, 0, 100]);

        const statsData = [
          ["Mean Molecular Weight", `${reactionResult.statistics.mean_mw?.toFixed(2) || "N/A"} g/mol`],
          ["Standard Deviation MW", `${reactionResult.statistics.std_mw?.toFixed(2) || "N/A"} g/mol`],
          ["Minimum Molecular Weight", `${reactionResult.statistics.min_mw?.toFixed(2) || "N/A"} g/mol`],
          ["Maximum Molecular Weight", `${reactionResult.statistics.max_mw?.toFixed(2) || "N/A"} g/mol`],
          ["Total Products Generated", reactionResult.statistics.total_products || 0],
          ["Successful Reactions", reactionResult.reactionResults?.length || 0],
          ["Failed Reactions", reactionResult.failedReactions?.length || 0],
        ];

        const headers = ["Statistic", "Value"];
        addTable(headers, statsData);
      }
    }

    // Footer
    checkPageBreak(20);
    doc.setFontSize(10);
    doc.setTextColor(100, 100, 100);
    doc.text("Generated by Drug Discovery Pipeline - Complete Analysis Report", margin, yPosition);
    doc.text(`Total Pages: ${doc.internal.getNumberOfPages()}`, margin, yPosition + 8);

    // Save the PDF
    const fileName = `Drug_Discovery_Analysis_${new Date().toISOString().split("T")[0]}.pdf`;
    doc.save(fileName);
    toast.success("Complete PDF report with all data, structures, and ADMET charts generated successfully!", {
      id: "pdf-success",
      duration: 4000,
      ...toastOptions,
    });
  };

  // Check if disease match is 100% to hide certain sections
  const is100PercentMatch = result?.disease?.predictedDiseases?.[0]?.DiseaseMatchness === "100%"
  const is0PercentMatch = result?.disease?.predictedDiseases?.[0]?.DiseaseMatchness === "0%"

  // Add this state at the top level of the RDkit component, alongside other useState declarations

  const loadImage = (base64) => {
    return new Promise((resolve, reject) => {
      const img = new Image();
      img.src = `data:image/png;base64,${base64}`;
      img.onload = () => resolve(img);
      img.onerror = () => reject(new Error('Image loading failed'));
    });
  };

  // Replace the previous useEffect with this corrected version
  useEffect(() => {
    if (
      !loading &&
      !error &&
      result &&
      reactionResult &&
      symptoms.length >= 3 &&
      !hasGeneratedPDF &&
      result.disease?.predictedDiseases?.[0]?.DiseaseMatchness !== "0%"
    ) {
      const timer = setTimeout(async () => {
        try {
          await generatePDF();
          setHasGeneratedPDF(true);
          toast.success("PDF report with molecular structures automatically generated and downloaded!", {
            id: "auto-pdf-success",
            duration: 4000,
            icon: "📄",
            ...toastOptions,
          });
        } catch (err) {
          console.error("Error generating PDF:", err);
          toast.error("Failed to generate PDF automatically.", {
            id: "auto-pdf-error",
            duration: 4000,
            ...toastOptions,
          });
        }
      }, 2000);

      return () => clearTimeout(timer);
    }
  }, [loading, error, result, reactionResult, symptoms, hasGeneratedPDF, smilesImages]);

  // Function to render section headers with dropdown functionality
  const renderSectionHeader = (title, key, count = null, color = "#5E81F4") => (
    <Box
      sx={{
        display: "flex",
        alignItems: "center",
        cursor: "pointer",
        p: 2,
        borderRadius: "4px",
        "&:hover": { bgcolor: "rgba(94, 129, 244, 0.1)" },
        transition: "background-color 0.2s ease",
      }}
      onClick={() => handleToggle(key)}
    >
      <Typography
        variant="h6"
        sx={{
          fontWeight: 600,
          flex: 1,
          color: "#E0E0E0",
          fontFamily: "'Inter', 'Barlow', sans-serif",
        }}
      >
        {title}
      </Typography>
      {count !== null && (
        <Chip
          label={count}
          size="small"
          sx={{
            mr: 1,
            bgcolor: color,
            color: "#0A192F",
            fontWeight: 600,
            fontFamily: "'Lato', sans-serif",
          }}
        />
      )}
      <IconButton size="small" sx={{ color: "#00F5D4" }}>
        {expanded[key] ? <ExpandLess /> : <ExpandMore />}
      </IconButton>
    </Box>
  )

  return (
    <Box sx={{ minHeight: "100vh", bgcolor: "#0A192F", py: 4, color: "#E0E0E0" }}>
      <Toaster position="top-right" />
      <Container maxWidth="xl">
        <Paper elevation={2} sx={{ p: 4, mb: 4, borderRadius: 2, bgcolor: "#172A45" }}>
          <Box sx={{ textAlign: "center", mb: 4 }}>
            <Typography
              variant="h4"
              sx={{ fontWeight: 700, color: "#00F5D4", mb: 1, fontFamily: "'Inter', 'Barlow', sans-serif" }}
            >
              Drug Discovery Pipeline
            </Typography>
            <Typography
              variant="subtitle1"
              sx={{ color: "#A0A0A0", maxWidth: 800, mx: "auto", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
            >
              Advanced New Drug Discovery Pipeline powered by <span style={{ color: "#00F5D4" ,fontWeight: 300 }}>Gemini, RDKit, Pubchem and many more Datasets</span> with  <span style={{ color: "#00F5D4" ,fontWeight: 600 }}>Accuracy of 46.7%</span>
            </Typography>
            {/* {(result || reactionResult || symptoms.length > 0) && (
              <Alert
                severity="info"
                sx={{
                  mt: 2,
                  bgcolor: "#0A192F",
                  border: "1px solid #00F5D4",
                  color: "#E0E0E0",
                  maxWidth: 600,
                  mx: "auto",
                }}
                icon={<Info sx={{ color: "#00F5D4" }} />}
              >
                <Typography variant="body2" sx={{ fontFamily: "'Roboto', 'Open Sans', sans-serif" }}>
                  📊 Data persisted from previous session. Results will be preserved on page refresh until you click
                  Reset.
                </Typography>
              </Alert>
            )} */}
          </Box>

          <Card variant="outlined" sx={{ mb: 4, borderRadius: 2, bgcolor: "#172A45", borderColor: "#00F5D4" }}>
            <CardContent sx={{ p: 3 }}>
              <Box component="form" onSubmit={handleSubmit}>
                <Grid container spacing={2} alignItems="center">
                  <Grid item xs={12} md={6}> {/* changed md={6} to md={12} */}
<TextField
  fullWidth
  label="Enter Symptoms"
  value={symptomInput}
  onChange={(e) => {
    const inputValue = e.target.value;
    
    // Prevent input if exceeds 500 characters
    if (inputValue.length > 500) {
      toast.error("Maximum 500 characters allowed", { id: "max-length" });
      return;
    }
    
    // Update the raw input value
    setSymptomInput(inputValue);
    
    // Process symptoms array for validation and storage
    const processedSymptoms = inputValue
      ? inputValue
          .split(",")
          .map((s) => s.trim())
          .filter((s) => s !== "")
      : [];
    
    setSymptoms(processedSymptoms);
  }}
  onBlur={() => {
    // Clean up the input on blur (optional)
    const cleanedInput = symptoms
      .filter((s) => s.trim() !== "")
      .join(", ");
    setSymptomInput(cleanedInput);
  }}
  placeholder="Enter symptoms separated by commas (e.g., fever, cough, fatigue)"
  type="text"
  name="symptoms"
  variant="outlined"
  className="custom-textfield"
  sx={{
    "& .MuiOutlinedInput-root": {
      bgcolor: "#0A192F",
      color: "#E0E0E0",
      "& fieldset": { borderColor: "#5E81F4" },
      "&:hover fieldset": { borderColor: "#00F5D4" },
      "&.Mui-focused fieldset": { borderColor: "#00F5D4" },
    },
    "& .MuiOutlinedInput-input": {
      cursor: "text !important",
    },
    "& .MuiInputLabel-root": {
      color: "#A0A0A0",
      fontFamily: "'Lato', sans-serif",
    },
    "& .MuiInputLabel-root.Mui-focused": {
      color: "#00F5D4",
    },
  }}
  InputProps={{
    startAdornment: <LocalHospital sx={{ mr: 1, color: "#A0A0A0" }} />,
  }}
  // helperText={`${symptomInput.length}/500 characters | ${symptoms.length} symptoms entered`}
/>
                  </Grid>

                  <Grid item xs={12} md={6}>
                    <Box sx={{ display: "flex", gap: 2, flexWrap: "wrap" }}>
                      <Button
                        type="submit"
                        variant="contained"
                        size="large"
                        disabled={loading || symptoms.join(", ").length > 500}
                        startIcon={loading ? <CircularProgress size={20} /> : <Biotech />}
                        sx={{
                          bgcolor: "#00F5D4",
                          color: "#0A192F",
                          "&:hover": { bgcolor: "#5E81F4" },
                          textTransform: "none",
                          fontFamily: "'Lato', sans-serif",
                          fontWeight: 600,
                        }}
                      >
                        Predict New Drug
                      </Button>
                      <Button
                        variant="outlined"
                        size="large"
                        onClick={handleReset}
                        disabled={loading}
                        sx={{
                          textTransform: "none",
                          borderColor: "#5E81F4",
                          color: "#E0E0E0",
                          "&:hover": { borderColor: "#00F5D4", color: "#00F5D4" },
                          fontFamily: "'Lato', sans-serif",
                          fontWeight: 600,
                        }}
                      >
                        Reset
                      </Button>
                      <Button
                        variant="outlined"
                        size="large"
                        onClick={() =>
                          handleExample([
                            "brittle toenails",
                            "cracking knuckles",
                            "ear blockage",
                            "fluttering heartbeat",
                            "sudden confusion",
                          ])
                        }
                        disabled={loading}
                        sx={{
                          textTransform: "none",
                          borderColor: "#5E81F4",
                          color: "#E0E0E0",
                          "&:hover": { borderColor: "#00F5D4", color: "#00F5D4" },
                          fontFamily: "'Lato', sans-serif",
                          fontWeight: 600,
                        }}
                      >
                        Example
                      </Button>
                    </Box>
                  </Grid>
                </Grid>
              </Box>
            </CardContent>
          </Card>

          {loading && (
            <Card
              variant="outlined"
              sx={{
                mb: 4,
                borderRadius: 4,
                bgcolor: "#172A45",
                borderColor: "#5E81F4",
                boxShadow: "0px 8px 24px rgba(0, 0, 0, 0.2)",
                overflow: "hidden",
                position: "relative",
                "&:before": {
                  content: '""',
                  position: "absolute",
                  top: 0,
                  left: 0,
                  right: 0,
                  height: 4,
                  background: "linear-gradient(90deg, #00F5D4 0%, #5E81F4 50%, #00F5D4 100%)",
                  animation: "shimmer 3s infinite linear",
                  "@keyframes shimmer": {
                    "0%": { transform: "translateX(-100%)" },
                    "100%": { transform: "translateX(100%)" },
                  },
                },
              }}
            >
              <CardHeader
                title={
                  <Typography
                    variant="h5"
                    sx={{
                      fontWeight: 700,
                      color: "#E0E0E0",
                      fontFamily: "'Inter', sans-serif",
                      letterSpacing: "0.5px",
                      position: "relative",
                      display: "inline-block",
                      "&:after": {
                        content: '""',
                        position: "absolute",
                        bottom: -8,
                        left: 0,
                        width: "100%",
                        height: 2,
                        background: "linear-gradient(90deg, #00F5D4 0%, transparent 100%)",
                        borderRadius: 2,
                      },
                    }}
                  >
                    Processing Pipeline
                    <Box
                      component="span"
                      sx={{
                        ml: 1,
                        display: "inline-block",
                        animation: "pulse 2s infinite",
                        "@keyframes pulse": {
                          "0%, 100%": { opacity: 1 },
                          "50%": { opacity: 0.5 },
                        },
                      }}
                    >
                      <Sync sx={{ fontSize: 24, color: "#00F5D4", verticalAlign: "middle" }} />
                    </Box>
                  </Typography>
                }
                subheader={
                  <Typography
                    variant="subtitle1"
                    sx={{
                      color: "#A0A0A0",
                      mt: 1,
                      fontFamily: "'Roboto', sans-serif",
                      display: "flex",
                      alignItems: "center",
                      gap: 1,
                    }}
                  >
                    <Box
                      sx={{
                        width: 8,
                        height: 8,
                        borderRadius: "50%",
                        bgcolor: "#00F5D4",
                        animation: "blink 1.5s infinite",
                        "@keyframes blink": {
                          "0%, 100%": { opacity: 1 },
                          "50%": { opacity: 0.3 },
                        },
                      }}
                    />
                    Real-time analysis progress tracking
                  </Typography>
                }
                sx={{
                  pb: 0,
                  position: "relative",
                  "&:after": {
                    content: '""',
                    position: "absolute",
                    bottom: 0,
                    left: "5%",
                    width: "90%",
                    height: 1,
                    background: "linear-gradient(90deg, transparent 0%, #5E81F4 50%, transparent 100%)",
                  },
                }}
              />

              <CardContent sx={{ pt: 3, pb: 4 }}>
                <LinearProgress
                  variant="determinate"
                  value={(currentStep / timelineSteps.length) * 100}
                  sx={{
                    mb: 6,
                    height: 8,
                    borderRadius: 4,
                    bgcolor: "#0A192F",
                    boxShadow: "inset 0 1px 3px rgba(0, 0, 0, 0.2)",
                    "& .MuiLinearProgress-bar": {
                      bgcolor: "#00F5D4",
                      borderRadius: 4,
                      position: "relative",
                      "&:after": {
                        content: '""',
                        position: "absolute",
                        top: 0,
                        right: 0,
                        bottom: 0,
                        width: 20,
                        background: "linear-gradient(90deg, rgba(0,245,212,0.8) 0%, rgba(0,245,212,0) 100%)",
                        animation: "shine 2s infinite",
                        "@keyframes shine": {
                          "0%": { opacity: 0 },
                          "50%": { opacity: 1 },
                          "100%": { opacity: 0 },
                        },
                      },
                    },
                  }}
                />

                <Grid container spacing={6}>
                  <Grid item xs={12} md={7}>
                    <Box
                      sx={{
                        position: "relative",
                        "&:before": {
                          content: '""',
                          position: "absolute",
                          left: 18,
                          top: 18,
                          bottom: 18,
                          width: 2,
                          bgcolor: "#5E81F4",
                          opacity: 0.2,
                          borderRadius: 2,
                        },
                      }}
                    >
                      {timelineSteps.map((step, index) => (
                        <Box
                          key={index}
                          sx={{
                            display: "flex",
                            alignItems: "center",
                            mb: 4,
                            position: "relative",
                            transition: "all 0.5s ease",
                            transform: index <= currentStep ? "translateX(0)" : "translateX(-10px)",
                            opacity: index <= currentStep ? 1 : 0.7,
                          }}
                        >
                          <Box
                            sx={{
                              width: 40,
                              height: 40,
                              borderRadius: "50%",
                              bgcolor: index <= currentStep ? "#00F5D4" : "#0A192F",
                              display: "flex",
                              alignItems: "center",
                              justifyContent: "center",
                              transition: "all 0.5s cubic-bezier(0.68, -0.55, 0.27, 1.55)",
                              zIndex: 2,
                              position: "relative",
                              border: `2px solid ${index <= currentStep ? "#5E81F4" : "#2E3A59"}`,
                              boxShadow: index <= currentStep ? "0 0 15px rgba(0, 245, 212, 0.5)" : "none",
                              transform: index === currentStep ? "scale(1.1)" : "scale(1)",
                              animation: index === currentStep ? "pulseGlow 2s infinite" : "none",
                              "@keyframes pulseGlow": {
                                "0%, 100%": { boxShadow: "0 0 15px rgba(0, 245, 212, 0.5)" },
                                "50%": { boxShadow: "0 0 25px rgba(0, 245, 212, 0.8)" },
                              },
                            }}
                          >
                            {index < currentStep ? (
                              <CheckCircle
                                sx={{
                                  color: "#0A192F",
                                  fontSize: 22,
                                  animation: "bounceIn 0.5s",
                                  "@keyframes bounceIn": {
                                    "0%": { transform: "scale(0.5)", opacity: 0 },
                                    "70%": { transform: "scale(1.2)", opacity: 1 },
                                    "100%": { transform: "scale(1)", opacity: 1 },
                                  },
                                }}
                              />
                            ) : index === currentStep ? (
                              <CircularProgress
                                size={22}
                                sx={{
                                  color: "#0A192F",
                                  "& .MuiCircularProgress-circle": {
                                    animationDuration: "1.5s",
                                  },
                                }}
                              />
                            ) : (
                              <Typography
                                sx={{
                                  color: "#A0A0A0",
                                  fontSize: 16,
                                  fontWeight: 700,
                                  fontFamily: "'Inter', sans-serif",
                                  transition: "all 0.3s ease",
                                }}
                              >
                                {index + 1}
                              </Typography>
                            )}
                          </Box>
                          <Box sx={{ ml: 3, flex: 1 }}>
                            <Typography
                              variant="body1"
                              sx={{
                                fontWeight: index === currentStep ? 700 : 500,
                                color: index === currentStep ? "#FFFFFF" : index < currentStep ? "#E0E0E0" : "#A0A0A0",
                                transition: "all 0.3s ease",
                                fontFamily: "'Inter', sans-serif",
                                fontSize: 16,
                                letterSpacing: "0.3px",
                                position: "relative",
                                display: "inline-block",
                                "&:after":
                                  index === currentStep
                                    ? {
                                      content: '""',
                                      position: "absolute",
                                      bottom: -4,
                                      left: 0,
                                      width: "100%",
                                      height: 2,
                                      background: "linear-gradient(90deg, #00F5D4 0%, transparent 100%)",
                                      borderRadius: 2,
                                      animation: "underlineExpand 0.5s forwards",
                                      "@keyframes underlineExpand": {
                                        "0%": { width: 0 },
                                        "100%": { width: "100%" },
                                      },
                                    }
                                    : {},
                              }}
                            >
                              {step}
                            </Typography>
                            {index === currentStep &&
                              (retryAttempts.disease > 0 ||
                                retryAttempts.protein > 0 ||
                                retryAttempts.reaction > 0 ||
                                retryAttempts.noReactions > 0) && (
                                <Typography
                                  variant="caption"
                                  sx={{
                                    color: "#FF4D6D",
                                    mt: 1,
                                    fontFamily: "'Roboto', sans-serif",
                                    display: "inline-flex",
                                    alignItems: "center",
                                    gap: 0.5,
                                    animation: "fadeIn 0.5s",
                                    "@keyframes fadeIn": {
                                      "0%": { opacity: 0, transform: "translateY(5px)" },
                                      "100%": { opacity: 1, transform: "translateY(0)" },
                                    },
                                  }}
                                >
                                  <Warning sx={{ fontSize: 14 }} />
                                  {retryAttempts.disease > 0 &&
                                    `Retrying Disease Prediction (Attempt ${retryAttempts.disease}/${maxRetries})...`}
                                  {retryAttempts.protein > 0 &&
                                    `Retrying Protein Prediction (Attempt ${retryAttempts.protein}/${maxRetries})...`}
                                  {retryAttempts.reaction > 0 &&
                                    `Retrying Reaction (Attempt ${retryAttempts.reaction}/${maxRetries})...`}
                                  {retryAttempts.noReactions > 0 &&
                                    `Retrying Drug Prediction (Attempt ${retryAttempts.noReactions}/${maxRetries})...`}
                                </Typography>
                              )}
                          </Box>
                          {index < timelineSteps.length - 1 && (
                            <Box
                              sx={{
                                position: "absolute",
                                left: 20,
                                top: 40,
                                width: "2px",
                                height: 52,
                                bgcolor: index < currentStep ? "#00F5D4" : "#0A192F",
                                transition: "all 0.5s ease",
                                zIndex: 1,
                                opacity: index < currentStep ? 0.8 : 0.3,
                              }}
                            />
                          )}
                        </Box>
                      ))}
                    </Box>
                  </Grid>

                  <Grid item xs={12} md={5}>
                    <Card
                      variant="outlined"
                      sx={{
                        height: "100%",
                        display: "flex",
                        alignItems: "center",
                        bgcolor: "#172A45",
                        borderColor: "#5E81F4",
                        borderRadius: 3,
                        boxShadow: "0 4px 20px rgba(0, 0, 0, 0.15)",
                        position: "relative",
                        overflow: "hidden",
                        "&:before": {
                          content: '""',
                          position: "absolute",
                          top: -50,
                          right: -50,
                          width: 100,
                          height: 100,
                          borderRadius: "50%",
                          bgcolor: "rgba(94, 129, 244, 0.1)",
                          filter: "blur(10px)",
                        },
                        "&:after": {
                          content: '""',
                          position: "absolute",
                          bottom: -30,
                          left: -30,
                          width: 80,
                          height: 80,
                          borderRadius: "50%",
                          bgcolor: "rgba(0, 245, 212, 0.1)",
                          filter: "blur(10px)",
                        },
                      }}
                    >
                      <CardContent sx={{ textAlign: "center", width: "100%", p: 4, position: "relative", zIndex: 1 }}>
                        <Box
                          sx={{
                            position: "relative",
                            height: "180px",
                            width: "100%",
                            mx: "auto",
                            mb: 3,
                            display: "flex",
                            justifyContent: "center",
                            alignItems: "flex-end",
                            gap: 4,
                          }}
                        >
                          {/* Main Flask */}
                          <Box
                            sx={{
                              position: "relative",
                              height: "120px",
                              width: "80px",
                              animation: "float 4s ease-in-out infinite",
                              "@keyframes float": {
                                "0%, 100%": { transform: "translateY(0)" },
                                "50%": { transform: "translateY(-10px)" },
                              },
                            }}
                          >
                            <Science
                              sx={{
                                fontSize: 70,
                                color: "#00F5D4",
                                filter: "drop-shadow(0 0 8px rgba(0, 245, 212, 0.5))",
                                position: "relative",
                                zIndex: 2,
                              }}
                            />
                            {/* Bubbles inside main flask */}
                            {[0, 1, 2, 3].map((i) => (
                              <Box
                                key={`bubble-${i}`}
                                sx={{
                                  position: "absolute",
                                  width: 6 + i * 2,
                                  height: 6 + i * 2,
                                  borderRadius: "50%",
                                  bgcolor: "rgba(0, 245, 212, 0.5)",
                                  bottom: 20 + i * 15,
                                  left: `calc(50% - ${3 + i}px)`,
                                  animation: `bubbleRise ${3 + i}s infinite ${i * 0.5}s`,
                                  "@keyframes bubbleRise": {
                                    "0%": { transform: "translateY(0) scale(0.5)", opacity: 0 },
                                    "10%": { opacity: 0.8 },
                                    "90%": { opacity: 0.6 },
                                    "100%": { transform: "translateY(-80px) scale(1.2)", opacity: 0 },
                                  },
                                }}
                              />
                            ))}
                            {/* Liquid in flask */}
                            <Box
                              sx={{
                                position: "absolute",
                                bottom: 0,
                                left: "50%",
                                transform: "translateX(-50%)",
                                width: "60%",
                                height: "40%",
                                bgcolor: "rgba(0, 245, 212, 0.2)",
                                borderTopLeftRadius: "40px 20px",
                                borderTopRightRadius: "40px 20px",
                                borderBottomLeftRadius: "10px",
                                borderBottomRightRadius: "10px",
                                animation: "liquidMove 8s infinite alternate",
                                "@keyframes liquidMove": {
                                  "0%, 100%": { height: "40%", transform: "translateX(-50%) skew(0deg)" },
                                  "50%": { height: "45%", transform: "translateX(-50%) skew(2deg)" },
                                },
                              }}
                            />
                          </Box>

                          {/* Small flasks */}
                          {[1, 2, 3].map((i) => (
                            <Box
                              key={`flask-${i}`}
                              sx={{
                                position: "relative",
                                height: "80px",
                                width: "50px",
                                animation: `floatSmall ${4 + i}s ease-in-out infinite ${i * 0.3}s`,
                                "@keyframes floatSmall": {
                                  "0%, 100%": { transform: "translateY(0)" },
                                  "50%": { transform: "translateY(-8px)" },
                                },
                              }}
                            >
                              <Science
                                sx={{
                                  fontSize: 40,
                                  color: i === 1 ? "#5E81F4" : i === 2 ? "#FF4D6D" : "#FFD166",
                                  filter: `drop-shadow(0 0 5px ${i === 1 ? "rgba(94, 129, 244, 0.5)" : i === 2 ? "rgba(255, 77, 109, 0.5)" : "rgba(255, 209, 102, 0.5)"})`,
                                  position: "relative",
                                  zIndex: 2,
                                }}
                              />
                              {/* Liquid in small flasks */}
                              <Box
                                sx={{
                                  position: "absolute",
                                  bottom: 0,
                                  left: "50%",
                                  transform: "translateX(-50%)",
                                  width: "60%",
                                  height: "30%",
                                  bgcolor:
                                    i === 1
                                      ? "rgba(94, 129, 244, 0.2)"
                                      : i === 2
                                        ? "rgba(255, 77, 109, 0.2)"
                                        : "rgba(255, 209, 102, 0.2)",
                                  borderTopLeftRadius: "30px 15px",
                                  borderTopRightRadius: "30px 15px",
                                  borderBottomLeftRadius: "8px",
                                  borderBottomRightRadius: "8px",
                                  animation: `liquidMoveSmall ${6 + i}s infinite alternate ${i * 0.5}s`,
                                  "@keyframes liquidMoveSmall": {
                                    "0%, 100%": { height: "30%", transform: "translateX(-50%) skew(0deg)" },
                                    "50%": { height: "35%", transform: "translateX(-50%) skew(1deg)" },
                                  },
                                }}
                              />
                            </Box>
                          ))}

                          {/* Connection lines between flasks */}
                          <Box
                            sx={{
                              position: "absolute",
                              top: "40%",
                              left: 0,
                              right: 0,
                              height: 2,
                              bgcolor: "rgba(94, 129, 244, 0.3)",
                              "&:before, &:after": {
                                content: '""',
                                position: "absolute",
                                top: -4,
                                width: 10,
                                height: 10,
                                borderRadius: "50%",
                                bgcolor: "rgba(0, 245, 212, 0.5)",
                                animation: "pulse 2s infinite",
                                "@keyframes pulse": {
                                  "0%, 100%": { transform: "scale(1)", opacity: 0.8 },
                                  "50%": { transform: "scale(1.5)", opacity: 0.3 },
                                },
                              },
                              "&:before": {
                                left: "25%",
                              },
                              "&:after": {
                                left: "75%",
                              },
                            }}
                          />

                          {/* Floating molecules */}
                          {[1, 2, 3, 4, 5].map((i) => (
                            <Box
                              key={`molecule-${i}`}
                              sx={{
                                position: "absolute",
                                width: 8,
                                height: 8,
                                borderRadius: "50%",
                                bgcolor: i % 2 === 0 ? "rgba(0, 245, 212, 0.6)" : "rgba(94, 129, 244, 0.6)",
                                top: `${Math.random() * 40 + 20}%`,
                                left: `${Math.random() * 80 + 10}%`,
                                animation: `moleculeFloat ${8 + i}s infinite ${i * 0.5}s`,
                                filter: "blur(1px)",
                                "@keyframes moleculeFloat": {
                                  "0%, 100%": { transform: "translate(0, 0)" },
                                  "25%": { transform: "translate(10px, -5px)" },
                                  "50%": { transform: "translate(5px, 5px)" },
                                  "75%": { transform: "translate(-5px, -3px)" },
                                },
                              }}
                            />
                          ))}

                          {/* Wave effects */}
                          {[0, 1, 2].map((i) => (
                            <Box
                              key={`wave-${i}`}
                              sx={{
                                position: "absolute",
                                top: 0,
                                left: 0,
                                right: 0,
                                bottom: 0,
                                borderRadius: "50%",
                                border: `2px solid ${i === 0 ? "#00F5D4" : i === 1 ? "#5E81F4" : "#FF4D6D"}`,
                                opacity: 0.2,
                                animation: `wave ${3 + i}s infinite ${i * 0.5}s`,
                                "@keyframes wave": {
                                  "0%": { transform: "scale(0.5)", opacity: 0.3 },
                                  "70%": { opacity: 0.1 },
                                  "100%": { transform: "scale(1.8)", opacity: 0 },
                                },
                              }}
                            />
                          ))}
                        </Box>

                        <Typography
                          variant="h6"
                          sx={{
                            fontWeight: 700,
                            color: "#FFFFFF",
                            mb: 2,
                            fontFamily: "'Inter', sans-serif",
                            textTransform: "uppercase",
                            letterSpacing: "1px",
                            position: "relative",
                            display: "inline-block",
                            "&:after": {
                              content: '""',
                              position: "absolute",
                              bottom: -8,
                              left: "25%",
                              width: "50%",
                              height: 2,
                              bgcolor: "#5E81F4",
                              borderRadius: 2,
                            },
                          }}
                        >
                          Current Process
                        </Typography>
                        <Typography
                          variant="body1"
                          sx={{
                            color: "#A0A0A0",
                            fontFamily: "'Roboto', sans-serif",
                            minHeight: 60,
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                            px: 2,
                          }}
                        >
                          <Box sx={{ animation: "fadeInOut 2s infinite" }}>{loadingMessages[currentStep]}</Box>
                        </Typography>
                        <Box sx={{ mt: 3 }}>
                          <Typography
                            variant="caption"
                            sx={{
                              color: "#5E81F4",
                              fontFamily: "'Roboto', sans-serif",
                              display: "inline-flex",
                              alignItems: "center",
                              gap: 1,
                            }}
                          >
                            <Box
                              sx={{
                                width: 10,
                                height: 10,
                                borderRadius: "50%",
                                bgcolor: "#5E81F4",
                                animation: "pulse 1.5s infinite",
                                "@keyframes pulse": {
                                  "0%, 100%": { transform: "scale(1)", opacity: 1 },
                                  "50%": { transform: "scale(1.3)", opacity: 0.7 },
                                },
                              }}
                            />
                            Processing in real-time
                          </Typography>
                        </Box>
                      </CardContent>
                    </Card>
                  </Grid>
                </Grid>
              </CardContent>
            </Card>
          )}

          {error && (
            <Alert
              severity="error"
              sx={{
                mb: 4,
                bgcolor: "#FF4D6D",
                color: "#0A192F",
                "& .MuiAlert-icon": { fontSize: 20, color: "#0A192F" },
                fontFamily: "'Roboto', 'Open Sans', sans-serif",
              }}
              onClose={() => setError("")}
              icon={<Warning />}
            >
              <Typography variant="body2" sx={{ fontWeight: 500 }}>
                {error}
              </Typography>
            </Alert>
          )}

          <Grid container spacing={4}>
            {/* Available Chemical Reactions Card */}
            <Grid item xs={12} md={6}>
              <Card
                variant="outlined"
                sx={{
                  height: "100%",
                  borderRadius: 2,
                  bgcolor: "#172A45",
                  borderColor: "#5E81F4",
                  display: "flex",
                  flexDirection: "column",
                }}
              >
                <CardHeader
                  title={
                    <Typography
                      variant="h6"
                      sx={{
                        fontWeight: 600,
                        color: "#E0E0E0",
                        fontFamily: "'Inter', 'Barlow', sans-serif",
                        letterSpacing: 0.5,
                      }}
                    >
                      Available Chemical Reactions
                    </Typography>
                  }
                  subheader={
                    <Typography
                      variant="body2"
                      sx={{
                        color: "#A0A0A0",
                        mt: 0.5,
                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                        lineHeight: 1.4,
                      }}
                    >
                      Browse supported reaction types and mechanisms
                    </Typography>
                  }
                  sx={{ pb: 1 }}
                />
                <Divider sx={{ bgcolor: "#5E81F4", opacity: 0.3, mx: 2 }} />
                <CardContent sx={{ flexGrow: 1, pt: 2 }}>
                  {reactionsLoading ? (
                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        height: "100%",
                        minHeight: 120,
                        gap: 2,
                      }}
                    >
                      <CircularProgress size={24} sx={{ color: "#00F5D4" }} />
                      <Typography
                        variant="subtitle2"
                        sx={{
                          color: "#E0E0E0",
                          fontFamily: "'Roboto', 'Open Sans', sans-serif",
                        }}
                      >
                        Loading available reactions...
                      </Typography>
                    </Box>
                  ) : availableReactions.length > 0 ? (
                    <Box
                      sx={{
                        display: "flex",
                        flexDirection: "column",
                        height: "100%",
                      }}
                    >
                      <Button
                        variant="contained"
                        onClick={handleMenuOpen}
                        startIcon={<Visibility />}
                        sx={{
                          bgcolor: "#00F5D4",
                          color: "#0A192F",
                          "&:hover": { bgcolor: "#5E81F4" },
                          textTransform: "none",
                          fontFamily: "'Lato', sans-serif",
                          fontWeight: 600,
                          alignSelf: "center",
                          mt: 9,
                          px: 4,
                          py: 1.5,
                        }}
                      >
                        View Available Reactions ({availableReactions.length})
                      </Button>
                      <Menu
                        anchorEl={anchorEl}
                        open={Boolean(anchorEl)}
                        onClose={handleMenuClose}
                        PaperProps={{
                          style: {
                            maxHeight: 500,
                            width: "700px",
                            backgroundColor: "#172A45",
                            color: "#E0E0E0",
                            border: "1px solid #5E81F4",
                            borderRadius: 8,
                          },
                        }}
                      >
                        {availableReactions.map((reaction, index) => (
                          <MenuItem
                            key={index}
                            onClick={handleMenuClose}
                            sx={{
                              py: 2,
                              bgcolor: "#172A45",
                              "&:hover": { bgcolor: "#0A192F" },
                              borderBottom:
                                index < availableReactions.length - 1 ? "1px solid rgba(94, 129, 244, 0.2)" : "none",
                            }}
                          >
                            <Box sx={{ width: "100%" }}>
                              <Box sx={{ display: "flex", alignItems: "center", mb: 1 }}>
                                <Chip
                                  label={reaction.type}
                                  size="small"
                                  sx={{
                                    mr: 2,
                                    bgcolor: "#00F5D4",
                                    color: "#0A192F",
                                    fontWeight: 600,
                                    fontFamily: "'Lato', sans-serif",
                                  }}
                                />
                                {reaction.priority && (
                                  <Chip
                                    label={reaction.priority.label}
                                    size="small"
                                    variant="outlined"
                                    sx={{
                                      borderColor: "#5E81F4",
                                      color: "#A0A0A0",
                                      fontFamily: "'Lato', sans-serif",
                                      "&:hover": {
                                        bgcolor: "rgba(94, 129, 244, 0.1)",
                                      },
                                    }}
                                  />
                                )}
                              </Box>
                              <Typography
                                variant="body2"
                                sx={{
                                  mb: 1,
                                  color: "#A0A0A0",
                                  fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  lineHeight: 1.6,
                                }}
                              >
                                {reaction.description}
                              </Typography>
                              <Typography
                                variant="body2"
                                sx={{
                                  fontFamily: "'Fira Code', monospace",
                                  bgcolor: "#0A192F",
                                  p: 1.5,
                                  borderRadius: 1,
                                  fontSize: "0.85rem",
                                  color: "#E0E0E0",
                                  borderLeft: "3px solid #00F5D4",
                                }}
                              >
                                <strong>SMARTS:</strong> {reaction.smarts}
                              </Typography>
                            </Box>
                          </MenuItem>
                        ))}
                      </Menu>
                    </Box>
                  ) : (
                    <Alert
                      severity="info"
                      icon={<Info sx={{ color: "#5E81F4" }} />}
                      sx={{
                        bgcolor: "#0A192F",
                        border: "1px solid #5E81F4",
                        color: "#E0E0E0",
                        mt: 2,
                      }}
                    >
                      <Typography
                        variant="body2"
                        sx={{
                          color: "#A0A0A0",
                          fontFamily: "'Roboto', 'Open Sans', sans-serif",
                          lineHeight: 1.6,
                        }}
                      >
                        No reactions currently available. Please check your connection or try again later.
                      </Typography>
                    </Alert>
                  )}
                </CardContent>
              </Card>
            </Grid>

            {/* Saved Symptom Groups Card */}
            <Grid item xs={12} md={6}>
              <Card
                variant="outlined"
                sx={{
                  height: "100%",
                  borderRadius: 2,
                  bgcolor: "#172A45",
                  borderColor: "#5E81F4",
                  display: "flex",
                  flexDirection: "column",
                }}
              >
                <CardHeader
                  title={
                    <Typography
                      variant="h6"
                      sx={{
                        fontWeight: 600,
                        color: "#E0E0E0",
                        fontFamily: "'Inter', 'Barlow', sans-serif",
                        letterSpacing: 0.5,
                      }}
                    >
                      Saved Symptom Groups
                    </Typography>
                  }
                  subheader={
                    <Typography
                      variant="body2"
                      sx={{
                        color: "#A0A0A0",
                        mt: 0.5,
                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                        lineHeight: 1.4,
                      }}
                    >
                      Select a group of previously used symptoms
                    </Typography>
                  }
                  sx={{ pb: 1 }}
                />
                <Divider sx={{ bgcolor: "#5E81F4", opacity: 0.3, mx: 2 }} />
                <CardContent sx={{ flexGrow: 1, pt: 2 }}>
                  {symptomsLoading ? (
                    <Box
                      sx={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        height: "100%",
                        minHeight: 120,
                        gap: 2,
                      }}
                    >
                      <CircularProgress size={24} sx={{ color: "#00F5D4" }} />
                      <Typography
                        variant="body2"
                        sx={{
                          color: "#E0E0E0",
                          fontFamily: "'Roboto', 'Open Sans', sans-serif",
                        }}
                      >
                        Loading saved symptoms...
                      </Typography>
                    </Box>
                  ) : savedSymptoms.length > 0 ? (
                    <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                      <FormControl fullWidth>
                        <InputLabel
                          id="saved-symptoms-label"
                          sx={{
                            color: "#A0A0A0",
                            fontFamily: "'Lato', sans-serif",
                            "&.Mui-focused": {
                              color: "#00F5D4",
                            },
                          }}
                        >
                          Select symptom group
                        </InputLabel>
                        <Select
                          labelId="saved-symptoms-label"
                          id="saved-symptoms-select"
                          value={
                            savedSymptoms.find((group) => JSON.stringify(group) === JSON.stringify(symptoms)) || []
                          }
                          label="Select symptom group"
                          onChange={(e) => handleSymptomGroupSelect(e.target.value)}
                          renderValue={(selected) =>
                            selected.length ? (
                              <Box
                                sx={{
                                  display: "flex",
                                  alignItems: "center",
                                  gap: 1,
                                }}
                              >
                                <MedicalServices sx={{ fontSize: 18, color: "#00F5D4" }} />
                                <Typography sx={{ fontFamily: "'Lato', sans-serif" }}>{selected.join(", ")}</Typography>
                              </Box>
                            ) : (
                              "Select a symptom group"
                            )
                          }
                          sx={{
                            bgcolor: "#0A192F",
                            color: "#E0E0E0",
                            "& .MuiOutlinedInput-notchedOutline": { borderColor: "#5E81F4" },
                            "&:hover .MuiOutlinedInput-notchedOutline": { borderColor: "#00F5D4" },
                            "&.Mui-focused .MuiOutlinedInput-notchedOutline": { borderColor: "#00F5D4" },
                            "& .MuiSvgIcon-root": { color: "#A0A0A0" },
                            height: "48px",
                          }}
                        >
                          <MenuItem value={[]} sx={{ bgcolor: "#172A45", "&:hover": { bgcolor: "#0A192F" } }}>
                            <em>Clear selection</em>
                          </MenuItem>
                          {savedSymptoms.map((symptomGroup, index) => (
                            <MenuItem
                              key={index}
                              value={symptomGroup}
                              sx={{
                                bgcolor: "#172A45",
                                "&:hover": { bgcolor: "#0A192F" },
                                borderBottom:
                                  index < savedSymptoms.length - 1 ? "1px solid rgba(94, 129, 244, 0.2)" : "none",
                              }}
                            >
                              <Box>
                                <Box sx={{ display: "flex", alignItems: "center", gap: 1, mb: 0.5 }}>
                                  <MedicalInformation sx={{ fontSize: 20, color: "#5E81F4" }} />
                                  <Typography
                                    variant="subtitle1"
                                    sx={{
                                      fontWeight: 500,
                                      color: "#E0E0E0",
                                      fontFamily: "'Inter', 'Barlow', sans-serif",
                                    }}
                                  >
                                    Symptom Group {index + 1}
                                  </Typography>
                                </Box>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                    ml: 3.5,
                                    lineHeight: 1.6,
                                  }}
                                >
                                  {symptomGroup.join(", ")}
                                </Typography>
                              </Box>
                            </MenuItem>
                          ))}
                        </Select>
                      </FormControl>
                      <TextField
                        label="Selected Symptoms"
                        value={symptoms.join(", ")}
                        onChange={(e) =>
                          setSymptoms(
                            e.target.value
                              ? e.target.value
                                .split(",")
                                .map((s) => s.trim())
                                .filter((s) => s)
                              : [],
                          )
                        }
                        fullWidth
                        margin="normal"
                        helperText="Selected symptoms will appear here. Edit to modify."
                        disabled={symptomsLoading}
                        multiline
                        rows={2}
                        sx={{
                          "& .MuiOutlinedInput-root": {
                            bgcolor: "#0A192F",
                            color: "#E0E0E0",
                            "& fieldset": { borderColor: "#5E81F4" },
                            "&:hover fieldset": { borderColor: "#00F5D4" },
                            "&.Mui-focused fieldset": { borderColor: "#00F5D4" },
                          },
                          "& .MuiInputLabel-root": {
                            color: "#A0A0A0",
                            fontFamily: "'Lato', sans-serif",
                            "&.Mui-focused": { color: "#00F5D4" },
                          },
                          "& .MuiFormHelperText-root": {
                            color: "#A0A0A0",
                            fontFamily: "'Roboto', 'Open Sans', sans-serif",
                            mx: 0,
                          },
                        }}
                      />
                    </Box>
                  ) : (
                    <Alert
                      severity="info"
                      icon={<Info sx={{ color: "#5E81F4" }} />}
                      sx={{
                        bgcolor: "#0A192F",
                        border: "1px solid #5E81F4",
                        color: "#E0E0E0",
                        mt: 2,
                      }}
                    >
                      <Typography
                        variant="body2"
                        sx={{
                          color: "#A0A0A0",
                          fontFamily: "'Roboto', 'Open Sans', sans-serif",
                          lineHeight: 1.6,
                        }}
                      >
                        No saved symptom groups found. Start by predicting a new drug with symptoms.
                      </Typography>
                    </Alert>
                  )}
                </CardContent>
              </Card>
            </Grid>
          </Grid>
        </Paper>

        {result && !loading && (
          <Box sx={{ display: "flex", flexDirection: "column", gap: 3 }}>
            {/* Download PDF Report Button - Moved outside all dropdowns */}
            {result && (result.disease || result.proteins || reactionResult) && (
              <Box sx={{ display: "flex", justifyContent: "flex-end", mb: 2 }}>
                <Button
                  variant="contained"
                  onClick={generatePDF}
                  startIcon={<GetApp />}
                  sx={{
                    bgcolor: "#00F5D4",
                    color: "#0A192F",
                    "&:hover": { bgcolor: "#5E81F4" },
                    textTransform: "none",
                    fontFamily: "'Lato', sans-serif",
                    fontWeight: 600,
                  }}
                >
                  Download PDF Report
                </Button>
              </Box>
            )}


            {/* Compare with Dataset Button - Add this after PDF button */}
{reactionResult && 
 reactionResult.reactionResults?.length > 0 && 
 !is100PercentMatch && 
 !is0PercentMatch && (
  <Button
    onClick={handleCompareWithDataset}
    disabled={comparisonLoading}
    startIcon={comparisonLoading ? <CircularProgress size={20} /> : <Sync />}
    sx={{
      bgcolor: "#5E81F4",
      color: "#E0E0E0",
      "&:hover": { bgcolor: "#00F5D4", color: "#0A192F" },
      textTransform: "none",
      fontFamily: "'Lato', sans-serif",
      fontWeight: 600,
      mb: 2,
    }}
  >
    {comparisonLoading ? "Analyzing..." : "Compare with Dataset"}
  </Button>
)}

{/* Comparison Results - Add this after the Compare button */}
{renderComparisonResults()}


            {/* Disease Prediction Results Section */}
            {result?.disease && (
              <Paper elevation={2} sx={{ borderRadius: 2, bgcolor: "#172A45" }}>
                {renderSectionHeader(
                  "Disease Prediction Results",
                  "disease-prediction",
                  result.disease.predictedDiseases?.length || 0,
                  "#F44336",
                )}
                <Collapse in={expanded["disease-prediction"]}>
                  <CardContent>
                    <Box sx={{ mb: 4 }}>
                      <Typography
                        variant="subtitle1"
                        sx={{ fontWeight: 600, mb: 2, color: "#E0E0E0", fontFamily: "'Inter', 'Barlow', sans-serif" }}
                      >
                        Predicted Diseases
                      </Typography>
                      {result.disease.predictedDiseases?.length > 0 ? (
                        result.disease.predictedDiseases.map((disease, index) => (
                          <Card
                            key={index}
                            variant="outlined"
                            sx={{ mb: 2, borderRadius: 2, bgcolor: "#0A192F", borderColor: "#5E81F4" }}
                          >
                            <CardContent sx={{ p: 2 }}>
                              <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                                <Typography
                                  variant="subtitle1"
                                  sx={{
                                    fontWeight: 600,
                                    flex: 1,
                                    color: "#E0E0E0",
                                    fontFamily: "'Inter', 'Barlow', sans-serif",
                                  }}
                                >
                                  {disease.diseaseName || "Unknown Disease"}
                                </Typography>
                                <Chip
                                  label={`Match: ${disease.DiseaseMatchness || "N/A"}`}
                                  size="small"
                                  sx={{
                                    bgcolor: "#70E000",
                                    color: "#0A192F",
                                    fontWeight: 500,
                                    fontFamily: "'Lato', sans-serif",
                                  }}
                                />
                              </Box>
                              {disease.diseaseCautions && disease.diseaseCautions.length > 0 && (
                                <Box>
                                  <Typography
                                    variant="body2"
                                    sx={{
                                      fontWeight: 600,
                                      mb: 1,
                                      color: "#A0A0A0",
                                      fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                    }}
                                  >
                                    Cautions:
                                  </Typography>
                                  <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                                    {disease.diseaseCautions.map((caution, idx) => (
                                      <Chip
                                        key={idx}
                                        label={caution}
                                        size="small"
                                        variant="outlined"
                                        sx={{
                                          borderColor: "#FF4D6D",
                                          color: "#FF4D6D",
                                          fontFamily: "'Lato', sans-serif",
                                        }}
                                      />
                                    ))}
                                  </Box>
                                </Box>
                              )}
                            </CardContent>
                          </Card>
                        ))
                      ) : (
                        <Typography
                          variant="body2"
                          sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                        >
                          No diseases predicted. Please try different symptoms or try again later.
                        </Typography>
                      )}
                    </Box>

                    {/* Detailed Disease Analysis Section */}
                    {result.disease.predictedDiseases[0]?.DiseaseMatchness !== "0%" && (
                      <>
                        <Divider sx={{ my: 3, bgcolor: "#5E81F4" }} />
                        {renderSectionHeader(
                          "Detailed Disease Analysis",
                          "detailed-disease-analysis",
                          Object.entries(result.disease?.DiseaseAnalysis || {}).length,
                          "#2196F3",
                        )}
                        <Collapse in={expanded["detailed-disease-analysis"]}>
                          <Box sx={{ mt: 2 }}>
                            {result.disease.predictedDiseases[0]?.DiseaseMatchness === "0%" ? (
                              <Typography
                                variant="body2"
                                sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                              >
                                No such disease is present in the history with all these symptoms.
                              </Typography>
                            ) : Object.entries(result.disease?.DiseaseAnalysis || {}).length > 0 ? (
                              Object.entries(result.disease.DiseaseAnalysis).map(([key, items]) => (
                                <Card
                                  key={key}
                                  variant="outlined"
                                  sx={{ mb: 2, borderRadius: 2, bgcolor: "#0A192F", borderColor: "#5E81F4" }}
                                >
                                  <CardContent sx={{ p: 2 }}>
                                    <Box
                                      sx={{
                                        display: "flex",
                                        alignItems: "center",
                                        cursor: "pointer",
                                      }}
                                      onClick={() => handleToggle(key)}
                                    >
                                      <Typography
                                        variant="body1"
                                        sx={{
                                          fontWeight: 600,
                                          flex: 1,
                                          textTransform: "capitalize",
                                          color: "#E0E0E0",
                                          fontFamily: "'Inter', 'Barlow', sans-serif",
                                        }}
                                      >
                                        {key.replace(/([A-Z])/g, " $1").trim()}
                                      </Typography>
                                      <Chip
                                        label={items.length}
                                        size="small"
                                        sx={{
                                          mr: 1,
                                          bgcolor: "#5E81F4",
                                          color: "#0A192F",
                                          fontFamily: "'Lato', sans-serif",
                                        }}
                                      />
                                      <IconButton size="small" sx={{ color: "#00F5D4" }}>
                                        {expanded[key] ? <ExpandLess /> : <ExpandMore />}
                                      </IconButton>
                                    </Box>
                                    <Collapse in={expanded[key]}>
                                      <Box sx={{ mt: 2 }}>
                                        {items.map((item, idx) => (
                                          <Paper
                                            key={idx}
                                            variant="outlined"
                                            sx={{ p: 2, mb: 2, bgcolor: "#172A45", borderColor: "#5E81F4" }}
                                          >
                                            <Typography
                                              variant="body2"
                                              sx={{
                                                fontWeight: 500,
                                                mb: 1,
                                                color: "#E0E0E0",
                                                fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                              }}
                                            >
                                              {item.summary || "No summary available"}
                                            </Typography>
                                            <Box sx={{ display: "flex", alignItems: "center", gap: 2, mt: 2 }}>
                                              <Chip
                                                label={`Source: ${item.source || "N/A"}`}
                                                size="small"
                                                variant="outlined"
                                                sx={{
                                                  borderColor: "#5E81F4",
                                                  color: "#A0A0A0",
                                                  fontFamily: "'Lato', sans-serif",
                                                }}
                                              />
                                              {item.url && (
                                                <Button
                                                  size="small"
                                                  href={item.url}
                                                  target="_blank"
                                                  rel="noopener noreferrer"
                                                  sx={{
                                                    textTransform: "none",
                                                    color: "#00F5D4",
                                                    fontFamily: "'Lato', sans-serif",
                                                  }}
                                                >
                                                  View Source
                                                </Button>
                                              )}
                                            </Box>
                                          </Paper>
                                        ))}
                                      </Box>
                                    </Collapse>
                                  </CardContent>
                                </Card>
                              ))
                            ) : (
                              <Typography
                                variant="body2"
                                sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                              >
                                No detailed analysis available.
                              </Typography>
                            )}
                          </Box>
                        </Collapse>
                      </>
                    )}
                  </CardContent>
                </Collapse>
              </Paper>
            )}

            {/* Target Proteins Section - Hidden when 100% match */}
            {result?.proteins?.TargetProteins && !is100PercentMatch && !is0PercentMatch && (
              <Paper elevation={2} sx={{ borderRadius: 2, bgcolor: "#172A45" }}>
                {renderSectionHeader(
                  "Target Proteins",
                  "target-proteins",
                  result.proteins.TargetProteins.length,
                  "#9C27B0",
                )}
                <Collapse in={expanded["target-proteins"]}>
                  <CardContent>
                    <Typography
                      variant="body2"
                      sx={{ color: "#A0A0A0", mb: 3, fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                    >
                      Identified protein targets for therapeutic intervention
                    </Typography>
                    {result.proteins.TargetProteins.length > 0 ? (
                      result.proteins.TargetProteins.map((protein, index) => (
                        <Card
                          key={index}
                          variant="outlined"
                          sx={{ mb: 2, borderRadius: 2, bgcolor: "#0A192F", borderColor: "#5E81F4" }}
                        >
                          <CardContent sx={{ p: 2 }}>
                            <Typography
                              variant="subtitle1"
                              sx={{
                                fontWeight: 600,
                                mb: 2,
                                color: "#E0E0E0",
                                fontFamily: "'Inter', 'Barlow', sans-serif",
                              }}
                            >
                              {protein.proteinName || "Unknown Protein"}
                            </Typography>
                            <Grid container spacing={2}>
                              <Grid item xs={12} sm={6}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Function:
                                </Typography>
                                <Typography
                                  variant="body2"
                                  sx={{ mb: 1, color: "#E0E0E0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                >
                                  {protein.proteinFunction || "N/A"}
                                </Typography>
                              </Grid>
                              <Grid item xs={12} sm={6}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  UniProt ID:
                                </Typography>
                                <Chip
                                  label={protein.proteinUniport || "N/A"}
                                  size="small"
                                  sx={{ fontFamily: "'Fira Code', monospace", bgcolor: "#5E81F4", color: "#0A192F" }}
                                />
                              </Grid>
                              <Grid item xs={12}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Description:
                                </Typography>
                                <Typography
                                  variant="body2"
                                  sx={{ mb: 1, color: "#E0E0E0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                >
                                  {protein.ProtienDiscription || "N/A"}
                                </Typography>
                              </Grid>
                              <Grid item xs={12}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Detailed Information:
                                </Typography>
                                <Typography variant="body2" sx={{ color: "#E0E0E0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}>
                                  {protein.proteinDetailedDiscription || "N/A"}
                                </Typography>
                              </Grid>
                            </Grid>
                          </CardContent>
                        </Card>
                      ))
                    ) : (
                      <Typography
                        variant="body2"
                        sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                      >
                        No target proteins identified.
                      </Typography>
                    )}
                  </CardContent>
                </Collapse>
              </Paper>
            )}

            {/* Target Ligands Section - Hidden when 100% match */}
            {result?.proteins?.TargetLigands && !is100PercentMatch && !is0PercentMatch && (
              <Paper elevation={2} sx={{ borderRadius: 2, bgcolor: "#172A45" }}>
                {renderSectionHeader(
                  "Target Ligands",
                  "target-ligands",
                  result.proteins.TargetLigands.length,
                  "#FF9800",
                )}
                <Collapse in={expanded["target-ligands"]}>
                  <CardContent>
                    <Typography
                      variant="body2"
                      sx={{ color: "#A0A0A0", mb: 3, fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                    >
                      Therapeutic compounds and molecular ligands
                    </Typography>
                    {result.proteins.TargetLigands.length > 0 ? (
                      result.proteins.TargetLigands.map((ligand, index) => (
                        <Card
                          key={index}
                          variant="outlined"
                          sx={{ mb: 2, borderRadius: 2, bgcolor: "#0A192F", borderColor: "#5E81F4" }}
                        >
                          <CardContent sx={{ p: 2 }}>
                            <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                              <Typography
                                variant="subtitle1"
                                sx={{
                                  fontWeight: 600,
                                  flex: 1,
                                  color: "#E0E0E0",
                                  fontFamily: "'Inter', 'Barlow', sans-serif",
                                }}
                              >
                                {ligand.ligandName || "Unknown Ligand"}
                              </Typography>
                              {ligand.LigandSmile === "Not applicable" && (
                                <Chip
                                  label="Protein"
                                  size="small"
                                  sx={{ bgcolor: "#FF4D6D", color: "#0A192F", fontFamily: "'Lato', sans-serif" }}
                                />
                              )}
                            </Box>

                            <Grid container spacing={2}>
                              <Grid item xs={12} sm={6}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Function:
                                </Typography>
                                <Typography
                                  variant="body2"
                                  sx={{ mb: 1, color: "#E0E0E0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                >
                                  {ligand.ligandFunction || "N/A"}
                                </Typography>
                              </Grid>
                              <Grid item xs={12} sm={6}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  DrugBank ID:
                                </Typography>
                                <Chip
                                  label={ligand.ligandDrugBankID || "N/A"}
                                  size="small"
                                  sx={{ fontFamily: "'Fira Code', monospace", bgcolor: "#5E81F4", color: "#0A192F" }}
                                />
                              </Grid>
                              <Grid item xs={12}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  SMILES:
                                </Typography>
                                <Paper
                                  variant="outlined"
                                  sx={{
                                    p: 1,
                                    bgcolor: "#0A192F",
                                    fontFamily: "'Fira Code', monospace",
                                    fontSize: "0.85rem",
                                    wordBreak: "break-all",
                                    color: "#E0E0E0",
                                    borderColor: "#5E81F4",
                                  }}
                                >
                                  {ligand.LigandSmile || "N/A"}
                                </Paper>
                              </Grid>
                              {ligand.LigandSmile &&
                                ligand.LigandSmile !== "Not applicable" &&
                                smilesImages[ligand.LigandSmile] ? (
                                <Grid item xs={12}>
                                  <Typography
                                    variant="body2"
                                    sx={{
                                      fontWeight: 600,
                                      color: "#A0A0A0",
                                      fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                    }}
                                  >
                                    Structure:
                                  </Typography>
                                  <img
                                    src={`data:image/png;base64,${smilesImages[ligand.LigandSmile]}`}
                                    alt={`Structure of ${ligand.ligandName}`}
                                    style={{ maxWidth: "200px", marginTop: "8px" }}
                                  />
                                </Grid>
                              ) : (
                                <Grid item xs={12}>
                                  <Typography
                                    variant="body2"
                                    sx={{
                                      fontWeight: 600,
                                      color: "#A0A0A0",
                                      fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                    }}
                                  >
                                    Structure:
                                  </Typography>
                                  <Typography
                                    variant="caption"
                                    sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                  >
                                    {ligand.LigandSmile === "Not applicable"
                                      ? "Structure not applicable (Protein)"
                                      : "Unable to load structure image"}
                                  </Typography>
                                </Grid>
                              )}
                              <Grid item xs={12}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Functional Groups:
                                </Typography>
                                <Typography
                                  variant="body2"
                                  sx={{ mb: 1, color: "#E0E0E0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                >
                                  {getFunctionalGroups(ligand.LigandSmile)}
                                </Typography>
                              </Grid>
                              <Grid item xs={12}>
                                <Typography
                                  variant="body2"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Description:
                                </Typography>
                                <Typography
                                  variant="body2"
                                  sx={{ color: "#E0E0E0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                >
                                  {ligand.LigandDiscription || "N/A"}
                                </Typography>
                              </Grid>
                              {ligand.LigandSmile !== "Not applicable" &&
                                reactionResult &&
                                reactionResult.reactants && (
                                  <Grid item xs={12}>
                                    <Typography
                                      variant="body2"
                                      sx={{
                                        fontWeight: 600,
                                        color: "#A0A0A0",
                                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                      }}
                                    >
                                      ADMET Properties:
                                    </Typography>
                                    {(() => {
                                      const reactant = reactionResult.reactants.find(
                                        (r) => r.smiles === ligand.LigandSmile,
                                      )
                                      if (
                                        !reactant ||
                                        !reactant.admet_properties ||
                                        reactant.admet_properties.length === 0
                                      ) {
                                        return (
                                          <Typography
                                            variant="body2"
                                            sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                          >
                                            Not available
                                          </Typography>
                                        )
                                      }
                                      return reactant.admet_properties.map((section, secIdx) => (
                                        <Box key={secIdx} sx={{ mt: 2 }}>
                                          <Typography
                                            variant="subtitle2"
                                            sx={{
                                              fontWeight: 600,
                                              mb: 1,
                                              color: "#E0E0E0",
                                              fontFamily: "'Inter', 'Barlow', sans-serif",
                                            }}
                                          >
                                            {section.section}
                                          </Typography>
                                          <TableContainer
                                            component={Paper}
                                            sx={{
                                              boxShadow: "0 4px 12px rgba(0,0,0,0.05)",
                                              borderRadius: 2,
                                              overflow: "hidden",
                                              bgcolor: "#0A192F",
                                            }}
                                          >
                                            <Table size="small">
                                              <TableHead>
                                                <TableRow sx={{ bgcolor: "#5E81F4" }}>
                                                  <TableCell
                                                    sx={{
                                                      fontWeight: 600,
                                                      color: "#0A192F",
                                                      py: 1.5,
                                                      px: 2,
                                                      borderBottom: "none",
                                                      fontFamily: "'Lato', sans-serif",
                                                    }}
                                                  >
                                                    Property
                                                  </TableCell>
                                                  <TableCell
                                                    sx={{
                                                      fontWeight: 600,
                                                      color: "#0A192F",
                                                      py: 1.5,
                                                      px: 2,
                                                      borderBottom: "none",
                                                      fontFamily: "'Lato', sans-serif",
                                                    }}
                                                  >
                                                    Prediction
                                                  </TableCell>
                                                  <TableCell
                                                    sx={{
                                                      fontWeight: 600,
                                                      color: "#0A192F",
                                                      py: 1.5,
                                                      px: 2,
                                                      borderBottom: "none",
                                                      fontFamily: "'Lato', sans-serif",
                                                    }}
                                                  >
                                                    Units
                                                  </TableCell>
                                                </TableRow>
                                              </TableHead>
                                              <TableBody>
                                                {section.properties.map((prop, propIdx) => (
                                                  <TableRow
                                                    key={propIdx}
                                                    sx={{
                                                      "&:nth-of-type(odd)": { bgcolor: "#172A45" },
                                                      "&:hover": { bgcolor: "#0A192F" },
                                                      transition: "background-color 0.3s ease",
                                                    }}
                                                  >
                                                    <TableCell
                                                      sx={{
                                                        py: 1.5,
                                                        px: 2,
                                                        borderBottom: "1px solid #5E81F4",
                                                        color: "#E0E0E0",
                                                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                      }}
                                                    >
                                                      {prop.name}
                                                    </TableCell>
                                                    <TableCell
                                                      sx={{
                                                        py: 1.5,
                                                        px: 2,
                                                        borderBottom: "1px solid #5E81F4",
                                                        color: "#E0E0E0",
                                                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                      }}
                                                    >
                                                      {prop.prediction}
                                                    </TableCell>
                                                    <TableCell
                                                      sx={{
                                                        py: 1.5,
                                                        px: 2,
                                                        borderBottom: "1px solid #5E81F4",
                                                        color: "#E0E0E0",
                                                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                      }}
                                                    >
                                                      {prop.units || "-"}
                                                    </TableCell>
                                                  </TableRow>
                                                ))}
                                              </TableBody>
                                            </Table>
                                          </TableContainer>
                                        </Box>
                                      ))
                                    })()}
                                  </Grid>
                                )}
                            </Grid>

                            {ligand.LigandSmile === "Not applicable" && (
                              <Alert
                                severity="warning"
                                sx={{
                                  mt: 2,
                                  bgcolor: "#FF4D6D",
                                  color: "#0A192F",
                                  fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                }}
                                icon={<Warning />}
                              >
                                <Typography variant="body2">
                                  This ligand is a protein and cannot be processed for chemical reactions.
                                </Typography>
                              </Alert>
                            )}
                          </CardContent>
                        </Card>
                      ))
                    ) : (
                      <Typography
                        variant="body2"
                        sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                      >
                        No target ligands identified.
                      </Typography>
                    )}
                  </CardContent>
                </Collapse>
              </Paper>
            )}

            {/* Chemical Reaction Results Section - Hidden when 100% match */}
            {reactionResult && !is100PercentMatch && !is0PercentMatch && (
              <Paper elevation={2} sx={{ borderRadius: 2, bgcolor: "#172A45" }}>
                {renderSectionHeader(
                  "Chemical Reaction Results",
                  "chemical-reactions",
                  `${reactionResult.reactionResults?.length || 0} Successful`,
                  "#70E000",
                )}
                <Collapse in={expanded["chemical-reactions"]}>
                  <CardContent>
                    <Typography
                      variant="body2"
                      sx={{ color: "#A0A0A0", mb: 3, fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                    >
                      Computational chemistry analysis and product generation
                    </Typography>
                    <Grid container spacing={3}>
                      <Grid item xs={12}>
                        <Typography
                          variant="subtitle1"
                          sx={{ fontWeight: 600, mb: 2, color: "#E0E0E0", fontFamily: "'Inter', 'Barlow', sans-serif" }}
                        >
                          Successful Reactions
                        </Typography>
                        {reactionResult.reactionResults.length === 0 ? (
                          <Alert severity="info" icon={<Info sx={{ color: "#5E81F4" }} />}>
                            <Typography
                              variant="body2"
                              sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                            >
                              No successful reactions found.{" "}
                              {retryAttempts.noReactions < maxRetries
                                ? "Retrying automatically..."
                                : "Please try different symptoms or check the backend."}
                            </Typography>
                          </Alert>
                        ) : (
                          reactionResult.reactionResults.map((reaction, index) => {
                            const mainProducts = deduplicateProducts(reaction.products.slice(0, 1))
                            return (
                              <Card
                                key={index}
                                variant="outlined"
                                sx={{ mb: 3, borderRadius: 2, bgcolor: "#0A192F", borderColor: "#5E81F4" }}
                              >
                                <CardContent sx={{ p: 2 }}>
                                  <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                                    <Typography
                                      variant="subtitle1"
                                      sx={{
                                        fontWeight: 600,
                                        flex: 1,
                                        color: "#E0E0E0",
                                        fontFamily: "'Inter', 'Barlow', sans-serif",
                                      }}
                                    >
                                      {reaction.reactionType || "Unknown Reaction"}
                                    </Typography>
                                    <Tooltip title={`Confidence: ${(reaction.confidence * 100).toFixed(1)}%`}>
                                      <Chip
                                        label={`${(reaction.confidence * 100).toFixed(1)}% Confidence`}
                                        size="small"
                                        sx={{
                                          bgcolor:
                                            reaction.confidence >= 0.7
                                              ? "#70E000" // Green for 70-100%
                                              : reaction.confidence >= 0.5
                                                ? "#FFD700" // Yellow for 50-70%
                                                : "#FF4D6D", // Red for 0-40%
                                          color: "#0A192F",
                                          fontWeight: 500,
                                          fontFamily: "'Lato', sans-serif",
                                          mr: 1,
                                        }}
                                      />
                                    </Tooltip>
                                    <IconButton
                                      size="small"
                                      onClick={() => handleToggle(`reaction-${index}`)}
                                      sx={{ color: "#00F5D4" }}
                                    >
                                      {expanded[`reaction-${index}`] ? <ExpandLess /> : <ExpandMore />}
                                    </IconButton>
                                  </Box>

                                  <Typography
                                    variant="body2"
                                    sx={{ mb: 2, color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                  >
                                    {reaction.description || "No description available"}
                                  </Typography>

                                  <Box sx={{ mb: 2 }}>
                                    <Typography
                                      variant="body2"
                                      sx={{
                                        fontWeight: 600,
                                        mb: 1,
                                        color: "#A0A0A0",
                                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                      }}
                                    >
                                      Reactants:
                                    </Typography>
                                    <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                                      {reaction.reactants.map((reactant, idx) => (
                                        <Box key={idx} sx={{ mb: 1 }}>
                                          <Chip
                                            label={reactant}
                                            size="small"
                                            variant="outlined"
                                            sx={{
                                              fontFamily: "'Fira Code', monospace",
                                              borderColor: "#5E81F4",
                                              color: "#E0E0E0",
                                            }}
                                          />
                                        </Box>
                                      ))}
                                    </Box>
                                  </Box>

                                  <Collapse in={expanded[`reaction-${index}`]}>
                                    <Divider sx={{ my: 2, bgcolor: "#5E81F4" }} />
                                    <Typography
                                      variant="subtitle2"
                                      sx={{
                                        fontWeight: 600,
                                        mb: 2,
                                        color: "#E0E0E0",
                                        fontFamily: "'Inter', 'Barlow', sans-serif",
                                      }}
                                    >
                                      Main Product
                                    </Typography>
                                    <Grid container spacing={2}>
                                      {mainProducts.map((product, prodIndex) => (
                                        <Grid item xs={12} key={prodIndex}>
                                          <Paper
                                            variant="outlined"
                                            sx={{ p: 2, borderRadius: 2, bgcolor: "#172A45", borderColor: "#5E81F4" }}
                                          >
                                            <Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
                                              <Typography
                                                variant="subtitle2"
                                                sx={{
                                                  fontWeight: 600,
                                                  flex: 1,
                                                  color: "#E0E0E0",
                                                  fontFamily: "'Inter', 'Barlow', sans-serif",
                                                }}
                                              >
                                                Product {prodIndex + 1}
                                              </Typography>
                                            </Box>
                                            <Grid container spacing={2}>
                                              <Grid item xs={12} md={6}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{
                                                    fontWeight: 600,
                                                    color: "#A0A0A0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  SMILES:
                                                </Typography>
                                                <Paper
                                                  variant="outlined"
                                                  sx={{
                                                    p: 1,
                                                    mt: 0.5,
                                                    bgcolor: "#0A192F",
                                                    fontFamily: "'Fira Code', monospace",
                                                    fontSize: "0.8rem",
                                                    wordBreak: "break-all",
                                                    color: "#E0E0E0",
                                                    borderColor: "#5E81F4",
                                                  }}
                                                >
                                                  {product.smiles || "N/A"}
                                                </Paper>
                                              </Grid>
                                              <Grid item xs={12} md={6}>
                                                {smilesImages[product.smiles] ? (
                                                  <Box sx={{ display: "flex", justifyContent: "center" }}>
                                                    <img
                                                      src={`data:image/png;base64,${smilesImages[product.smiles]}`}
                                                      alt={`Product ${prodIndex + 1}`}
                                                      style={{ maxWidth: "100%", height: "auto", maxHeight: "200px" }}
                                                    />
                                                  </Box>
                                                ) : (
                                                  <Typography
                                                    variant="caption"
                                                    sx={{
                                                      color: "#A0A0A0",
                                                      fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                    }}
                                                  >
                                                    Unable to load structure image
                                                  </Typography>
                                                )}
                                              </Grid>
                                              <Grid item xs={6} md={3}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{
                                                    fontWeight: 600,
                                                    color: "#A0A0A0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  Molecular Weight:
                                                </Typography>
                                                <Typography
                                                  variant="body2"
                                                  sx={{
                                                    color: "#E0E0E0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  {product.molecular_weight?.toFixed(2) || "N/A"} g/mol
                                                </Typography>
                                              </Grid>
                                              <Grid item xs={6} md={3}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{
                                                    fontWeight: 600,
                                                    color: "#A0A0A0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  LogP:
                                                </Typography>
                                                <Typography
                                                  variant="body2"
                                                  sx={{
                                                    color: "#E0E0E0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  {product.logP?.toFixed(2) || "N/A"}
                                                </Typography>
                                              </Grid>
                                              <Grid item xs={6} md={3}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{
                                                    fontWeight: 600,
                                                    color: "#A0A0A0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  TPSA:
                                                </Typography>
                                                <Typography
                                                  variant="body2"
                                                  sx={{
                                                    color: "#E0E0E0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  {product.tpsa?.toFixed(2) || "N/A"} Å²
                                                </Typography>
                                              </Grid>
                                              <Grid item xs={6} md={3}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{
                                                    fontWeight: 600,
                                                    color: "#A0A0A0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  H-Donors/Acceptors:
                                                </Typography>
                                                <Typography
                                                  variant="body2"
                                                  sx={{
                                                    color: "#E0E0E0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  {product.num_h_donors || 0}/{product.num_h_acceptors || 0}
                                                </Typography>
                                              </Grid>
                                              <Grid item xs={12}>
                                                <Typography
                                                  variant="caption"
                                                  sx={{
                                                    fontWeight: 600,
                                                    color: "#A0A0A0",
                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                  }}
                                                >
                                                  Functional Groups:
                                                </Typography>
                                                <Box sx={{ display: "flex", flexWrap: "wrap", gap: 0.5, mt: 0.5 }}>
                                                  {product.functional_groups?.length > 0 ? (
                                                    product.functional_groups.map((group, groupIdx) => (
                                                      <Chip
                                                        key={groupIdx}
                                                        label={group}
                                                        size="small"
                                                        variant="outlined"
                                                        sx={{
                                                          fontSize: "0.7rem",
                                                          borderColor: "#5E81F4",
                                                          color: "#E0E0E0",
                                                          fontFamily: "'Lato', sans-serif",
                                                        }}
                                                      />
                                                    ))
                                                  ) : (
                                                    <Typography
                                                      variant="body2"
                                                      sx={{
                                                        color: "#A0A0A0",
                                                        fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                      }}
                                                    >
                                                      None identified
                                                    </Typography>
                                                  )}
                                                </Box>
                                              </Grid>
                                              {product.admet_properties && product.admet_properties.length > 0 && (
                                                <Grid item xs={12}>
                                                  <Typography
                                                    variant="caption"
                                                    sx={{
                                                      fontWeight: 600,
                                                      color: "#A0A0A0",
                                                      mt: 2,
                                                      fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                    }}
                                                  >
                                                    ADMET Properties:
                                                  </Typography>
                                                  {product.admet_properties.map((section, secIdx) => (
                                                    <Box key={secIdx} sx={{ mt: 2 }}>
                                                      <Typography
                                                        variant="subtitle2"
                                                        sx={{
                                                          fontWeight: 600,
                                                          mb: 1,
                                                          color: "#E0E0E0",
                                                          fontFamily: "'Inter', 'Barlow', sans-serif",
                                                        }}
                                                      >
                                                        {section.section}
                                                      </Typography>
                                                      <TableContainer
                                                        component={Paper}
                                                        sx={{
                                                          boxShadow: "0 4px 12px rgba(0,0,0,0.05)",
                                                          borderRadius: 2,
                                                          overflow: "hidden",
                                                          bgcolor: "#0A192F",
                                                        }}
                                                      >
                                                        <Table size="small">
                                                          <TableHead>
                                                            <TableRow sx={{ bgcolor: "#5E81F4" }}>
                                                              <TableCell
                                                                sx={{
                                                                  fontWeight: 600,
                                                                  color: "#0A192F",
                                                                  py: 1.5,
                                                                  px: 2,
                                                                  borderBottom: "none",
                                                                  fontFamily: "'Lato', sans-serif",
                                                                }}
                                                              >
                                                                Property
                                                              </TableCell>
                                                              <TableCell
                                                                sx={{
                                                                  fontWeight: 600,
                                                                  color: "#0A192F",
                                                                  py: 1.5,
                                                                  px: 2,
                                                                  borderBottom: "none",
                                                                  fontFamily: "'Lato', sans-serif",
                                                                }}
                                                              >
                                                                Prediction
                                                              </TableCell>
                                                              <TableCell
                                                                sx={{
                                                                  fontWeight: 600,
                                                                  color: "#0A192F",
                                                                  py: 1.5,
                                                                  px: 2,
                                                                  borderBottom: "none",
                                                                  fontFamily: "'Lato', sans-serif",
                                                                }}
                                                              >
                                                                Units
                                                              </TableCell>
                                                            </TableRow>
                                                          </TableHead>
                                                          <TableBody>
                                                            {section.properties.map((prop, propIdx) => (
                                                              <TableRow
                                                                key={propIdx}
                                                                sx={{
                                                                  "&:nth-of-type(odd)": { bgcolor: "#172A45" },
                                                                  "&:hover": { bgcolor: "#0A192F" },
                                                                  transition: "background-color 0.3s ease",
                                                                }}
                                                              >
                                                                <TableCell
                                                                  sx={{
                                                                    py: 1.5,
                                                                    px: 2,
                                                                    borderBottom: "1px solid #5E81F4",
                                                                    color: "#E0E0E0",
                                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                                  }}
                                                                >
                                                                  {prop.name}
                                                                </TableCell>
                                                                <TableCell
                                                                  sx={{
                                                                    py: 1.5,
                                                                    px: 2,
                                                                    borderBottom: "1px solid #5E81F4",
                                                                    color: "#E0E0E0",
                                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                                  }}
                                                                >
                                                                  {prop.prediction}
                                                                </TableCell>
                                                                <TableCell
                                                                  sx={{
                                                                    py: 1.5,
                                                                    px: 2,
                                                                    borderBottom: "1px solid #5E81F4",
                                                                    color: "#E0E0E0",
                                                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                                                  }}
                                                                >
                                                                  {prop.units || "-"}
                                                                </TableCell>
                                                              </TableRow>
                                                            ))}
                                                          </TableBody>
                                                        </Table>
                                                      </TableContainer>
                                                    </Box>
                                                  ))}
                                                </Grid>
                                              )}
                                            </Grid>
                                          </Paper>
                                        </Grid>
                                      ))}
                                    </Grid>

                                    <Box sx={{ mt: 4 }}>
                                      <Typography
                                        variant="subtitle2"
                                        sx={{
                                          fontWeight: 600,
                                          mb: 2,
                                          color: "#E0E0E0",
                                          fontFamily: "'Inter', 'Barlow', sans-serif",
                                        }}
                                      >
                                        ADMET Properties Comparison
                                      </Typography>
                                      {prepareAdmetGraphData(index) ? (
                                        <Box sx={{ height: 400, width: "100%" }}>
                                          <Radar
                                            key={`radar-${index}`}
                                            data={prepareAdmetGraphData(index)}
                                            options={{
                                              responsive: true,
                                              maintainAspectRatio: false,
                                              plugins: {
                                                legend: {
                                                  position: "top",
                                                  labels: { color: "#E0E0E0", font: { family: "'Lato', sans-serif" } },
                                                },
                                                title: {
                                                  display: true,
                                                  text: "ADMET Properties Analysis",
                                                  color: "#E0E0E0",
                                                  font: { family: "'Inter', 'Barlow', sans-serif", size: 16 },
                                                },
                                                tooltip: {
                                                  callbacks: {
                                                    label: (tooltipItem) => {
                                                      return `${tooltipItem.dataset.label}: ${tooltipItem.raw}`
                                                    },
                                                  },
                                                  backgroundColor: "#172A45",
                                                  titleColor: "#E0E0E0",
                                                  bodyColor: "#E0E0E0",
                                                  bodyFont: { family: "'Roboto', 'Open Sans', sans-serif" },
                                                },
                                              },
                                              scales: {
                                                r: {
                                                  beginAtZero: true,
                                                  min: 0,
                                                  max: 100,
                                                  ticks: {
                                                    stepSize: 25,
                                                    color: "#000000",
                                                    font: { family: "'Lato', sans-serif" },
                                                  },
                                                  pointLabels: {
                                                    font: {
                                                      size: 12,
                                                      family: "'Lato', sans-serif",
                                                    },
                                                    color: "#E0E0E0",
                                                  },
                                                  grid: {
                                                    color: "#5E81F4",
                                                  },
                                                },
                                              },
                                            }}
                                          />
                                        </Box>
                                      ) : (
                                        <Typography
                                          variant="body2"
                                          sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                        >
                                          No ADMET data available for comparison.
                                        </Typography>
                                      )}
                                    </Box>
                                  </Collapse>
                                </CardContent>
                              </Card>
                            )
                          })
                        )}
                      </Grid>

                      <Grid item xs={12}>
                        <Card variant="outlined" sx={{ borderRadius: 2, bgcolor: "#172A45", borderColor: "#5E81F4" }}>
                          <CardHeader
                            title={
                              <Typography
                                variant="subtitle1"
                                sx={{ fontWeight: 600, color: "#E0E0E0", fontFamily: "'Inter', 'Barlow', sans-serif" }}
                              >
                                Overall Statistics
                              </Typography>
                            }
                          />
                          <CardContent sx={{ pt: 0 }}>
                            <Grid container spacing={2}>
                              <Grid item xs={3}>
                                <Typography
                                  variant="caption"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Mean MW:
                                </Typography>
                                <Typography
                                  variant="h6"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#E0E0E0",
                                    fontFamily: "'Inter', 'Barlow', sans-serif",
                                  }}
                                >
                                  {reactionResult.statistics?.mean_mw?.toFixed(2) || "N/A"}
                                </Typography>
                              </Grid>
                              <Grid item xs={3}>
                                <Typography
                                  variant="caption"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Std MW:
                                </Typography>
                                <Typography
                                  variant="h6"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#E0E0E0",
                                    fontFamily: "'Inter', 'Barlow', sans-serif",
                                  }}
                                >
                                  {reactionResult.statistics?.std_mw?.toFixed(2) || "N/A"}
                                </Typography>
                              </Grid>
                              <Grid item xs={3}>
                                <Typography
                                  variant="caption"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Min MW:
                                </Typography>
                                <Typography
                                  variant="h6"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#E0E0E0",
                                    fontFamily: "'Inter', 'Barlow', sans-serif",
                                  }}
                                >
                                  {reactionResult.statistics?.min_mw?.toFixed(2) || "N/A"}
                                </Typography>
                              </Grid>
                              <Grid item xs={3}>
                                <Typography
                                  variant="caption"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#A0A0A0",
                                    fontFamily: "'Roboto', 'Open Sans', sans-serif",
                                  }}
                                >
                                  Max MW:
                                </Typography>
                                <Typography
                                  variant="h6"
                                  sx={{
                                    fontWeight: 600,
                                    color: "#E0E0E0",
                                    fontFamily: "'Inter', 'Barlow', sans-serif",
                                  }}
                                >
                                  {reactionResult.statistics?.max_mw?.toFixed(2) || "N/A"}
                                </Typography>
                              </Grid>
                            </Grid>
                          </CardContent>
                        </Card>
                      </Grid>

                      {/* Failed Reactions Section - Hidden when 100% match */}
                      <Grid item xs={12}>
                        {renderSectionHeader(
                          "Failed Reactions",
                          "failed-reactions",
                          reactionResult.failedReactions?.length || 0,
                          "#FF4D6D",
                        )}
                        <Collapse in={expanded["failed-reactions"]}>
                          <Box sx={{ mt: 2 }}>
                            {reactionResult.failedReactions?.length === 0 ? (
                              <Alert severity="success" icon={<CheckCircle sx={{ color: "#70E000" }} />}>
                                <Typography
                                  variant="body2"
                                  sx={{ color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                >
                                  All reactions completed successfully!
                                </Typography>
                              </Alert>
                            ) : (
                              reactionResult.failedReactions?.map((fail, index) => (
                                <Card
                                  key={index}
                                  variant="outlined"
                                  sx={{ mb: 2, borderRadius: 2, bgcolor: "#0A192F", borderColor: "#5E81F4" }}
                                >
                                  <CardContent sx={{ p: 2 }}>
                                    <Typography
                                      variant="subtitle2"
                                      sx={{
                                        fontWeight: 600,
                                        mb: 1,
                                        color: "#E0E0E0",
                                        fontFamily: "'Inter', 'Barlow', sans-serif",
                                      }}
                                    >
                                      {fail.reactionType || "Unknown Reaction"}
                                    </Typography>
                                    <Typography
                                      variant="body2"
                                      sx={{ mb: 1, color: "#A0A0A0", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                    >
                                      Reactants: {fail.reactants?.join(", ") || "N/A"}
                                    </Typography>
                                    <Typography
                                      variant="body2"
                                      sx={{ color: "#FF4D6D", fontFamily: "'Roboto', 'Open Sans', sans-serif" }}
                                    >
                                      Reason: {fail.reason || "Unknown error"}
                                    </Typography>
                                  </CardContent>
                                </Card>
                              ))
                            )}
                          </Box>
                        </Collapse>
                      </Grid>
                    </Grid>
                  </CardContent>
                </Collapse>
              </Paper>
            )}
          </Box>
        )}
      </Container>
    </Box>
  )
}

export default RDkit;


