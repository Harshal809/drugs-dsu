import React, { useState, useEffect, useCallback, useRef } from 'react';
import axios from 'axios';
import { RefreshCw, ChevronDown, ChevronUp, AlertCircle, Info, AlertTriangle, Beaker, Activity, Shield, Eye } from 'lucide-react';
import { useAuthStore } from '../../Store/auth.store';
import { useNavigate } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode === 'development' ? API_BASE_URL : '/api',
  withCredentials: true,
});

const ToxicityPrediction = () => {
  const isMounted = useRef(true);
  const { user, checkAuth, logout, checkingAuth } = useAuthStore();
  const navigate = useNavigate();

  const [compound, setCompound] = useState('');
  const [symptomGroupIndex, setSymptomGroupIndex] = useState('');
  const [symptomGroups, setSymptomGroups] = useState([]);
  const [productCompoundGroups, setProductCompoundGroups] = useState([]);
  const [result, setResult] = useState(null);
  const [history, setHistory] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [notification, setNotification] = useState(null);
  const [isResultOpen, setIsResultOpen] = useState(false);
  const [isSideEffectsOpen, setIsSideEffectsOpen] = useState(false);
  const [geminiAnalysis, setGeminiAnalysis] = useState(null);
  const [geminiLoading, setGeminiLoading] = useState(false);
  const [isGeminiAnalysisOpen, setIsGeminiAnalysisOpen] = useState(false);
  const [activeTab, setActiveTab] = useState('predict');
  const [selectedHistoryItemId, setSelectedHistoryItemId] = useState(null);
  const [isHistoryResultOpen, setIsHistoryResultOpen] = useState(false);
  const [isHistorySideEffectsOpen, setIsHistorySideEffectsOpen] = useState(false);
  const [isHistoryGeminiOpen, setIsHistoryGeminiOpen] = useState(false);

  useEffect(() => {
    return () => {
      isMounted.current = false;
    };
  }, []);

  const handleAuthError = useCallback((err) => {
    if (!isMounted.current) return;
    if (err.response?.status === 401) {
      const message = err.response.data.message || 'Not authorized';
      setError(message);
      setNotification(err.response.data.suggestion || 'Please log in to continue.');
      if (message.includes('login')) {
        logout();
        navigate('/login');
      }
    } else {
      setError(err.response?.data?.message || 'An error occurred');
      setNotification(err.response?.data?.suggestion || 'Please check your inputs and try again.');
    }
  }, [logout, navigate]);

  const fetchSymptomsAndProducts = async () => {
    if (!user?._id) return;
    try {
      const response = await axiosInstance.get(`/getdata/getsymptoms-product/${user._id}`);
      const { symptoms, productSmiles } = response.data;
      setSymptomGroups(symptoms || []);
      setProductCompoundGroups(productSmiles || []);
      if (symptoms?.length > 0) {
        setSymptomGroupIndex('0');
        if (productSmiles?.[0]?.length > 0) {
          setCompound(productSmiles[0][0]);
        }
      }
    } catch (err) {
      console.error("Error fetching symptoms and products:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch symptoms and products");
      setNotification(err.response?.data?.suggestion || 'Please try again or contact support.');
      setSymptomGroups([]);
      setProductCompoundGroups([]);
    }
  };

  const fetchHistory = async () => {
    if (!user?._id) return;
    try {
      const response = await axiosInstance.get('/toxicity/history');
      const fetchedHistory = response.data.history || [];
      const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "{}");
      const updatedHistory = fetchedHistory.map(item => ({
        ...item,
        symptomsGrp: storedSymptomGroups[item._id] || item.symptoms || "Not available"
      }));
      setHistory(updatedHistory);
    } catch (err) {
      console.error('Error fetching history:', err.response?.data || err.message);
      handleAuthError(err);
    }
  };

  const getGeminiAnalysis = async (compoundString, symptomsString) => {
    console.log("Querying Gemini for toxicity analysis:", compoundString, symptomsString);
    setGeminiLoading(true);
    try {
      const response = await axiosInstance.post('/toxicity/gemini-analysis', {
        compound: compoundString,
        symptoms: symptomsString
      });
      const analysis = response.data.analysis;
      await saveGeminiAnalysis(compoundString, symptomsString, analysis);
      return analysis;
    } catch (err) {
      console.error("Error in getGeminiAnalysis:", err.response?.data || err.message);
      throw new Error(err.response?.data?.message || "Failed to generate Gemini analysis: " + err.message);
    } finally {
      setGeminiLoading(false);
    }
  };

  const saveGeminiAnalysis = async (compoundString, symptomsString, analysis) => {
    try {
      const response = await axiosInstance.post('/toxicity/save-analysis', {
        compoundName: compoundString,
        symptoms: symptomsString,
        geminiAnalysis: analysis
      });
      console.log("Successfully saved Gemini analysis:", response.data);
      return response.data.historyId;
    } catch (err) {
      console.error("Failed to save Gemini analysis:", err.response?.data || err.message);
      throw err;
    }
  };

  useEffect(() => {
    let isActive = true;
    const initialize = async () => {
      console.log('Starting initialization - checkingAuth:', checkingAuth, 'user:', user);
      setLoading(true);
      try {
        if (user && user._id) {
          console.log('User already authenticated, skipping checkAuth');
        } else {
          console.log('Checking authentication...');
          await checkAuth();
        }
        const currentUser = useAuthStore.getState().user;
        console.log('Current user after auth check:', currentUser);
        if (!isActive) return;
        if (!currentUser || !currentUser._id) {
          setError('Authentication failed. Please log in.');
          setNotification('Please log in to access toxicity prediction.');
          navigate('/login');
          return;
        }
        console.log('Fetching symptoms, products, and history...');
        await fetchSymptomsAndProducts();
        await fetchHistory();
        console.log('Initial data fetched successfully');
      } catch (err) {
        if (!isActive) return;
        console.error('Initialization error:', err);
        setError('Failed to verify authentication. Please try refreshing the page or logging in again.');
        setNotification('Please try again or contact support.');
      } finally {
        if (isActive) {
          setLoading(false);
        }
      }
    };
    initialize();
    return () => {
      isActive = false;
    };
  }, [checkAuth]);

  const predictToxicity = async (compoundString, symptomsString) => {
    console.log("Running toxicity prediction with payload:", { compound: compoundString, symptoms: symptomsString });
    try {
      const response = await axiosInstance.post('/toxicity/predict', { compound: compoundString, symptoms: symptomsString });
      console.log("Backend toxicity prediction response (Gemini-based):", response.data);
      return response.data.result;
    } catch (err) {
      console.error("Error calling backend toxicity API:", err.response?.data || err.message);
      throw new Error(err.response?.data?.message || "Failed to predict toxicity: " + err.message);
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!compound.trim()) {
      setError('Please select a valid molecule name');
      setNotification('Choose a molecule from the dropdown menu.');
      return;
    }
    if (symptomGroupIndex === '' || !symptomGroups[symptomGroupIndex]?.length) {
      setError('Please select a valid symptom group');
      setNotification('Choose a symptom group from the dropdown menu.');
      return;
    }
    if (!user?._id) {
      setError('User not authenticated');
      setNotification('Please log in to continue.');
      return;
    }
    setLoading(true);
    setError(null);
    setNotification(null);
    setGeminiAnalysis(null);
    setIsGeminiAnalysisOpen(false);
    setIsResultOpen(false);
    setIsSideEffectsOpen(false);

    try {
      const symptomsString = symptomGroups[symptomGroupIndex]?.join(", ") || "";
      const normalizeSymptoms = (symptoms) =>
        symptoms
          .toLowerCase()
          .split(/[,;]+/)
          .map(s => s.trim())
          .filter(s => s.length > 0)
          .sort()
          .join(", ");

      const normalizedCompound = compound.trim().toLowerCase();
      const normalizedSymptoms = normalizeSymptoms(symptomsString);
      const matchingHistoryItem = history.find(item => {
        const historyCompound = item.compoundName?.trim().toLowerCase() || item.smiles?.trim().toLowerCase();
        const historySymptoms = normalizeSymptoms(item.symptomsGrp || "");
        return historyCompound === normalizedCompound && historySymptoms === normalizedSymptoms;
      });

      if (matchingHistoryItem) {
        setActiveTab('saved');
        setSelectedHistoryItemId(matchingHistoryItem._id);
        setIsHistoryResultOpen(true);
        setIsHistorySideEffectsOpen(false);
        setIsHistoryGeminiOpen(!!matchingHistoryItem.geminiAnalysis);
        setNotification(`Result for "${compound}" with symptoms "${symptomsString}" already generated. View it in the Saved Results tab below.`);
        setTimeout(() => {
          const element = document.getElementById(`history-item-${matchingHistoryItem._id}`);
          if (element) {
            element.scrollIntoView({ behavior: 'smooth', block: 'start' });
          }
        }, 100);
        setLoading(false);
        return;
      }

      const toxicityData = await predictToxicity(compound, symptomsString);
      console.log("Toxicity Data:", toxicityData); // Debug
      if (!toxicityData) {
        throw new Error("No toxicity data returned from backend.");
      }

      const updatedResult = { ...toxicityData, symptoms: symptomsString };
      setResult(updatedResult);
      setIsResultOpen(true);

      setGeminiLoading(true);
      const analysisText = await getGeminiAnalysis(compound, symptomsString);
      setGeminiAnalysis(analysisText);

      const historyResponse = await axiosInstance.get('/toxicity/history');
      const latestHistory = historyResponse.data.history[0];
      if (latestHistory) {
        const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "{}");
        storedSymptomGroups[latestHistory._id] = symptomsString;
        localStorage.setItem(`symptomGroups_${user._id}`, JSON.stringify(storedSymptomGroups));
      }

      setIsGeminiAnalysisOpen(true);
      await fetchHistory();
    } catch (err) {
      console.error('Error in handleSubmit:', err);
      setError(err.message || 'Error predicting toxicity analysis');
      setNotification(err.response?.data?.suggestion || 'Please verify the molecule name and error, then try again.');
    } finally {
      setLoading(false);
      setGeminiLoading(false);
    }
  };

  const renderToxicityDetails = (result) => {
    if (!result) {
      return (
        <div className="text-center py-4">
          <h3 className="text-gray-400">No details available</h3>
        </div>
      );
    }

    return (
      <div className="space-y-4">
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div className="bg-gradient-to-r from-blue-500/10 to-purple-500/10 p-4 rounded-lg border border-blue-300/20">
            <h4 className="font-semibold text-blue-300 mb-2 flex items-center">
              <Activity className="w-4 h-4 mr-2" />
              Basic Information
            </h4>
            <p className="text-gray-300 text-sm-1"><span className="font-medium">Symptoms:</span> {result?.symptoms || 'Not available'}</p>
            <p className="text-gray-300 text-sm"><span className="font-medium">Compound:</span></p>
            <div className="bg-gray-800/50 p-2 rounded mt-1 overflow-x-auto">
              <code className="text-green-400 text-xs whitespace-nowrap">{result.smiles || 'Unknown'}</code>
            </div>
            {result.structuralAnalysis ? (
              <>
                <p className="text-gray-300 text-sm mt-2"><span class="font-medium">Molecular Weight:</span> {result.structuralAnalysis?.molecularWeight || 'Unknown'}</p>
                <p className="text-gray-300 text-sm"><span class="font-medium">LogP:</span> {result.structuralAnalysis?.logP || 'Unknown'}</p>
                <p className="text-gray-300 text-sm"><span class="font-medium">Polar Surface Area:</span> {result.structuralAnalysis?.polarSurfaceArea || 'Unknown'}</p>
              </>
            ) : (
              <p className="text-gray-400 text-sm mt-2">No structural analysis available.</p>
            )}
          </div>
          <div className="bg-gradient-to-r from-red-500/20 to-orange-500/20 p-4 rounded-lg border border-red-300/20">
            <h4 className="font-semibold text-red-300 mb-6 flex items-center">
              <Shield className="w-4 h-4 mr-2" />
              Toxicity Data
            </h4>
            <p className="text-gray-300 text-sm-1"><span class="font-medium">LD50:</span> {result?.acuteToxicity?.data?.LD50 || 'Unknown'}</p>
            <p className="text-gray-300 text-sm mb-1"><span class="font-medium">Toxicity Class:</span> {result?.acuteToxicity?.toxicityClass || 'Unknown'}</p>
            <p className="text-gray-300 text-sm mb-1"><span class="font-medium">Hepatotoxicity:</span> {result?.endpoints?.hepatotoxicity || 'Unknown'}</p>
            <p className="text-gray-300 text-sm"><span class="font-medium">Carcinogenicity:</span> {result?.endpoints?.carcinogenicity || 'Unknown'}</p>
            {result?.acuteToxicity?.supplemental && (
              <p className="text-gray-400 text-xs mt-2">italic mt {result?.acuteToxicity?.supplemental}</p>
            )}
          </div>
        </div>
      </div>
    );
  };

  const renderSideEffectsDetails = (sideEffects) => {
    if (!sideEffects || sideEffects.length === 0) {
      return (
        <div className="text-center py-8">
          <AlertTriangle className="w-12 h-12 text-yellow-400/50 mx-auto mb-3" />
          <p className="text-gray-400">
            No specific side effects data available. Inferred effects have been provided based on symptoms.
          </p>
        </div>
      );
    }

    return (
      <div className="space-y-3">
  {sideEffects.map((effect, index) => (
    <motion.div
      key={index}
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.1 }}
      className="bg-gradient-to-r from-yellow-500/10 to-orange-500/10 border-l-4 border-yellow-400 p-4 rounded-r-lg"
    >
      <p className="font-semibold text-yellow-300 flex items-center mb-2">
        <AlertTriangle className="w-4 h-4 mr-2" />
        {effect.name}
      </p>
      <p className="text-gray-300 text-sm leading-relaxed">{effect.description}</p>
      <p className="text-gray-400 text-xs mt-1">Severity: {effect.severity}</p>
      <p className="text-gray-400 text-xs">Organ System: {effect.organSystem}</p>
      <p className="text-gray-400 text-xs">Likelihood: {effect.likelihood}</p>
    </motion.div>
  ))}
</div>

    );
  };

  const renderGeminiAnalysis = (analysis) => {
    if (!analysis) {
      return (
        <div className="text-center py-8">
          <Info className="w-12 h-12 text-blue-400/50 mx-auto mb-3" />
          <p className="text-gray-400">No toxicity analysis available for this prediction.</p>
        </div>
      );
    }

    const {
      isInputValid,
      inputError,
      smiles,
      moleculeName,
      isNovelCompound,
      structuralAnalysis,
      acuteToxicity,
      endpoints,
      sideEffects,
      mechanisms,
      structureConcerns,
      safetyRecommendations,
      qsarAnalysis,
      toxicophoreAnalysis,
      developmentRecommendations
    } = analysis;

    return (
      <div className="space-y-6">
        {!isInputValid && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-yellow-500/10 p-6 rounded-lg border border-yellow-300/20"
          >
            <h4 className="font-semibold text-yellow-300 mb-3 flex items-center">
              <AlertTriangle className="w-5 h-5 mr-2" />
              Input Validation Error
            </h4>
            <p className="text-gray-300 text-sm mb-2">
              <span className="font-medium">Compound:</span> {moleculeName || 'Unknown'}
            </p>
            <p className="text-gray-300 text-sm mb-2">
              <span className="font-medium">Error:</span> {inputError || 'Invalid input'}
            </p>
          </motion.div>
        )}

        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-gradient-to-r from-blue-500/10 to-purple-500/10 p-6 rounded-lg border border-blue-300/20"
        >
          <h4 className="font-semibold text-blue-300 mb-3 flex items-center">
            <Info className="w-5 h-5 mr-2" />
            Basic Information
          </h4>
          <p className="text-gray-300 text-sm mb-2">
            <span className="font-medium">Compound:</span> {moleculeName || 'N/A'}
          </p>
          <p className="text-gray-300 text-sm mb-2">
            <span className="font-medium">SMILES:</span> {smiles || 'N/A'}
          </p>
          <p className="text-gray-300 text-sm mb-2">
            <span className="font-medium">Novel Compound:</span> {isNovelCompound ? 'Yes' : 'No'}
          </p>
          {structuralAnalysis && (
            <>
              <p className="text-gray-300 text-sm mb-2">
                <span className="font-medium">Molecular Weight:</span> {structuralAnalysis.molecularWeight || 'Unknown'}
              </p>
              <p className="text-gray-300 text-sm mb-2">
                <span className="font-medium">LogP:</span> {structuralAnalysis.logP || 'Unknown'}
              </p>
              <p className="text-gray-300 text-sm mb-2">
                <span className="font-medium">Polar Surface Area:</span> {structuralAnalysis.polarSurfaceArea || 'Unknown'}
              </p>
              <p className="text-gray-300 text-sm mb-2">
                <span className="font-medium">Functional Groups:</span> {structuralAnalysis.functionalGroups?.join(', ') || 'None identified'}
              </p>
            </>
          )}
        </motion.div>

        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-gradient-to-r from-red-500/10 to-orange-500/10 p-6 rounded-lg border border-red-300/20"
        >
          <h4 className="font-semibold text-red-300 mb-3 flex items-center">
            <Shield className="w-5 h-5 mr-2" />
            Acute Toxicity
          </h4>
          <p className="text-gray-300 text-sm mb-2">
            <span className="font-medium">LD50:</span> {acuteToxicity?.LD50 || 'Unknown'}
          </p>
          <p className="text-gray-300 text-sm mb-2">
            <span className="font-medium">Toxicity Class:</span> {acuteToxicity?.toxicityClass || 'Unknown'}
          </p>
          <p className="text-gray-300 text-sm mb-2">
            <span className="font-medium">Rationale:</span> {acuteToxicity?.rationale || 'Not provided'}
          </p>
          {acuteToxicity?.supplemental && (
            <p className="text-gray-400 text-xs italic">{acuteToxicity.supplemental}</p>
          )}
        </motion.div>

        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-gradient-to-r from-indigo-500/10 to-purple-500/10 p-6 rounded-lg border border-indigo-300/20"
        >
          <h4 className="font-semibold text-indigo-300 mb-3 flex items-center">
            <Beaker className="w-5 h-5 mr-2" />
            Toxicological Endpoints
          </h4>
          <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
            {Object.entries(endpoints || {}).map(([key, value]) => (
              key !== 'supplemental' && (
                <li key={key}>
                  {key.replace(/([A-Z])/g, ' $1').replace(/^./, str => str.toUpperCase())}: {value || 'Unknown'}
                </li>
              )
            ))}
          </ul>
          {endpoints?.supplemental && (
            <p className="text-gray-400 text-xs italic mt-2">{endpoints.supplemental}</p>
          )}
        </motion.div>

        {mechanisms?.length > 0 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-gradient-to-r from-green-500/10 to-teal-500/10 p-6 rounded-lg border border-green-300/20"
          >
            <h4 className="font-semibold text-green-300 mb-3 flex items-center">
              <Activity className="w-5 h-5 mr-2" />
              Mechanisms of Toxicity
            </h4>
            <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
              {mechanisms.map((mechanism, index) => (
                <li key={index}>
                  <span className="font-medium">{mechanism.feature}:</span> {mechanism.pathway} ({mechanism.toxicityType})
                </li>
              ))}
            </ul>
          </motion.div>
        )}

        {structureConcerns?.length > 0 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-gradient-to-r from-orange-500/10 to-red-500/10 p-6 rounded-lg border border-orange-300/20"
          >
            <h4 className="font-semibold text-orange-300 mb-3 flex items-center">
              <AlertCircle className="w-5 h-5 mr-2" />
              Structural Concerns
            </h4>
            <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
              {structureConcerns.map((concern, index) => (
                <li key={index}>
                  <span className="font-medium">{concern.substructure}:</span> {concern.concern} (Risk: {concern.riskLevel})
                </li>
              ))}
            </ul>
          </motion.div>
        )}

        {safetyRecommendations?.length > 0 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-gradient-to-r from-emerald-500/10 to-blue-500/10 p-6 rounded-lg border border-emerald-300/20"
          >
            <h4 className="font-semibold text-emerald-300 mb-3 flex items-center">
              <Shield className="w-5 h-5 mr-2" />
              Safety Recommendations
            </h4>
            <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
              {safetyRecommendations.map((recommendation, index) => (
                <li key={index}>
                  <span className="font-medium">Test:</span> {recommendation.test} <br />
                  <span className="font-medium">Priority:</span> {recommendation.priority} <br />
                  <span className="font-medium">Rationale:</span> {recommendation.rationale}
                </li>
              ))}
            </ul>
          </motion.div>
        )}

        {(qsarAnalysis?.properties?.length > 0 || toxicophoreAnalysis?.length > 0) && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-gradient-to-r from-cyan-500/10 to-teal-500/10 p-6 rounded-lg border border-cyan-300/20"
          >
            <h4 className="font-semibold text-cyan-300 mb-3 flex items-center">
              <Info className="w-5 h-5 mr-2" />
              QSAR and Toxicophore Analysis
            </h4>
            {qsarAnalysis?.properties?.length > 0 && (
              <div className="mb-4">
                <p className="text-gray-300 text-sm font-medium mb-2">Physicochemical Properties:</p>
                <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
                  {qsarAnalysis.properties.map((prop, index) => (
                    <li key={index}>
                      <span className="font-medium">{prop.name}:</span> {prop.predictedValue} ({prop.toxicologicalImplication})
                    </li>
                  ))}
                </ul>
              </div>
            )}
            {qsarAnalysis?.overallRisk && (
              <p className="text-gray-300 text-sm mb-2">
                <span className="font-medium">Overall Risk:</span> {qsarAnalysis.overallRisk}
              </p>
            )}
            {qsarAnalysis?.symptomContext?.length > 0 && (
              <div className="mb-4">
                <p className="text-gray-300 text-sm font-medium mb-2">Symptom Contextual Analysis:</p>
                <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
                  {qsarAnalysis.symptomContext.map((context, index) => (
                    <li key={index}>
                      <span className="font-medium">{context.symptom}:</span> {context.structuralInsight} (Risk: {context.riskAssessment})
                    </li>
                  ))}
                </ul>
              </div>
            )}
            {toxicophoreAnalysis?.length > 0 && (
              <div>
                <p className="text-gray-300 text-sm font-medium mb-2">Toxicophore Analysis:</p>
                <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
                  {toxicophoreAnalysis.map((toxicophore, index) => (
                    <li key={index}>
                      <span className="font-medium">{toxicophore.substructure}:</span> {toxicophore.concern} (Prevalence: {toxicophore.prevalence || 'Unknown'}, Mitigation: {toxicophore.mitigation || 'None'})
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </motion.div>
        )}

        {developmentRecommendations?.length > 0 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="bg-gradient-to-r from-gray-500/10 to-gray-600/10 p-6 rounded-lg border border-gray-300/20"
          >
            <h4 className="font-semibold text-gray-300 mb-3 flex items-center">
              <Info className="w-5 h-5 mr-2" />
              Development Recommendations
            </h4>
            <ul className="list-disc list-inside text-gray-300 text-sm space-y-1">
              {developmentRecommendations.map((recommendation, index) => (
                <li key={index}>
                  <span className="font-medium">Phase:</span> {recommendation.phase} <br />
                  <span className="font-medium">Recommendation:</span> {recommendation.recommendation} <br />
                  <span className="font-medium">Importance:</span> {recommendation.importance}
                </li>
              ))}
            </ul>
          </motion.div>
        )}
      </div>
    );
  };

  const handleViewHistoryItem = (item, e) => {
    e.stopPropagation();
    console.log("Viewing history item:", item);
    if (selectedHistoryItemId === item._id) {
      setSelectedHistoryItemId(null);
      setIsHistoryResultOpen(false);
      setIsHistorySideEffectsOpen(false);
      setIsHistoryGeminiOpen(false);
    } else {
      setSelectedHistoryItemId(item._id);
      setIsHistoryResultOpen(true);
      setIsHistorySideEffectsOpen(false);
      setIsHistoryGeminiOpen(!!item.geminiAnalysis);
    }
  };

  return (
    <div className="min-h-screen">
      <div className="container mx-auto px-4 py-12 max-w-6xl">
        <motion.div
          initial={{ opacity: 0, y: -30 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.6 }}
          className="text-center mb-12"
        >
          <h3 className="text-5xl font-bold bg-gradient-to-r from-emerald-400 to-blue-400 bg-clip-text text-transparent mb-4">
            Toxicity & Side Effects Prediction
          </h3>
          <div className="flex flex-wrap justify-center items-center gap-6 text-gray-400">
            <span className="flex items-center">
              {/* <Shield className="w-4 h-4 mr-2" /> */}
              <p className="text-sm text-accent-secondary font-mono">POWERED BY GEMINI</p>
            </span>
          </div>
        </motion.div>

        <div className="flex gap-4 justify-center mb-6">
          <motion.button
            className={`flex items-center gap-2 py-3 px-6 rounded-lg transition-all duration-300 text-sm font-medium ${activeTab === 'predict'
                ? 'bg-gradient-to-r from-emerald-500 to-blue-500 text-white shadow-md'
                : 'bg-gray-800/50 text-gray-400 hover:text-white hover:bg-gray-700/50 border border-gray-700/50'
              }`}
            onClick={() => setActiveTab('predict')}
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
            initial={{ opacity: 0, x: -10 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.2, duration: 0.3 }}
          >
            <Beaker className="w-5 h-5" />
            <span>Predict Toxicity & Side Effects</span>
          </motion.button>

          <motion.button
            className={`flex items-center gap-2 py-3 px-6 rounded-lg transition-all duration-300 text-sm font-medium ${activeTab === 'saved'
                ? 'bg-gradient-to-r from-emerald-500 to-blue-500 text-white shadow-md'
                : 'bg-gray-800/50 text-gray-400 hover:text-white hover:bg-gray-700/50 border border-gray-700/50'
              }`}
            onClick={() => setActiveTab('saved')}
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
            initial={{ opacity: 0, x: 10 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.2, duration: 0.3 }}
          >
            <Activity className="w-5 h-5" />
            <span>Saved Results</span>
          </motion.button>
        </div>

        <AnimatePresence>
          {activeTab === 'predict' && (
            <motion.div
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -20 }}
              transition={{ duration: 0.4 }}
              className="bg-gray-800/50 backdrop-blur-lg p-8 rounded-xl shadow-2xl border border-gray-700/50"
            >
              <form onSubmit={handleSubmit} className="space-y-6">
                <div className="space-y-2">
                  <label className="block text-gray-300 font-medium">Molecule Name</label>
                  <select
                    value={compound}
                    onChange={(e) => setCompound(e.target.value)}
                    className="w-full p-3 rounded-lg bg-gray-900/50 border border-gray-600 text-gray-300 focus:outline-none focus:ring-2 focus:ring-emerald-500 transition-all"
                    disabled={loading || geminiLoading}
                  >
                    <option value="">Select a molecule</option>
                    {productCompoundGroups[symptomGroupIndex]?.map((comp, index) => (
                      <option key={index} value={comp}>
                        {comp}
                      </option>
                    ))}
                  </select>
                </div>

                <div className="space-y-2">
                  <label className="block text-gray-300 font-medium">Symptom Group</label>
                  <select
                    value={symptomGroupIndex}
                    onChange={(e) => {
                      setSymptomGroupIndex(e.target.value);
                      setCompound(productCompoundGroups[e.target.value]?.[0] || '');
                    }}
                    className="w-full p-3 rounded-lg bg-gray-900/50 border border-gray-600 text-gray-300 focus:outline-none focus:ring-2 focus:ring-emerald-500 transition-all"
                    disabled={loading || geminiLoading}
                  >
                    <option value="">Select a symptom group</option>
                    {symptomGroups.map((group, index) => (
                      <option key={index} value={index}>
                        {group.join(', ') || `Group ${index + 1}`}
                      </option>
                    ))}
                  </select>
                </div>

                {error && (
                  <motion.div
                    initial={{ opacity: 0, y: 10 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="bg-red-500/10 border-l-4 border-red-500 p-4 rounded-r-lg flex items-start"
                  >
                    <AlertCircle className="w-5 h-5 text-red-400 mr-2 mt-1" />
                    <div>
                      <p className="text-red-300 font-medium">{error}</p>
                      {notification && <p className="text-gray-400 text-sm">{notification}</p>}
                    </div>
                  </motion.div>
                )}

                <motion.button
                  type="submit"
                  disabled={loading || geminiLoading}
                  className={`w-full py-3 px-6 rounded-lg font-medium transition-all duration-300 flex items-center justify-center gap-2 ${loading || geminiLoading
                      ? 'bg-gray-600 cursor-not-allowed'
                      : 'bg-gradient-to-r from-emerald-500 to-blue-500 hover:from-emerald-600 hover:to-blue-600 text-white'
                    }`}
                  whileHover={{ scale: loading || geminiLoading ? 1 : 1.02 }}
                  whileTap={{ scale: loading || geminiLoading ? 1 : 0.98 }}
                >
                  {loading || geminiLoading ? (
                    <>
                      <RefreshCw className="w-5 h-5 animate-spin" />
                      Processing...
                    </>
                  ) : (
                    <>
                      <Beaker className="w-5 h-5" />
                      Predict Toxicity
                    </>
                  )}
                </motion.button>
              </form>

              {result && (
                <div className="mt-8 space-y-6">
                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="bg-gray-900/50 p-6 rounded-lg border border-gray-700/50"
                  >
                    <button
                      onClick={() => setIsResultOpen(!isResultOpen)}
                      className="w-full flex justify-between items-center text-gray-300 font-semibold"
                    >
                      <span className="flex items-center">
                        <Info className="w-5 h-5 mr-2" />
                        Toxicity Details
                      </span>
                      {isResultOpen ? <ChevronUp className="w-5 h-5" /> : <ChevronDown className="w-5 h-5" />}
                    </button>
                    <AnimatePresence>
                      {isResultOpen && (
                        <motion.div
                          initial={{ height: 0, opacity: 0 }}
                          animate={{ height: 'auto', opacity: 1 }}
                          exit={{ height: 0, opacity: 0 }}
                          transition={{ duration: 0.3 }}
                          className="mt-4 overflow-hidden"
                        >
                          {renderToxicityDetails(result)}
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </motion.div>

                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="bg-gray-900/50 p-6 rounded-lg border border-gray-700/50"
                  >
                    <button
                      onClick={() => setIsSideEffectsOpen(!isSideEffectsOpen)}
                      className="w-full flex justify-between items-center text-gray-300 font-semibold"
                    >
                      <span className="flex items-center">
                        <AlertTriangle className="w-5 h-5 mr-2" />
                        Side Effects
                      </span>
                      {isSideEffectsOpen ? <ChevronUp className="w-5 h-5" /> : <ChevronDown className="w-5 h-5" />}
                    </button>
                    <AnimatePresence>
                      {isSideEffectsOpen && (
                        <motion.div
                          initial={{ height: 0, opacity: 0 }}
                          animate={{ height: 'auto', opacity: 1 }}
                          exit={{ height: 0, opacity: 0 }}
                          transition={{ duration: 0.3 }}
                          className="mt-4 overflow-hidden"
                        >
                          {renderSideEffectsDetails(result.sideEffects)}
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </motion.div>

                  <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="bg-gray-900/50 p-6 rounded-lg border border-gray-700/50"
                  >
                    <button
                      onClick={() => setIsGeminiAnalysisOpen(!isGeminiAnalysisOpen)}
                      className="w-full flex justify-between items-center text-gray-300 font-semibold"
                    >
                      <span className="flex items-center">
                        <Shield className="w-5 h-5 mr-2" />
                        Gemini Toxicity Analysis
                      </span>
                      {isGeminiAnalysisOpen ? <ChevronUp className="w-5 h-5" /> : <ChevronDown className="w-5 h-5" />}
                    </button>
                    <AnimatePresence>
                      {isGeminiAnalysisOpen && (
                        <motion.div
                          initial={{ height: 0, opacity: 0 }}
                          animate={{ height: 'auto', opacity: 1 }}
                          exit={{ height: 0, opacity: 0 }}
                          transition={{ duration: 0.3 }}
                          className="mt-4 overflow-hidden"
                        >
                          {geminiLoading ? (
                            <div className="text-center py-8">
                              <RefreshCw className="w-8 h-8 text-emerald-400 animate-spin mx-auto mb-3" />
                              <p className="text-gray-400">Generating Gemini analysis...</p>
                            </div>
                          ) : (
                            renderGeminiAnalysis(geminiAnalysis)
                          )}
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </motion.div>
                </div>
              )}
            </motion.div>
          )}

          {activeTab === 'saved' && (
            <motion.div
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -20 }}
              transition={{ duration: 0.4 }}
              className="bg-gray-800/50 backdrop-blur-lg p-8 rounded-xl shadow-2xl border border-gray-700/50"
            >
              <h3 className="text-2xl font-semibold text-gray-300 mb-6">Saved Toxicity Predictions</h3>
              {history.length === 0 ? (
                <div className="text-center py-8">
                  <Info className="w-12 h-12 text-blue-400/50 mx-auto mb-3" />
                  <p className="text-gray-400">No saved predictions found. Make a new prediction to get started.</p>
                </div>
              ) : (
                <div className="space-y-4">
                  {history.map((item) => (
                    <motion.div
                      key={item._id}
                      id={`history-item-${item._id}`}
                      initial={{ opacity: 0, y: 10 }}
                      animate={{ opacity: 1, y: 0 }}
                      className="bg-gray-900/50 p-6 rounded-lg border border-gray-700/50"
                    >
                      <div className="flex justify-between items-center">
                        <div>
                          <h4 className="text-gray-300 font-semibold">
                            {item.compoundName || item.smiles || 'Unknown Compound'}
                          </h4>
                          <p className="text-gray-400 text-sm">{item.symptomsGrp}</p>
                        </div>
                        <motion.button
                          onClick={(e) => handleViewHistoryItem(item, e)}
                          className="flex items-center gap-2 px-4 py-2 rounded-lg bg-emerald-500/20 text-emerald-300 hover:bg-emerald-500/30 transition-all"
                          whileHover={{ scale: 1.02 }}
                          whileTap={{ scale: 0.98 }}
                        >
                          <Eye className="w-4 h-4" />
                          {selectedHistoryItemId === item._id ? 'Hide' : 'View'}
                        </motion.button>
                      </div>

                      {selectedHistoryItemId === item._id && (
                        <div className="mt-4 space-y-6">
                          <div className="bg-gray-900/70 p-6 rounded-lg">
                            <button
                              onClick={() => setIsHistoryResultOpen(!isHistoryResultOpen)}
                              className="w-full flex justify-between items-center text-gray-300 font-semibold"
                            >
                              <span className="flex items-center">
                                <Info className="w-5 h-5 mr-2" />
                                Toxicity Details
                              </span>
                              {isHistoryResultOpen ? <ChevronUp className="w-5 h-5" /> : <ChevronDown className="w-5 h-5" />}
                            </button>
                            <AnimatePresence>
                              {isHistoryResultOpen && (
                                <motion.div
                                  initial={{ height: 0, opacity: 0 }}
                                  animate={{ height: 'auto', opacity: 1 }}
                                  exit={{ height: 0, opacity: 0 }}
                                  transition={{ duration: 0.3 }}
                                  className="mt-4 overflow-hidden"
                                >
                                  {renderToxicityDetails(item.toxicityResult)}
                                </motion.div>
                              )}
                            </AnimatePresence>
                          </div>

                          <div className="bg-gray-900/70 p-6 rounded-lg">
                            <button
                              onClick={() => setIsHistorySideEffectsOpen(!isHistorySideEffectsOpen)}
                              className="w-full flex justify-between items-center text-gray-300 font-semibold"
                            >
                              <span className="flex items-center">
                                <AlertTriangle className="w-5 h-5 mr-2" />
                                Side Effects
                              </span>
                              {isHistorySideEffectsOpen ? <ChevronUp className="w-5 h-5" /> : <ChevronDown className="w-5 h-5" />}
                            </button>
                            <AnimatePresence>
                              {isHistorySideEffectsOpen && (
                                <motion.div
                                  initial={{ height: 0, opacity: 0 }}
                                  animate={{ height: 'auto', opacity: 1 }}
                                  exit={{ height: 0, opacity: 0 }}
                                  transition={{ duration: 0.3 }}
                                  className="mt-4 overflow-hidden"
                                >
                                  {renderSideEffectsDetails(item.toxicityResult?.sideEffects)}
                                </motion.div>
                              )}
                            </AnimatePresence>
                          </div>

                          <div className="bg-gray-900/70 p-6 rounded-lg">
                            <button
                              onClick={() => setIsHistoryGeminiOpen(!isHistoryGeminiOpen)}
                              className="w-full flex justify-between items-center text-gray-300 font-semibold"
                            >
                              <span className="flex items-center">
                                <Shield className="w-5 h-5 mr-2" />
                                Gemini Toxicity Analysis
                              </span>
                              {isHistoryGeminiOpen ? <ChevronUp className="w-5 h-5" /> : <ChevronDown className="w-5 h-5" />}
                            </button>
                            <AnimatePresence>
                              {isHistoryGeminiOpen && (
                                <motion.div
                                  initial={{ height: 0, opacity: 0 }}
                                  animate={{ height: 'auto', opacity: 1 }}
                                  exit={{ height: 0, opacity: 0 }}
                                  transition={{ duration: 0.3 }}
                                  className="mt-4 overflow-hidden"
                                >
                                  {renderGeminiAnalysis(item.geminiAnalysis)}
                                </motion.div>
                              )}
                            </AnimatePresence>
                          </div>
                        </div>
                      )}
                    </motion.div>
                  ))}
                </div>
              )}
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </div>
  );
};

export default ToxicityPrediction;

