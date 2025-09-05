import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import * as NGL from 'ngl';
import { Atom, Database, History, Loader2, AlertCircle, CheckCircle, Eye, EyeOff, RotateCcw, Zap } from 'lucide-react';
import toast from 'react-hot-toast';
const AlphaFoldExplorer = () => {
  const [uniprotId, setUniprotId] = useState('');
  const [predictionResult, setPredictionResult] = useState(null);
  const [jobStatus, setJobStatus] = useState(null);
  const [summary, setSummary] = useState(null);
  const [annotations, setAnnotations] = useState(null);
  const [previousJobs, setPreviousJobs] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [activeTab, setActiveTab] = useState('predict');
  const [structureRendered, setStructureRendered] = useState(false);
  const [nglComponents, setNglComponents] = useState(null);
  const [viewportLoading, setViewportLoading] = useState(false);
  const stageRef = useRef(null);

  const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
  const axiosInstance = axios.create({
    baseURL: import.meta.mode === "development" ? API_BASE_URL : '/api',
    withCredentials: true,
  });

  // Handle prediction submission
  const handlePredictionSubmit = async (e) => {
    e.preventDefault();
    if (!uniprotId || !/^[A-Z0-9]{5,10}$/i.test(uniprotId)) {
      setError('Please enter a valid UniProt ID (5-10 alphanumeric characters)');
      return;
    }
    setLoading(true);
    setError('');
    setPredictionResult(null);
    setJobStatus(null);
    setStructureRendered(false);
    
    try {
      const response = await axiosInstance.post(`${API_BASE_URL}/alphafold/predict`, { uniprot_id: uniprotId }, {
        timeout: 30000,
        headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
      });
      setPredictionResult(response.data);
    } catch (err) {
      setError(err.response?.data?.error || err.message || 'Prediction submission failed');
      setLoading(false);
    }
  };
  // Custom toast theme
  const toastOptions = {
    style: {
      background: '#172A45', // secondary
      color: '#E0E0E0', // text-primary
      border: '1px solid #5E81F4', // accent-secondary
      borderRadius: '8px',
      padding: '12px',
      fontFamily: 'Roboto, Open Sans, sans-serif', // body font
    },
    success: {
      style: {
        borderColor: '#70E000', // success
      },
      iconTheme: {
        primary: '#70E000', // success
        secondary: '#E0E0E0', // text-primary
      },
    },
  };
  // Handle UniProt data fetch
  const handleUniProtSubmit = async (e) => {
    e.preventDefault();
    if (!uniprotId || !/^[A-Z0-9]{5,10}$/i.test(uniprotId)) {
      setError('Please enter a valid UniProt ID (5-10 alphanumeric characters)');
      toast.error("Please enter a valid UniProt ID (5-10 alphanumeric characters)", toastOptions);
      return;
    }
    setLoading(true);
    setError('');
    setSummary(null);
    setAnnotations(null);

    try {
      const [summaryRes, annotationsRes] = await Promise.all([
        axiosInstance.get(`${API_BASE_URL}/alphafold/uniprot/summary/${uniprotId}`, { timeout: 30000 }),
        axiosInstance.get(`${API_BASE_URL}/alphafold/uniprot/annotations/${uniprotId}`, { timeout: 30000 }),
      ]);
      setSummary(summaryRes.data);
      setAnnotations(annotationsRes.data);
    } catch (err) {
      setError(err.response?.data?.error || err.message || 'Failed to fetch UniProt data');
      toast.error('Failed to fetch UniProt data', toastOptions);
    } finally {
      setLoading(false);
    }
  };

  // Fetch previous jobs
  const fetchPreviousJobs = async () => {
    try {
      const response = await axiosInstance.get(`${API_BASE_URL}/alphafold/previous-jobs`, {
        headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
      });
      setPreviousJobs(response.data);
    } catch (err) {
      setError('Failed to fetch previous jobs');
    }
  };

  // Poll job status
  useEffect(() => {
    let intervalId;
    if (predictionResult && predictionResult.jobId && !jobStatus?.pdbUrl && !structureRendered) {
      setLoading(true);
      const pollStatus = async () => {
        try {
          const response = await axiosInstance.get(`${API_BASE_URL}/alphafold/status/${predictionResult.jobId}`, {
            headers: { Authorization: `Bearer ${localStorage.getItem('token')}` }
          });
          setJobStatus(response.data);
          if (response.data.status === 'completed' || response.data.status === 'failed') {
            clearInterval(intervalId);
            fetchPreviousJobs();
          }
        } catch (err) {
          setError(err.response?.data?.error || 'Failed to fetch job status');
          clearInterval(intervalId);
          setLoading(false);
        }
      };

      pollStatus();
      intervalId = setInterval(pollStatus, 5000);
    }

    return () => {
      if (intervalId) clearInterval(intervalId);
    };
  }, [predictionResult, structureRendered]);

  // Initialize NGL stage and load structure
  useEffect(() => {
    if (jobStatus?.pdbUrl && !structureRendered) {
      if (!stageRef.current) {
        stageRef.current = new NGL.Stage('viewport', {
          backgroundColor: '#0A192F',
          quality: 'high',
          antialias: true,
          clipNear: 0,
          clipFar: 100,
          clipDist: 10
        });
      } else {
        stageRef.current.removeAllComponents();
      }

      const handleResize = () => {
        stageRef.current.handleResize();
      };

      window.addEventListener('resize', handleResize);

      setViewportLoading(true);
      setError(null);

      stageRef.current.loadFile(jobStatus.pdbUrl, { ext: 'pdb' }).then((structure) => {
        structure.addRepresentation('cartoon', {
          colorScheme: 'residueindex',
          colorScale: 'rainbow',
          opacity: 0.85,
          side: 'front'
        });

        structure.addRepresentation('cartoon', {
          colorScheme: 'sstruc',
          opacity: 0.7,
          visible: false,
          name: 'secondary_structure'
        });

        structure.addRepresentation('ball+stick', {
          sele: 'hetero and not water',
          colorScheme: 'element',
          radius: 0.3,
          multipleBond: 'symmetric'
        });

        structure.addRepresentation('surface', {
          sele: 'protein',
          colorScheme: 'electrostatic',
          opacity: 0.5,
          visible: false,
          name: 'surface_view'
        });

        structure.addRepresentation('licorice', {
          sele: 'CYS or HIS or ASP or GLU or LYS or ARG',
          colorValue: 'red',
          radius: 0.25,
          name: 'active_site'
        });

        structure.autoView(500);

        setNglComponents({
          main: structure,
          secondary: stageRef.current.getRepresentationsByName('secondary_structure')[0],
          surface: stageRef.current.getRepresentationsByName('surface_view')[0],
          activeSite: stageRef.current.getRepresentationsByName('active_site')[0]
        });

        setStructureRendered(true);
        setViewportLoading(false);
        setLoading(false);
      }).catch((err) => {
        console.error('Structure loading failed:', err);
        setError(`Failed to load structure: ${err.message}`);
        setViewportLoading(false);
        setLoading(false);
        setStructureRendered(true);
      });

      return () => {
        window.removeEventListener('resize', handleResize);
      };
    }
  }, [jobStatus?.pdbUrl, structureRendered]);

  // Reset structure and prediction data when switching tabs
  useEffect(() => {
    setPredictionResult(null);
    setJobStatus(null);
    setStructureRendered(false);
    setNglComponents(null);
    setError('');
    setLoading(false);
    setViewportLoading(false);

    if (stageRef.current) {
      stageRef.current.removeAllComponents();
    }
  }, [activeTab]);

  // Initial fetch of previous jobs
  useEffect(() => {
    fetchPreviousJobs();
  }, []);

  // Cleanup NGL stage on component unmount
  useEffect(() => {
    return () => {
      if (stageRef.current) {
        stageRef.current.dispose();
        stageRef.current = null;
      }
    };
  }, []);

  const loadPreviousStructure = (pdbUrl) => {
    setJobStatus({ pdbUrl });
    setStructureRendered(false);
  };

  const toggleRepresentation = (repName) => {
    if (nglComponents && nglComponents[repName]) {
      nglComponents[repName].toggleVisibility();
    }
  };

  const StatusIcon = ({ status }) => {
    switch (status) {
      case 'completed':
        return <CheckCircle className="w-5 h-5 text-success" />;
      case 'failed':
        return <AlertCircle className="w-5 h-5 text-error" />;
      case 'processing':
        return <Loader2 className="w-5 h-5 text-accent animate-spin" />;
      default:
        return <Loader2 className="w-5 h-5 text-accent-secondary animate-spin" />;
    }
  };

  const TabButton = ({ id, icon: Icon, label, isActive, onClick }) => (
    <button
      onClick={onClick}
      className={`flex items-center gap-2 px-4 py-3 rounded-lg font-label font-semibold text-sm transition-all duration-300 relative overflow-hidden group ${
        isActive 
          ? 'bg-gradient-to-r from-accent to-accent-secondary text-primary shadow-lg shadow-accent/25' 
          : 'bg-secondary/50 text-text-secondary hover:bg-secondary hover:text-text-primary border border-secondary'
      }`}
    >
      <Icon className="w-4 h-4" />
      {label}
      {isActive && (
        <div className="absolute inset-0 bg-gradient-to-r from-accent/20 to-accent-secondary/20 animate-pulse-glow"></div>
      )}
    </button>
  );

  return (
    <div className="min-h-screen  font-body">
      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="text-center mb-12 animate-fade-in">
          <div className="flex items-center justify-center gap-3 mb-4">
            <div className="p-3 bg-gradient-to-br from-accent to-accent-secondary rounded-xl shadow-lg shadow-accent/25">
              <Atom className="w-8 h-8 text-primary" />
            </div>
            <h1 className="text-4xl md:text-5xl font-heading font-bold bg-gradient-to-r from-accent to-accent-secondary bg-clip-text text-transparent">
              AlphaFold Explorer
            </h1>
          </div>
          <p className="text-text-secondary font-label text-lg max-w-2xl mx-auto">
            Advanced protein structure 3D prediction with GOOGLE Deepmind AlphaFold
          </p>
        </div>

        {/* Navigation Tabs */}
        <div className="flex flex-wrap justify-center gap-3 mb-8 animate-slide-up">
          <TabButton
            id="predict"
            icon={Zap}
            label="Predict Structure"
            isActive={activeTab === 'predict'}
            onClick={() => setActiveTab('predict')}
          />
          <TabButton
            id="uniprot"
            icon={Database}
            label="UniProt Analysis"
            isActive={activeTab === 'uniprot'}
            onClick={() => setActiveTab('uniprot')}
          />
          <TabButton
            id="previous"
            icon={History}
            label="Previous Jobs"
            isActive={activeTab === 'previous'}
            onClick={() => setActiveTab('previous')}
          />
        </div>

        {/* Main Content */}
        <div className="max-w-6xl mx-auto">
          <div className="bg-gradient-to-br from-secondary via-secondary to-secondary/80 rounded-2xl shadow-2xl border border-accent/20 p-6 md:p-8 backdrop-blur-sm animate-slide-up">
            
            {/* Input Form */}
            {activeTab !== 'previous' && (
              <form onSubmit={activeTab === 'predict' ? handlePredictionSubmit : handleUniProtSubmit} className="mb-8">
                <div className="mb-6">
                  <label className="block text-text-primary font-label font-semibold mb-3 text-lg" htmlFor="uniprotId">
                    UniProt Identifier
                  </label>
                  <div className="relative group">
                    <input
                      type="text"
                      id="uniprotId"
                      value={uniprotId}
                      onChange={(e) => setUniprotId(e.target.value.toUpperCase())}
                      className="w-full p-4 bg-primary border border-accent/30 rounded-xl focus:outline-none focus:ring-2 focus:ring-accent focus:border-accent text-text-primary placeholder-text-secondary font-code text-lg transition-all duration-300 group-hover:border-accent/50"
                      placeholder="Enter UniProt ID (e.g., P00520)"
                      required
                    />
                    <div className="absolute inset-0 rounded-xl bg-gradient-to-r from-accent/5 to-accent-secondary/5 opacity-0 group-hover:opacity-100 transition-opacity duration-300 pointer-events-none"></div>
                  </div>
                  <div className="mt-3 flex flex-wrap gap-2">
                    <span className="text-text-secondary font-label text-sm">Examples:</span>
                    {['P00520', 'P69905', 'P0DTD1'].map(id => (
                      <button
                        key={id}
                        type="button"
                        onClick={() => setUniprotId(id)}
                        className="px-3 py-1 bg-accent/10 text-accent border border-accent/30 rounded-md hover:bg-accent/20 transition-colors duration-200 font-code text-sm"
                      >
                        {id}
                      </button>
                    ))}
                  </div>
                </div>
                
                <button
                  type="submit"
                  disabled={loading}
                  className="w-full bg-gradient-to-r from-accent to-accent-secondary text-primary px-6 py-4 rounded-xl font-label font-bold text-lg hover:shadow-lg hover:shadow-accent/25 disabled:opacity-50 disabled:cursor-not-allowed transition-all duration-300 transform hover:scale-[1.02] relative overflow-hidden group"
                >
                  <div className="flex items-center justify-center gap-3">
                    {loading ? (
                      <>
                        <Loader2 className="w-5 h-5 animate-spin" />
                        Processing...
                      </>
                    ) : (
                      <>
                        {activeTab === 'predict' ? <Zap className="w-5 h-5" /> : <Database className="w-5 h-5" />}
                        {activeTab === 'predict' ? 'Submit Prediction' : 'Analyze Protein'}
                      </>
                    )}
                  </div>
                  <div className="absolute inset-0 bg-gradient-to-r from-white/20 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-300"></div>
                </button>
              </form>
            )}

            {/* Previous Jobs Tab */}
            {activeTab === 'previous' && (
              <div className="animate-slide-up">
                <div className="flex items-center gap-3 mb-6">
                  <History className="w-6 h-6 text-accent" />
                  <h2 className="text-2xl font-heading font-bold text-text-primary">Previous Predictions</h2>
                </div>
                {previousJobs.length === 0 ? (
                  <div className="text-center py-12">
                    <div className="w-16 h-16 bg-secondary/50 rounded-full flex items-center justify-center mx-auto mb-4">
                      <History className="w-8 h-8 text-text-secondary" />
                    </div>
                    <p className="text-text-secondary font-label text-lg">No previous predictions found</p>
                  </div>
                ) : (
                  <div className="grid gap-4">
                    {previousJobs.map(job => (
                      <div
                        key={job.jobId}
                        className="bg-primary/50 border border-accent/20 p-6 rounded-xl hover:border-accent/40 transition-all duration-300 cursor-pointer group transform hover:scale-[1.01]"
                        onClick={() => loadPreviousStructure(job.pdbUrl)}
                      >
                        <div className="flex items-start justify-between mb-4">
                          <div className="flex items-center gap-3">
                            <StatusIcon status={job.status} />
                            <div>
                              <p className="text-text-primary font-label font-semibold text-lg">
                                UniProt ID: <span className="font-code text-accent">{job.uniprotId}</span>
                              </p>
                              <p className="text-text-secondary font-label text-sm">
                                Job ID: {job.jobId}
                              </p>
                            </div>
                          </div>
                          <span className={`px-3 py-1 rounded-full text-xs font-label font-semibold ${
                            job.status === 'completed' ? 'bg-success/20 text-success border border-success/30' :
                            job.status === 'failed' ? 'bg-error/20 text-error border border-error/30' : 
                            'bg-accent/20 text-accent border border-accent/30'
                          }`}>
                            {job.status}
                          </span>
                        </div>
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm">
                          <div>
                            <span className="text-text-secondary font-label">Created:</span>
                            <span className="text-text-primary font-code ml-2">{new Date(job.createdAt).toLocaleString()}</span>
                          </div>
                          {job.completedAt && (
                            <div>
                              <span className="text-text-secondary font-label">Completed:</span>
                              <span className="text-text-primary font-code ml-2">{new Date(job.completedAt).toLocaleString()}</span>
                            </div>
                          )}
                        </div>
                        {job.error && (
                          <div className="mt-4 p-3 bg-error/10 border border-error/30 rounded-lg">
                            <p className="text-error font-label text-sm">Error: {job.error}</p>
                          </div>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            )}

            {/* Error Display */}
            {error && (
              <div className="mb-6 p-4 bg-error/10 border border-error/30 text-error rounded-xl animate-slide-up">
                <div className="flex items-center gap-3">
                  <AlertCircle className="w-5 h-5" />
                  <span className="font-label">{error}</span>
                </div>
              </div>
            )}

            {/* Prediction Results */}
            {activeTab === 'predict' && (predictionResult || jobStatus) && (
              <div className="bg-gradient-to-br from-accent/5 to-accent-secondary/5 rounded-xl p-6 border border-accent/20 animate-slide-up">
                <div className="flex items-center gap-3 mb-6">
                  <Zap className="w-6 h-6 text-accent" />
                  <h2 className="text-2xl font-heading font-bold text-text-primary">Prediction Status</h2>
                </div>
                
                {predictionResult && (
                  <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                    <div className="bg-secondary/50 p-4 rounded-lg border border-accent/20">
                      <p className="text-text-secondary font-label text-sm mb-1">UniProt ID</p>
                      <p className="text-text-primary font-code text-lg">{predictionResult.uniprotId}</p>
                    </div>
                    <div className="bg-secondary/50 p-4 rounded-lg border border-accent/20">
                      <p className="text-text-secondary font-label text-sm mb-1">Status</p>
                      <div className="flex items-center gap-2">
                        <StatusIcon status={jobStatus?.status || predictionResult.status} />
                        <span className="text-text-primary font-label font-semibold">
                          {jobStatus?.status || predictionResult.status}
                        </span>
                      </div>
                    </div>
                    <div className="bg-secondary/50 p-4 rounded-lg border border-accent/20">
                      <p className="text-text-secondary font-label text-sm mb-1">Submitted</p>
                      <p className="text-text-primary font-code text-sm">{new Date(predictionResult.createdAt).toLocaleString()}</p>
                    </div>
                  </div>
                )}

                {/* Structure Viewer */}
                {jobStatus?.pdbUrl && (
                  <div className="mt-6">
                    <div className="flex flex-wrap gap-2 mb-4">
                      {[
                        { name: 'secondary', label: 'Secondary Structure', icon: Eye },
                        { name: 'surface', label: 'Surface View', icon: EyeOff },
                        { name: 'activeSite', label: 'Active Sites', icon: Zap },
                      ].map(({ name, label, icon: Icon }) => (
                        <button
                          key={name}
                          onClick={() => toggleRepresentation(name)}
                          className="flex items-center gap-2 px-3 py-2 bg-accent/10 text-accent border border-accent/30 rounded-lg hover:bg-accent/20 transition-colors duration-200 font-label text-sm"
                        >
                          <Icon className="w-4 h-4" />
                          {label}
                        </button>
                      ))}
                      <button
                        onClick={() => stageRef.current?.autoView()}
                        className="flex items-center gap-2 px-3 py-2 bg-accent-secondary/10 text-accent-secondary border border-accent-secondary/30 rounded-lg hover:bg-accent-secondary/20 transition-colors duration-200 font-label text-sm"
                      >
                        <RotateCcw className="w-4 h-4" />
                        Reset View
                      </button>
                    </div>
                    
                    <div className="relative">
                      <div
                        id="viewport"
                        className="w-full h-96 md:h-[500px] bg-primary rounded-xl border border-accent/20 relative overflow-hidden"
                      >
                        {/* {viewportLoading ? (
                          <div className="absolute inset-0 flex items-center justify-center bg-primary/90">
                            <div className="text-center">
                              <Loader2 className="w-12 h-12 text-accent animate-spin mx-auto mb-4" />
                              <p className="text-text-primary font-label text-lg">Loading 3D Structure...</p>
                              <p className="text-text-secondary font-label text-sm">Rendering molecular visualization</p>
                            </div>
                          </div>
                        ) : structureRendered ? (
                          <div className="absolute inset-0 flex items-center justify-center">
                            <div className="text-center">
                              <div className="w-32 h-32 bg-gradient-to-br from-accent/20 to-accent-secondary/20 rounded-full flex items-center justify-center mx-auto mb-4 animate-pulse-glow">
                                <Atom className="w-16 h-16 text-accent" />
                              </div>
                              <p className="text-text-primary font-label text-lg">3D Structure Rendered</p>
                              <p className="text-text-secondary font-label text-sm">Interactive molecular view ready</p>
                            </div>
                          </div>
                        ) : null} */}
                      </div>
                    </div>
                  </div>
                )}

                {jobStatus?.error && (
                  <div className="mt-4 p-3 bg-error/10 border border-error/30 rounded-lg">
                    <p className="text-error font-label">Error: {jobStatus.error}</p>
                  </div>
                )}
              </div>
            )}

            {/* UniProt Results */}
            {activeTab === 'uniprot' && summary && (
              <div className="bg-gradient-to-br from-accent-secondary/5 to-accent/5 rounded-xl p-6 border border-accent-secondary/20 animate-slide-up">
                <div className="flex items-center gap-3 mb-6">
                  <Database className="w-6 h-6 text-accent-secondary" />
                  <h2 className="text-2xl font-heading font-bold text-text-primary">Protein Information</h2>
                </div>

                {/* Basic Info Grid */}
                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 mb-8">
                  {[
                    { label: 'Protein Name', value: summary.protein?.recommendedName?.fullName?.value || 'Not available' },
                    { label: 'Gene', value: summary.gene?.[0]?.name?.value || 'Not available' },
                    { label: 'Organism', value: summary.organism?.scientificName || 'Not available' },
                    { label: 'Length', value: summary.sequence?.length ? `${summary.sequence.length} amino acids` : 'Not available' }
                  ].map(({ label, value }) => (
                    <div key={label} className="bg-secondary/50 p-4 rounded-lg border border-accent-secondary/20">
                      <p className="text-text-secondary font-label text-sm mb-1">{label}</p>
                      <p className="text-text-primary font-body text-sm font-medium break-words">{value}</p>
                    </div>
                  ))}
                </div>

                {/* Functional Analysis */}
                <div className="bg-gradient-to-br from-secondary/30 to-secondary/50 rounded-xl p-6 border border-accent-secondary/20 mb-6">
                  <h3 className="text-xl font-heading font-bold text-text-primary mb-6 flex items-center gap-2">
                    <Zap className="w-5 h-5 text-accent-secondary" />
                    Functional Analysis
                  </h3>
                  
                  <div className="space-y-6">
                    {/* Function */}
                    <div>
                      <h4 className="font-label font-semibold text-text-primary mb-3 flex items-center gap-2">
                        <CheckCircle className="w-4 h-4 text-success" />
                        Protein Function
                      </h4>
                      <div className="bg-primary/50 p-4 rounded-lg border border-accent-secondary/20">
                        <p className="text-text-primary font-body text-sm leading-relaxed">
                          {summary.comments?.find(c => c.type === 'FUNCTION')?.text?.[0]?.value || 'Function not specified in database'}
                        </p>
                      </div>
                    </div>

                    {/* Features */}
                    {annotations?.features?.length > 0 && (
                      <div>
                        <h4 className="font-label font-semibold text-text-primary mb-3 flex items-center gap-2">
                          <Atom className="w-4 h-4 text-accent" />
                          Structural Features
                        </h4>
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                          {annotations.features.slice(0, 6).map((feat, i) => (
                            <div key={i} className="bg-primary/50 p-3 rounded-lg border border-accent-secondary/20 flex items-center gap-3">
                              <span className="inline-block bg-accent/20 text-accent text-xs px-2 py-1 rounded-full font-label font-semibold">
                                {feat.type}
                              </span>
                              <span className="text-text-primary font-body text-sm flex-1">{feat.description || 'No description'}</span>
                            </div>
                          ))}
                        </div>
                      </div>
                    )}

                    {/* Properties Grid */}
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      <div>
                        <h4 className="font-label font-semibold text-text-primary mb-2">Subcellular Location</h4>
                        <div className="bg-primary/50 p-3 rounded-lg border border-accent-secondary/20">
                          <p className="text-text-primary font-body text-sm">
                            {summary.comments?.find(c => c.type === 'SUBCELLULAR_LOCATION')?.text?.[0]?.value || 'Location unknown'}
                          </p>
                        </div>
                      </div>
                      <div>
                        <h4 className="font-label font-semibold text-text-primary mb-2">Catalytic Activity</h4>
                        <div className="bg-primary/50 p-3 rounded-lg border border-accent-secondary/20">
                          <p className="text-text-primary font-body text-sm">
                            {summary.comments?.find(c => c.type === 'CATALYTIC_ACTIVITY')?.text?.[0]?.value || 'Non-enzymatic protein'}
                          </p>
                        </div>
                      </div>
                    </div>
                  </div>
                </div>

                {/* Sequence Display */}
                {summary?.sequence?.sequence && (
                  <div>
                    <h3 className="text-xl font-heading font-bold text-text-primary mb-4 flex items-center gap-2">
                      <Database className="w-5 h-5 text-accent-secondary" />
                      Amino Acid Sequence
                    </h3>
                    <div className="bg-primary/50 p-4 rounded-lg border border-accent-secondary/20 overflow-x-auto">
                      <pre className="font-code text-xs text-text-primary whitespace-pre-wrap break-all leading-relaxed">
                        {summary.sequence.sequence.match(/.{1,10}/g)?.map((chunk, i) => (
                          <div key={i} className="flex">
                            <span className="text-text-secondary mr-4 w-12 text-right">
                              {(i * 10 + 1).toString().padStart(4, ' ')}
                            </span>
                            <span className="text-accent font-medium tracking-wider">
                              {chunk.split('').join(' ')}
                            </span>
                          </div>
                        )) || 'No sequence data available'}
                      </pre>
                    </div>
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default AlphaFoldExplorer;