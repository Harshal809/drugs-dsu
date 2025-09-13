import React, { useState, useEffect, useRef } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";
import jsPDF from "jspdf";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const PY_API_BASE_URL = import.meta.env.VITE_PY_API_BASE_URL || "http://10.1.40.155:5001";

const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
});

const pythonAxios = axios.create({
  baseURL: PY_API_BASE_URL,
  withCredentials: false,
  timeout: 30000, // 30 second timeout
});

// Utility functions
const cleanAbstract = (abstract) => {
  if (!abstract || abstract === "Abstract not available.") return abstract;
  return abstract
    .replace(/<[^>]+>/g, "")
    .replace(/\s+/g, " ")
    .trim();
};

const isValidPaper = (paper) => {
  return (
    paper &&
    typeof paper === "object" &&
    paper.title &&
    paper.authors &&
    paper.abstract &&
    paper.introduction &&
    paper.methodology &&
    paper.resultsAndDiscussion &&
    paper.conclusion &&
    Array.isArray(paper.keywords) &&
    Array.isArray(paper.references)
  );
};

// Clean Gemini response to extract valid JSON
const cleanGeminiResponse = (response) => {
  let cleaned = response.replace(/```json|```/g, "").trim();
  cleaned = cleaned.replace(/^```+|```$/g, "").trim();
  const jsonMatch = cleaned.match(/({.*?})/s);
  return jsonMatch ? jsonMatch[0] : cleaned;
};

// Enhanced parsing function for LLM output with markdown-like formatting
const parseCitationsText = (text) => {
  if (!text || typeof text !== "string") return [];
  
  // Split by --- separated papers
  const paperBlocks = text.split(/\n?-{3,}\n?/g)
    .map(block => block.trim())
    .filter(Boolean);
  
  const papers = [];
  
  for (const block of paperBlocks) {
    if (!block.trim()) continue;
    
    const paper = {
      title: "",
      authors: "",
      year: "",
      abstract: "",
      doi: "",
      url: "",
      rawContent: block.trim(),
      number: "",
      relatedTo: ""
    };
    
    // Extract header **Paper N: Related to ...**
    const headerMatch = block.match(/^\*\*Paper (\d+): Related to (.+?)\*\*/);
    if (headerMatch) {
      paper.number = headerMatch[1];
      paper.relatedTo = headerMatch[2].trim();
    }
    
    // Extract title: after "Title\n"
    const titleMatch = block.match(/Title\s*\n(.+?)(?=\nAuthors:|\n-{3,}|$)/s);
    if (titleMatch) {
      paper.title = titleMatch[1].trim();
    }
    
    // Extract authors
    const authorsMatch = block.match(/Authors:\s*(.+?)(?=\nPublished:|\nAbstract:|\n-{3,}|$)/s);
    if (authorsMatch) {
      paper.authors = authorsMatch[1].trim();
    }
    
    // Extract publication year
    const yearMatch = block.match(/Published:\s*(.+?)(?=\nAbstract:|\nLink|\n-{3,}|$)/s);
    if (yearMatch) {
      paper.year = yearMatch[1].trim();
    }
    
    // Extract abstract
    const abstractMatch = block.match(/Abstract:\s*(.+?)(?=\nLink|\n-{3,}|$)/s);
    if (abstractMatch) {
      paper.abstract = cleanAbstract(abstractMatch[1].trim());
    }
    
    // Extract URL/Link
    const urlMatch = block.match(/Link of the research paper:\s*(https?:\/\/[^\s\n]+)/i);
    if (urlMatch) {
      paper.url = urlMatch[1].trim();
    } else {
      // Fallback: Look for any URL in the text
      const fallbackUrlMatch = block.match(/(https?:\/\/[^\s\n]+)/);
      if (fallbackUrlMatch) {
        paper.url = fallbackUrlMatch[1].trim();
      } else {
        paper.url = "No URL available";
      }
    }
    
    // Only add paper if it has at least a title
    if (paper.title) {
      papers.push(paper);
    }
  }
  
  return papers;
};

// Component to render formatted text with markdown-like styling
const FormattedText = ({ text, className = "" }) => {
  if (!text) return null;
  
  const processText = (inputText) => {
    const parts = inputText.split(/(\*\*[^*]+\*\*|https?:\/\/[^\s]+|\n|•|\*\s+)/g);
    
    return parts.map((part, index) => {
      if (!part) return null;
      
      // Handle bold text
      if (part.match(/^\*\*(.+)\*\*$/)) {
        const boldText = part.replace(/^\*\*|\*\*$/g, '');
        return (
          <strong key={`bold-${index}`} className="font-bold text-text-primary">
            {boldText}
          </strong>
        );
      }
      
      // Handle URLs
      if (part.match(/^https?:\/\//)) {
        return (
          <a
            key={`link-${index}`}
            href={part}
            target="_blank"
            rel="noopener noreferrer"
            className="text-accent hover:text-accent-secondary underline transition-colors duration-300"
          >
            {part}
          </a>
        );
      }
      
      // Handle line breaks
      if (part === '\n') {
        return <br key={`br-${index}`} />;
      }
      
      // Handle bullet points
      if (part === '•' || part.match(/^\*\s+/)) {
        return (
          <span key={`bullet-${index}`} className="text-accent font-bold">
            • 
          </span>
        );
      }
      
      // Regular text
      return <span key={`text-${index}`}>{part}</span>;
    });
  };
  
  return (
    <div className={`formatted-text ${className}`}>
      {processText(text)}
    </div>
  );
};

// Component to render a research paper card
const ResearchPaperCard = ({ paper, index }) => {
  return (
    <div
      className="border-l-4 border-accent pl-4 transform transition-all duration-500 animate-slide-up bg-white/5 p-4 rounded-lg"
      style={{ animationDelay: `${index * 100}ms` }}
    >
      {/* Title */}
      {paper.title && (
        <h5 className="text-lg sm:text-xl font-heading font-bold text-text-primary mb-3">
          <FormattedText text={paper.title} />
        </h5>
      )}
      
      {/* Authors */}
      {paper.authors && (
        <div className="mb-2">
          <span className="font-bold text-text-primary">Authors: </span>
          <FormattedText text={paper.authors} className="text-text-secondary" />
        </div>
      )}
      
      {/* Publication Year */}
      {paper.year && (
        <div className="mb-2">
          <span className="font-bold text-text-primary">Published: </span>
          <FormattedText text={paper.year} className="text-text-secondary" />
        </div>
      )}
      
      {/* Abstract */}
      {paper.abstract && (
        <div className="mb-3">
          <span className="font-bold text-text-primary">Abstract: </span>
          <FormattedText text={paper.abstract} className="text-text-secondary leading-relaxed" />
        </div>
      )}
      
      {/* DOI */}
      {paper.doi && (
        <div className="mb-2">
          <span className="font-bold text-text-primary">DOI: </span>
          <FormattedText text={paper.doi} className="text-text-secondary" />
        </div>
      )}
      
      {/* URL/Link */}
      {paper.url && paper.url !== "No URL available" && (
        <div className="mb-2">
          <span className="font-bold text-text-primary">Link: </span>
          <FormattedText text={paper.url} />
        </div>
      )}
    </div>
  );
};

function Airesearchgenerator() {
  const [activeTab, setActiveTab] = useState("related");
  const [symptomGroups, setSymptomGroups] = useState([]);
  const [productSmilesGroups, setProductSmilesGroups] = useState([]);
  const [selectedSymptomGroupIndex, setSelectedSymptomGroupIndex] = useState("");
  const [selectedSmiles, setSelectedSmiles] = useState("");
  const [researchPapers, setResearchPapers] = useState([]);
  const [savedPapers, setSavedPapers] = useState([]);
  const [savedGeneratedPapers, setSavedGeneratedPapers] = useState([]);
  const [researchSummary, setResearchSummary] = useState("");
  const [generatedPaper, setGeneratedPaper] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const toastShown = useRef(false);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    if (!toastShown.current) {
      toast(
        "First, select a symptom group and SMILES string",
        {
          position: "top-right",
          duration: 8000,
          style: {
            background: "#fefcbf",
            color: "#92400e",
            border: "1px solid #f59e0b",
          },
          icon: "⚠️",
        }
      );
      toastShown.current = true;
    }
  }, []);

  useEffect(() => {
    const initializeApp = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError("Authentication failed. Please log in.");
        return;
      }
      await fetchSymptomsAndProducts();
      await fetchSavedPapers();
      await fetchSavedGeneratedPapers();
    };
    initializeApp();
  }, [checkAuth]);

  const fetchSymptomsAndProducts = async () => {
    if (!user?._id) return;
    setLoading(true);
    try {
      const response = await axiosInstance.get(`/getdata/getsymptoms-product/${user._id}`);
      const { symptoms, productSmiles } = response.data;
      setSymptomGroups(symptoms || []);
      setProductSmilesGroups(productSmiles || []);
      if (symptoms?.length > 0) {
        setSelectedSymptomGroupIndex("0");
        if (productSmiles?.[0]?.length > 0) {
          setSelectedSmiles(productSmiles[0][0]);
          localStorage.setItem(`symptomGroups_${user._id}`, JSON.stringify(symptoms));
        }
      }
    } catch (err) {
      console.error("Error fetching symptoms and products:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch symptoms and products");
    } finally {
      setLoading(false);
    }
  };

  const fetchSavedPapers = async () => {
    if (!user?._id) return;
    try {
      const response = await axiosInstance.get("/researchPaper/saved-research-papers");
      setSavedPapers(response.data.papers || []);
    } catch (err) {
      console.error("Error fetching saved papers:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch saved papers");
    }
  };

  const fetchSavedGeneratedPapers = async () => {
    if (!user?._id) return;
    try {
      const response = await axiosInstance.get("/researchPaper/saved-generated-research-papers");
      const papers = response.data.papers || [];
      const sortedPapers = papers.sort((a, b) => {
        const dateA = a.createdAt ? new Date(a.createdAt) : new Date(0);
        const dateB = b.createdAt ? new Date(b.createdAt) : new Date(0);
        return dateB - dateA;
      });
      setSavedGeneratedPapers(sortedPapers);
    } catch (err) {
      console.error("Error fetching saved generated papers:", err.response?.data || err.message);
      setError(err.response?.data?.message || "Failed to fetch saved generated papers");
    }
  };

  const checkIfPapersExist = async (symptoms, smiles) => {
    if (!user?._id || !symptoms || !smiles) return false;
    try {
      const response = await axiosInstance.get("/researchPaper/check-saved-papers", {
        params: { symptoms, smiles },
      });
      return response.data.exists;
    } catch (err) {
      console.error("Error checking saved papers:", err.response?.data || err.message);
      return false;
    }
  };

  const checkIfGeneratedPaperExists = async (symptoms, smiles) => {
    if (!user?._id || !symptoms || !smiles) return false;
    try {
      const response = await axiosInstance.get("/researchPaper/check-saved-generated-papers", {
        params: { symptoms, smiles },
      });
      console.log("Check generated paper response:", response.data);
      return response.data.exists;
    } catch (err) {
      console.error("Error checking saved generated papers:", err.response?.data || err.message);
      return false;
    }
  };

  const savePapers = async (symptoms, smiles, papers) => {
    if (!user?._id || !papers.length) return;
    const payload = {
      userId: user._id,
      molecule: { symptoms, smiles },
      papers,
    };
    try {
      await axiosInstance.post("/researchPaper/save-research-papers", payload);
      await fetchSavedPapers();
      toast.success("Research papers saved successfully!");
    } catch (err) {
      console.error("Error saving papers:", err.response?.data || err.message);
      toast.error("Failed to save research papers");
    }
  };

  const saveGeneratedPaper = async (symptoms, smiles, paper) => {
    if (!user?._id || !paper) return;
    const payload = {
      userId: user._id,
      molecule: { symptoms, smiles },
      paper,
    };
    try {
      await axiosInstance.post("/researchPaper/save-generated-research-paper", payload);
      await fetchSavedGeneratedPapers();
      toast.success("Generated research paper saved successfully!");
    } catch (err) {
      console.error("Error saving generated paper:", err.response?.data || err.message);
      toast.error("Failed to save generated research paper");
    }
  };

const fetchResearchPapers = async (symptoms, smiles) => {
  try {
    console.log("Making request to Python API:", {
      url: `${pythonAxios.defaults.baseURL}/api/research-papers`,
      data: { symptom_group: symptoms, smiles_string: smiles }
    });

    const response = await pythonAxios.post("/api/research-papers", {
      symptom_group: symptoms,
      smiles_string: smiles,
    });

    console.log("Python API Response:", response.data);
    
    const data = response.data;
    
    // Use parsed_papers if available (from new backend)
    if (data.parsed_papers && Array.isArray(data.parsed_papers) && data.parsed_papers.length > 0) {
      return {
        papers: data.parsed_papers,
        message: data.message
      };
    }
    // Fallback to old parsing method
    else if (data.papers && typeof data.papers === 'string' && data.papers.length > 0) {
      const parsedPapers = parseCitationsText(data.papers);
      return {
        papers: parsedPapers,
        message: data.message
      };
    } else if (Array.isArray(data.papers)) {
      return {
        papers: data.papers,
        message: data.message
      };
    } else {
      return {
        papers: [],
        message: data.message || "No research papers found"
      };
    }
  } catch (err) {
    console.error("Error calling Python API:", err);
    
    if (err.response) {
      const errorMsg = err.response.data?.message || err.response.statusText || 'Server error';
      throw new Error(`API Error (${err.response.status}): ${errorMsg}`);
    } else if (err.request) {
      throw new Error("No response from server. Please check if the Python backend is running on http://10.1.40.155:5001");
    } else {
      throw new Error(`Request setup error: ${err.message}`);
    }
  }
};

  const generateResearchPaper = async (symptoms, smiles) => {
    try {
      const prompt = `You are an expert in chemical informatics and academic writing. Given the SMILES string "${smiles}" and associated symptoms "${symptoms}", generate a high-quality research paper in IEEE format. The paper should be structured with the following sections: Title, Authors, Abstract, Keywords, I. Introduction, II. Methodology, III. Results and Discussion, IV. Conclusion, References. Ensure the content is scientifically accurate, relevant to the compound's therapeutic applications for the symptoms, and includes realistic data and references. Return the paper in clean JSON format (no Markdown or code fences) with fields: title (string), authors (string), abstract (string), keywords (array of strings), introduction (string), methodology (string), resultsAndDiscussion (string), conclusion (string), references (array of strings). Example output: { "title": "Therapeutic Applications of Compound X", "authors": "Test", "abstract": "This paper investigates...", "keywords": ["compound X", "therapeutics", "symptoms"], "introduction": "Introduction text...", "methodology": "Methodology text...", "resultsAndDiscussion": "Results text...", "conclusion": "Conclusion text...", "references": ["Ref 1", "Ref 2"] }`;
      const response = await axiosInstance.post("/researchPaper/proxy/gemini", { prompt });
      const content = response.data.content;
      try {
        const cleanedContent = cleanGeminiResponse(content);
        const paper = JSON.parse(cleanedContent);
        if (!isValidPaper(paper)) throw new Error("Invalid paper format");
        return paper;
      } catch (parseError) {
        console.error("Error parsing Gemini response:", parseError);
        return null;
      }
    } catch (err) {
      console.error("Error generating research paper:", err.response?.data || err.message);
      return null;
    }
  };

  const handleResearchClick = async () => {
    if (selectedSymptomGroupIndex === "" || !selectedSmiles) {
      toast.error("Please select both a symptom group and SMILES string");
      return;
    }

    const selectedSymptoms = symptomGroups[selectedSymptomGroupIndex]?.join(", ") || "";

    const papersExist = await checkIfPapersExist(selectedSymptoms, selectedSmiles);
    if (papersExist) {
      toast("Research papers already saved. Redirecting to Saved Research Papers.", {
        type: "info",
      });
      setActiveTab("saved");
      await fetchSavedPapers();
      return;
    }

    setLoading(true);
    setError(null);
    setResearchPapers([]);
    setResearchSummary("");

    try {
      const result = await fetchResearchPapers(selectedSymptoms, selectedSmiles);

      if (result.papers.length === 0) {
        setError(
          result.message ||
          "No research papers found for the selected combination. Please try a different selection."
        );
        return;
      }

      setResearchPapers(result.papers);
      setResearchSummary(
        `Found ${result.papers.length} research papers related to the compound with SMILES "${selectedSmiles}" for treating symptoms: ${selectedSymptoms}.`
      );

      await savePapers(selectedSymptoms, selectedSmiles, result.papers);
    } catch (err) {
      setError(err.message || "Failed to fetch research papers from AI service. Please try again.");
      toast.error(err.message || "Failed to fetch research papers");
    } finally {
      setLoading(false);
    }
  };

  const handleGeneratePaperClick = async () => {
    if (selectedSymptomGroupIndex === "" || !selectedSmiles) {
      toast.error("Please select both a symptom group and SMILES string");
      return;
    }
    const selectedSymptoms = symptomGroups[selectedSymptomGroupIndex]?.join(", ") || "";
    setLoading(true);
    setError(null);
    setGeneratedPaper(null);

    try {
      const paperExists = await checkIfGeneratedPaperExists(selectedSymptoms, selectedSmiles);
      if (paperExists) {
        const confirmRedirect = window.confirm(
          "A generated research paper for this SMILES and symptoms already exists. Do you want to view it in Saved Generated Papers?"
        );
        if (confirmRedirect) {
          toast("Redirecting to Saved Generated Papers.", { type: "info" });
          setActiveTab("savedGenerated");
          await fetchSavedGeneratedPapers();
        }
        return;
      }

      const paper = await generateResearchPaper(selectedSymptoms, selectedSmiles);
      if (paper) {
        setGeneratedPaper(paper);
        await saveGeneratedPaper(selectedSymptoms, selectedSmiles, paper);
        toast.success("Research paper generated and saved successfully!");
      } else {
        setError("Failed to generate research paper. Please try again.");
        toast.error("Failed to generate research paper");
      }
    } catch (err) {
      console.error("Error in generating research paper:", err.response?.data || err.message);
      setError("Failed to generate research paper. Please try again.");
      toast.error("Failed to generate research paper");
    } finally {
      setLoading(false);
    }
  };

  const exportToPDF = (paper) => {
    if (!paper) {
      toast.error("No paper to download");
      return;
    }
    const pdf = new jsPDF("p", "mm", "a4");
    const pageWidth = pdf.internal.pageSize.getWidth();
    const pageHeight = pdf.internal.pageSize.getHeight();
    const margin = 15;
    let yPosition = margin;

    const checkPageBreak = (requiredHeight) => {
      if (yPosition + requiredHeight > pageHeight - margin) {
        pdf.addPage();
        yPosition = margin;
      }
    };

    const addTextWithPagination = (text, fontSize, x, y, maxWidth) => {
      pdf.setFontSize(fontSize);
      pdf.setFont("times", "normal");
      const lines = pdf.splitTextToSize(text, maxWidth);
      const lineHeight = fontSize * 0.4;
      lines.forEach((line) => {
        checkPageBreak(lineHeight);
        pdf.text(line, x, yPosition);
        yPosition += lineHeight;
      });
      return yPosition;
    };

    pdf.setFontSize(16);
    pdf.setFont("times", "bold");
    const titleLines = pdf.splitTextToSize(paper.title, pageWidth - 2 * margin);
    const titleLineHeight = 16 * 0.4;
    titleLines.forEach((line) => {
      checkPageBreak(titleLineHeight);
      const lineWidth = pdf.getTextWidth(line);
      pdf.text(line, (pageWidth - lineWidth) / 2, yPosition);
      yPosition += titleLineHeight;
    });
    yPosition += 10;

    pdf.setFontSize(12);
    pdf.setFont("times", "normal");
    const authorsWidth = pdf.getTextWidth(paper.authors);
    checkPageBreak(8);
    pdf.text(paper.authors, (pageWidth - authorsWidth) / 2, yPosition);
    yPosition += 15;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("Abstract", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.abstract, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 10;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("Keywords", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.keywords.join(", "), 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("I. Introduction", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.introduction, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("II. Methodology", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.methodology, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("III. Results and Discussion", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.resultsAndDiscussion, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("IV. Conclusion", margin, yPosition);
    yPosition += 6;
    yPosition = addTextWithPagination(paper.conclusion, 10, margin, yPosition, pageWidth - 2 * margin);
    yPosition += 15;

    pdf.setFontSize(12);
    pdf.setFont("times", "bold");
    checkPageBreak(8);
    pdf.text("References", margin, yPosition);
    yPosition += 6;
    pdf.setFontSize(10);
    pdf.setFont("times", "normal");
    paper.references.forEach((ref) => {
      checkPageBreak(6);
      const lines = pdf.splitTextToSize(ref, pageWidth - 2 * margin);
      lines.forEach((line) => {
        checkPageBreak(5);
        pdf.text(line, margin, yPosition);
        yPosition += 5;
      });
      yPosition += 2;
    });

    pdf.save(`${paper.title.replace(/\s+/g, "_")}.pdf`);
    toast.success("PDF downloaded successfully!");
  };

  const handleTabChange = (tab) => {
    setActiveTab(tab);
    setResearchPapers([]);
    setResearchSummary("");
    setGeneratedPaper(null);
    setError(null);
  };

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-primary text-text-primary animate-pulse">
        <p className="text-lg font-body">Verifying authentication...</p>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-primary">
        <div className="text-center p-6 bg-secondary rounded-lg shadow-lg max-w-md w-full transform transition-all duration-500 ease-out animate-slide-up">
          <p className="text-text-secondary font-body mb-4">Please log in to access AI Research Generator</p>
          <button
            className="px-4 py-2 bg-accent text-primary rounded-lg hover:bg-accent-secondary transition-all duration-300 transform hover:scale-105"
            onClick={() => (window.location.href = "/login")}
          >
            Go to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen py-6 px-4 sm:px-6 lg:px-8 bg-primary">
      <div className="max-w-7xl mx-auto">
        <h1 className="text-3xl sm:text-4xl font-heading font-bold text-accent mb-6 text-center transform transition-all duration-500 ease-out animate-slide-down">
          AI Research Generator
          <p className="text-sm text-accent-secondary font-mono">(POWERED BY GEMINI)</p>
        </h1>

        <div className="flex flex-col sm:flex-row justify-center mb-6 space-y-2 sm:space-y-0 sm:space-x-4">
          {["related", "saved", "generate", "savedGenerated"].map((tab) => (
            <button
              key={tab}
              className={`px-4 sm:px-6 py-2 rounded-lg font-label text-sm sm:text-base transition-all duration-300 transform hover:scale-105 ${
                activeTab === tab
                  ? "bg-accent text-primary shadow-lg"
                  : "bg-secondary text-text-secondary hover:bg-accent-secondary"
              }`}
              onClick={() => handleTabChange(tab)}
            >
              {tab === "related" && "Related Research Papers"}
              {tab === "saved" && "Saved Research Papers"}
              {tab === "generate" && "Generate Research Paper"}
              {tab === "savedGenerated" && "Saved Generated Papers"}
            </button>
          ))}
        </div>

        <div className="bg-secondary p-4 sm:p-6 rounded-xl shadow-lg border border-accent-secondary transform transition-all duration-500 ease-out animate-fade-in">
          {activeTab === "related" && (
            <>
              <h2 className="text-xl sm:text-2xl font-heading font-semibold text-accent mb-4 sm:mb-6">
                Related Research Papers
              </h2>

              <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 sm:gap-6 mb-4 sm:mb-6">
                <div className="relative">
                  <label className="block text-sm sm:text-base font-label text-text-primary mb-2">
                    Select Symptom Group
                  </label>
                  <select
                    value={selectedSymptomGroupIndex}
                    onChange={(e) => {
                      setSelectedSymptomGroupIndex(e.target.value);
                      setSelectedSmiles(productSmilesGroups[e.target.value]?.[0] || "");
                    }}
                    className="w-full p-2 sm:p-3 border border-accent-secondary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent transition-all duration-300 text-text-primary bg-primary font-body"
                    disabled={loading || symptomGroups.length === 0}
                  >
                    {symptomGroups.length === 0 ? (
                      <option value="">No symptom groups available</option>
                    ) : (
                      symptomGroups.map((group, index) => (
                        <option key={index} value={index}>
                          {group.join(", ")}
                        </option>
                      ))
                    )}
                  </select>
                </div>
                <div className="relative">
                  <label className="block text-sm sm:text-base font-label text-text-primary mb-2">
                    Select SMILES String
                  </label>
                  <select
                    value={selectedSmiles}
                    onChange={(e) => setSelectedSmiles(e.target.value)}
                    className="w-full p-2 sm:p-3 border border-accent-secondary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent transition-all duration-300 text-text-primary bg-primary font-body"
                    disabled={loading || !selectedSymptomGroupIndex || productSmilesGroups[selectedSymptomGroupIndex]?.length === 0}
                  >
                    {productSmilesGroups[selectedSymptomGroupIndex]?.length === 0 ? (
                      <option value="">No SMILES available</option>
                    ) : (
                      productSmilesGroups[selectedSymptomGroupIndex]?.map((smiles, index) => (
                        <option key={index} value={smiles}>
                          {smiles}
                        </option>
                      ))
                    )}
                  </select>
                </div>
              </div>

              <button
                onClick={handleResearchClick}
                disabled={loading || selectedSymptomGroupIndex === "" || !selectedSmiles}
                className="w-full py-2 sm:py-3 px-4 bg-accent text-primary rounded-lg hover:bg-accent-secondary disabled:bg-gray-400 disabled:cursor-not-allowed transition-all duration-300 transform hover:scale-95 relative overflow-hidden"
              >
                <span className="relative z-10">
                  {loading ? (
                    <span className="flex items-center justify-center">
                      <svg
                        className="animate-spin h-5 w-5 mr-2 text-primary"
                        viewBox="0 0 24 24"
                      >
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                        <path
                          className="opacity-75"
                          fill="currentColor"
                          d="M4 12a8 8 0 018-8V0l4 4-4 4V4a4 4 0 00-4 4h-4z"
                        />
                      </svg>
                      Fetching Research...
                    </span>
                  ) : (
                    "Fetch Related Research"
                  )}
                </span>
                <span className="absolute inset-0 bg-accent-secondary opacity-0 hover:opacity-20 transition-opacity duration-300" />
              </button>

              {error && (
                <div className="mt-4 sm:mt-6 bg-error bg-opacity-10 border border-error text-error px-4 py-3 rounded-lg flex justify-between items-center transform transition-all duration-500 animate-slide-up">
                  <p className="text-sm sm:text-base font-body">{error}</p>
                  <button
                    className="text-error underline hover:text-error/80 font-label text-sm sm:text-base"
                    onClick={() => setError(null)}
                  >
                    Dismiss
                  </button>
                </div>
              )}

              {(researchSummary || researchPapers.length > 0) && (
                <div className="mt-6 sm:mt-8 bg-primary p-4 sm:p-6 rounded-xl border border-accent-secondary transform transition-all duration-500 animate-slide-up">
                  <h3 className="text-lg sm:text-xl font-heading font-semibold text-accent mb-4">
                    Related Research Information
                  </h3>

                  {researchPapers.length > 0 && (
                    <div className="mb-6 sm:mb-8">
                      <h4 className="text-md sm:text-lg font-heading font-semibold text-text-primary mb-4">
                        Research Papers ({researchPapers.length} found)
                      </h4>
                      <div className="overflow-x-auto">
                        <table className="min-w-full bg-primary border border-accent-secondary rounded-lg overflow-hidden">
                          <thead>
                            <tr className="bg-accent text-primary text-left">
                              <th className="px-4 py-3 font-label">#</th>
                              <th className="px-4 py-3 font-label">Title</th>
                              <th className="px-4 py-3 font-label">Authors</th>
                              <th className="px-4 py-3 font-label">Year</th>
                              <th className="px-4 py-3 font-label">Link</th>
                              <th className="px-4 py-3 font-label">Details</th>
                            </tr>
                          </thead>
                          <tbody>
                            {researchPapers.map((paper, index) => (
                              <tr key={index} className="border-t border-accent-secondary hover:bg-secondary/50 transition-colors duration-200">
                                <td className="px-4 py-3 text-text-primary font-body">
                                  {paper.number || (index + 1)}
                                </td>
                                <td className="px-4 py-3 text-text-primary font-body max-w-xs">
                                  <div className="truncate" title={paper.title}>
                                    {paper.title}
                                  </div>
                                </td>
                                <td className="px-4 py-3 text-text-secondary font-body max-w-xs">
                                  <div className="truncate" title={paper.authors}>
                                    {paper.authors}
                                  </div>
                                </td>
                                <td className="px-4 py-3 text-text-secondary font-body">
                                  {paper.year}
                                </td>
                                <td className="px-4 py-3">
                                  {paper.url && paper.url !== "No URL available" ? (
                                    <a
                                      href={paper.url}
                                      target="_blank"
                                      rel="noopener noreferrer"
                                      className="text-accent hover:text-accent-secondary underline transition-colors duration-300 font-label text-sm"
                                    >
                                      View Paper
                                    </a>
                                  ) : (
                                    <span className="text-text-secondary font-body text-sm">No Link</span>
                                  )}
                                </td>
                                <td className="px-4 py-3">
                                  <details className="group cursor-pointer">
                                    <summary className="flex items-center gap-2 text-accent hover:text-accent-secondary font-label text-sm transition-all duration-300 hover:scale-105">
                                      <svg 
                                        className="w-4 h-4 transition-transform duration-300 group-open:rotate-180" 
                                        fill="none" 
                                        stroke="currentColor" 
                                        viewBox="0 0 24 24"
                                      >
                                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                                      </svg>
                                      <span className="font-semibold">View Details</span>
                                    </summary>
                                    <div className="mt-4 p-6 bg-gradient-to-br from-primary to-secondary rounded-xl border border-accent-secondary shadow-lg backdrop-blur-sm max-w-2xl transform transition-all duration-500 ease-out animate-slide-up">
                                      <div className="space-y-4">
                                        {paper.abstract && (
                                          <div className="bg-white/5 p-4 rounded-lg border border-accent-secondary/30">
                                            <div className="flex items-center gap-2 mb-3">
                                              <svg className="w-5 h-5 text-accent" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                                              </svg>
                                              <h4 className="text-accent font-heading font-semibold text-lg">Abstract</h4>
                                            </div>
                                            <p className="text-text-secondary leading-relaxed text-sm bg-white/5 p-3 rounded-lg border-l-4 border-accent">
                                              {paper.abstract}
                                            </p>
                                          </div>
                                        )}
                                        
                                        {paper.rawContent && (
                                          <div className="bg-white/5 p-4 rounded-lg border border-accent-secondary/30">
                                            <div className="flex items-center gap-2 mb-3">
                                              <svg className="w-5 h-5 text-accent" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
                                              </svg>
                                              <h4 className="text-accent font-heading font-semibold text-lg">Full Research Details</h4>
                                            </div>
                                            <div className="text-sm text-text-secondary max-h-60 overflow-y-auto bg-white/5 p-4 rounded-lg border border-accent-secondary/20 research-details-scroll">
                                              <FormattedText text={paper.rawContent} />
                                            </div>
                                          </div>
                                        )}

                                        <div className="flex items-center justify-between pt-2 border-t border-accent-secondary/30">
                                          <div className="flex items-center gap-2 text-xs text-text-secondary">
                                            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
                                            </svg>
                                            <span>Research Paper #{paper.number || (index + 1)}</span>
                                          </div>
                                          <div className="flex items-center gap-1">
                                            <div className="w-2 h-2 bg-accent rounded-full animate-pulse"></div>
                                            <span className="text-xs text-accent font-medium">Active Research</span>
                                          </div>
                                        </div>
                                      </div>
                                    </div>
                                  </details>
                                </td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  )}
                </div>
              )}
            </>
          )}

          {activeTab === "saved" && (
            <>
              <h2 className="text-xl sm:text-2xl font-heading font-semibold text-accent mb-4 sm:mb-6">
                Saved Research Papers
              </h2>
              {savedPapers.length > 0 ? (
                <div className="space-y-6 sm:space-y-10">
                  {savedPapers.map((entry, index) => {
                    if (!entry.molecule || !entry.molecule.symptoms || !entry.molecule.smiles || !entry.papers) {
                      console.warn(`Skipping invalid saved research entry at index ${index}:`, entry);
                      return null;
                    }
                    const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "[]");
                    const symptomGroupIndex = storedSymptomGroups.findIndex(group => group.join(", ") === entry.molecule.symptoms);
                    const displaySymptoms = symptomGroupIndex !== -1 ? `Group ${symptomGroupIndex + 1}: ${entry.molecule.symptoms}` : entry.molecule.symptoms;

                    return (
                      <div
                        key={index}
                        className="bg-primary p-4 sm:p-6 rounded-xl border border-accent-secondary transform transition-all duration-500 animate-slide-up"
                        style={{ animationDelay: `${index * 100}ms` }}
                      >
                        <h3 className="text-lg sm:text-xl font-heading font-semibold text-text-primary mb-4">
                          Symptoms: {displaySymptoms} (SMILES: {entry.molecule.smiles})
                        </h3>
                        <div className="space-y-6 sm:space-y-8">
                          {entry.papers.map((paper, paperIndex) => (
                            <ResearchPaperCard key={paperIndex} paper={paper} index={paperIndex} />
                          ))}
                        </div>
                      </div>
                    );
                  })}
                </div>
              ) : (
                <p className="text-text-secondary font-body text-center text-sm sm:text-base transform transition-all duration-500 animate-fade-in">
                  No saved research papers found.
                </p>
              )}
            </>
          )}

          {activeTab === "generate" && (
            <>
              <h2 className="text-xl sm:text-2xl font-heading font-semibold text-accent mb-4 sm:mb-6">
                Generate Research Paper
              </h2>

              <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 sm:gap-6 mb-4 sm:mb-6">
                <div className="relative">
                  <label className="block text-sm sm:text-base font-label text-text-primary mb-2">
                    Select Symptom Group
                  </label>
                  <select
                    value={selectedSymptomGroupIndex}
                    onChange={(e) => {
                      setSelectedSymptomGroupIndex(e.target.value);
                      setSelectedSmiles(productSmilesGroups[e.target.value]?.[0] || "");
                    }}
                    className="w-full p-2 sm:p-3 border border-accent-secondary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent transition-all duration-300 text-text-primary bg-primary font-body"
                    disabled={loading || symptomGroups.length === 0}
                  >
                    {symptomGroups.length === 0 ? (
                      <option value="">No symptom groups available</option>
                    ) : (
                      symptomGroups.map((group, index) => (
                        <option key={index} value={index}>
                          {group.join(", ")}
                        </option>
                      ))
                    )}
                  </select>
                </div>
                <div className="relative">
                  <label className="block text-sm sm:text-base font-label text-text-primary mb-2">
                    Select SMILES String
                  </label>
                  <select
                    value={selectedSmiles}
                    onChange={(e) => setSelectedSmiles(e.target.value)}
                    className="w-full p-2 sm:p-3 border border-accent-secondary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent transition-all duration-300 text-text-primary bg-primary font-body"
                    disabled={loading || !selectedSymptomGroupIndex || productSmilesGroups[selectedSymptomGroupIndex]?.length === 0}
                  >
                    {productSmilesGroups[selectedSymptomGroupIndex]?.length === 0 ? (
                      <option value="">No SMILES available</option>
                    ) : (
                      productSmilesGroups[selectedSymptomGroupIndex]?.map((smiles, index) => (
                        <option key={index} value={smiles}>
                          {smiles}
                        </option>
                      ))
                    )}
                  </select>
                </div>
              </div>

              <button
                onClick={handleGeneratePaperClick}
                disabled={loading || selectedSymptomGroupIndex === "" || !selectedSmiles}
                className="w-full py-2 sm:py-3 px-4 bg-accent text-primary rounded-lg hover:bg-accent-secondary disabled:bg-gray-400 disabled:cursor-not-allowed transition-all duration-300 transform hover:scale-105 relative overflow-hidden"
              >
                <span className="relative z-10">
                  {loading ? (
                    <span className="flex items-center justify-center">
                      <svg
                        className="animate-spin h-5 w-5 mr-2 text-primary"
                        viewBox="0 0 24 24"
                      >
                        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                        <path
                          className="opacity-75"
                          fill="currentColor"
                          d="M4 12a8 8 0 018-8V0l4 4-4 4V4a4 4 0 00-4 4h-4z"
                        />
                      </svg>
                      Generating Paper...
                    </span>
                  ) : (
                    "Generate Research Paper"
                  )}
                </span>
                <span className="absolute inset-0 bg-accent-secondary opacity-0 hover:opacity-20 transition-opacity duration-300" />
              </button>

              {error && (
                <div className="mt-4 sm:mt-6 bg-error bg-opacity-10 border border-error text-error px-4 py-3 rounded-lg flex justify-between items-center transform transition-all duration-500 animate-slide-up">
                  <p className="text-sm sm:text-base font-body">{error}</p>
                  <button
                    className="text-error underline hover:text-error/80 font-label text-sm sm:text-base"
                    onClick={() => setError(null)}
                  >
                    Dismiss
                  </button>
                </div>
              )}

              {generatedPaper && (
                <div className="mt-6 sm:mt-8 p-4 sm:p-6 rounded-xl border border-accent-secondary bg-primary transform transition-all duration-500 animate-slide-up">
                  <h3 className="text-lg sm:text-xl font-heading font-semibold text-accent mb-4">
                    Generated Research Paper (IEEE Format)
                  </h3>

                  <div className="flex justify-end mb-4">
                    <button
                      onClick={() => exportToPDF(generatedPaper)}
                      className="px-3 sm:px-4 py-2 bg-success text-primary rounded-lg hover:bg-success/80 transition-all duration-300 transform hover:scale-105"
                    >
                      Download as PDF
                    </button>
                  </div>

                  <div className="space-y-6 sm:space-y-8">
                    <h4 className="text-xl sm:text-2xl font-heading font-bold text-center text-text-primary">
                      {generatedPaper.title || "Untitled"}
                    </h4>
                    <p className="text-center text-text-secondary font-body text-sm sm:text-base">
                      {generatedPaper.authors || "Unknown Authors"}
                    </p>

                    <div className="transform transition-all duration-500 animate-slide-up">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        Abstract
                      </h5>
                      <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                        {generatedPaper.abstract || "Abstract not available"}
                      </p>
                    </div>

                    <div className="transform transition-all duration-500 animate-slide-up delay-100">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        Keywords
                      </h5>
                      <p className="text-sm sm:text-base text-text-secondary font-body">
                        {Array.isArray(generatedPaper.keywords) ? generatedPaper.keywords.join(", ") : "No keywords available"}
                      </p>
                    </div>

                    <div className="transform transition-all duration-500 animate-slide-up delay-200">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        I. Introduction
                      </h5>
                      <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                        {generatedPaper.introduction || "Introduction not available"}
                      </p>
                    </div>

                    <div className="transform transition-all duration-500 animate-slide-up delay-300">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        II. Methodology
                      </h5>
                      <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                        {generatedPaper.methodology || "Methodology not available"}
                      </p>
                    </div>

                    <div className="transform transition-all duration-500 animate-slide-up delay-400">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        III. Results and Discussion
                      </h5>
                      <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                        {generatedPaper.resultsAndDiscussion || "Results and Discussion not available"}
                      </p>
                    </div>

                    <div className="transform transition-all duration-500 animate-slide-up delay-500">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        IV. Conclusion
                      </h5>
                      <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                        {generatedPaper.conclusion || "Conclusion not available"}
                      </p>
                    </div>

                    <div className="transform transition-all duration-500 animate-slide-up delay-600">
                      <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                        References
                      </h5>
                      <ul className="list-none text-sm sm:text-base space-y-2 text-text-secondary font-body">
                        {Array.isArray(generatedPaper.references) && generatedPaper.references.length > 0 ? (
                          generatedPaper.references.map((ref, index) => (
                            <li key={index}>{ref}</li>
                          ))
                        ) : (
                          <li>No references available</li>
                        )}
                      </ul>
                    </div>
                  </div>
                </div>
              )}
            </>
          )}

          {activeTab === "savedGenerated" && (
            <>
              <h2 className="text-xl sm:text-2xl font-heading font-semibold text-accent mb-4 sm:mb-6">
                Saved Generated Research Papers
              </h2>
              {savedGeneratedPapers.length > 0 ? (
                <div className="space-y-6 sm:space-y-10">
                  {savedGeneratedPapers.map((entry, index) => {
                    if (
                      !entry ||
                      !entry.molecule ||
                      !entry.molecule.symptoms ||
                      !entry.molecule.smiles ||
                      !entry.paper
                    ) {
                      console.warn(`Skipping invalid saved generated research entry at index ${index}:`, entry);
                      return null;
                    }
                    const storedSymptomGroups = JSON.parse(localStorage.getItem(`symptomGroups_${user._id}`) || "[]");
                    const symptomGroupIndex = storedSymptomGroups.findIndex(group => group.join(", ") === entry.molecule.symptoms);
                    const displaySymptoms = symptomGroupIndex !== -1 ? `Group ${symptomGroupIndex + 1}: ${entry.molecule.symptoms}` : entry.molecule.symptoms;
                    const paper = entry.paper;
                    if (!isValidPaper(paper)) {
                      console.warn(`Skipping invalid paper in entry at index ${index}:`, paper);
                      return null;
                    }

                    return (
                      <div
                        key={index}
                        className="p-4 sm:p-6 rounded-xl border border-accent-secondary bg-primary transform transition-all duration-500 animate-slide-up"
                        style={{ animationDelay: `${index * 100}ms` }}
                      >
                        <h3 className="text-lg sm:text-xl font-heading font-semibold text-text-primary mb-4">
                          Symptoms: {displaySymptoms} (SMILES: {entry.molecule.smiles})
                        </h3>
                        <div className="flex justify-end mb-4">
                          <button
                            onClick={() => exportToPDF(paper)}
                            className="px-3 sm:px-4 py-2 bg-success text-primary rounded-lg hover:bg-success/80 transition-all duration-300 transform hover:scale-105"
                          >
                            Download as PDF
                          </button>
                        </div>
                        <div className="space-y-6 sm:space-y-8">
                          <h4 className="text-xl sm:text-2xl font-heading font-bold text-center text-text-primary">
                            {paper.title || "Untitled"}
                          </h4>
                          <p className="text-center text-text-secondary font-body text-sm sm:text-base">
                            {paper.authors || "Unknown Authors"}
                          </p>

                          <div className="transform transition-all duration-500 animate-slide-up">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              Abstract
                            </h5>
                            <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                              {paper.abstract || "Abstract not available"}
                            </p>
                          </div>

                          <div className="transform transition-all duration-500 animate-slide-up delay-100">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              Keywords
                            </h5>
                            <p className="text-sm sm:text-base text-text-secondary font-body">
                              {Array.isArray(paper.keywords) ? paper.keywords.join(", ") : "No keywords available"}
                            </p>
                          </div>

                          <div className="transform transition-all duration-500 animate-slide-up delay-200">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              I. Introduction
                            </h5>
                            <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                              {paper.introduction || "Introduction not available"}
                            </p>
                          </div>

                          <div className="transform transition-all duration-500 animate-slide-up delay-300">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              II. Methodology
                            </h5>
                            <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                              {paper.methodology || "Methodology not available"}
                            </p>
                          </div>

                          <div className="transform transition-all duration-500 animate-slide-up delay-400">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              III. Results and Discussion
                            </h5>
                            <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                              {paper.resultsAndDiscussion || "Results and Discussion not available"}
                            </p>
                          </div>

                          <div className="transform transition-all duration-500 animate-slide-up delay-500">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              IV. Conclusion
                            </h5>
                            <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                              {paper.conclusion || "Conclusion not available"}
                            </p>
                          </div>

                          <div className="transform transition-all duration-500 animate-slide-up delay-600">
                            <h5 className="text-lg sm:text-xl font-heading font-semibold mb-2 text-text-primary">
                              References
                            </h5>
                            <ul className="list-none text-sm sm:text-base space-y-2 text-text-secondary font-body">
                              {Array.isArray(paper.references) && paper.references.length > 0 ? (
                                paper.references.map((ref, index) => (
                                  <li key={index}>{ref}</li>
                                ))
                              ) : (
                                <li>No references available</li>
                              )}
                            </ul>
                          </div>
                        </div>
                      </div>
                    );
                  })}
                </div>
              ) : (
                <p className="text-text-secondary font-body text-center text-sm sm:text-base transform transition-all duration-500 animate-fade-in">
                  No saved generated research papers found.
                </p>
              )}
            </>
          )}
        </div>
      </div>
    </div>
  );
}

export default Airesearchgenerator;