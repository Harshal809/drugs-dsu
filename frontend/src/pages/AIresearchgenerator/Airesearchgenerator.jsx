import React, { useState, useEffect, useRef } from "react";
import axios from "axios";
import { toast } from "react-hot-toast";
import { useAuthStore } from "../../Store/auth.store.js";
import jsPDF from "jspdf";



const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
});





// Utility functions
const cleanAbstract = (abstract) => {
  if (!abstract || abstract === "Abstract not available.") return abstract;
  return abstract
    .replace(/<\/?[^>]+(>|$)/g, "")
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
  cleaned = cleaned.replace(/^`+|`+$/g, "").trim();
  const jsonMatch = cleaned.match(/(\[.*?\]|\{.*?\})/s);
  return jsonMatch ? jsonMatch[0] : cleaned;
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
    const prompt = `
You are a biomedical research assistant. Given the symptoms: "${symptoms}", your tasks are:

1. Generate intelligent academic search queries relevant to these symptoms and potential SMILES-derived compounds.
2. Search for highly relevant, peer-reviewed research articles using official APIs or scraping of trusted sources such as IEEE Xplore, PubMed, and CrossRef.
3. Extract and return real-time research data in strict JSON format with the following keys:

[
  {
    "title": "Research paper title",
    "authors": "Author1, Author2, ...",
    "year": "Publication year (YYYY)",
    "abstract": "Full abstract content of the paper",
    "doi": "Digital Object Identifier (DOI)",
    "url": "Direct link to the research paper",
    "is_simulated": false
  }
]

⚠️ Important Guidelines:
- Do NOT return any hardcoded, static, or example content.
- If no actual papers are found after querying all sources, only then return 3–5 realistic simulated entries with "is_simulated": true.
- Ensure each result is verifiable and highly relevant to the symptoms.

Make sure all metadata is accurate and based on **real-time** web search or API/scraping.
`;

    const response = await axiosInstance.post("/researchPaper/proxy/gemini", { prompt });
    const content = response.data.content;

    try {
      const cleanedContent = cleanGeminiResponse(content);
      const papers = JSON.parse(cleanedContent);

      if (!Array.isArray(papers)) {
        throw new Error("Invalid papers format");
      }

      return papers.map(paper => ({
        ...paper,
        is_simulated: paper.is_simulated || false,
        abstract: cleanAbstract(paper.abstract),
      }));
    } catch (parseError) {
      console.error("Error parsing Gemini response:", parseError);
      return generateFallbackPapers(smiles, symptoms);
    }
  } catch (err) {
    console.error("Error fetching research papers:", err.response?.data || err.message);
    return generateFallbackPapers(smiles, symptoms);
  }
};


  const generateResearchPaper = async (symptoms, smiles) => {
    try {
      const prompt = `
        You are an expert in chemical informatics and academic writing. Given the SMILES string "${smiles}" and associated symptoms "${symptoms}", generate a high-quality research paper in IEEE format. The paper should be structured with the following sections: Title, Authors, Abstract, Keywords, I. Introduction, II. Methodology, III. Results and Discussion, IV. Conclusion, References. Ensure the content is scientifically accurate, relevant to the compound's therapeutic applications for the symptoms, and includes realistic data and references. Return the paper in clean JSON format (no Markdown or code fences) with fields: title (string), authors (string), abstract (string), keywords (array of strings), introduction (string), methodology (string), resultsAndDiscussion (string), conclusion (string), references (array of strings). Example output:
        {
          "title": "Therapeutic Applications of Compound X",
          "authors": "Test",
          "abstract": "This paper investigates...",
          "keywords": ["compound X", "therapeutics", "symptoms"],
          "introduction": "Introduction text...",
          "methodology": "Methodology text...",
          "resultsAndDiscussion": "Results text...",
          "conclusion": "Conclusion text...",
          "references": ["Ref 1", "Ref 2"]
        }
      `;
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
      const papers = await fetchResearchPapers(selectedSymptoms, selectedSmiles);
      setResearchPapers(papers);
      setResearchSummary(`Found ${papers.length} research papers related to the compound with SMILES "${selectedSmiles}" for treating symptoms: ${selectedSymptoms}.`);
      await savePapers(selectedSymptoms, selectedSmiles, papers);
    } catch (err) {
      setError("Failed to fetch research papers. Please try again.");
      toast.error("Failed to fetch research papers");
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

  const constructIeeeUrl = (doi) => {
    if (!doi) return `https://ieeexplore.ieee.org/document/${generateUniqueDoi()}`;
    const doiSuffix = doi.includes("10.") ? doi.split("/").pop() : doi.replace(/[^0-9]/g, "");
    return `https://ieeexplore.ieee.org/document/${doiSuffix}`;
  };

  const generateUniqueDoi = (year = 2025, index = 0) => {
    const randomId = Math.floor(Math.random() * 1000000) + index;
    return `10.1109/TBME.${year}.${randomId}`;
  };

  const generateFallbackPapers = (smiles, symptoms) => {
    const inferredFeatures = smiles.includes("c") ? "aromatic rings" : "aliphatic chains";
    return [
      {
        title: `Computational Analysis of ${inferredFeatures} for ${symptoms} Treatment`,
        authors: "R. Johnson, S. Lee, T. Brown",
        year: "2022",
        abstract: `This paper explores the role of ${inferredFeatures} derived from SMILES structures like ${smiles.substring(0, 10)}... in predicting drug efficacy for symptoms including ${symptoms}.`,
        doi: "10.1109/TBME.2022.987654",
        url: "https://ieeexplore.ieee.org/document/987654",
        is_simulated: true,
      },
      {
        title: `Synthesis and Properties of Novel ${inferredFeatures} Compounds for ${symptoms}`,
        authors: "M. Davis, P. Kim, Q. Zhang",
        year: "2020",
        abstract: `Investigates synthesis pathways for compounds with ${inferredFeatures}, leveraging SMILES-based modeling for therapeutic applications targeting ${symptoms}.`,
        doi: "10.1109/TBME.2020.876543",
        url: "https://ieeexplore.ieee.org/document/876543",
        is_simulated: true,
      },
    ];
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

                  {researchSummary && (
                    <div className="mb-6 sm:mb-8 transform transition-all duration-500 delay-100 animate-slide-up">
                      <h4 className="text-md sm:text-lg font-heading font-semibold text-text-primary mb-2">
                        Research Context
                      </h4>
                      <p className="text-sm sm:text-base leading-relaxed text-text-secondary font-body">
                        {researchSummary}
                      </p>
                    </div>
                  )}

                  {researchPapers.length > 0 && (
                    <div className="mb-6 sm:mb-8">
                      <h4 className="text-md sm:text-lg font-heading font-semibold text-text-primary mb-4">
                        Newly Fetched Research Papers
                      </h4>
                      <div className="space-y-6 sm:space-y-8">
                        {researchPapers.map((paper, index) => (
                          <div
                            key={index}
                            className="border-l-4 border-accent pl-4 transform transition-all duration-500 animate-slide-up"
                            style={{ animationDelay: `${index * 100}ms` }}
                          >
                            <h5 className="text-lg sm:text-xl font-heading font-bold text-text-primary mb-2 uppercase">
                              {paper.title}
                            </h5>
                            <p className="text-text-secondary font-body text-sm sm:text-base mb-1">
                              <span className="font-bold">Authors:</span> {paper.authors}
                            </p>
                            <p className="text-text-secondary font-body text-sm sm:text-base mb-1">
                              <span className="font-bold">Published:</span> {paper.year}
                            </p>
                            <p className="text-text-secondary font-body text-sm sm:text-base mb-1">
                              <span className="font-bold">Abstract:</span> {paper.abstract}
                            </p>
                            {paper.doi && paper.url !== "No URL available" && (
                              <p className="text-accent font-body text-sm sm:text-base">
                                <span className="font-bold">DOI:</span>{" "}
                                <a
                                  href={paper.url}
                                  target="_blank"
                                  rel="noopener noreferrer"
                                  className="underline hover:text-accent-secondary transition-colors duration-300"
                                >
                                  {paper.doi}
                                </a>
                              </p>
                            )}
                            {paper.is_simulated && (
                              <p className="text-text-secondary font-body text-sm sm:text-base italic">
                                Note: This is a simulated paper generated for illustrative purposes.
                              </p>
                            )}
                          </div>
                        ))}
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
                            <div
                              key={paperIndex}
                              className="border-l-4 border-accent pl-4 transform transition-all duration-500 animate-slide-up"
                              style={{ animationDelay: `${paperIndex * 100}ms` }}
                            >
                              <h5 className="text-lg sm:text-xl font-heading font-bold text-text-primary mb-2 uppercase">
                                {paper.title || "Untitled"}
                              </h5>
                              <p className="text-text-secondary font-body text-sm sm:text-base mb-1">
                                <span className="font-bold">Authors:</span> {paper.authors || "Unknown Authors"}
                              </p>
                              <p className="text-text-secondary font-body text-sm sm:text-base mb-1">
                                <span className="font-bold">Published:</span> {paper.year || "Unknown Year"}
                              </p>
                              <p className="text-text-secondary font-body text-sm sm:text-base mb-1">
                                <span className="font-bold">Abstract:</span> {paper.abstract || "Abstract not available"}
                              </p>
                              {paper.doi && paper.url && paper.url !== "No URL available" && (
                                <p className="text-accent font-body text-sm sm:text-base">
                                  <span className="font-bold">DOI:</span>{" "}
                                  <a
                                    href={paper.url}
                                    target="_blank"
                                    rel="noopener noreferrer"
                                    className="underline hover:text-accent-secondary transition-colors duration-300"
                                  >
                                    {paper.doi}
                                  </a>
                                </p>
                              )}
                              {paper.is_simulated && (
                                <p className="text-text-secondary font-body text-sm sm:text-base italic">
                                  Note: This is a simulated paper generated for illustrative purposes.
                                </p>
                              )}
                            </div>
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