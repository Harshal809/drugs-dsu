"use client"

import { useState, useEffect } from "react"
import axios from "axios"
import { toast } from "react-hot-toast"
import { useAuthStore } from "../../Store/auth.store.js"
import { jsPDF } from "jspdf"
import { Copy, FileDown, Info, ChevronDown, ChevronUp, Loader2, AlertCircle, X, Beaker, Dna, Search } from "lucide-react"

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api"
const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
})

const ProteinStructureEvolution = () => {
  const [formData, setFormData] = useState({ smilesoffirst: "", smilesofsecond: "", newmoleculetitle: "" })
  const [molecules, setMolecules] = useState([])
  const [realTimeOutput, setRealTimeOutput] = useState("")
  const [formattedOutput, setFormattedOutput] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const [expandedInfoId, setExpandedInfoId] = useState(null)
  const [fetchingMolecules, setFetchingMolecules] = useState(false)
  const [searchFirst, setSearchFirst] = useState("")
  const [searchSecond, setSearchSecond] = useState("")
  const [suggestionsFirst, setSuggestionsFirst] = useState([])
  const [suggestionsSecond, setSuggestionsSecond] = useState([])
  const [loadingSuggestions, setLoadingSuggestions] = useState({ first: false, second: false })

  const { user, checkAuth, checkingAuth } = useAuthStore()

  useEffect(() => {
    const initializeApp = async () => {
      await checkAuth()
      if (!useAuthStore.getState().user) {
        setError("Authentication failed. Please log in.")
        return
      }
      await fetchAllMolecules()
    }
    initializeApp()
  }, [checkAuth])

  const fetchAllMolecules = async () => {
    if (!user?._id) return
    setFetchingMolecules(true)
    try {
      const { data } = await axiosInstance.get("/protein/generatednewmolecule")
      setMolecules(data.molecules || [])
    } catch (err) {
      console.error("Fetch error:", err)
      setMolecules([])
      setError(err.response?.data?.message || "Failed to fetch molecules.")
    } finally {
      setFetchingMolecules(false)
    }
  }

  const fetchPubChemSuggestions = async (query, type) => {
    if (!query) {
      type === "first" ? setSuggestionsFirst([]) : setSuggestionsSecond([])
      return
    }
    setLoadingSuggestions(prev => ({ ...prev, [type]: true }))
    try {
      const response = await axios.get(
        `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(query)}?limit=5`
      )
      const suggestions = response.data.dictionary_terms?.compound || []
      type === "first" ? setSuggestionsFirst(suggestions) : setSuggestionsSecond(suggestions)
    } catch (err) {
      console.error("PubChem autocomplete error:", err)
      toast.error("Failed to fetch suggestions from PubChem")
    } finally {
      setLoadingSuggestions(prev => ({ ...prev, [type]: false }))
    }
  }

  const fetchSmilesFromPubChem = async (name, type) => {
    try {
      const response = await axios.get(
        `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/property/CanonicalSMILES/JSON`
      )
      const smiles = response.data.PropertyTable.Properties[0]?.CanonicalSMILES
      if (smiles) {
        setFormData(prev => ({
          ...prev,
          [type === "first" ? "smilesoffirst" : "smilesofsecond"]: smiles
        }))
        toast.success(`SMILES fetched for ${name}`)
      } else {
        toast.error("No SMILES found for this compound")
      }
    } catch (err) {
      console.error("PubChem SMILES fetch error:", err)
      toast.error("Failed to fetch SMILES from PubChem")
    }
  }

  const formatOutput = (output) => {
    if (!output) return null;

    try {
      const result = JSON.parse(output);
      const formatted = {
        sections: [
          {
            title: "New SMILES",
            points: [result.newSmiles || "Not available"]
          },
          {
            title: "IUPAC Name",
            points: [result.newIupacName || "Not available"]
          },
          {
            title: "Conversion Details",
            points: typeof result.conversionDetails === "string"
              ? result.conversionDetails
                  .split(".")
                  .filter(sentence => sentence.trim())
                  .map(sentence => sentence.trim() + ".")
              : Array.isArray(result.conversionDetails)
              ? result.conversionDetails
              : ["Not available"]
          },
          {
            title: "Potential Diseases",
            points: typeof result.potentialDiseases === "string"
              ? result.potentialDiseases
                  .split(".")
                  .filter(sentence => sentence.trim())
                  .map(sentence => sentence.trim() + ".")
              : Array.isArray(result.potentialDiseases)
              ? result.potentialDiseases
              : ["Not available"]
          }
        ]
      };
      return formatted;
    } catch (error) {
      console.error("Error parsing output:", error);
      return {
        sections: [
          {
            title: "Error",
            points: ["Failed to parse the output. Please check the input or try again."]
          }
        ]
      };
    }
  };

  const formatMoleculeInfo = (info) => {
    if (!info) return null;

    try {
      const parsedInfo = JSON.parse(info);
      return {
        sections: [
          {
            title: "Input SMILES",
            points: [parsedInfo.inputSmiles || "Not available"]
          },
          {
            title: "New SMILES",
            points: [parsedInfo.newSmiles || "Not available"]
          },
          {
            title: "IUPAC Name",
            points: [parsedInfo.newIupacName || "Not available"]
          },
          {
            title: "Conversion Details",
            points: Array.isArray(parsedInfo.conversionDetails)
              ? parsedInfo.conversionDetails
              : typeof parsedInfo.conversionDetails === "string"
              ? parsedInfo.conversionDetails.split(".").filter(s => s.trim()).map(s => s.trim() + ".")
              : ["Not available"]
          },
          {
            title: "Potential Diseases",
            points: Array.isArray(parsedInfo.potentialDiseases)
              ? parsedInfo.potentialDiseases
              : typeof parsedInfo.potentialDiseases === "string"
              ? parsedInfo.potentialDiseases.split(".").filter(s => s.trim()).map(s => s.trim() + ".")
              : ["Not available"]
          }
        ]
      };
    } catch (error) {
      console.error("Error parsing molecule info:", error);
      return {
        sections: [
          {
            title: "Error",
            points: ["Failed to parse molecule information. Data may be corrupted."]
          }
        ]
      };
    }
  };

  const handleInputChange = (e) => {
    setFormData((prev) => ({ ...prev, [e.target.name]: e.target.value }));
  };

  const handleSearchChange = (e, type) => {
    const value = e.target.value;
    type === "first" ? setSearchFirst(value) : setSearchSecond(value);
    fetchPubChemSuggestions(value, type);
  };

  const handleSuggestionSelect = (suggestion, type) => {
    type === "first" ? setSearchFirst(suggestion) : setSearchSecond(suggestion);
    fetchSmilesFromPubChem(suggestion, type);
    type === "first" ? setSuggestionsFirst([]) : setSuggestionsSecond([]);
  };

  const handleGenerate = async (e) => {
    e.preventDefault();
    const { smilesoffirst, smilesofsecond, newmoleculetitle } = formData;
    if (!smilesoffirst || !smilesofsecond || !newmoleculetitle) {
      setError("All fields are required.");
      return;
    }
    if (!user?._id) {
      setError("Please log in to generate a molecule.");
      return;
    }

    setLoading(true);
    setError(null);
    setRealTimeOutput("");
    setFormattedOutput(null);

    try {
      const response = await fetch(`${API_BASE_URL}/protein/generatenewmolecule/${user._id}`, {
        method: "POST",
        headers: { "Content-Type": "application/json", Accept: "text/event-stream" },
        credentials: "include",
        body: JSON.stringify(formData),
      });

      if (!response.ok) throw new Error((await response.json().catch(() => ({}))).message || "Server error");

      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let fullResponse = "";

      while (true) {
        const { done, value } = await reader.read();
        if (done) {
          setLoading(false);
          toast.success("Molecule generated successfully!");
          await fetchAllMolecules();
          setFormData({ smilesoffirst: "", smilesofsecond: "", newmoleculetitle: "" });
          setSearchFirst("");
          setSearchSecond("");
          setFormattedOutput(formatOutput(fullResponse));
          break;
        }
        const chunk = decoder
          .decode(value)
          .split("\n")
          .map((line) => line.replace(/^data:\s*/, "").trim())
          .filter((line) => line && line !== "[DONE]")
          .join(" ");
        fullResponse += chunk + " ";
        setRealTimeOutput(fullResponse);
      }
    } catch (err) {
      console.error("Generate error:", err);
      setError(err.message || "Failed to generate molecule.");
      setLoading(false);
    }
  };

  const handleCopySmiles = (smiles) => {
    if (!smiles || smiles === "Not available") {
      toast.error("No SMILES available to copy.");
      return;
    }
    navigator.clipboard
      .writeText(smiles)
      .then(() => toast.success("SMILES copied!"))
      .catch((err) => toast.error("Copy failed: " + err.message));
  };

  const toggleInfo = (id) => setExpandedInfoId(expandedInfoId === id ? null : id);

  const exportToPDF = (molecule) => {
    const doc = new jsPDF();
    const margin = 20;
    const maxWidth = doc.internal.pageSize.width - 2 * margin;
    let y = margin;

    const addText = (text, size, bold = false) => {
      doc.setFontSize(size);
      doc.setFont(undefined, bold ? "bold" : "normal");
      const lines = doc.splitTextToSize(text, maxWidth);
      lines.forEach((line) => {
        if (y + size / 2 > doc.internal.pageSize.height - margin) {
          doc.addPage();
          y = margin;
        }
        doc.text(line, margin, y);
        y += size / 2 + 2;
      });
    };

    addText("Molecule Details", 16, true);
    addText(`Title: ${molecule.newmoleculetitle}`, 12);
    addText(`SMILES: ${molecule.newSmiles || "Not available"}`, 12);
    addText(`IUPAC Name: ${molecule.newIupacName || "Not available"}`, 12);
    addText(`Created: ${new Date(molecule.created).toLocaleString()}`, 12);
    addText("Information:", 12, true);

    const formattedInfo = formatMoleculeInfo(molecule.information);
    formattedInfo.sections.forEach(section => {
      addText(section.title + ":", 10, true);
      section.points.forEach(point => {
        addText(`- ${point}`, 10);
      });
    });

    doc.save(`${molecule.newmoleculetitle}_details.pdf`);
  };

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gradient-to-b from-slate-50 to-slate-100">
        <div className="flex items-center gap-2 text-lg font-semibold text-slate-700">
          <Loader2 className="h-5 w-5 animate-spin" />
          Verifying Authentication...
        </div>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen bg-gradient-to-b from-slate-50 to-slate-100">
        <div className="bg-white p-8 rounded-xl shadow-lg text-center transform transition hover:scale-105 max-w-md w-full">
          <h2 className="text-2xl font-bold text-slate-800 mb-4">Access Denied</h2>
          <p className="text-slate-600 mb-6">Please log in to explore Protein Structure Evolution.</p>
          <a
            href="/login"
            className="inline-block px-6 py-3 bg-gradient-to-r from-blue-600 to-blue-500 text-white rounded-lg hover:from-blue-700 hover:to-blue-600 transition w-full"
          >
            Login Now
          </a>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gradient-to-b from-slate-50 to-slate-100 py-12 px-4 sm:px-6">
      <div className="max-w-5xl mx-auto space-y-8">
        <div className="text-center space-y-2">
          <div className="inline-flex items-center justify-center p-2 bg-slate-100 rounded-full mb-2">
            <Dna className="h-8 w-8 text-blue-600" />
          </div>
          <h1 className="text-4xl font-extrabold tracking-tight bg-gradient-to-r from-blue-600 to-blue-500 text-transparent bg-clip-text">
            New Drug Molecule Generation
          </h1>
          <p className="text-sm text-slate-600 font-medium">Powered by Gemini AI</p>
        </div>

        {error && (
          <div className="mb-8 p-4 bg-red-50 border-l-4 border-red-500 text-red-700 rounded-lg flex justify-between items-center shadow-md animate-in fade-in slide-in-from-top-5 duration-300">
            <div className="flex items-start">
              <AlertCircle className="h-5 w-5 mt-0.5 mr-2 flex-shrink-0" />
              <span>{error}</span>
            </div>
            <button
              onClick={() => setError(null)}
              className="text-red-900 hover:text-red-600 font-medium p-1 rounded-full hover:bg-red-100"
            >
              <X className="h-5 w-5" />
            </button>
          </div>
        )}

        {/* Generate New Molecule */}
        <div className="bg-white p-6 rounded-xl shadow-md hover:shadow-lg transition-shadow duration-300">
          <div className="mb-4">
            <h2 className="text-2xl font-semibold text-slate-800 flex items-center gap-2">
              <Beaker className="h-6 w-6 text-blue-500" />
              Generate New Molecule
            </h2>
            <p className="text-slate-500 text-sm mt-1">
              Search for molecules by name to fetch SMILES from PubChem or enter SMILES manually
            </p>
          </div>

          <form onSubmit={handleGenerate} className="space-y-5">
            <div className="space-y-4">
              <div className="grid gap-2">
                <label htmlFor="newmoleculetitle" className="text-sm font-medium text-slate-700">
                  New Molecule Title
                </label>
                <input
                  id="newmoleculetitle"
                  name="newmoleculetitle"
                  type="text"
                  value={formData.newmoleculetitle}
                  onChange={handleInputChange}
                  placeholder="Enter a descriptive title for the new molecule"
                  disabled={loading}
                  required
                  className="w-full p-3 border border-slate-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition bg-white"
                />
              </div>

              {/* First Molecule Search */}
              <div className="grid gap-2 relative">
                <label htmlFor="searchFirst" className="text-sm font-medium text-slate-700">
                  Search First Molecule
                </label>
                <div className="relative">
                  <input
                    id="searchFirst"
                    type="text"
                    value={searchFirst}
                    onChange={(e) => handleSearchChange(e, "first")}
                    placeholder="Search molecule by name (e.g., Aspirin)"
                    disabled={loading}
                    className="w-full p-3 pl-10 border border-slate-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition bg-white"
                  />
                  <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-slate-400" />
                </div>
                {loadingSuggestions.first && (
                  <div className="absolute top-full left-0 w-full bg-white border border-slate-200 rounded-b-lg shadow-md p-2 z-10">
                    <Loader2 className="h-4 w-4 animate-spin mx-auto" />
                  </div>
                )}
                {suggestionsFirst.length > 0 && !loadingSuggestions.first && (
                  <ul className="absolute top-full left-0 w-full bg-white border border-slate-200 rounded-b-lg shadow-md max-h-40 overflow-auto z-10">
                    {suggestionsFirst.map((suggestion, index) => (
                      <li
                        key={index}
                        onClick={() => handleSuggestionSelect(suggestion, "first")}
                        className="px-3 py-2 text-sm text-slate-700 hover:bg-blue-50 cursor-pointer"
                      >
                        {suggestion}
                      </li>
                    ))}
                  </ul>
                )}
                <input
                  id="smilesoffirst"
                  name="smilesoffirst"
                  type="text"
                  value={formData.smilesoffirst}
                  onChange={handleInputChange}
                  placeholder="SMILES will be auto-filled or enter manually"
                  disabled={loading}
                  required
                  className="w-full p-3 border border-slate-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition bg-white font-mono text-sm"
                />
              </div>

              {/* Second Molecule Search */}
              <div className="grid gap-2 relative">
                <label htmlFor="searchSecond" className="text-sm font-medium text-slate-700">
                  Search Second Molecule
                </label>
                <div className="relative">
                  <input
                    id="searchSecond"
                    type="text"
                    value={searchSecond}
                    onChange={(e) => handleSearchChange(e, "second")}
                    placeholder="Search molecule by name (e.g., Ibuprofen)"
                    disabled={loading}
                    className="w-full p-3 pl-10 border border-slate-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition bg-white"
                  />
                  <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-slate-400" />
                </div>
                {loadingSuggestions.second && (
                  <div className="absolute top-full left-0 w-full bg-white border border-slate-200 rounded-b-lg shadow-md p-2 z-10">
                    <Loader2 className="h-4 w-4 animate-spin mx-auto" />
                  </div>
                )}
                {suggestionsSecond.length > 0 && !loadingSuggestions.second && (
                  <ul className="absolute top-full left-0 w-full bg-white border border-slate-200 rounded-b-lg shadow-md max-h-40 overflow-auto z-10">
                    {suggestionsSecond.map((suggestion, index) => (
                      <li
                        key={index}
                        onClick={() => handleSuggestionSelect(suggestion, "second")}
                        className="px-3 py-2 text-sm text-slate-700 hover:bg-blue-50 cursor-pointer"
                      >
                        {suggestion}
                      </li>
                    ))}
                  </ul>
                )}
                <input
                  id="smilesofsecond"
                  name="smilesofsecond"
                  type="text"
                  value={formData.smilesofsecond}
                  onChange={handleInputChange}
                  placeholder="SMILES will be auto-filled or enter manually"
                  disabled={loading}
                  required
                  className="w-full p-3 border border-slate-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition bg-white font-mono text-sm"
                />
              </div>
            </div>

            <button
              type="submit"
              disabled={loading}
              className="w-full py-3 bg-gradient-to-r from-blue-600 to-blue-500 text-white rounded-lg hover:from-blue-700 hover:to-blue-600 disabled:opacity-50 disabled:cursor-not-allowed transition flex items-center justify-center font-medium"
            >
              {loading ? (
                <>
                  <Loader2 className="mr-2 h-5 w-5 animate-spin" />
                  Generating...
                </>
              ) : (
                "Generate Molecule"
              )}
            </button>
          </form>

          {realTimeOutput && (
            <div className="mt-6">
              <div className="flex items-center gap-2 mb-2">
                <div className="h-2 w-2 rounded-full bg-blue-500 animate-pulse" />
                <h3 className="text-sm font-medium text-slate-700">Live Output</h3>
              </div>
              <div className="h-96 w-full overflow-auto rounded-md border border-slate-200 bg-slate-50 p-4">
                {formattedOutput ? (
                  <div className="text-sm text-slate-700">
                    {formattedOutput.sections.map((section, index) => (
                      <div key={index} className="mb-4">
                        <h4 className="text-base font-semibold text-slate-800 mb-2">{section.title}</h4>
                        <ul className="list-disc pl-5 space-y-1">
                          {section.points.map((point, idx) => (
                            <li key={idx} className="text-sm text-slate-600">
                              {point}
                            </li>
                          ))}
                        </ul>
                      </div>
                    ))}
                  </div>
                ) : (
                  <div className="text-sm text-slate-700 whitespace-pre-wrap font-mono">
                    {realTimeOutput}
                  </div>
                )}
              </div>
            </div>
          )}
        </div>

        {/* Your Molecules */}
        <div className="bg-white p-6 rounded-xl shadow-md hover:shadow-lg transition-shadow duration-300">
          <div className="flex items-center justify-between mb-4">
            <h2 className="text-2xl font-semibold text-slate-800 flex items-center gap-2">
              <Dna className="h-6 w-6 text-blue-500" />
              Your Molecules
            </h2>
            <span className="inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-semibold bg-slate-100 text-slate-800 border border-slate-200">
              {molecules.length} total
            </span>
          </div>
          <p className="text-slate-500 text-sm mb-4">View and manage your previously generated molecules</p>

          {fetchingMolecules && (
            <div className="flex items-center justify-center py-8">
              <Loader2 className="h-8 w-8 animate-spin text-blue-500" />
            </div>
          )}

          {!fetchingMolecules && molecules.length === 0 && (
            <div className="text-center py-12 text-slate-500 border border-dashed border-slate-200 rounded-lg">
              <Beaker className="h-12 w-12 mx-auto mb-3 text-slate-300" />
              <p className="text-lg font-medium">No molecules generated yet</p>
              <p className="text-sm mt-1">Use the form above to create your first molecule</p>
            </div>
          )}

          {!fetchingMolecules && molecules.length > 0 && (
            <div className="h-[32rem] overflow-auto pr-1">
              <ul className="space-y-4">
                {molecules.map((molecule) => (
                  <li key={molecule.id} className="group">
                    <div className="border border-slate-200 rounded-lg overflow-hidden transition-all duration-200 group-hover:border-blue-200 group-hover:shadow-sm">
                      <div className="p-4">
                        <div className="flex flex-col space-y-3">
                          <div className="space-y-1">
                            <div className="flex items-center justify-between">
                              <h3 className="font-medium text-slate-900">{molecule.newmoleculetitle}</h3>
                              <span className="inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium bg-slate-100 text-slate-800">
                                {new Date(molecule.created).toLocaleDateString()}
                              </span>
                            </div>

                            <div className="flex items-center space-x-2 bg-slate-50 p-2 rounded-md">
                              <div className="flex-1 overflow-hidden">
                                <p className="text-xs text-slate-500">SMILES:</p>
                                <p className="text-sm font-mono text-slate-700 truncate">
                                  {molecule.newSmiles || "Not available"}
                                </p>
                              </div>
                              <button
                                onClick={() => handleCopySmiles(molecule.newSmiles)}
                                className="p-1.5 text-slate-500 hover:text-blue-600 hover:bg-slate-100 rounded-md"
                              >
                                <Copy className="h-4 w-4" />
                              </button>
                            </div>
                          </div>

                          <div className="flex space-x-2">
                            <button
                              onClick={() => toggleInfo(molecule.id)}
                              className="flex-1 flex items-center justify-center gap-1 px-3 py-2 border border-slate-200 bg-white hover:bg-slate-50 text-slate-700 rounded-md text-sm font-medium"
                            >
                              <Info className="h-4 w-4" />
                              {expandedInfoId === molecule.id ? "Hide Info" : "Show Info"}
                              {expandedInfoId === molecule.id ? (
                                <ChevronUp className="h-3 w-3 ml-auto" />
                              ) : (
                                <ChevronDown className="h-3 w-3 ml-auto" />
                              )}
                            </button>
                            <button
                              onClick={() => exportToPDF(molecule)}
                              className="flex-1 flex items-center justify-center gap-1 px-3 py-2 border border-slate-200 bg-white hover:bg-slate-50 text-slate-700 rounded-md text-sm font-medium"
                            >
                              <FileDown className="h-4 w-4" />
                              Export PDF
                            </button>
                          </div>
                        </div>

                        {expandedInfoId === molecule.id && (
                          <>
                            <div className="h-[1px] w-full bg-slate-200 my-3" />
                            <div className="h-96 w-full overflow-auto rounded-md bg-slate-50 p-3">
                              <div className="text-sm text-slate-700">
                                {formatMoleculeInfo(molecule.information) ? (
                                  formatMoleculeInfo(molecule.information).sections.map((section, index) => (
                                    <div key={index} className="mb-4">
                                      <h4 className="text-base font-semibold text-slate-800 mb-2">{section.title}</h4>
                                      <ul className="list-disc pl-5 space-y-1">
                                        {section.points.map((point, idx) => (
                                          <li key={idx} className="text-sm text-slate-600">
                                            {point}
                                          </li>
                                        ))}
                                      </ul>
                                    </div>
                                  ))
                                ) : (
                                  <p>No information available for this molecule.</p>
                                )}
                              </div>
                            </div>
                          </>
                        )}
                      </div>
                    </div>
                  </li>
                ))}
              </ul>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default ProteinStructureEvolution;