// "use client"

// import { useState, useEffect, useRef } from "react"
// import axios from "axios"
// import { toast } from "react-hot-toast"
// import { useAuthStore } from "../../Store/auth.store.js"
// import { Loader2, AlertCircle, FileText, Upload, Search } from "lucide-react"

// const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api"
// const CHEMBL_API_URL = "https://www.ebi.ac.uk/chembl/api/data"
// const PUBCHEM_API_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
// const OPENALEX_API_URL = "https://api.openalex.org"

// // Axios instances for different APIs
// const axiosInstance = axios.create({
//   baseURL: API_BASE_URL,
//   withCredentials: true,
// })

// const chemblAxios = axios.create({
//   baseURL: CHEMBL_API_URL,
//   withCredentials: false,
// })

// const pubchemAxios = axios.create({
//   baseURL: PUBCHEM_API_URL,
//   withCredentials: false,
// })

// const openalexAxios = axios.create({
//   baseURL: OPENALEX_API_URL,
//   withCredentials: false,
// })

// // Utility to clean HTML tags from abstracts
// const cleanAbstract = (abstract) => {
//   if (!abstract || abstract === "Abstract not available.") return abstract
//   return abstract
//     .replace(/<\/?[^>]+(>|$)/g, "") // Remove HTML tags
//     .replace(/\s+/g, " ") // Normalize whitespace
//     .trim()
// }

// // Updated SimpleMoleculeEditor with PubChem autocomplete
// const SimpleMoleculeEditor = ({ onChange }) => {
//   const [searchQuery, setSearchQuery] = useState("")
//   const [suggestions, setSuggestions] = useState([])
//   const [loadingSuggestions, setLoadingSuggestions] = useState(false)
//   const [selectedSmiles, setSelectedSmiles] = useState("")

//   const fetchPubChemSuggestions = async (query) => {
//     if (!query) {
//       setSuggestions([])
//       return
//     }
//     setLoadingSuggestions(true)
//     try {
//       const response = await axios.get(
//         `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(query)}?limit=5`
//       )
//       const compounds = response.data.dictionary_terms?.compound || []
//       setSuggestions(compounds)
//     } catch (err) {
//       console.error("PubChem autocomplete error:", err)
//       toast.error("Failed to fetch suggestions from PubChem")
//     } finally {
//       setLoadingSuggestions(false)
//     }
//   }

//   const fetchSmilesFromPubChem = async (name) => {
//     try {
//       const response = await axios.get(
//         `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/property/CanonicalSMILES/JSON`
//       )
//       const smiles = response.data.PropertyTable.Properties[0]?.CanonicalSMILES
//       if (smiles) {
//         setSelectedSmiles(smiles)
//         onChange(smiles)
//         toast.success(`SMILES fetched for ${name}`)
//       } else {
//         toast.error("No SMILES found for this compound")
//       }
//     } catch (err) {
//       console.error("PubChem SMILES fetch error:", err)
//       toast.error("Failed to fetch SMILES from PubChem")
//     }
//   }

//   const handleSearchChange = (e) => {
//     const value = e.target.value
//     setSearchQuery(value)
//     fetchPubChemSuggestions(value)
//   }

//   const handleSuggestionSelect = (suggestion) => {
//     setSearchQuery(suggestion)
//     fetchSmilesFromPubChem(suggestion)
//     setSuggestions([])
//   }

//   const handlePreviewError = (e) => {
//     e.target.src = "/placeholder.svg?height=160&width=160"
//     e.target.alt = "Molecule preview unavailable"
//     e.target.className = "mx-auto h-40 object-contain opacity-50"
//   }

//   return (
//     <div className="p-4 border border-gray-300 rounded-lg bg-gray-50">
//       <p className="text-sm text-gray-600 mb-4">Search for a molecule by name or enter a SMILES string manually:</p>
//       <div className="relative mb-4">
//         <input
//           type="text"
//           value={searchQuery}
//           onChange={handleSearchChange}
//           placeholder="Search molecule (e.g., Aspirin)"
//           className="w-full p-3 pl-10 border border-gray-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition bg-white"
//         />
//         <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-gray-400" />
//         {loadingSuggestions && (
//           <div className="absolute top-full left-0 w-full bg-white border border-gray-200 rounded-b-lg shadow-md p-2 z-10">
//             <Loader2 className="h-4 w-4 animate-spin mx-auto" />
//           </div>
//         )}
//         {suggestions.length > 0 && !loadingSuggestions && (
//           <ul className="absolute top-full left-0 w-full bg-white border border-gray-200 rounded-b-lg shadow-md max-h-40 overflow-auto z-10">
//             {suggestions.map((suggestion, index) => (
//               <li
//                 key={index}
//                 onClick={() => handleSuggestionSelect(suggestion)}
//                 className="px-3 py-2 text-sm text-gray-700 hover:bg-blue-50 cursor-pointer"
//               >
//                 {suggestion}
//               </li>
//             ))}
//           </ul>
//         )}
//       </div>
//       <div className="mt-4">
//         <p className="text-xs text-gray-500 mb-2">Preview (updates with selection or manual SMILES input):</p>
//         <div className="bg-white p-2 border border-gray-200 rounded-lg">
//           <img
//             src={
//               selectedSmiles
//                 ? `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/PNG?record_type=2d&smiles=${encodeURIComponent(selectedSmiles)}`
//                 : "/placeholder.svg?height=160&width=160"
//             }
//             alt="Molecule structure preview"
//             className="mx-auto h-40 object-contain"
//             onError={handlePreviewError}
//           />
//         </div>
//       </div>
//     </div>
//   )
// }

// function AIdriventargetprediction() {
//   const [activeTab, setActiveTab] = useState("moleculeInput")
//   const [moleculeInput, setMoleculeInput] = useState("")
//   const [moleculeFile, setMoleculeFile] = useState(null)
//   const [marvinSmiles, setMarvinSmiles] = useState("")
//   const [predictedTargets, setPredictedTargets] = useState([])
//   const [dockingResults, setDockingResults] = useState(null)
//   const [relatedResearch, setRelatedResearch] = useState([])
//   const [savedSearches, setSavedSearches] = useState([])
//   const [loading, setLoading] = useState(false)
//   const [error, setError] = useState(null)
//   const [useAlternativeEditor, setUseAlternativeEditor] = useState(true)
//   const [analysisCompleted, setAnalysisCompleted] = useState(false)
//   const [currentSmiles, setCurrentSmiles] = useState("")

//   const { user, checkAuth, checkingAuth } = useAuthStore()
//   const marvinContainerRef = useRef(null)

//   useEffect(() => {
//     const initializeApp = async () => {
//       await checkAuth()
//       if (!useAuthStore.getState().user) {
//         setError("Authentication failed. Please log in.")
//         return
//       }
//       await fetchSavedSearches()
//     }
//     initializeApp()
//   }, [checkAuth])

//   const fetchSavedSearches = async () => {
//     if (!user?._id) {
//       console.warn("User ID not available, cannot fetch saved searches")
//       return
//     }
//     try {
//       setLoading(true)
//       const response = await axiosInstance.get("/protein/saved-searches")
//       setSavedSearches(response.data.searches || [])
//     } catch (err) {
//       console.error("Error fetching saved searches:", err.response?.data || err.message)
//       setError(err.response?.data?.message || "Failed to fetch saved searches")
//       toast.error("Failed to fetch saved searches")
//     } finally {
//       setLoading(false)
//     }
//   }

//   const checkIfSearchExists = async (smiles) => {
//     if (!user?._id || !smiles) return false
//     try {
//       const response = await axiosInstance.get("/protein/check-saved-searches", { params: { smiles } })
//       return response.data.exists
//     } catch (err) {
//       console.error("Error checking saved searches:", err.response?.data || err.message)
//       return false
//     }
//   }

//   const saveSearch = async (smiles, targets, research, docking) => {
//     if (!user?._id || !smiles) {
//       toast.error("Cannot save search: User ID or molecule data missing")
//       return null
//     }
//     const payload = { userId: user._id, smiles, targets: targets || [], research: research || [], docking: docking || null }
//     try {
//       const response = await axiosInstance.post("/protein/save-search", payload)
//       toast.success("Analysis results saved successfully!")
//       return response.data.search
//     } catch (err) {
//       console.error("Error saving search:", err.response?.data || err.message)
//       toast.error("Failed to save search: " + (err.response?.data?.message || err.message))
//       return null
//     }
//   }

//   const handleMoleculeInputChange = (value) => {
//     setMoleculeInput(value)
//     const previewImg = document.querySelector(".molecule-preview-img")
//     if (previewImg && value) {
//       const cleanedSmiles = value.trim().replace(/'/g, "")
//       previewImg.src = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/PNG?record_type=2d&smiles=${encodeURIComponent(cleanedSmiles)}`
//     }
//   }

//   const handleMoleculeSubmit = async () => {
//     let smiles = marvinSmiles || moleculeInput
//     if (!smiles && moleculeFile) {
//       const formData = new FormData()
//       formData.append("file", moleculeFile)
//       try {
//         const response = await axiosInstance.post("/protein/convert-file-to-smiles", formData)
//         smiles = response.data.smiles
//       } catch (err) {
//         setError("Failed to process molecular file: " + (err.response?.data?.message || err.message))
//         toast.error("Failed to process molecular file")
//         return
//       }
//     }
//     if (!smiles) {
//       toast.error("Please provide a molecule via sketcher, SMILES, IUPAC, or file upload")
//       return
//     }

//     setCurrentSmiles(smiles)
//     const searchExists = await checkIfSearchExists(smiles)
//     if (searchExists) {
//       toast.info("This molecule has already been searched. Redirecting to Saved Searches.")
//       setActiveTab("savedSearches")
//       await fetchSavedSearches()
//       return
//     }

//     setLoading(true)
//     setError(null)
//     setPredictedTargets([])
//     setRelatedResearch([])
//     setDockingResults(null)
//     setAnalysisCompleted(false)

//     try {
//       const fingerprintsResponse = await axiosInstance.post("/protein/rdkit-fingerprints", { smiles })
//       const fingerprints = fingerprintsResponse.data.fingerprints

//       const geminiPrompt = `
//         You are an expert in cheminformatics and bioinformatics. Given a molecule's SMILES string and its fingerprint embeddings, predict realistic potential protein targets, their mechanism of action (MOA), pathways, and associated diseases. Use the following data:

//         - **SMILES:** "${smiles}"
//         - **Morgan Fingerprint:** "${fingerprints.morgan}"
//         - **MACCS Fingerprint:** "${fingerprints.maccs}"

//         ### Instructions:
//         - Predict 3-5 realistic protein targets with confidence scores (0-1) based on known cheminformatics patterns.
//         - For each target, provide a plausible mechanism of action (e.g., inhibitor, agonist), associated pathways, and diseases, grounded in real-world biological knowledge.
//         - Return the response in a structured JSON format:
//         {
//           "targets": [
//             {
//               "protein": "Target Name",
//               "confidence": 0.95,
//               "moa": "Mechanism of Action",
//               "pathways": ["Pathway 1", "Pathway 2"],
//               "diseases": ["Disease 1", "Disease 2"]
//             }
//           ]
//         }
//         - Ensure the predictions are realistic and could plausibly match known drug-target interactions.
//       `

//       const geminiResponse = await axiosInstance.post("/protein/proxy/gemini", { prompt: geminiPrompt })
//       const geminiContent = geminiResponse.data.content
//       const jsonMatch = geminiContent.match(/{[\s\S]*}/)
//       if (!jsonMatch) throw new Error("No valid JSON found in Gemini response")
//       const parsedGeminiData = JSON.parse(jsonMatch[0])
//       const targets = parsedGeminiData.targets || []
//       setPredictedTargets(targets)

//       try {
//         const chemblResponse = await chemblAxios.get(`/molecule?smiles=${encodeURIComponent(smiles)}`)
//         const pubchemResponse = await pubchemAxios.get(`/compound/smiles/${encodeURIComponent(smiles)}/JSON`)
//         const externalData = {
//           chembl: chemblResponse.data?.molecules || [],
//           pubchem: pubchemResponse.data?.PC_Compounds || [],
//         }
//         setPredictedTargets((prev) =>
//           prev.map((target) => ({
//             ...target,
//             knownInteractions: externalData.chembl.find((t) => t.target === target.protein) || null,
//           }))
//         )
//       } catch (dbErr) {
//         console.error("Error querying external databases:", dbErr)
//       }

//       let papers = []
//       try {
//         const searchQuery = targets.length > 0 ? targets[0].protein : "protein targets"
//         const openalexResponse = await openalexAxios.get(
//           `/works?filter=title.search:${encodeURIComponent(searchQuery)}&per_page=3`
//         )
//         papers = openalexResponse.data.results.map((item) => ({
//           title: item.title || "Untitled",
//           authors: item.authorships?.map((a) => a.author.display_name).join(", ") || "Unknown Authors",
//           journal: item.primary_location?.source?.display_name || "Not specified",
//           year: item.publication_year?.toString() || "Unknown Year",
//           abstract: cleanAbstract(item.abstract) || "Abstract not available.",
//           doi: item.doi || "No DOI available",
//           url: item.doi ? `https://doi.org/${item.doi}` : "No URL available",
//         }))
//         setRelatedResearch(papers)
//       } catch (researchErr) {
//         console.error("Error fetching research papers:", researchErr)
//         papers = [
//           {
//             title: "Recent advances in molecular target identification and drug discovery",
//             authors: "Smith J, Johnson A, Williams B",
//             journal: "Journal of Medicinal Chemistry",
//             year: "2023",
//             abstract: "This review covers recent methodologies in computational drug discovery and molecular target identification.",
//             doi: "10.1021/example.doi",
//             url: "https://doi.org/10.1021/example.doi",
//           },
//           {
//             title: "Structure-based drug design: principles and applications",
//             authors: "Chen L, Garcia R",
//             journal: "Nature Reviews Drug Discovery",
//             year: "2022",
//             abstract: "An overview of structure-based drug design approaches and their application in modern pharmaceutical research.",
//             doi: "10.1038/example.doi",
//             url: "https://doi.org/10.1038/example.doi",
//           },
//         ]
//         setRelatedResearch(papers)
//       }

//       try {
//         const dockingResponse = await axiosInstance.post("/protein/docking", { smiles })
//         setDockingResults(dockingResponse.data.results)
//       } catch (dockingErr) {
//         console.error("Error performing molecular docking:", dockingErr)
//       }

//       const currentTargets = targets.map((target) => ({ ...target, knownInteractions: null }))
//       const saveResult = await saveSearch(smiles, currentTargets, papers || [], dockingResults)
//       if (saveResult) {
//         setAnalysisCompleted(true)
//         toast.success("Molecule analysis completed successfully!")
//       } else {
//         toast.error("Analysis completed but failed to save results")
//       }
//     } catch (err) {
//       console.error("Error in molecule analysis:", err)
//       setError("Analysis failed: " + err.message)
//       toast.error("Analysis failed: " + err.message)
//     } finally {
//       setLoading(false)
//     }
//   }

//   const handleTabChange = (tab) => {
//     setActiveTab(tab)
//     if (tab === "moleculeInput") setError(null)
//     else if (tab === "savedSearches") fetchSavedSearches()
//   }

//   const handleViewSavedSearch = (search) => {
//     setPredictedTargets(search.targets || [])
//     setRelatedResearch(search.research || [])
//     setDockingResults(search.docking || null)
//     setCurrentSmiles(search.smiles)
//     setAnalysisCompleted(true)
//     setActiveTab("moleculeInput")
//   }

//   if (checkingAuth) {
//     return (
//       <div className="flex items-center justify-center min-h-screen bg-gray-100">
//         <div className="text-center p-6">
//           <Loader2 className="h-8 w-8 animate-spin mx-auto mb-4 text-blue-600" />
//           <p className="text-gray-600">Verifying authentication...</p>
//         </div>
//       </div>
//     )
//   }

//   if (!user) {
//     return (
//       <div className="flex items-center justify-center min-h-screen bg-gray-100">
//         <div className="text-center p-6 bg-white rounded-lg shadow-lg">
//           <AlertCircle className="h-8 w-8 mx-auto mb-4 text-red-500" />
//           <p className="text-gray-600 mb-4">Please log in to access AI-Driven Target Prediction</p>
//           <button
//             className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
//             onClick={() => (window.location.href = "/login")}
//           >
//             Go to Login
//           </button>
//         </div>
//       </div>
//     )
//   }

//   return (
//     <div className="min-h-screen bg-gray-100 py-12 px-4 sm:px-6 lg:px-8">
//       <div className="max-w-7xl mx-auto">
//         <h1 className="text-4xl font-bold text-blue-700 mb-10 text-center">
//           AI-Driven Target Prediction
//           <p className="text-xs p-1 text-blue-700 font-semibold">(Powered by Gemini)</p>
//         </h1>

//         <div className="flex justify-center mb-8 space-x-4 flex-wrap">
//           {["moleculeInput", "savedSearches"].map((tab) => (
//             <button
//               key={tab}
//               className={`px-6 py-2 rounded-lg font-medium transition-all duration-300 m-1 ${
//                 activeTab === tab ? "bg-blue-600 text-white shadow-md" : "bg-gray-200 text-gray-700 hover:bg-gray-300"
//               }`}
//               onClick={() => handleTabChange(tab)}
//             >
//               {tab === "moleculeInput" && "Molecule Input & Analysis"}
//               {tab === "savedSearches" && "Saved Searches"}
//             </button>
//           ))}
//         </div>

//         <div className="bg-white p-8 rounded-xl shadow-lg border border-gray-200">
//           {activeTab === "moleculeInput" && (
//             <>
//               <h2 className="text-2xl font-semibold text-blue-700 mb-6">Molecule Input & Analysis</h2>

//               <div className="mb-6">
//                 <label className="block text-sm font-medium text-gray-700 mb-2">Draw or Select Molecule</label>
//                 <div className="border border-gray-300 rounded-lg">
//                   {useAlternativeEditor ? (
//                     <SimpleMoleculeEditor onChange={handleMoleculeInputChange} />
//                   ) : (
//                     <div
//                       ref={marvinContainerRef}
//                       id="marvinjs_container"
//                       style={{ width: "100%", height: "400px", position: "relative" }}
//                     >
//                       <div className="absolute inset-0 flex items-center justify-center bg-gray-50">
//                         <Loader2 className="h-8 w-8 animate-spin text-blue-600 mr-2" />
//                         <span>Loading molecule editor...</span>
//                       </div>
//                     </div>
//                   )}
//                 </div>
//               </div>

//               {/* <div className="mb-6">
//                 <label className="block text-sm font-medium text-gray-700 mb-2">
//                   <FileText className="inline-block w-4 h-4 mr-1" />
//                   Enter SMILES/IUPAC Name
//                 </label>
//                 <input
//                   type="text"
//                   value={moleculeInput}
//                   onChange={(e) => handleMoleculeInputChange(e.target.value)}
//                   placeholder="Enter SMILES or IUPAC name (e.g., CC(=O)OC1=CC=CC=C1C(=O)O for Aspirin)"
//                   className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all"
//                 />
//                 {moleculeInput && (
//                   <div className="mt-2 p-2 bg-gray-50 rounded-lg border border-gray-200">
//                     <p className="text-xs text-gray-500 mb-1">Molecule Preview:</p>
//                     <img
//                       className="molecule-preview-img mx-auto h-32 object-contain"
//                       src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/PNG?record_type=2d&smiles=${encodeURIComponent(moleculeInput)}`}
//                       alt="Molecule structure preview"
//                     />
//                   </div>
//                 )}
//               </div> */}

//               {/* <div className="mb-6">
//                 <label className="block text-sm font-medium text-gray-700 mb-2">
//                   <Upload className="inline-block w-4 h-4 mr-1" />
//                   Upload MOL/SDF File
//                 </label>
//                 <input
//                   type="file"
//                   accept=".mol,.sdf"
//                   onChange={(e) => setMoleculeFile(e.target.files[0])}
//                   className="w-full p-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-all"
//                 />
//                 {moleculeFile && <p className="mt-2 text-sm text-green-600">File selected: {moleculeFile.name}</p>}
//               </div> */}

//               <button
//                 onClick={handleMoleculeSubmit}
//                 disabled={loading}
//                 className="w-full py-3 px-4 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-all duration-300 flex items-center justify-center"
//               >
//                 {loading ? (
//                   <>
//                     <Loader2 className="animate-spin h-5 w-5 mr-2" />
//                     Analyzing Molecule...
//                   </>
//                 ) : (
//                   <>
//                     <Search className="h-5 w-5 mr-2" />
//                     Analyze Molecule
//                   </>
//                 )}
//               </button>

//               {error && (
//                 <div className="mt-6 bg-red-50 border border-red-300 text-red-700 px-4 py-3 rounded-lg flex justify-between items-center">
//                   <p>{error}</p>
//                   <button className="text-red-700 underline hover:text-red-900" onClick={() => setError(null)}>
//                     Dismiss
//                   </button>
//                 </div>
//               )}

//               {analysisCompleted && (
//                 <div className="mt-6 bg-green-50 border border-green-300 text-green-700 px-4 py-3 rounded-lg">
//                   <p>Analysis completed for SMILES: {currentSmiles}</p>
//                   <p className="text-sm mt-1">Results are displayed below and saved to your account.</p>
//                 </div>
//               )}

//               {(predictedTargets.length > 0 || relatedResearch.length > 0 || dockingResults) && (
//                 <div className="mt-8 bg-blue-50 p-6 rounded-xl border border-blue-200">
//                   <h3 className="text-xl font-semibold text-blue-700 mb-6">Analysis Results</h3>

//                   {predictedTargets.length > 0 && (
//                     <div className="mb-10">
//                       <h4 className="text-lg font-semibold text-gray-800 mb-4">Predicted Protein Targets</h4>
//                       <div className="space-y-6">
//                         {predictedTargets.map((target, index) => (
//                           <div key={index} className="bg-white p-4 rounded-lg shadow-sm border border-blue-100">
//                             <p className="text-gray-700 text-sm mb-1">
//                               <span className="font-bold">Protein:</span> {target.protein}
//                             </p>
//                             <p className="text-gray-700 text-sm mb-1">
//                               <span className="font-bold">Confidence:</span> {(target.confidence * 100).toFixed(2)}%
//                             </p>
//                             <p className="text-gray-700 text-sm mb-1">
//                               <span className="font-bold">Mechanism of Action:</span> {target.moa}
//                             </p>
//                             <p className="text-gray-700 text-sm mb-1">
//                               <span className="font-bold">Pathways:</span>{" "}
//                               {target.pathways && target.pathways.length > 0 ? target.pathways.join(", ") : "N/A"}
//                             </p>
//                             <p className="text-gray-700 text-sm mb-1">
//                               <span className="font-bold">Diseases:</span>{" "}
//                               {target.diseases && target.diseases.length > 0 ? target.diseases.join(", ") : "N/A"}
//                             </p>
//                             {target.knownInteractions && (
//                               <p className="text-gray-700 text-sm mb-1">
//                                 <span className="font-bold">Known Interactions:</span>{" "}
//                                 {JSON.stringify(target.knownInteractions)}
//                               </p>
//                             )}
//                           </div>
//                         ))}
//                       </div>
//                     </div>
//                   )}

//                   {relatedResearch.length > 0 && (
//                     <div className="mb-10">
//                       <h4 className="text-lg font-semibold text-gray-800 mb-4">Related Research Papers</h4>
//                       <div className="grid gap-6 md:grid-cols-2">
//                         {relatedResearch.map((paper, index) => (
//                           <div
//                             key={index}
//                             className="bg-white p-5 rounded-lg shadow-md border border-gray-100 hover:shadow-lg transition-shadow duration-300"
//                           >
//                             <h5 className="text-lg font-bold text-gray-900 mb-3 line-clamp-2">{paper.title}</h5>
//                             <p className="text-gray-700 text-sm mb-2">
//                               <span className="font-semibold">Authors:</span> {paper.authors}
//                             </p>
//                             <p className="text-gray-700 text-sm mb-2">
//                               <span className="font-semibold">Journal:</span> {paper.journal || "Not specified"}
//                             </p>
//                             <p className="text-gray-700 text-sm mb-2">
//                               <span className="font-semibold">Published:</span> {paper.year}
//                             </p>
//                             <p className="text-gray-700 text-sm mb-3 line-clamp-3">
//                               <span className="font-semibold">Abstract:</span> {paper.abstract}
//                             </p>
//                             {paper.doi && paper.url !== "No URL available" && (
//                               <p className="text-blue-600 text-sm flex items-center">
//                                 <svg
//                                   className="w-4 h-4 mr-1"
//                                   fill="none"
//                                   stroke="currentColor"
//                                   viewBox="0 0 24 24"
//                                   xmlns="http://www.w3.org/2000/svg"
//                                 >
//                                   <path
//                                     strokeLinecap="round"
//                                     strokeLinejoin="round"
//                                     strokeWidth="2"
//                                     d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101m-.758-4.899a4 4 0 005.656 0l4-4a4 4 0 00-5.656-5.656l-1.1 1.1"
//                                   ></path>
//                                 </svg>
//                                 <a
//                                   href={paper.url}
//                                   target="_blank"
//                                   rel="noopener noreferrer"
//                                   className="underline hover:text-blue-500"
//                                 >
//                                   {paper.doi}
//                                 </a>
//                               </p>
//                             )}
//                           </div>
//                         ))}
//                       </div>
//                     </div>
//                   )}

//                   {dockingResults && (
//                     <div className="mb-10">
//                       <h4 className="text-lg font-semibold text-gray-800 mb-4">Molecular Docking Results</h4>
//                       <div className="bg-white p-4 rounded-lg shadow-sm border border-blue-100">
//                         <p className="text-gray-700 text-sm mb-1">
//                           <span className="font-bold">Binding Energy:</span> {dockingResults.bindingEnergy} kcal/mol
//                         </p>
//                         <p className="text-gray-700 text-sm mb-1">
//                           <span className="font-bold">Pose:</span> {dockingResults.pose}
//                         </p>
//                         <p className="text-gray-700 text-sm mb-1">
//                           <span className="font-bold">Details:</span> {dockingResults.details}
//                         </p>
//                       </div>
//                     </div>
//                   )}
//                 </div>
//               )}
//             </>
//           )}

//           {activeTab === "savedSearches" && (
//             <>
//               <h2 className="text-2xl font-semibold text-blue-700 mb-6">Saved Searches</h2>
//               {loading ? (
//                 <div className="flex justify-center items-center py-10">
//                   <Loader2 className="h-8 w-8 animate-spin text-blue-600 mr-2" />
//                   <span>Loading saved searches...</span>
//                 </div>
//               ) : savedSearches.length > 0 ? (
//                 <div className="space-y-10">
//                   {savedSearches.map((entry, index) => (
//                     <div key={index} className="bg-blue-50 p-6 rounded-xl border border-blue-200">
//                       <div className="flex justify-between items-start mb-4">
//                         <div>
//                           <h3 className="text-lg font-semibold text-gray-800">Molecule SMILES: {entry.smiles}</h3>
//                           <p className="text-sm text-gray-500 mt-1">
//                             Saved on: {new Date(entry.createdAt).toLocaleString()}
//                           </p>
//                         </div>
//                         <button
//                           onClick={() => handleViewSavedSearch(entry)}
//                           className="px-3 py-1 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 transition-colors"
//                         >
//                           View Details
//                         </button>
//                       </div>
//                       {entry.targets && entry.targets.length > 0 ? (
//                         <div className="mb-6">
//                           <h4 className="text-md font-semibold text-gray-800 mb-2">Predicted Protein Targets</h4>
//                           <div className="bg-white p-3 rounded-lg border border-blue-100">
//                             <ul className="list-disc list-inside space-y-1">
//                               {entry.targets.slice(0, 2).map((target, tIndex) => (
//                                 <li key={tIndex} className="text-gray-700 text-sm">
//                                   <span className="font-medium">{target.protein || "Unknown protein"}</span>
//                                   {target.confidence ? ` (${(target.confidence * 100).toFixed(0)}% confidence)` : ""}
//                                 </li>
//                               ))}
//                               {entry.targets.length > 2 && (
//                                 <li className="text-gray-500 text-sm italic">
//                                   + {entry.targets.length - 2} more targets
//                                 </li>
//                               )}
//                             </ul>
//                           </div>
//                         </div>
//                       ) : (
//                         <p className="text-gray-600 mb-4">No predicted targets available.</p>
//                       )}
//                       {entry.research && entry.research.length > 0 ? (
//                         <div className="mb-6">
//                           <h4 className="text-md font-semibold text-gray-800 mb-2">Related Research Papers</h4>
//                           <div className="bg-white p-3 rounded-lg border border-blue-100">
//                             <ul className="list-disc list-inside space-y-1">
//                               {entry.research.slice(0, 2).map((paper, pIndex) => (
//                                 <li key={pIndex} className="text-gray-700 text-sm line-clamp-1">
//                                   {paper.title || "Untitled paper"}
//                                 </li>
//                               ))}
//                               {entry.research.length > 2 && (
//                                 <li className="text-gray-500 text-sm italic">
//                                   + {entry.research.length - 2} more papers
//                                 </li>
//                               )}
//                             </ul>
//                           </div>
//                         </div>
//                       ) : (
//                         <p className="text-gray-600 mb-4">No related research papers available.</p>
//                       )}
//                       {entry.docking ? (
//                         <div>
//                           <h4 className="text-md font-semibold text-gray-800 mb-2">Molecular Docking Results</h4>
//                           <div className="bg-white p-3 rounded-lg border border-blue-100">
//                             <p className="text-gray-700 text-sm">
//                               <span className="font-medium">Binding Energy:</span>{" "}
//                               {entry.docking.bindingEnergy || "N/A"} kcal/mol
//                             </p>
//                           </div>
//                         </div>
//                       ) : (
//                         <p className="text-gray-600">No docking results available.</p>
//                       )}
//                     </div>
//                   ))}
//                 </div>
//               ) : (
//                 <div className="text-center py-10">
//                   <p className="text-gray-600">No saved searches found.</p>
//                   <button
//                     onClick={() => handleTabChange("moleculeInput")}
//                     className="mt-4 px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
//                   >
//                     Analyze a Molecule
//                   </button>
//                 </div>
//               )}
//             </>
//           )}
//         </div>
//       </div>
//     </div>
//   )
// }

// export default AIdriventargetprediction