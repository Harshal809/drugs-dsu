import { useState, useEffect, useRef } from "react"
import { useNavigate, useLocation } from "react-router-dom"
import {
  Bot,
  Microscope,
  Atom,
  Brain,
  X,
  Mic,
  MicOff,
  Dna,
  Pill,
  ChevronRight,
  ChevronLeft,
  CheckCircle,
  HelpCircle,
  ArrowLeft,
  Home,
  FlaskConical,
  DollarSign,
  FileText,
  AlertTriangle,
  FileBox,
  Newspaper,
  BrainCog,
  Sparkles,
  Zap,
} from "lucide-react"
import { useAuthStore } from "../../Store/auth.store.js"

export default function Chatbot() {
  const [isPanelOpen, setIsPanelOpen] = useState(false)
  const [response, setResponse] = useState("")
  const [isLoading, setIsLoading] = useState(false)
  const [isSpeaking, setIsSpeaking] = useState(false)
  const [isListening, setIsListening] = useState(false)
  const [isRecognizing, setIsRecognizing] = useState(false)
  const [transcript, setTranscript] = useState("")
  const [conversationHistory, setConversationHistory] = useState([])
  const [tooltipVisible, setTooltipVisible] = useState("")
  const [mode, setMode] = useState(null) // null, "doubt", or "beginner"
  const [subMode, setSubMode] = useState(null) // null, "waiting_for_selection", "waiting_for_tool_selection", "waiting_for_step_selection", "explaining_tab", "waiting_for_question", "resolving_doubt", "tour_complete", "process_ended"
  const [targetRoute, setTargetRoute] = useState(null)
  const [currentToolIndex, setCurrentToolIndex] = useState(0)
  const [selectedOption, setSelectedOption] = useState(null) // "dashboard" or "drugDiscovery"
  const [selectedTool, setSelectedTool] = useState(null) // Track manually selected tool or step
  const [completedSteps, setCompletedSteps] = useState(
    JSON.parse(localStorage.getItem("drugDiscoveryProgress") || "[]"),
  ) // Persist progress
  const recognitionRef = useRef(null)
  const isCleaningUp = useRef(false)
  const recognitionState = useRef("idle")
  const conversationRef = useRef(null)
  const lastProcessedRoute = useRef(null)
  const recognitionLock = useRef(false)
  const { user } = useAuthStore()
  const navigate = useNavigate()
  const location = useLocation()
  const [isMinimized, setIsMinimized] = useState(false)

  // Updated dashboard routes with detailed descriptions
 const dashboardRoutes = [
  {
    path: "/dashboard/protein-structure-mutation",
    name: "New Drug Discovery",
    icon: "FlaskConical",
    description: "Explore and generate new drug compounds",
    details:
      "Leverage AI and mutation-based models to explore novel drug compounds and simulate molecular interactions for drug discovery.",
  },
  {
    path: "/dashboard/ai-naming",
    name: "AI Naming Suggestions",
    icon: "BrainCog",
    description: "Generate creative and compliant drug names",
    details:
      "Use AI to suggest effective and regulatory-friendly pharmaceutical names that reflect the drug’s purpose and chemistry.",
  },
  {
    path: "/dashboard/cost-estimation",
    name: "Cost Estimation",
    icon: "DollarSign",
    description: "Estimate drug synthesis and production costs",
    details:
      "Predict cost parameters for synthesis, scale-up, and manufacturing to streamline budgeting and financial planning.",
  },
  {
    path: "/dashboard/protein-structure",
    name: "Protein Structure Generation",
    icon: "Dna",
    description: "Generate 3D protein structures",
    details:
      "Design protein models, analyze structural binding sites, and visualize protein folding to support structure-based drug design.",
  },
  {
    path: "/dashboard/getalphafoldstrcture",
    name: "AlphaFold 3D Predictions",
    icon: "Atom",
    description: "Access accurate protein structure predictions",
    details:
      "Use AlphaFold’s powerful predictions to study protein folding, binding pockets, and conformations for target validation.",
  },
  {
    path: "/dashboard/ai-research-paper-generator",
    name: "Research Paper Generation",
    icon: "FileText",
    description: "AI-based scientific writing",
    details:
      "Generate documentation, experimental reports, and research articles using AI based on your drug discovery inputs.",
  },
  {
    path: "/dashboard/toxicityPrediction",
    name: "Toxicity & Side Effects",
    icon: "AlertTriangle",
    description: "Predict toxicology and adverse reactions",
    details:
      "Analyze compound structures to identify possible toxic effects and interactions early in the development cycle.",
  },
  {
    path: "/dashboard/voice-text-notes",
    name: "Audio Note Capture",
    icon: "Mic",
    description: "Convert voice input to structured notes",
    details:
      "Capture audio inputs, convert them to searchable text, and organize your lab notes using voice-to-text transcription.",
  },
  {
    path: "/dashboard/summary",
    name: "Summarization",
    icon: "FileBox",
    description: "Summarize research insights",
    details:
      "Generate comprehensive summaries from your research journey including results, analysis, and proposed next steps.",
  },
  {
    path: "/dashboard/live-news",
    name: "NewsFeed",
    icon: "Newspaper",
    description: "Get live pharmaceutical updates",
    details:
      "Stay informed on the latest trends, FDA updates, and global research news related to biotechnology and pharmaceuticals.",
  },
];


  const drugDiscoverySteps = [
  {
    path: "/dashboard/protein-structure-mutation",
    name: "Step 1: New Drug Discovery",
    description: "Use AI and mutation modeling for new compounds",
    icon: "FlaskConical",
    details:
      "Start your discovery process by simulating compound mutations, predicting interactions, and exploring drug-like properties.",
  },
  {
    path: "/dashboard/ai-naming",
    name: "Step 2: AI Naming Suggestions",
    description: "Generate suitable names for your drug candidates",
    icon: "BrainCog",
    details:
      "Use AI to produce brandable, memorable names that reflect mechanism of action and comply with naming regulations.",
  },
  {
    path: "/dashboard/cost-estimation",
    name: "Step 3: Cost Estimation",
    description: "Financial insights for drug production",
    icon: "DollarSign",
    details:
      "Estimate costs including synthesis, formulation, trials, and marketing to ensure financial feasibility.",
  },
  {
    path: "/dashboard/protein-structure",
    name: "Step 4: Protein Structure Generation",
    description: "3D modeling of protein targets",
    icon: "Dna",
    details:
      "Understand your target protein by visualizing its 3D structure and preparing it for molecular docking or AI analysis.",
  },
  {
    path: "/dashboard/getalphafoldstrcture",
    name: "Step 5: AlphaFold 3D Predictions",
    description: "Leverage AlphaFold’s database",
    icon: "Atom",
    details:
      "Access AlphaFold’s structural predictions for thousands of proteins, aiding in target identification and design.",
  },
  {
    path: "/dashboard/ai-research-paper-generator",
    name: "Step 6: Research Paper Generation",
    description: "Compile reports and papers",
    icon: "FileText",
    details:
      "Automatically generate structured scientific documents that capture your experimental and AI-generated findings.",
  },
  {
    path: "/dashboard/toxicityPrediction",
    name: "Step 7: Toxicity & Side Effects",
    description: "Analyze safety and side effects",
    icon: "AlertTriangle",
    details:
      "Predict toxicity levels and adverse reactions based on AI models trained on large molecular and clinical datasets.",
  },
  {
    path: "/dashboard/voice-text-notes",
    name: "Step 8: Audio Note Capture",
    description: "Voice-based note-taking",
    icon: "Mic",
    details:
      "Use your voice to create structured, searchable notes—perfect for recording lab ideas or experiment outcomes.",
  },
  {
    path: "/dashboard/summary",
    name: "Step 9: Summarization",
    description: "Project-wide AI summary",
    icon: "FileBox",
    details:
      "Wrap up with a detailed AI-driven overview including methodology, results, costs, risks, and future directions.",
  },
];


  // Enhanced suggested questions per tool/step
const suggestedQuestions = {
  "/dashboard/protein-structure-mutation": [
    "How does mutation-based drug discovery work?",
    "What inputs are needed for structure mutation?",
    "Can I compare mutated compounds?",
  ],
  "/dashboard/ai-naming": [
    "How does AI generate drug names?",
    "Can I customize the naming criteria?",
    "Are these names compliant with international standards?",
  ],
  "/dashboard/cost-estimation": [
    "What cost components are analyzed?",
    "How accurate are the cost models?",
    "Can I adjust for custom manufacturing setups?",
  ],
  "/dashboard/protein-structure": [
    "How do I generate a 3D protein structure?",
    "What visualization tools are supported?",
    "Can I simulate binding with other compounds?",
  ],
  "/dashboard/getalphafoldstrcture": [
    "How do I access AlphaFold predictions?",
    "What is the confidence score for predictions?",
    "Can I visualize the folding pathways?",
  ],
  "/dashboard/ai-research-paper-generator": [
    "What structure does the research paper follow?",
    "Can I upload datasets to include?",
    "Is citation management automated?",
  ],
  "/dashboard/toxicityPrediction": [
    "How is toxicity assessed?",
    "Which side effects are predicted?",
    "Can I view molecular reasons for side effects?",
  ],
  "/dashboard/voice-text-notes": [
    "Can I record multiple audio notes?",
    "How accurate is the transcription?",
    "Is language translation available?",
  ],
  "/dashboard/summary": [
    "What is included in the AI summary?",
    "Can I export as PDF or DOCX?",
    "Does it include graphs and charts?",
  ],
  "/dashboard/live-news": [
    "Which sources are considered reliable?",
    "How frequently is data updated?",
    "Can I subscribe to specific topics?",
  ],
};


  // Map icon names to components
  const iconMap = {
    Home: <Home className="h-4 w-4" />,
    FlaskConical: <FlaskConical className="h-4 w-4" />,
    DollarSign: <DollarSign className="h-4 w-4" />,
    Dna: <Dna className="h-4 w-4" />,
    Atom: <Atom className="h-4 w-4" />,
    BrainCog: <BrainCog className="h-4 w-4" />,
    FileText: <FileText className="h-4 w-4" />,
    AlertTriangle: <AlertTriangle className="h-4 w-4" />,
    Mic: <Mic className="h-4 w-4" />,
    FileBox: <FileBox className="h-4 w-4" />,
    Newspaper: <Newspaper className="h-4 w-4" />,
    Brain: <Brain className="h-4 w-4" />,
    Microscope: <Microscope className="h-4 w-4" />,
    Bot: <Bot className="h-4 w-4" />,
    Pill: <Pill className="h-4 w-4" />,
  }

  // Normalize path for route comparison
  const normalizePath = (path) => {
    if (!path || typeof path !== "string") {
      console.warn(`Invalid path: ${path}, returning default '/dashboard'`)
      return "/dashboard"
    }
    return path.replace(/\/+$/, "").toLowerCase()
  }

  // Initialize Web Speech API for recognition (doubt mode only)
  const createRecognition = () => {
    const SpeechRecognition = window.SpeechRecognition || window.webkitSpeechRecognition
    if (!SpeechRecognition) {
      console.error("SpeechRecognition API not supported.")
      return null
    }
    const recognition = new SpeechRecognition()
    recognition.lang = "en-US"
    recognition.interimResults = false
    recognition.maxAlternatives = 1
    recognition.continuous = true
    return recognition
  }

  // Configure recognition handlers (doubt mode only)
  const setupRecognition = (recognition) => {
    recognition.onstart = () => {
      console.log("Recognition started")
      setIsRecognizing(true)
      recognitionState.current = "active"
    }

    recognition.onresult = async (event) => {
      const speechResult = event.results[event.results.length - 1][0].transcript.trim()
      console.log(`Captured transcript: "${speechResult}"`)
      setTranscript(speechResult)
      setConversationHistory((prev) => [...prev, { type: "user", text: speechResult }])

      if (speechResult.toLowerCase().includes("stop") || speechResult.toLowerCase().includes("thank you")) {
        console.log("Stop command detected, stopping conversation")
        await stopConversation()
        return
      }

      await fetchResponse(speechResult)
    }

    recognition.onend = async () => {
      console.log(`Recognition ended - isListening: ${isListening}, isSpeaking: ${isSpeaking}`)
      setIsRecognizing(false)
      recognitionState.current = "idle"
      if (!isCleaningUp.current && !isSpeaking && isListening && mode === "doubt") {
        console.log("Restarting recognition")
        await startRecognition()
      }
    }

    recognition.onerror = async (event) => {
      console.error(`Speech recognition error: ${event.error}`)
      setIsRecognizing(false)
      recognitionState.current = "idle"
      if (event.error === "no-speech" && !isCleaningUp.current && isListening && mode === "doubt") {
        const message = "I didn't hear anything. Please try again."
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: message }])
        await speakResponse(message)
        await startRecognition()
      } else {
        recognitionRef.current = null
        if (!isCleaningUp.current && isListening && mode === "doubt") {
          await startRecognition()
        }
      }
    }
  }

  // Start speech recognition (doubt mode only)
  const startRecognition = async (retryCount = 0) => {
    if (recognitionLock.current || isCleaningUp.current || isSpeaking) {
      console.log("Recognition blocked")
      if (isSpeaking && retryCount < 3) {
        await new Promise((resolve) => setTimeout(resolve, 500))
        await startRecognition(retryCount + 1)
      }
      return
    }

    recognitionLock.current = true

    try {
      if (recognitionRef.current && (isRecognizing || recognitionState.current === "active")) {
        recognitionRef.current.stop()
        await new Promise((resolve) => setTimeout(resolve, 200))
      }

      if (!recognitionRef.current) {
        recognitionRef.current = createRecognition()
        if (!recognitionRef.current) {
          recognitionLock.current = false
          const errorMessage = "Speech recognition is not supported on this device."
          setConversationHistory((prev) => [...prev, { type: "jarvis", text: errorMessage }])
          await speakResponse(errorMessage)
          setIsListening(false)
          return
        }
        setupRecognition(recognitionRef.current)
      }

      recognitionRef.current.start()
      setIsListening(true)
      setIsRecognizing(true)
    } catch (error) {
      console.error(`Error starting recognition: ${error}`)
      setIsRecognizing(false)
      recognitionState.current = "idle"

      if (retryCount < 2) {
        recognitionRef.current = null
        await new Promise((resolve) => setTimeout(resolve, 200))
        await startRecognition(retryCount + 1)
      } else {
        const message = "Sorry, I couldn't start speech recognition. Please try again."
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: message }])
        await speakResponse(message)
        setIsListening(false)
      }
    } finally {
      recognitionLock.current = false
    }
  }

  // Speak response using Web Speech API
  const speakResponse = async (text) => {
    return new Promise((resolve) => {
      if (!text || typeof text !== "string") {
        console.error("Invalid text for speech synthesis:", text)
        setIsSpeaking(false)
        resolve()
        return
      }

      if (recognitionRef.current && isRecognizing) {
        try {
          recognitionRef.current.stop()
          setIsRecognizing(false)
          recognitionState.current = "idle"
        } catch (error) {
          console.error(`Error stopping recognition: ${error}`)
        }
      }

      window.speechSynthesis.cancel()

      const utterance = new SpeechSynthesisUtterance(text)
      utterance.lang = "en-US"
      utterance.rate = 0.9
      utterance.onstart = () => {
        setIsSpeaking(true)
      }
      utterance.onend = async () => {
        setIsSpeaking(false)
        if (mode === "doubt" && isListening && !isCleaningUp.current) {
          await startRecognition()
        }
        resolve()
      }
      utterance.onerror = (event) => {
        console.error(`Speech synthesis error: ${event.error}`)
        setIsSpeaking(false)
        resolve()
      }
      window.speechSynthesis.speak(utterance)
    })
  }

  // Start conversation
  const startConversation = async (selectedMode) => {
    isCleaningUp.current = false
    setIsPanelOpen(true)
    setIsListening(false)
    setIsRecognizing(false)
    setTranscript("")
    setResponse("")
    setConversationHistory([])
    setMode(selectedMode)
    setSubMode(null)
    setTargetRoute(null)
    setCurrentToolIndex(0)
    setSelectedTool(null)
    setSelectedOption(null)
    setCompletedSteps(JSON.parse(localStorage.getItem("drugDiscoveryProgress") || "[]"))

    if (recognitionRef.current) {
      try {
        recognitionRef.current.stop()
      } catch (error) {
        console.error(`Error stopping recognition: ${error}`)
      }
      recognitionRef.current = null
    }

    let welcomeMessage = ""
    if (selectedMode === "doubt") {
      welcomeMessage =
        "I am Jarvis, here to assist with your medical and pharmaceutical questions. Please speak your question."
      setConversationHistory([{ type: "jarvis", text: welcomeMessage }])
      await speakResponse(welcomeMessage)
      setIsListening(true)
      await startRecognition()
    } else if (selectedMode === "beginner") {
      welcomeMessage =
        "Welcome to your drug discovery journey! I'm Jarvis, and I'll guide you through the dashboard tools or the drug discovery process. Please select an option."
      setConversationHistory([{ type: "jarvis", text: welcomeMessage }])
      await speakResponse(welcomeMessage)
      setSubMode("waiting_for_selection")
    }
  }

  // Stop conversation
  const stopConversation = async () => {
    // Stop listening and recognition first
    setIsListening(false)
    setIsRecognizing(false)

    if (recognitionRef.current) {
      try {
        recognitionRef.current.stop()
      } catch (error) {
        console.error(`Error stopping recognition: ${error}`)
      }
      recognitionRef.current = null
    }

    // Add goodbye message to conversation history
    const goodbyeMessage =
      "Thank you for using Jarvis! It was a pleasure assisting you with your drug discovery journey. Have a great day!"
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: goodbyeMessage }])

    // Speak the goodbye message and wait for it to complete
    try {
      await speakResponse(goodbyeMessage)
    } catch (error) {
      console.error("Error speaking goodbye message:", error)
    }

    // Only after speech is complete, start cleanup
    isCleaningUp.current = true

    // Cancel any remaining speech
    window.speechSynthesis.cancel()

    // Reset states
    setMode(null)
    setSubMode(null)
    setTargetRoute(null)
    setCurrentToolIndex(0)
    setSelectedTool(null)
    setSelectedOption(null)
    setCompletedSteps([])
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]))

    // Close the panel after everything is complete
    setIsPanelOpen(false)
    isCleaningUp.current = false
  }

  // Fetch response from Gemini API
  const fetchResponse = async (query) => {
    setIsLoading(true)
    const contextInfo = selectedTool ? `Current tool context: ${selectedTool.name}. ` : ""
    const problemStatement = `You are Jarvis, an AI assistant specialized in medical and pharmaceutical domains. ${contextInfo}Answer the user's query comprehensively and accurately. Keep your response concise but informative. Query: ${query}`

    try {
      const response = await fetch(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=AIzaSyDyujm50dHMYvn1V50dDDqcAhgUqCOuUGU`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            contents: [{ parts: [{ text: problemStatement }] }],
          }),
        },
      )

      if (!response.ok) {
        throw new Error(`HTTP error! Status: ${response.status}`)
      }

      const data = await response.json()
      const text = data.candidates[0].content.parts[0].text
      setResponse(text)
      setConversationHistory((prev) => [...prev, { type: "jarvis", text }])
      await speakResponse(text)

      // Store the Q&A in local storage
      const jarvisEntry = {
        title: `Query: ${query.slice(0, 50)}`,
        questions: query,
        answers: text,
        toolContext: selectedTool?.name || "General",
        userId: user?._id || "66f172e9c7e124f5c4b3c2d1",
        createdAt: new Date().toISOString(),
      }
      const storedEntries = JSON.parse(localStorage.getItem("jarvisEntries") || "[]")
      storedEntries.push(jarvisEntry)
      localStorage.setItem("jarvisEntries", JSON.stringify(storedEntries))

      if (mode === "beginner" && (subMode === "resolving_doubt" || subMode === "waiting_for_question")) {
        setSubMode("waiting_for_question")
        const followUpPrompt = "Would you like to ask another question or select a different tool/step?"
        setConversationHistory((prev) => [...prev, { type: "jarvis", text: followUpPrompt }])
        await speakResponse(followUpPrompt)
      }
    } catch (error) {
      console.error(`Error fetching Gemini API: ${error}`)
      const errorMessage =
        "Sorry, I encountered an error while processing your question. Would you like to try again or select another tool/step?"
      setResponse(errorMessage)
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: errorMessage }])
      await speakResponse(errorMessage)

      if (mode === "beginner") {
        setSubMode("waiting_for_question")
      }
    } finally {
      setIsLoading(false)
    }
  }

  // Handle dashboard tools selection
  const handleDashboardSelection = async () => {
    setSelectedOption("dashboard")
    setSubMode("waiting_for_tool_selection")
    const prompt = "You chose to explore dashboard tools. Please select a tool to learn about."
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
    await speakResponse(prompt)
  }

  // Handle drug discovery selection
  const handleDrugDiscoverySelection = async () => {
    setSelectedOption("drugDiscovery")
    setSubMode("waiting_for_step_selection")
    const prompt = "You chose the drug discovery process. Please select the current step to begin."
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
    await speakResponse(prompt)
  }

  // Handle manual tool selection with delay
  const handleToolSelection = async (tool) => {
    setSelectedTool({ ...tool, icon: tool.icon })
    const prompt = `You selected ${tool.name}. Now navigating to that page for more information.`
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
    await speakResponse(prompt)
    await new Promise((resolve) => setTimeout(resolve, 1000)) // 1-second delay
    navigate(tool.path, { replace: true, state: { fromJarvis: true } })
    setTargetRoute(tool.path)
    setSubMode("explaining_tab")
    setCurrentToolIndex(dashboardRoutes.findIndex((t) => t.path === tool.path))
    lastProcessedRoute.current = tool.path
    await speakFieldsForRoute(tool.path)
  }

  // Handle manual step selection with delay
  const handleStepSelection = async (step) => {
    setSelectedTool({ ...step, icon: step.icon })
    const prompt = `You selected ${step.name}. Now navigating to that page for more information.`
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
    await speakResponse(prompt)
    await new Promise((resolve) => setTimeout(resolve, 1000)) // 1-second delay
    navigate(step.path, { replace: true, state: { fromJarvis: true } })
    setTargetRoute(step.path)
    setSubMode("explaining_tab")
    setCurrentToolIndex(drugDiscoverySteps.findIndex((s) => s.path === step.path))
    setCompletedSteps((prev) => {
      const newSteps = [...prev, step.path]
      localStorage.setItem("drugDiscoveryProgress", JSON.stringify(newSteps))
      return newSteps
    })
    await speakFieldsForRoute(step.path)
  }

  // Proceed to next step in the process
  const proceedToNextStep = async () => {
    setCompletedSteps((prev) => {
      const newSteps = [...prev, drugDiscoverySteps[currentToolIndex].path]
      localStorage.setItem("drugDiscoveryProgress", JSON.stringify(newSteps))
      return newSteps
    })
    const completionPrompt = `Step ${currentToolIndex + 1} completed!`
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: completionPrompt }])
    await speakResponse(completionPrompt)

    if (currentToolIndex < drugDiscoverySteps.length - 1) {
      const nextStep = drugDiscoverySteps[currentToolIndex + 1]
      const prompt = `Moving to ${nextStep.name}.`
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
      await speakResponse(prompt)
      await new Promise((resolve) => setTimeout(resolve, 1000)) // 1-second delay
      navigate(nextStep.path, { replace: true, state: { fromJarvis: true } })
      setTargetRoute(nextStep.path)
      setSubMode("explaining_tab")
      setCurrentToolIndex(currentToolIndex + 1)
      lastProcessedRoute.current = nextStep.path
      await speakFieldsForRoute(nextStep.path)
    } else {
      const prompt =
        "You've reached the end of the drug discovery process. Click 'End Process' to complete or explore other options."
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
      await speakResponse(prompt)
      setSubMode("tour_complete")
    }
  }

  // End drug discovery process
  const endProcess = async () => {
    const congratsPrompt =
      "Congratulations! You've successfully completed the drug discovery process! Would you like to restart the process or explore dashboard tools?"
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: congratsPrompt }])
    await speakResponse(congratsPrompt)
    setSubMode("process_ended")
    setCurrentToolIndex(0)
    setCompletedSteps([])
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]))
  }

  // Restart the current process
  const restartProcess = async () => {
    setCompletedSteps([])
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]))
    const routes = selectedOption === "drugDiscovery" ? drugDiscoverySteps : dashboardRoutes
    const prompt = `Starting over with ${routes[0].name}.`
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
    await speakResponse(prompt)
    await new Promise((resolve) => setTimeout(resolve, 1000)) // 1-second delay
    navigate(routes[0].path, { replace: true, state: { fromJarvis: true } })
    setTargetRoute(routes[0].path)
    setSubMode("explaining_tab")
    setCurrentToolIndex(0)
    lastProcessedRoute.current = routes[0].path
    await speakFieldsForRoute(routes[0].path)
  }

  // Enhanced speak form fields and tab info based on route
  const speakFieldsForRoute = async (route) => {
    const normalizedRoute = normalizePath(route)
    const foundRoute =
      dashboardRoutes.find((r) => normalizePath(r.path) === normalizedRoute) ||
      drugDiscoverySteps.find((s) => normalizePath(s.path) === normalizedRoute)

    if (!foundRoute) {
      const message = "It seems you're on an unrecognized page. Let's return to the dashboard home."
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: message }])
      await speakResponse(message)
      navigate("/dashboard", { replace: true, state: { fromJarvis: true } })
      lastProcessedRoute.current = "/dashboard"
      return
    }

    const fullMessage = `Now at ${foundRoute.name}. ${foundRoute.details || foundRoute.description}`
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: fullMessage }])
    await speakResponse(fullMessage)

    if (mode === "beginner") {
      setSubMode("waiting_for_question")
      const questionPrompt =
        "Do you have any questions about this tool? If not, you can select another step or proceed to the next step."
      setConversationHistory((prev) => [...prev, { type: "jarvis", text: questionPrompt }])
      await speakResponse(questionPrompt)
    }
  }

  // Get route name for display
  const getRouteName = (route) => {
    const foundRoute =
      dashboardRoutes.find((r) => normalizePath(r.path) === normalizePath(route)) ||
      drugDiscoverySteps.find((s) => normalizePath(s.path) === normalizePath(route))
    return foundRoute ? foundRoute.name : "Unknown"
  }

  // Scroll to bottom of conversation
  useEffect(() => {
    if (conversationRef.current) {
      conversationRef.current.scrollTop = conversationRef.current.scrollHeight
    }
  }, [conversationHistory])

  // Clean up on unmount
  useEffect(() => {
    return () => {
      isCleaningUp.current = true
      window.speechSynthesis.cancel()
      if (recognitionRef.current) {
        try {
          recognitionRef.current.stop()
        } catch (error) {
          console.error(`Error stopping recognition: ${error}`)
        }
      }
      recognitionRef.current = null
    }
  }, [])

  // Save progress to localStorage when completedSteps changes
  useEffect(() => {
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify(completedSteps))
  }, [completedSteps])

  const resetToInitialSelection = async () => {
    const prompt = "Returning to the initial selection screen. How can I help you today?"
    setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
    await speakResponse(prompt)
    setMode(null)
    setSubMode(null)
    setSelectedOption(null)
    setTargetRoute(null)
    setSelectedTool(null)
    setCurrentToolIndex(0)
    setCompletedSteps([])
    localStorage.setItem("drugDiscoveryProgress", JSON.stringify([]))
  }

  return (
    <div className="font-body">
      {/* Trigger Button with enhanced design */}
      <button
        onClick={() => setIsPanelOpen(true)}
        disabled={isPanelOpen}
        onMouseEnter={() => setTooltipVisible("trigger")}
        onMouseLeave={() => setTooltipVisible("")}
        aria-label={isPanelOpen ? "Close Jarvis Assistant" : "Open Jarvis Assistant"}
        className={`fixed bottom-6 right-6 h-16 w-16 rounded-full shadow-2xl transition-all duration-300 transform hover:scale-110 flex items-center justify-center z-50 group animate-float ${isPanelOpen
            ? "bg-gradient-to-br from-accent to-accent-secondary hover:shadow-accent/30"
            : "bg-gradient-to-br from-accent to-accent-secondary hover:shadow-accent/30 animate-pulse-glow"
          }`}
      >
        <Bot className="h-7 w-7 text-primary transition-transform duration-300 group-hover:rotate-12" />
        <div className="absolute inset-0 rounded-full bg-gradient-to-br from-accent/20 to-accent-secondary/20 animate-gradient-shift"></div>
      </button>

      {/* Enhanced Tooltip
      {tooltipVisible === "trigger" && (
        <div className="fixed bottom-24 right-4 px-3 py-2 bg-secondary/95 border border-accent/20 text-text-primary text-sm rounded-xl shadow-xl z-50 backdrop-blur-sm animate-slide-up">
          <div className="flex items-center space-x-2">
            <Sparkles className="h-4 w-4 text-accent" />
            <span className="font-label">{isPanelOpen ? "Close Jarvis" : "Ask Jarvis"}</span>
          </div>
        </div>
      )} */}

      {/* Enhanced Pointing Message */}
      {!isPanelOpen && (
        <div className="fixed bottom-24 right-6 z-40 animate-bounce-subtle">
          <div className="relative bg-gradient-to-br from-secondary to-secondary/80 backdrop-blur-sm rounded-2xl shadow-2xl p-4 mb-2 border border-accent/20 animate-scale-in">
            <div className="absolute right-6 bottom-0 transform translate-y-full w-4 h-4 bg-secondary  border-r border-b border-accent/20"></div>
            <div className="flex items-center space-x-3">
              <div className="bg-gradient-to-br from-accent/30 to-accent-secondary/30 p-2 rounded-full animate-pulse-glow">
                <Bot className="h-6 w-6 text-accent" />
              </div>
              <div>
                <div className="text-base font-heading font-semibold text-text-primary">Jarvis AI Assistant</div>
                <div className="text-sm text-text-secondary font-body">Your pharmaceutical research companion</div>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Enhanced Right Side Panel */}
      {isPanelOpen && (
        <div className="fixed top-16 right-4 h-[92vh] w-[32rem]  bg-gradient-to-br from-secondary to-secondary/95 backdrop-blur-xl shadow-2xl z-50 flex flex-col border border-accent/20 rounded-2xl animate-slide-down">
          {/* Enhanced Header */}
          <div className="flex items-center justify-between p-5 border-b border-accent/20 bg-gradient-to-r from-accent/10 to-accent-secondary/10 rounded-t-2xl">
            <div className="flex items-center space-x-3">
              {mode && (
                <button
                  onClick={resetToInitialSelection}
                  aria-label="Back to initial selection"
                  className="rounded-full h-10 w-10 bg-primary/80 hover:bg-primary text-accent flex items-center justify-center shadow-lg transition-all duration-300 hover:scale-110 border border-accent/20"
                >
                  <ArrowLeft className="h-5 w-5" />
                </button>
              )}
              <div className="bg-gradient-to-br from-accent/30 to-accent-secondary/30 p-2 rounded-full animate-pulse-glow">
                <Dna className="h-7 w-7 text-accent" />
              </div>
              <div>
                <h2 className="text-2xl font-heading font-bold bg-gradient-to-r from-accent to-accent-secondary bg-clip-text text-transparent">
                  JARVIS
                </h2>
                <div className="text-xs text-text-secondary font-code">AI Research Assistant</div>
              </div>
            </div>
            <div className="flex items-center space-x-2">
              <button
                onClick={() => {
                  const helpPrompt =
                    "In beginner mode, you can type questions or follow guided steps. In doubt mode, use voice input to ask questions. Use the buttons to navigate or say 'stop' to end the conversation."
                  setConversationHistory((prev) => [...prev, { type: "jarvis", text: helpPrompt }])
                  speakResponse(helpPrompt)
                }}
                aria-label="Help"
                className="rounded-full h-10 w-10 bg-primary/80 hover:bg-primary text-accent-secondary flex items-center justify-center shadow-lg transition-all duration-300 hover:scale-110 border border-accent/20"
              >
                <HelpCircle className="h-5 w-5" />
              </button>
              <button
                onClick={stopConversation}
                aria-label="Close panel"
                className="rounded-full h-10 w-10 bg-primary/80 hover:bg-primary text-error flex items-center justify-center shadow-lg transition-all duration-300 hover:scale-110 border border-error/20"
              >
                <X className="h-5 w-5" />
              </button>
            </div>
          </div>

          {!isMinimized && (
            <div className="flex flex-col h-full overflow-hidden">
              {/* Enhanced Welcome Screen */}
              {!mode && (
                <div className="flex-1 p-6 flex flex-col items-center justify-center space-y-6 overflow-y-auto bg-gradient-to-br from-primary/20 to-accent/5">
                  <div className="text-center animate-scale-in">
                    <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-6 rounded-2xl w-20 h-20 mx-auto mb-6 flex items-center justify-center  border border-accent/30">
                      <Bot className="h-10 w-10 text-accent" />
                    </div>
                    <h3 className="text-2xl font-heading font-bold text-text-primary mb-3">Welcome to Jarvis</h3>
                  </div>
                    <p className="text-base text-text-secondary font-body leading-relaxed">Choose an option to start your drug discovery journey.</p>
                    <p className="text-base text-text-secondary font-body leading-relaxed">Please be clear & load enough in conversation
                    </p>

                  <div className="w-full space-y-4 max-w-sm">
                    <button
                      onClick={() => startConversation("beginner")}
                      className="w-full py-4 px-6 bg-gradient-to-r from-accent to-accent-secondary text-primary rounded-2xl hover:from-accent/90 hover:to-accent-secondary/90 transition-all duration-300 text-base font-heading font-semibold shadow-2xl transform hover:scale-105 hover:shadow-accent/30 border border-accent/30"
                      aria-label="Start as a beginner"
                    >
                      <div className="flex items-center justify-center space-x-3">
                        <Brain className="h-6 w-6" />
                        <span>Beginner Guide</span>
                        <Sparkles className="h-5 w-5 animate-pulse" />
                      </div>
                    </button>

                    <button
                      onClick={() => startConversation("doubt")}
                      className="w-full py-4 px-6 bg-gradient-to-r from-accent-secondary to-accent text-primary rounded-2xl hover:from-accent-secondary/90 hover:to-accent/90 transition-all duration-300 text-base font-heading font-semibold shadow-2xl transform hover:scale-105 hover:shadow-accent-secondary/30 border border-accent-secondary/30"
                      aria-label="Ask a question"
                    >
                      <div className="flex items-center justify-center space-x-3">
                        <Mic className="h-6 w-6" />
                        <span>Voice Questions</span>
                        <Zap className="h-5 w-5 animate-pulse" />
                      </div>
                    </button>
                  </div>

                  <div className="w-full pt-6 border-t border-accent/20 max-w-sm">
                    <h4 className="text-sm font-label font-medium text-text-secondary mb-4 text-center">Quick Actions</h4>
                    <div className="grid grid-cols-2 gap-3">
                      <button
                        onClick={() => {
                          const helpPrompt =
                            "Here are some example questions you can ask: 'How do I analyze protein structures?', 'What is drug discovery?', 'Explain molecular docking'"
                          setConversationHistory([{ type: "jarvis", text: helpPrompt }])
                          speakResponse(helpPrompt)
                          setMode("doubt")
                        }}
                        className="p-3 bg-primary/60 rounded-xl text-sm text-text-secondary hover:bg-primary/80 hover:text-text-primary transition-all duration-300 font-body border border-accent/10 hover:border-accent/30"
                      >
                        Examples
                      </button>
                      <button
                        onClick={() => {
                          const tutorialPrompt =
                            "I'll guide you through our drug discovery tools step by step. Let's start!"
                          setConversationHistory([{ type: "jarvis", text: tutorialPrompt }])
                          speakResponse(tutorialPrompt)
                          setMode("beginner")
                          setSubMode("waiting_for_selection")
                        }}
                        className="p-3 bg-primary/60 rounded-xl text-sm text-text-secondary hover:bg-primary/80 hover:text-text-primary transition-all duration-300 font-body border border-accent/10 hover:border-accent/30"
                      >
                        Tutorial
                      </button>
                    </div>
                  </div>
                </div>
              )}

              {mode && (
                <div className="flex flex-col h-full overflow-hidden">
                  {/* Enhanced Progress Header for drug discovery */}
                  {mode === "beginner" &&
                    selectedOption === "drugDiscovery" &&
                    subMode !== "waiting_for_selection" &&
                    subMode !== "waiting_for_step_selection" &&
                    subMode !== "process_ended" && (
                      <div className="p-4 bg-gradient-to-r from-accent/10 to-accent-secondary/10 border-b border-accent/20 flex-shrink-0 animate-slide-down">
                        <div className="flex items-center justify-between mb-3">
                          <h3 className="text-base font-heading font-medium text-accent">Drug Discovery Process</h3>
                          <div className="text-sm text-accent-secondary font-code">
                            Step {currentToolIndex + 1} of {drugDiscoverySteps.length}
                          </div>
                        </div>
                        <div className="w-full bg-primary/40 border border-accent/20 rounded-full h-3 overflow-hidden">
                          <div
                            className="bg-gradient-to-r from-accent to-accent-secondary h-3 rounded-full transition-all duration-700 animate-gradient-shift"
                            style={{ width: `${((currentToolIndex + 1) / drugDiscoverySteps.length) * 100}%` }}
                          ></div>
                        </div>
                        <div className="flex items-center space-x-2 overflow-x-auto pb-2 mt-3">
                          {drugDiscoverySteps.map((step, index) => (
                            <div
                              key={index}
                              className={`flex-shrink-0 flex items-center transition-all duration-300 ${completedSteps.includes(step.path)
                                  ? "text-success"
                                  : index === currentToolIndex
                                    ? "text-accent"
                                    : "text-text-secondary"
                                }`}
                            >
                              <div
                                className={`h-8 w-8 rounded-full flex items-center justify-center text-sm font-bold transition-all duration-300 ${completedSteps.includes(step.path)
                                    ? "bg-success/20 border-2 border-success animate-pulse-glow"
                                    : index === currentToolIndex
                                      ? "bg-accent/20 border-2 border-accent animate-pulse-glow"
                                      : "bg-primary/40 border border-text-secondary/30"
                                  }`}
                              >
                                {completedSteps.includes(step.path) ? (
                                  <CheckCircle className="h-4 w-4" />
                                ) : (
                                  <span className="font-code">{index + 1}</span>
                                )}
                              </div>
                              {index < drugDiscoverySteps.length - 1 && (
                                <ChevronRight className="h-4 w-4 mx-2 text-text-secondary" />
                              )}
                            </div>
                          ))}
                        </div>
                      </div>
                    )}

                  {/* Enhanced Current step details */}
                  {mode === "beginner" &&
                    selectedOption === "drugDiscovery" &&
                    subMode !== "waiting_for_selection" &&
                    subMode !== "waiting_for_step_selection" &&
                    subMode !== "process_ended" && (
                      <div className="p-4 bg-gradient-to-r from-secondary/80 to-primary/20 border-b border-accent/20 flex-shrink-0 animate-slide-down">
                        <div className="flex items-start space-x-4">
                          <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-3 rounded-xl border border-accent/30 animate-pulse-glow">
                            {iconMap[drugDiscoverySteps[currentToolIndex]?.icon]}
                          </div>
                          <div className="min-w-0 flex-1">
                            <h4 className="text-lg font-heading font-medium text-text-primary truncate mb-1">
                              {drugDiscoverySteps[currentToolIndex]?.name}
                            </h4>
                            <p className="text-sm text-text-secondary line-clamp-2 font-body leading-relaxed">
                              {drugDiscoverySteps[currentToolIndex]?.description}
                            </p>
                          </div>
                        </div>
                      </div>
                    )}

                  {/* Enhanced Manual navigation buttons */}
                  {mode === "beginner" &&
                    selectedOption === "drugDiscovery" &&
                    subMode !== "waiting_for_selection" &&
                    subMode !== "waiting_for_step_selection" &&
                    subMode !== "process_ended" && (
                      <div className="flex justify-between p-4 bg-secondary/60 border-b border-accent/20 flex-shrink-0">
                        <button
                          disabled={true}
                          className="flex items-center space-x-2 px-4 py-2 rounded-xl text-sm font-heading font-medium text-text-secondary/50 cursor-not-allowed bg-primary/20 border border-text-secondary/20"
                          aria-label="Previous step (disabled)"
                        >
                          <ChevronLeft className="h-4 w-4" />
                          <span>Previous</span>
                        </button>
                        {currentToolIndex === drugDiscoverySteps.length - 1 ? (
                          <button
                            onClick={endProcess}
                            disabled={isSpeaking}
                            className={`flex items-center space-x-2 px-4 py-2 rounded-xl text-sm font-heading font-medium transition-all duration-300 ${isSpeaking
                                ? "text-text-secondary/50 cursor-not-allowed bg-primary/20"
                                : "text-primary bg-gradient-to-r from-success to-accent hover:from-success/90 hover:to-accent/90 transform hover:scale-105 shadow-lg hover:shadow-success/30"
                              } border border-success/30`}
                            aria-label="End drug discovery process"
                          >
                            <span>Complete</span>
                            <Sparkles className="h-4 w-4" />
                          </button>
                        ) : (
                          <button
                            onClick={proceedToNextStep}
                            disabled={isSpeaking}
                            className={`flex items-center space-x-2 px-4 py-2 rounded-xl text-sm font-heading font-medium transition-all duration-300 ${isSpeaking
                                ? "text-text-secondary/50 cursor-not-allowed bg-primary/20"
                                : "text-primary bg-gradient-to-r from-accent to-accent-secondary hover:from-accent/90 hover:to-accent-secondary/90 transform hover:scale-105 shadow-lg hover:shadow-accent/30"
                              } border border-accent/30`}
                            aria-label="Next step"
                          >
                            <span>Next</span>
                            <ChevronRight className="h-4 w-4" />
                          </button>
                        )}
                      </div>
                    )}

                  {/* Enhanced Manual selection for beginner mode */}
                  {mode === "beginner" && subMode === "waiting_for_selection" && (
                    <div className="p-5 flex-shrink-0 overflow-y-auto bg-gradient-to-br from-primary/10 to-accent/5 animate-slide-up">
                      <h3 className="text-lg font-heading font-medium text-text-primary mb-4">Select an option to continue:</h3>
                      <div className="grid grid-cols-1 gap-4">
                        <button
                          onClick={() => handleDashboardSelection()}
                          className="flex items-center space-x-4 p-4 bg-gradient-to-r from-secondary/80 to-secondary/60 rounded-2xl border border-accent/20 hover:border-accent/40 hover:bg-secondary/90 transition-all duration-300 transform hover:scale-105 shadow-lg hover:shadow-accent/20"
                          aria-label="Explore dashboard tools"
                        >
                          <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-3 rounded-xl border border-accent/30 animate-pulse-glow">
                            {iconMap["Home"]}
                          </div>
                          <div className="text-left min-w-0 flex-1">
                            <div className="text-base font-heading font-medium text-text-primary">Dashboard Tools</div>
                            <div className="text-sm text-text-secondary font-body">Explore all available tools and features</div>
                          </div>
                        </button>
                        <button
                          onClick={() => handleDrugDiscoverySelection()}
                          className="flex items-center space-x-4 p-4 bg-gradient-to-r from-secondary/80 to-secondary/60 rounded-2xl border border-accent/20 hover:border-accent/40 hover:bg-secondary/90 transition-all duration-300 transform hover:scale-105 shadow-lg hover:shadow-accent/20"
                        >
                          <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-3 rounded-xl border border-accent/30 animate-pulse-glow">
                            {iconMap["FlaskConical"]}
                          </div>
                          <div className="text-left min-w-0 flex-1">
                            <div className="text-base font-heading font-medium text-text-primary">Drug Discovery Process</div>
                            <div className="text-sm text-text-secondary font-body">Guided step-by-step drug development journey</div>
                          </div>
                        </button>
                      </div>
                      <div className="flex justify-between gap-3 mt-6">
                        <button
                          onClick={() => {
                            stopConversation()
                            setIsPanelOpen(false)
                          }}
                          className="flex-1 px-5 py-3 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-error to-error/80 hover:from-error/90 hover:to-error/70 transition-all duration-300 shadow-lg hover:shadow-error/30 transform hover:scale-105 border border-error/30"
                          aria-label="End conversation"
                        >
                          End Conversation
                        </button>
                      </div>
                    </div>
                  )}

                  {/* Enhanced Tool selection for dashboard */}
                  {mode === "beginner" && subMode === "waiting_for_tool_selection" && (
                    <div className="p-5 bg-gradient-to-br from-primary/10 to-accent/5 border-b border-accent/20 flex-shrink-0 overflow-y-auto max-h-72 animate-slide-up">
                      <h3 className="text-lg font-heading font-medium text-text-primary mb-4">Select a tool to explore:</h3>
                      <div className="grid grid-cols-1 gap-3 max-h-48 overflow-y-auto pb-2 space-y-2">
                        {dashboardRoutes.map((tool, index) => (
                          <button
                            key={index}
                            onClick={() => handleToolSelection(tool)}
                            className="flex items-center space-x-4 p-4 bg-gradient-to-r from-secondary/70 to-secondary/50 rounded-xl border border-accent/20 hover:border-accent/40 hover:bg-secondary/80 transition-all duration-300 transform hover:scale-105 shadow-lg hover:shadow-accent/20 group"
                            aria-label={`Select ${tool.name} tool`}
                            title={tool.description}
                          >
                            <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-2 rounded-lg border border-accent/30 group-hover:animate-pulse-glow transition-all duration-300">
                              {iconMap[tool.icon]}
                            </div>
                            <div className="text-left min-w-0 flex-1">
                              <div className="text-sm font-heading font-medium truncate text-text-primary">{tool.name}</div>
                              <div className="text-xs text-text-secondary line-clamp-1 font-body">{tool.description}</div>
                            </div>
                          </button>
                        ))}
                      </div>
                      <div className="flex justify-between gap-3 mt-4">
                        <button
                          onClick={async () => {
                            const prompt = "Returning to option selection."
                            setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
                            await speakResponse(prompt)
                            setSubMode("waiting_for_selection")
                            setSelectedOption(null)
                          }}
                          className="flex-1 px-4 py-2 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-accent to-accent-secondary hover:from-accent/90 hover:to-accent-secondary/90 transition-all duration-300 shadow-lg hover:shadow-accent/30 transform hover:scale-105 border border-accent/30"
                          aria-label="Go back to option selection"
                          title="Go back to option selection"
                        >
                          Back
                        </button>
                        <button
                          onClick={() => {
                            stopConversation()
                            setIsPanelOpen(false)
                          }}
                          className="flex-1 px-4 py-2 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-error to-error/80 hover:from-error/90 hover:to-error/70 transition-all duration-300 shadow-lg hover:shadow-error/30 transform hover:scale-105 border border-error/30"
                          aria-label="End conversation for tool selection"
                        >
                          End Conversation
                        </button>
                      </div>
                    </div>
                  )}

                  {/* Enhanced Step selection for drug discovery */}
                  {mode === "beginner" && subMode === "waiting_for_step_selection" && (
                    <div className="p-5 bg-gradient-to-br from-primary/10 to-accent/5 border-b border-accent/20 flex-shrink-0 overflow-y-auto max-h-96 animate-slide-up">
                      <h3 className="text-lg font-heading font-medium text-text-primary mb-4">Select a step to explore:</h3>
                      <div className="grid grid-cols-1 gap-3">
                        {drugDiscoverySteps.map((step, index) => (
                          <button
                            key={index}
                            onClick={() => handleStepSelection(step)}
                            disabled={index !== currentToolIndex || completedSteps.includes(step.path)}
                            className={`flex items-center space-x-4 p-4 bg-gradient-to-r rounded-xl border transition-all duration-300 transform ${index !== currentToolIndex || completedSteps.includes(step.path)
                                ? "from-primary/30 to-primary/20 border-text-secondary/20 opacity-50 cursor-not-allowed"
                                : "from-secondary/70 to-secondary/50 border-accent/20 hover:border-accent/40 hover:bg-secondary/80 hover:scale-105 shadow-lg hover:shadow-accent/20 group"
                              }`}
                            aria-label={`${index !== currentToolIndex || completedSteps.includes(step.path)
                                ? `${step.name} (locked)`
                                : `Select ${step.name}`
                              }`}
                            title={step.description}
                          >
                            <div className={`p-2 rounded-lg border transition-all duration-300 ${completedSteps.includes(step.path)
                                ? "bg-success/20 border-success/30"
                                : index === currentToolIndex
                                  ? "bg-gradient-to-br from-accent/20 to-accent-secondary/20 border-accent/30 group-hover:animate-pulse-glow"
                                  : "bg-primary/20 border-text-secondary/20"
                              }`}>
                              {completedSteps.includes(step.path) ? (
                                <CheckCircle className="h-4 w-4 text-success" />
                              ) : (
                                <div className="text-accent">{iconMap[step.icon]}</div>
                              )}
                            </div>
                            <div className="text-left min-w-0 flex-1">
                              <div className="text-sm font-heading font-medium truncate text-text-primary">{step.name}</div>
                              <div className="text-xs text-text-secondary line-clamp-2 font-body">{step.description}</div>
                            </div>
                          </button>
                        ))}
                      </div>
                      <div className="flex justify-between gap-3 mt-4">
                        <button
                          onClick={async () => {
                            const prompt = "Returning to option selection."
                            setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
                            await speakResponse(prompt)
                            setSubMode("waiting_for_selection")
                            setSelectedOption(null)
                          }}
                          className="flex-1 px-4 py-2 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-accent to-accent-secondary hover:from-accent/90 hover:to-accent-secondary/90 transition-all duration-300 shadow-lg hover:shadow-accent/30 transform hover:scale-105 border border-accent/30"
                          aria-label="Go back to option selection"
                          title="Go back to option selection"
                        >
                          Back
                        </button>
                        <button
                          onClick={() => {
                            stopConversation()
                            setIsPanelOpen(false)
                          }}
                          className="flex-1 px-4 py-2 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-error to-error/80 hover:from-error/90 hover:to-error/70 transition-all duration-300 shadow-lg hover:shadow-error/30 transform hover:scale-105 border border-error/30"
                          aria-label="End conversation for step selection"
                        >
                          End Conversation
                        </button>
                      </div>
                    </div>
                  )}

                  {/* Enhanced Process ended options */}
                  {mode === "beginner" && subMode === "process_ended" && (
                    <div className="p-5 bg-gradient-to-br from-success/10 to-accent/5 border-b border-success/20 flex-shrink-0 overflow-y-auto animate-slide-up">
                      <div className="text-center mb-4">
                        <div className="bg-gradient-to-br from-success/20 to-accent/20 p-4 rounded-2xl w-16 h-16 mx-auto mb-3 flex items-center justify-center animate-pulse-glow border border-success/30">
                          <Sparkles className="h-8 w-8 text-success" />
                        </div>
                        <h3 className="text-xl font-heading font-bold text-success mb-2">Process Completed!</h3>
                        <p className="text-sm text-text-secondary font-body">Congratulations on completing the drug discovery journey</p>
                      </div>
                      <div className="grid grid-cols-1 gap-3">
                        <button
                          onClick={() => restartProcess()}
                          className="flex items-center space-x-3 p-4 bg-gradient-to-r from-secondary/70 to-secondary/50 rounded-xl border border-accent/20 hover:border-accent/40 hover:bg-secondary/80 transition-all duration-300 transform hover:scale-105 shadow-lg hover:shadow-accent/20"
                          aria-label="Restart drug discovery process"
                        >
                          <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-2 rounded-lg border border-accent/30 animate-pulse-glow">
                            {iconMap["FlaskConical"]}
                          </div>
                          <div className="text-left min-w-0 flex-1">
                            <div className="text-sm font-heading font-medium text-text-primary">Restart Drug Discovery</div>
                            <div className="text-xs text-text-secondary font-body">Start the process again from the beginning</div>
                          </div>
                        </button>
                        <button
                          onClick={() => handleDashboardSelection()}
                          className="flex items-center space-x-3 p-4 bg-gradient-to-r from-secondary/70 to-secondary/50 rounded-xl border border-accent/20 hover:border-accent/40 hover:bg-secondary/80 transition-all duration-300 transform hover:scale-105 shadow-lg hover:shadow-accent/20"
                          aria-label="Explore dashboard tools"
                        >
                          <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-2 rounded-lg border border-accent/30 animate-pulse-glow">
                            {iconMap["Home"]}
                          </div>
                          <div className="text-left min-w-0 flex-1">
                            <div className="text-sm font-heading font-medium text-text-primary">Explore Dashboard Tools</div>
                            <div className="text-xs text-text-secondary font-body">Try other available tools and features</div>
                          </div>
                        </button>
                        <button
                          onClick={() => {
                            setSubMode("waiting_for_selection")
                            const prompt = "Returning to option selection."
                            setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
                            speakResponse(prompt)
                          }}
                          className="flex items-center space-x-3 p-4 bg-gradient-to-r from-secondary/70 to-secondary/50 rounded-xl border border-accent/20 hover:border-accent/40 hover:bg-secondary/80 transition-all duration-300 transform hover:scale-105 shadow-lg hover:shadow-accent/20"
                          aria-label="Return to option selection"
                        >
                          <div className="bg-gradient-to-br from-accent/20 to-accent-secondary/20 p-2 rounded-lg border border-accent/30 animate-pulse-glow">
                            {iconMap["Bot"]}
                          </div>
                          <div className="text-left min-w-0 flex-1">
                            <div className="text-sm font-heading font-medium text-text-primary">Back to Options</div>
                            <div className="text-xs text-text-secondary font-body">Choose between dashboard or drug discovery</div>
                          </div>
                        </button>
                      </div>
                      <div className="flex justify-between gap-3 mt-4">
                        <button
                          onClick={() => {
                            stopConversation()
                            setIsPanelOpen(false)
                          }}
                          className="flex-1 px-4 py-2 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-error to-error/80 hover:from-error/90 hover:to-error/70 transition-all duration-300 shadow-lg hover:shadow-error/30 transform hover:scale-105 border border-error/30"
                          aria-label="End conversation for step selection"
                        >
                          End Conversation
                        </button>
                      </div>
                    </div>
                  )}

                  {/* Enhanced Conversation History */}
                  <div
                    ref={conversationRef}
                    className="flex-1 overflow-y-auto p-5 space-y-4 min-h-[200px] max-h-[50vh] bg-gradient-to-br from-primary/5 to-accent/5"
                  >
                    {conversationHistory.map((msg, index) => (
                      <div key={index} className={`flex ${msg.type === "user" ? "justify-end" : "justify-start"} animate-slide-up`}>
                        <div
                          className={`max-w-[85%] ${msg.type === "user"
                              ? "bg-gradient-to-br from-accent to-accent-secondary text-primary rounded-tr-none shadow-lg border border-accent/30"
                              : "bg-gradient-to-br from-secondary/90 to-secondary/70 text-text-primary rounded-tl-none border border-accent/20 shadow-lg backdrop-blur-sm"
                            } rounded-2xl px-5 py-4 mb-2 transition-all duration-300 hover:shadow-xl`}
                        >
                          {msg.type === "jarvis" && (
                            <div className="flex items-center space-x-2 mb-2 pb-2 border-b border-accent/20">
                              <Bot className="h-4 w-4 text-accent animate-pulse" />
                              <span className="text-xs font-code font-bold text-accent">JARVIS AI</span>
                              <div className="flex space-x-1">
                                <div className="w-1 h-1 bg-accent rounded-full animate-pulse"></div>
                                <div className="w-1 h-1 bg-accent-secondary rounded-full animate-pulse" style={{ animationDelay: '0.2s' }}></div>
                                <div className="w-1 h-1 bg-accent rounded-full animate-pulse" style={{ animationDelay: '0.4s' }}></div>
                              </div>
                            </div>
                          )}
                          <p className="text-sm whitespace-pre-wrap break-words font-body leading-relaxed">{msg.text}</p>
                        </div>
                      </div>
                    ))}
                    {(isLoading || isSpeaking) && (
                      <div className="flex justify-start animate-slide-up">
                        <div className="max-w-[85%] w-full py-4 px-5 bg-gradient-to-br from-secondary/90 to-secondary/70 border border-accent/20 rounded-2xl shadow-lg backdrop-blur-sm">
                          <div className="flex items-center space-x-2 pb-2 border-b border-accent/20 mb-3">
                            <Bot className="h-4 w-4 text-accent animate-pulse" />
                            <span className="text-xs font-code font-bold text-accent">JARVIS AI</span>
                            <div className="text-xs text-text-secondary font-body">Processing...</div>
                          </div>
                          <div className="flex space-x-2 items-center">
                            <div className="h-3 w-2 bg-accent rounded-full animate-wave"></div>
                            <div className="h-4 w-2 bg-accent-secondary rounded-full animate-wave" style={{ animationDelay: '0.1s' }}></div>
                            <div className="h-5 w-2 bg-accent rounded-full animate-wave" style={{ animationDelay: '0.2s' }}></div>
                            <div className="h-4 w-2 bg-accent-secondary rounded-full animate-wave" style={{ animationDelay: '0.3s' }}></div>
                            <div className="h-3 w-2 bg-accent rounded-full animate-wave" style={{ animationDelay: '0.4s' }}></div>
                          </div>
                        </div>
                      </div>
                    )}
                  </div>

                  {mode && (
                    // Wrap the entire block in a fragment so the comment is inside valid JSX
                    <>
                      {/* Enhanced Status Bar */}
                      <div className="px-5 py-3 mt-2 bg-gradient-to-r from-secondary/60 to-primary/40 border-t border-accent/20 flex-shrink-0 backdrop-blur-sm">
                        <div className="flex items-center justify-between text-sm">
                          <div className="flex items-center space-x-3">
                            <div
                              className={`w-3 h-3 rounded-full transition-all duration-300 ${isListening
                                  ? "bg-success animate-pulse shadow-lg shadow-success/50"
                                  : isSpeaking
                                    ? "bg-accent animate-pulse shadow-lg shadow-accent/50"
                                    : "bg-text-secondary/50"
                                }`}
                            ></div>
                            <span
                              className={`font-label font-medium ${isListening
                                  ? "text-success"
                                  : isSpeaking
                                    ? "text-accent"
                                    : "text-text-secondary"
                                }`}
                            >
                              {isListening ? "Listening..." : isSpeaking ? "Speaking..." : "Ready"}
                            </span>
                          </div>
                          <div className="flex items-center space-x-2 text-text-secondary">
                            <span className="font-code text-xs">
                              {conversationHistory.length} messages
                            </span>
                            {mode === "doubt" && (
                              <div className="flex items-center space-x-1">
                                <Mic className="h-3 w-3" />
                                <span className="text-xs font-code">Voice Mode</span>
                              </div>
                            )}
                          </div>
                        </div>
                      </div>
                    </>
                  )}


                  {/* Enhanced Input Area */}
                  <div className="p-5 border-t border-accent/20 flex-shrink-0 bg-gradient-to-r from-secondary/40 to-primary/20 backdrop-blur-sm">
                    {mode === "doubt" && (
                      <div className="flex items-center justify-between px-5 py-3 rounded-2xl bg-gradient-to-r from-secondary/80 to-secondary/60 border border-accent/20 shadow-lg backdrop-blur-sm">
                        <div className="flex items-center space-x-3 text-sm">
                          {isListening && mode === "doubt" ? (
                            <>
                              <div className="relative">
                                <div className="absolute inset-0 bg-success/20 rounded-full animate-ping"></div>
                                <Mic className="h-6 w-6 text-success animate-pulse" />
                              </div>
                              <div>
                                <span className="text-base font-heading font-medium text-success">Listening...</span>
                                <div className="text-xs text-text-secondary font-body">Speak your question now</div>
                              </div>
                            </>
                          ) : (
                            <>
                              <MicOff className="h-6 w-6 text-error animate-pulse" />
                              <div>
                                <span className="text-base font-heading font-medium text-error">Voice Paused</span>
                                <div className="text-xs text-text-secondary font-body">Click Start to begin</div>
                              </div>
                            </>
                          )}
                        </div>
                        <button
                          onClick={() => (isListening ? stopConversation() : startConversation("doubt"))}
                          className={`flex items-center px-5 py-2 rounded-xl text-sm font-heading font-semibold text-primary transition-all duration-300 transform hover:scale-105 shadow-lg border ${isListening
                              ? "bg-gradient-to-r from-error to-error/80 hover:from-error/90 hover:to-error/70 hover:shadow-error/30 border-error/30"
                              : "bg-gradient-to-r from-success to-accent hover:from-success/90 hover:to-accent/90 hover:shadow-success/30 border-success/30"
                            }`}
                          aria-label="Start or stop voice conversation"
                        >
                          {isListening ? "Stop" : "Start"}
                        </button>
                      </div>
                    )}
                    {mode === "beginner" && subMode === "waiting_for_question" && (
                      <div className="flex flex-col space-y-4">
                        <input
                          type="text"
                          placeholder="Ask a question about this tool..."
                          className="w-full px-5 py-3 rounded-2xl border border-accent/20 focus:outline-none focus:ring-2 focus:ring-accent focus:border-accent/40 text-sm bg-secondary/60 text-text-primary placeholder-text-secondary/70 font-body backdrop-blur-sm transition-all duration-300"
                          onKeyDown={async (e) => {
                            if (e.key === "Enter" && e.target.value.trim()) {
                              setSubMode("resolving_doubt")
                              const query = e.target.value.trim()
                              setConversationHistory((prev) => [...prev, { type: "user", text: query }])
                              await fetchResponse(query)
                              e.target.value = ""
                            }
                          }}
                          aria-label="Ask a question about this tool"
                        />
                        {suggestedQuestions[targetRoute]?.length && (
                          <div className="flex flex-wrap gap-2">
                            {suggestedQuestions[targetRoute].map((q, index) => (
                              <button
                                key={index}
                                onClick={() => {
                                  setSubMode("resolving_doubt")
                                  setConversationHistory((prev) => [...prev, { type: "user", text: q }])
                                  fetchResponse(q)
                                }}
                                className="px-4 py-2 bg-gradient-to-r from-accent/20 to-accent-secondary/20 text-accent-secondary rounded-xl text-xs hover:from-accent/30 hover:to-accent-secondary/30 hover:text-accent transition-all duration-300 font-body border border-accent/20 hover:border-accent/40 transform hover:scale-105"
                                aria-label={`Ask suggested question: ${q}`}
                                title={`Ask suggested question: ${q}`}
                              >
                                {q}
                              </button>
                            ))}
                          </div>
                        )}
                        <div className="flex justify-end gap-3">
                          {selectedOption === "drugDiscovery" && (
                            <button
                              onClick={() => proceedToNextStep()}
                              disabled={isSpeaking}
                              className="flex-1 px-5 py-3 rounded-xl bg-gradient-to-r from-accent to-accent-secondary text-sm font-heading font-semibold text-primary hover:from-accent/90 hover:to-accent-secondary/90 transition-all duration-300 shadow-lg hover:shadow-accent/30 transform hover:scale-105 border border-accent/30"
                              aria-label="Proceed to next step in drug discovery"
                            >
                              Next Step
                            </button>
                          )}
                          <button
                            onClick={async () => {
                              const prompt = `Returning to ${selectedOption === "drugDiscovery" ? "step" : "tool"} selection.`
                              setConversationHistory((prev) => [...prev, { type: "jarvis", text: prompt }])
                              await speakResponse(prompt)
                              setSubMode(
                                selectedOption === "drugDiscovery"
                                  ? "waiting_for_step_selection"
                                  : "waiting_for_tool_selection",
                              )
                            }}
                            disabled={isSpeaking}
                            className="flex-1 px-5 py-3 rounded-xl text-sm font-heading font-semibold text-primary bg-gradient-to-r from-accent-secondary to-accent hover:from-accent-secondary/90 hover:to-accent/90 transition-all duration-300 shadow-lg hover:shadow-accent-secondary/30 transform hover:scale-105 border border-accent-secondary/30"
                            aria-label={`Select another ${selectedOption === "drugDiscovery" ? "step" : "tool"}`}
                          >
                            Select Another {selectedOption === "drugDiscovery" ? "Step" : "Tool"}
                          </button>
                        </div>
                      </div>
                    )}

                    {/* Enhanced Help text */}
                    <div className="mt-4 py-3 px-5 bg-gradient-to-r from-accent/10 to-accent-secondary/10 rounded-xl border border-accent/20 backdrop-blur-sm">
                      <div className="text-xs text-center text-text-secondary leading-relaxed font-body flex items-center justify-center space-x-2">
                        <Sparkles className="h-4 w-4 text-accent animate-pulse" />
                        <span>Jarvis specializes in molecular structures, drug discovery, and pharmaceutical research</span>
                        <Zap className="h-4 w-4 text-accent-secondary animate-pulse" />
                      </div>
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  )
}