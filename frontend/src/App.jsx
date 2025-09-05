// App.jsx
import React from 'react';
import { Routes, Route, Navigate } from 'react-router-dom';
import { Toaster } from 'react-hot-toast';
import { useAuthStore } from './Store/auth.store.js';

import Navbar from './components/Navbar.jsx';
import Homepage from './pages/Homepage/Homepage.jsx';
import Signup from './pages/SIgnup/Signuppage.jsx';
import Login from './pages/Login/Loginpage.jsx';
import VerifyPhone from './pages/VerifyPhone/VerifyPhone.jsx'
import Dashboard from './components/Dashboard.jsx';
import DashboardHome from './components/Dashboardhome.jsx';
import ProteinStructureEvolution from './pages/Protienstructureevolution/Protienstructureevolution.jsx';
import ProteinStructureApp from './pages/Proteinstructureapp/ProteinStructureApp.jsx';
import AIResearchPaperGenerator from "./pages/AIresearchgenerator/Airesearchgenerator.jsx";
// import AIDrivenTargetPrediction from "./pages/AIdriventargetprediction/AIdriventargetprediction.jsx";
import Costestimation from "./pages/Costestimination/Costestimination.jsx";
import DrugDiscoveryRecommendation from "./pages/ToxicityPrediction/ToxicityPrediction.jsx";
import LiveNews from "./pages/Livenews/Livenews.jsx";
// import Message from "./pages/Message/Message.jsx";
import GetAlphaFoldStructure from "./pages/Alphafold/Alphafold.jsx"
import AINamingSuggestion from './pages/AINamingSuggestion/AINamingSuggestion.jsx';
import ToxicityPrediction from './pages/ToxicityPrediction/ToxicityPrediction.jsx';
import Summary from './pages/Summary/Summary.jsx';
import VoiceToTextNotes from './pages/VoiceToTextNotes/VoicetoTextNotes.jsx';
import Chatbot from './pages/Chatbot/Chatbot.jsx';
import SymptomChecker from './pages/SymptomChecker/SymptomChecker.jsx';
import RDkit from './pages/RDkit/RDkit.jsx';
import DashboardPage from './components/Dashboard.jsx';
// import Toxicityprediction from './pages/Toxicityprediction/Toxicityprediction.jsx';

const ProtectedRoute = ({ children }) => {
  const { user } = useAuthStore();
  console.log(user)
  return user ? children : <Navigate to="/login" replace />;
};

function App() {


  return (
    <div className="min-h-screen bg-red-400">
      {/* <Navbar /> */}
      <Routes>
        <Route path="/" element={<Homepage />} />
        <Route path="/signup" element={<Signup />} />
        <Route path="/login" element={<Login />} />
        <Route path="/verifyotp" element={<VerifyPhone />} />
        <Route path='/dashboard' element={
          <ProtectedRoute>
            <DashboardPage />
          </ProtectedRoute>
        } >
          <Route index element={<DashboardHome />} /> {/* Default route for /dashboard */}
          <Route
            path="protein-structure"
            element={
              <ProtectedRoute>
                <ProteinStructureApp />
              </ProtectedRoute>
            }
          />
          <Route
            path="protein-structure-mutation"
            element={
              <ProtectedRoute>
                <ProteinStructureEvolution />
              </ProtectedRoute>
            }
          />
          <Route
            path="cost-estimation"
            element={
              <ProtectedRoute>
                <Costestimation />
              </ProtectedRoute>
            }
          />
          <Route
            path="ai-research-paper-generator"
            element={
              <ProtectedRoute>
                <AIResearchPaperGenerator />
              </ProtectedRoute>
            }
          />
          {/* <Route
            path="ai-driven-target-prediction"
            element={
              <ProtectedRoute>
                <AIDrivenTargetPrediction />
              </ProtectedRoute>
            }
          /> */}
          <Route
            path="ai-naming"
            element={
              <ProtectedRoute>
                <AINamingSuggestion />
              </ProtectedRoute>
            }
          />
          <Route
            path="voice-text-notes"
            element={
              <ProtectedRoute>
                <VoiceToTextNotes />
              </ProtectedRoute>
            }
          />
          {/* <Route
            path="drug-discovery-recommendation"
            element={
              <ProtectedRoute>
                <DrugDiscoveryRecommendation />
              </ProtectedRoute>
            }
          /> */}
          <Route
            path="live-news"
            element={
              <ProtectedRoute>
                <LiveNews />
              </ProtectedRoute>
            }
          />

          <Route
            path="getalphafoldstrcture"
            element={
              <ProtectedRoute>
                <GetAlphaFoldStructure />
              </ProtectedRoute>
            }
          />

          <Route
            path="sideeffect-prediction"
            element={
              <ProtectedRoute>
                <ToxicityPrediction />
              </ProtectedRoute>
            }
          />
          <Route
            path="chatbot"
            element={
              <ProtectedRoute>
                <Chatbot />
              </ProtectedRoute>
            }
          />

          <Route
            path="summary"
            element={
              <ProtectedRoute>
                <Summary />
              </ProtectedRoute>
            }
          />
          <Route
            path="symptom-checker"
            element={
              <ProtectedRoute>
                <SymptomChecker />
              </ProtectedRoute>
            }
          />


          <Route
            path="newdrug-discovery"
            element={
              <ProtectedRoute>
                <RDkit />
              </ProtectedRoute>
            }
          />
        </Route>



      </Routes>
      {/* <Toaster /> */}
    </div>
  );
}

export default App;