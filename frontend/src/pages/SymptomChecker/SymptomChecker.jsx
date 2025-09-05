import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { useAuthStore } from '../../Store/auth.store.js';
import { Search, AlertCircle, Loader2 } from 'lucide-react';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';

const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

const SymptomChecker = () => {
  const [symptoms, setSymptoms] = useState([]);
  const [selectedSymptoms, setSelectedSymptoms] = useState([]);
  const [searchTerm, setSearchTerm] = useState('');
  const [diagnosis, setDiagnosis] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [age, setAge] = useState('');
  const [sex, setSex] = useState('male');
  const { user, checkAuth } = useAuthStore();

  useEffect(() => {
    const initializeSymptoms = async () => {
      try {
        setLoading(true);
        const response = await axiosInstance.get('/symptoms/list');
        setSymptoms(response.data);
      } catch (err) {
        setError('Failed to load symptoms. Please try again later.');
        console.error('Error loading symptoms:', err);
      } finally {
        setLoading(false);
      }
    };

    initializeSymptoms();
  }, []);

  const handleSymptomSelect = (symptom) => {
    if (!selectedSymptoms.find((s) => s.ID === symptom.ID)) {
      setSelectedSymptoms([...selectedSymptoms, symptom]);
    }
    setSearchTerm('');
  };

  const handleRemoveSymptom = (symptomId) => {
    setSelectedSymptoms(selectedSymptoms.filter((s) => s.ID !== symptomId));
  };

  const handleDiagnosis = async () => {
    if (selectedSymptoms.length === 0) {
      setError('Please select at least one symptom');
      return;
    }
    if (!age || isNaN(age) || age < 1 || age > 120) {
      setError('Please enter a valid age between 1 and 120');
      return;
    }

    try {
      setLoading(true);
      setError(null);
      const response = await axiosInstance.post('/symptoms/diagnosis', {
        symptoms: selectedSymptoms.map((s) => s.ID),
        age: parseInt(age),
        sex,
      });
      setDiagnosis(response.data);
    } catch (err) {
      setError('Failed to get diagnosis. Please try again.');
      console.error('Diagnosis error:', err);
    } finally {
      setLoading(false);
    }
  };

  const filteredSymptoms = symptoms.filter((symptom) =>
    symptom.Name.toLowerCase().includes(searchTerm.toLowerCase())
  );

  if (!user) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-gray-50">
        <div className="text-center p-6 bg-white rounded-lg shadow-lg max-w-md">
          <AlertCircle className="h-12 w-12 text-red-500 mx-auto mb-4" />
          <h2 className="text-2xl font-bold text-gray-800 mb-2">Access Denied</h2>
          <p className="text-gray-600 mb-4">Please log in to use the Symptom Checker</p>
          <button
            onClick={() => (window.location.href = '/login')}
            className="bg-blue-600 text-white px-6 py-2 rounded-lg hover:bg-blue-700 transition-colors"
          >
            Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 py-12">
      <div className="max-w-4xl mx-auto px-4">
        <div className="bg-white rounded-xl shadow-lg p-6">
          <h1 className="text-3xl font-bold text-gray-800 mb-8 text-center">
            Symptom Checker
            <p className="text-sm text-blue-600 mt-1">Powered by Open Disease API</p>
          </h1>

          {/* Age and Sex Inputs */}
          <div className="mb-8 flex gap-4">
            <div className="flex-1">
              <label className="block text-sm font-medium text-gray-700 mb-1">Age</label>
              <input
                type="number"
                value={age}
                onChange={(e) => setAge(e.target.value)}
                placeholder="Enter age"
                className="w-full px-4 py-3 border rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
                required
              />
            </div>
            <div className="flex-1">
              <label className="block text-sm font-medium text-gray-700 mb-1">Sex</label>
              <select
                value={sex}
                onChange={(e) => setSex(e.target.value)}
                className="w-full px-4 py-3 border rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                <option value="male">Male</option>
                <option value="female">Female</option>
              </select>
            </div>
          </div>

          {/* Symptom Search */}
          <div className="mb-8">
            <div className="relative">
              <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 text-gray-400" />
              <input
                type="text"
                value={searchTerm}
                onChange={(e) => setSearchTerm(e.target.value)}
                placeholder="Search symptoms..."
                className="w-full pl-10 pr-4 py-3 border rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
              />
            </div>

            {/* Symptom Suggestions */}
            {searchTerm && (
              <div className="mt-2 max-h-60 overflow-y-auto border rounded-lg bg-white shadow-sm">
                {filteredSymptoms.length > 0 ? (
                  filteredSymptoms.map((symptom) => (
                    <button
                      key={symptom.ID}
                      onClick={() => handleSymptomSelect(symptom)}
                      className="w-full text-left px-4 py-2 hover:bg-gray-50 transition-colors"
                    >
                      {symptom.Name}
                    </button>
                  ))
                ) : (
                  <p className="px-4 py-2 text-gray-500">No matching symptoms found</p>
                )}
              </div>
            )}
          </div>

          {/* Selected Symptoms */}
          <div className="mb-8">
            <h2 className="text-lg font-semibold mb-4">Selected Symptoms:</h2>
            <div className="flex flex-wrap gap-2">
              {selectedSymptoms.length > 0 ? (
                selectedSymptoms.map((symptom) => (
                  <div
                    key={symptom.ID}
                    className="bg-blue-50 px-3 py-1 rounded-full flex items-center"
                  >
                    <span className="text-blue-700">{symptom.Name}</span>
                    <button
                      onClick={() => handleRemoveSymptom(symptom.ID)}
                      className="ml-2 text-blue-500 hover:text-blue-700"
                    >
                      ×
                    </button>
                  </div>
                ))
              ) : (
                <p className="text-gray-500">No symptoms selected</p>
              )}
            </div>
          </div>

          {/* Get Diagnosis Button */}
          <button
            onClick={handleDiagnosis}
            disabled={loading || selectedSymptoms.length === 0 || !age}
            className={`w-full py-3 rounded-lg font-medium transition-colors ${
              loading || selectedSymptoms.length === 0 || !age
                ? 'bg-gray-300 cursor-not-allowed'
                : 'bg-blue-600 hover:bg-blue-700 text-white'
            }`}
          >
            {loading ? (
              <div className="flex items-center justify-center">
                <Loader2 className="animate-spin mr-2" />
                Getting Diagnosis...
              </div>
            ) : (
              'Get Diagnosis'
            )}
          </button>

          {/* Error Message */}
          {error && (
            <div className="mt-4 p-4 bg-red-50 text-red-700 rounded-lg flex items-center">
              <AlertCircle className="mr-2" />
              {error}
            </div>
          )}

          {/* Diagnosis Results */}
          {diagnosis && (
            <div className="mt-8 bg-blue-50 p-6 rounded-lg">
              <h3 className="text-xl font-semibold mb-4">Possible Conditions:</h3>
              <div className="space-y-4">
                {diagnosis.length > 0 ? (
                  diagnosis.map((condition, index) => (
                    <div key={index} className="bg-white p-4 rounded-lg shadow-sm">
                      <h4 className="font-medium text-lg">{condition.Issue.Name}</h4>
                      <p className="text-gray-600 mt-1">
                        Probability: {condition.Issue.Accuracy}%
                      </p>
                      <p className="text-gray-600 mt-1">
                        Condition ID: {condition.Issue.IcdName}
                      </p>
                    </div>
                  ))
                ) : (
                  <p className="text-gray-600">No conditions found based on selected symptoms.</p>
                )}
              </div>
              <p className="mt-4 text-sm text-gray-600">
                Note: This is not a medical diagnosis. Please consult with a healthcare professional for proper medical advice.
              </p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default SymptomChecker;