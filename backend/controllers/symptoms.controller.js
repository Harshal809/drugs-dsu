import axios from 'axios';

const RAPIDAPI_URL = 'https://open-disease-api.p.rapidapi.com';
const RAPIDAPI_KEY = '40ca9e335cmsh545e7151d01d6e9p1f5d32jsn4e1365dfd96e';

export const getSymptoms = async (req, res) => {
  try {
    const response = await axios.get(`${RAPIDAPI_URL}/symptoms`, {
      headers: {
        'X-RapidAPI-Key': RAPIDAPI_KEY,
        'X-RapidAPI-Host': 'open-disease-api.p.rapidapi.com',
      },
    });

    // Validate response data
    if (!response.data || !Array.isArray(response.data)) {
      return res.status(500).json({ error: 'Invalid symptoms data from API' });
    }

    // Transform response to match frontend expectations
    const symptoms = response.data.map((symptom) => ({
      ID: symptom.id || symptom.ID, // Handle variations in API response
      Name: symptom.name || symptom.Name,
    }));

    res.json(symptoms);
  } catch (error) {
    console.error('Symptoms error:', error.response?.data || error.message);
    res.status(500).json({ error: 'Failed to fetch symptoms' });
  }
};

// Get diagnosis based on symptoms
export const getDiagnosis = async (req, res) => {
  try {
    const { symptoms, age = 30, sex = 'male' } = req.body;
    if (!symptoms || !Array.isArray(symptoms) || symptoms.length === 0) {
      return res.status(400).json({ error: 'Invalid or no symptoms provided' });
    }

    const response = await axios.post(
      `${RAPIDAPI_URL}/predict`,
      { symptoms, age, gender: sex }, // Include age and gender as per API requirements
      {
        headers: {
          'X-RapidAPI-Key': RAPIDAPI_KEY,
          'X-RapidAPI-Host': 'open-disease-api.p.rapidapi.com',
          'Content-Type': 'application/json',
        },
      }
    );

    // Validate response data
    if (!response.data || !Array.isArray(response.data)) {
      return res.status(500).json({ error: 'Invalid diagnosis data from API' });
    }

    // Transform response to match frontend expectations
    const conditions = response.data.map((condition) => ({
      Issue: {
        Name: condition.disease || condition.name,
        Accuracy: (condition.probability || condition.accuracy || 0) * 100, // Convert to percentage if needed
        IcdName: condition.id || condition.icd || 'Unknown',
      },
    }));

    res.json(conditions);
  } catch (error) {
    console.error('Diagnosis error:', error.response?.data || error.message);
    res.status(500).json({ error: 'Failed to get diagnosis' });
  }
};