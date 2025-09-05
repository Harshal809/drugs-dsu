// src/api/costEstimation.js
import axios from 'axios';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode==="development" ? API_BASE_URL : '/api',
  withCredentials: true,
});

export const postCostEstimation = async (smiles) => {
  const response = await axiosInstance.post('/costestimation/cost-estimation', { smiles });
  return response.data;
};

export const getCostEstimations = async (userId) => {
  const response = await axiosInstance.get(`/costestimation/getcostestimation/${userId}`);
  return response.data;
};