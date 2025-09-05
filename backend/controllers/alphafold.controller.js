import axios from 'axios';
import { Job } from '../models/alphafold.model.js';
import { v4 as uuidv4 } from 'uuid';

const UNIPROT_API = 'https://www.ebi.ac.uk/proteins/api';
const ALPHAFOLD_API = 'https://alphafold.ebi.ac.uk/files';
const TIMEOUT = 30000;

const validateUniprotId = (uniprotId) => {
  if (!uniprotId || typeof uniprotId !== 'string') {
    return { valid: false, message: 'UniProt ID must be a non-empty string' };
  }
  const validUniprot = /^[A-Z0-9]{6,10}$/i.test(uniprotId);
  if (!validUniprot) {
    return { valid: false, message: 'Invalid UniProt ID format' };
  }
  return { valid: true };
};

const fetchAlphaFoldPDB = async (uniprotId) => {
  try {
    const pdbUrl = `${ALPHAFOLD_API}/AF-${uniprotId}-F1-model_v4.pdb`;
    const response = await axios.get(pdbUrl, { responseType: 'text', timeout: TIMEOUT });
    return { data: response.data, url: pdbUrl };
  } catch (error) {
    console.error(`Failed to fetch AlphaFold PDB for ${uniprotId}: ${error.message}`);
    throw new Error(`AlphaFold PDB fetch failed: ${error.message}`);
  }
};

const processAlphaFold = async (job) => {
  try {
    const { data, url } = await fetchAlphaFoldPDB(job.uniprot_id);
    job.pdbUrl = url;
    job.status = 'completed';
    job.completedAt = new Date();
    job.updatedAt = new Date();
    await job.save();
    console.log(`Job ${job.pythonJobId} for ${job.uniprot_id} completed successfully`);
  } catch (error) {
    console.error(`Error processing job ${job.pythonJobId}: ${error.message}`);
    job.status = 'failed';
    job.error = error.message;
    job.updatedAt = new Date();
    await job.save();
  }
};

export const submitPrediction = async (req, res) => {
  let job;
  try {
    const { uniprot_id } = req.body;
    console.log('Received prediction request:', uniprot_id);

    const validation = validateUniprotId(uniprot_id);
    if (!validation.valid) {
      return res.status(400).json({ error: validation.message });
    }

    // Check if job already exists for this user and uniprot_id
    job = await Job.findOne({ uniprot_id, userId: req.user?._id });
    if (job && job.status === 'completed') {
      return res.status(200).json({
        jobId: job._id,
        pythonJobId: job.pythonJobId,
        status: job.status,
        uniprotId: job.uniprot_id,
        pdbUrl: job.pdbUrl,
        createdAt: job.createdAt,
        updatedAt: job.updatedAt,
        completedAt: job.completedAt,
      });
    }

    // Create new job with a unique pythonJobId
    job = new Job({
      uniprot_id,
      pythonJobId: uuidv4(),
      status: 'pending',
      userId: req.user?._id,
      createdAt: new Date(),
      updatedAt: new Date(),
    });
    await job.save();

    job.status = 'processing';
    job.updatedAt = new Date();
    await job.save();

    // Process in background
    processAlphaFold(job).catch(async (err) => {
      console.error(`Background AlphaFold error for ${job.pythonJobId}: ${err.message}`);
      job.status = 'failed';
      job.error = err.message;
      job.updatedAt = new Date();
      await job.save();
    });

    res.status(202).json({
      jobId: job._id,
      pythonJobId: job.pythonJobId,
      status: job.status,
      uniprotId: uniprot_id,
      createdAt: job.createdAt,
      updatedAt: job.updatedAt,
    });
  } catch (error) {
    console.error('Submission error:', error);
    if (job) {
      job.status = 'failed';
      job.error = error.message;
      job.updatedAt = new Date();
      await job.save();
    }
    res.status(500).json({ error: error.message || 'Failed to submit prediction' });
  }
};

export const getPreviousJobs = async (req, res) => {
  try {
    const jobs = await Job.find({ userId: req.user?._id })
      .sort({ createdAt: -1 })
      .limit(10);
    res.json(jobs.map(job => ({
      jobId: job._id,
      pythonJobId: job.pythonJobId,
      uniprotId: job.uniprot_id,
      status: job.status,
      pdbUrl: job.pdbUrl,
      createdAt: job.createdAt,
      updatedAt: job.updatedAt,
      completedAt: job.completedAt,
      error: job.error,
    })));
  } catch (error) {
    console.error('Previous jobs error:', error);
    res.status(500).json({ error: 'Failed to fetch previous jobs' });
  }
};

export const getUniprotSummary = async (req, res) => {
  try {
    const { uniprotId } = req.params;
    const response = await axios.get(`${UNIPROT_API}/proteins/${uniprotId}`, { timeout: TIMEOUT });
    res.json(response.data);
  } catch (error) {
    console.error(`Failed to fetch UniProt summary: ${error.message}`);
    res.status(error.response?.status || 500).json({
      error: error.message || 'Failed to fetch UniProt summary',
    });
  }
};

export const getUniprotAnnotations = async (req, res) => {
  try {
    const { uniprotId } = req.params;
    const response = await axios.get(`${UNIPROT_API}/features/${uniprotId}`, { timeout: TIMEOUT });
    res.json(response.data);
  } catch (error) {
    console.error(`Failed to fetch UniProt annotations: ${error.message}`);
    res.status(error.response?.status || 500).json({
      error: error.message || 'Failed to fetch annotations',
    });
  }
};

export const getJobStatus = async (req, res) => {
  try {
    const { jobId } = req.params;
    const job = await Job.findById(jobId);
    if (!job) {
      return res.status(404).json({ error: 'Job not found' });
    }
    res.json({
      jobId: job._id,
      pythonJobId: job.pythonJobId,
      status: job.status,
      uniprotId: job.uniprot_id,
      pdbUrl: job.pdbUrl,
      createdAt: job.createdAt,
      updatedAt: job.updatedAt,
      completedAt: job.completedAt,
      error: job.error,
    });
  } catch (error) {
    console.error('Job status error:', error);
    res.status(500).json({ error: 'Failed to fetch job status' });
  }
};

