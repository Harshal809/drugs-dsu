import { Redis } from '@upstash/redis'
import dotenv from "dotenv";

dotenv.config();

// Create Redis instance with better error handling
export const redis = new Redis({
  url: process.env.UPSTASH_REDIS_REST_URL,
  token: process.env.UPSTASH_REDIS_REST_TOKEN,
  retry: {
    retries: 3,
    backoff: (retryCount) => Math.exp(retryCount) * 50,
  },
});