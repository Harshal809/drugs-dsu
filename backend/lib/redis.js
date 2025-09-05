import { Redis } from '@upstash/redis'
import dotenv from "dotenv";

dotenv.config();


export const redis = new Redis({
  url: 'https://immense-goldfish-16178.upstash.io',
  token: 'AT8yAAIjcDE5YTI1YWJlN2Q0YmE0NTQ3OWFjNTZmYWMwZDViMjdhMHAxMA',
})

