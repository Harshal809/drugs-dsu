import bcrypt from 'bcrypt';
import jwt from 'jsonwebtoken';
import { User } from '../models/auth.model.js';
import { generateTokenAndCookie } from '../utils/generateTokenAndCookie.js';
import twilio from 'twilio';
import { redis } from '../lib/redis.js';

// Twilio setup
const accountSid =process.env.ACCOUNT_SID;

const authToken = process.env.AUTH_TOKEN;
const client = twilio(accountSid, authToken);

// Generate OTP
const generateOTP = () => Math.floor(100000 + Math.random() * 900000);

// Send OTP via Twilio
const sendOTP = async (phoneNumber, otp) => {
  const fromNumber =  process.env.TWILIO_PHONE_NUMBER;
  if (!fromNumber) {
    console.error('TWILIO_PHONE_NUMBER is not defined in environment variables');
    throw new Error('Server configuration error: Missing Twilio phone number');
  }

  try {
    console.log(`Sending OTP to ${phoneNumber} from ${fromNumber}`);
    const message = await client.messages.create({
      body: `Your OTP is: ${otp}`,
      from: fromNumber,
      to: phoneNumber,
    });
    console.log(`OTP ${otp} sent to ${phoneNumber}, SID: ${message.sid}`);
  } catch (error) {
    console.error('Failed to send OTP via Twilio:', {
      message: error.message,
      code: error.code,
      moreInfo: error.moreInfo,
    });
    throw new Error(`Failed to send OTP: ${error.message} (Code: ${error.code})`);
  }
};

// Store OTP in Redis
const storeOTP = async (phoneNumber, otp) => {
  try {
    await redis.set(`otp:${phoneNumber}`, otp, 'EX', 300); // 5 minutes expiry
    console.log(`OTP stored for ${phoneNumber}`);
  } catch (error) {
    console.error('Failed to store OTP in Redis:', error.message);
    throw new Error('Temporary storage unavailable');
  }
};

// Verify OTP
const verifyOTP = async (phoneNumber, otp) => {
  try {
    const storedOTP = await redis.get(`otp:${phoneNumber}`);
    console.log(`Verifying OTP for ${phoneNumber}: stored=${storedOTP}, provided=${otp}`);
    return storedOTP === otp.toString(); // Ensure type consistency
  } catch (error) {
    console.error('Failed to verify OTP from Redis:', error.message);
    throw new Error('OTP verification unavailable');
  }
};

export const signup = async (req, res) => {
  const { firstName, lastName, username, email, phoneNumber, password, confirmPassword } = req.body;

  try {
    if (!firstName || !lastName || !username || !email || !phoneNumber || !password || !confirmPassword) {
      return res.status(400).json({ success: false, message: 'All fields are required' });
    }

    const userAlreadyExists = await User.findOne({ $or: [{ email }, { username }] });
    if (userAlreadyExists) {
      return res.status(400).json({
        success: false,
        message: userAlreadyExists.email === email ? 'Email already exists' : 'Username already exists',
      });
    }

    if (password !== confirmPassword) {
      return res.status(400).json({ success: false, message: 'Passwords do not match' });
    }

    const hashedPassword = await bcrypt.hash(password, 10);
    const otp = generateOTP();

    try {
      await redis.set(
        `signup:${phoneNumber}`,
        JSON.stringify({ firstName, lastName, username, email, password: hashedPassword,confirmPassword: hashedPassword, phoneNumber }),
        'EX',
        300
      );
      console.log(`Signup data stored for ${phoneNumber}`);
    } catch (redisError) {
      console.error('Redis set error in signup:', redisError.message);
      return res.status(503).json({ success: false, message: 'Temporary storage unavailable, please try again later' });
    }

    await storeOTP(phoneNumber, otp);
    await sendOTP(phoneNumber, otp);

    return res.status(200).json({
      success: true,
      message: 'OTP sent to your phone number',
      phoneNumber,
    });
  } catch (error) {
    console.error('Signup error:', error.message);
    return res.status(400).json({ success: false, message: error.message });
  }
};

export const verifyPhone = async (req, res) => {
  const { phoneNumber, otp } = req.body;

  try {
    if (!phoneNumber || !otp) {
      return res.status(400).json({ success: false, message: 'Phone number and OTP are required' });
    }

    const isValidOTP = await verifyOTP(phoneNumber, otp);
    if (!isValidOTP) {
      console.log(`OTP verification failed for ${phoneNumber}: Invalid or expired OTP`);
      return res.status(400).json({ success: false, message: 'Invalid or expired OTP' });
    }

    let userDetailsStr;
    try {
      userDetailsStr = await redis.get(`signup:${phoneNumber}`);
      console.log(`Retrieved signup data for ${phoneNumber}: ${userDetailsStr}`);
    } catch (redisError) {
      console.error('Redis get error in verifyPhone:', redisError.message);
      return res.status(503).json({ success: false, message: 'Temporary storage unavailable, please try again' });
    }

    if (!userDetailsStr) {
      console.log(`No signup data found for ${phoneNumber}`);
      return res.status(400).json({ success: false, message: 'Registration session expired. Please sign up again' });
    }

    const userDetails = JSON.parse(userDetailsStr);
    const user = new User({
      ...userDetails,
      otp: parseInt(otp),
    });

    await user.save();
    generateTokenAndCookie(res, user._id);

    try {
      await redis.del(`signup:${phoneNumber}`);
      await redis.del(`otp:${phoneNumber}`);
      console.log(`Cleaned up Redis data for ${phoneNumber}`);
    } catch (redisError) {
      console.error('Redis cleanup error:', redisError.message);
    }

    return res.status(201).json({
      success: true,
      message: 'User created successfully',
      user: { ...user._doc, password: undefined },
    });
  } catch (error) {
    console.error('Verify phone error:', error.message);
    return res.status(400).json({ success: false, message: error.message });
  }
};

export const login = async (req, res) => {
  const { email, password } = req.body;

  try {
    const user = await User.findOne({ email });
    if (!user) {
      return res.status(400).json({ success: false, message: 'Invalid Credentials' });
    }

    const isPasswordValid = await bcrypt.compare(password, user.password);
    if (!isPasswordValid) {
      return res.status(400).json({ success: false, message: 'Invalid Credentials' });
    }

    generateTokenAndCookie(res, user._id);
    user.lastLogin = new Date();
    await user.save();

    return res.status(200).json({
      success: true,
      message: 'Logged in Successfully',
      user: { ...user._doc, password: undefined },
    });
  } catch (error) {
    console.error('Login error:', error.message);
    return res.status(400).json({ success: false, message: `Login Failed: ${error.message}` });
  }
};

export const logout = async (req, res) => {
  try {
    res.clearCookie('token');
    return res.status(200).json({ success: true, message: 'Logged out successfully' });
  } catch (error) {
    console.error('Logout error:', error.message);
    return res.status(500).json({ success: false, message: 'Server error' });
  }
};

export const getProfile = async (req, res) => {
  try {
    if (!req.user) {
      return res.status(401).json({ success: false, message: 'Not authenticated' });
    }
    return res.status(200).json({
      success: true,
      user: { ...req.user._doc, password: undefined },
    });
  } catch (error) {
    console.error('Get profile error:', error.message);
    return res.status(500).json({ success: false, message: 'Server error', error: error.message });
  }
};