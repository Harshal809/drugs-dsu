import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAuthStore } from '../../Store/auth.store.js';
import { toast } from 'react-hot-toast';
import { ShieldCheck } from 'lucide-react';

// A component for the floating shapes in the background
const BackgroundShape = ({ className, children }) => (
  <div className={`absolute rounded-full filter blur-sm opacity-20 ${className}`}>
    {children}
  </div>
);

function VerifyPhone() {
  const [otp, setOtp] = useState('');
  const { verifyPhone, loading, user } = useAuthStore(); // Assuming user object contains phoneNumber
  const navigate = useNavigate();

  const handleSubmit = async (e) => {
    e.preventDefault();

    if (!otp || otp.length !== 6 || !/^\d+$/.test(otp)) {
      return toast.error('Please enter a valid 6-digit OTP');
    }

    if (!user?.phoneNumber) {
      toast.error('Phone number not found. Please sign up again.');
      return navigate('/signup');
    }

    try {
      await verifyPhone({ phoneNumber: user.phoneNumber, otp });
      toast.success('Phone number verified successfully!');
      setOtp('');
      navigate('/dashboard');
    } catch (error) {
      const errorMessage = error.response?.data?.message || 'OTP verification failed';
      toast.error(errorMessage);
    }
  };

  return (
    <div className="min-h-screen flex items-center justify-center bg-primary text-text-primary p-4 overflow-hidden relative">
      {/* Animated Background Shapes */}
      <BackgroundShape className="w-48 h-48 bg-accent/50 top-1/4 left-1/4 animate-pulse" />
      <BackgroundShape className="w-60 h-60 bg-accent-secondary/50 bottom-1/4 right-1/4 animate-pulse delay-1000" />
      <BackgroundShape className="w-32 h-32 bg-accent/30 bottom-10 left-10 animate-bounce" />
      <BackgroundShape className="w-40 h-40 bg-accent-secondary/40 top-20 right-20 animate-bounce delay-500" />
      
      {/* You can add more complex shapes or icons here if you like */}
      {/* Example with a simple text-based DNA icon */}
      <div className="absolute top-1/2 left-20 text-accent/20 text-8xl font-code animate-pulse select-none">§</div>
      <div className="absolute bottom-1/3 right-32 text-accent-secondary/20 text-6xl font-code animate-pulse delay-700 select-none">{"{...}"}</div>


      {/* OTP Verification Card */}
      <div className="w-full max-w-md p-8 bg-secondary rounded-xl shadow-lg border border-accent/20 z-10">
        <div className="text-center mb-8">
          <div className="flex justify-center items-center mx-auto w-16 h-16 bg-accent/10 rounded-full border-2 border-accent/30 mb-4">
            <ShieldCheck className="w-8 h-8 text-accent" />
          </div>
          <h2 className="text-3xl font-bold text-text-primary font-heading mb-2">
            Verify Your Account
          </h2>
          <p className="text-text-secondary font-body">
            An OTP has been sent to {user?.phoneNumber || 'your phone number'}.
          </p>
        </div>

        <form className="space-y-6" onSubmit={handleSubmit}>
          <div>
            <label
              htmlFor="otp"
              className="block text-sm font-medium text-text-primary font-label mb-2 text-center"
            >
              Enter 6-Digit Code
            </label>
            <input
              id="otp"
              type="text"
              value={otp}
              onChange={(e) => setOtp(e.target.value)}
              className="block w-full text-center tracking-[1em] text-2xl px-3 py-4 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent font-code"
              placeholder="••••••"
              maxLength={6}
              disabled={loading}
              autoComplete="one-time-code"
            />
          </div>

          <button
            type="submit"
            disabled={loading}
            className="w-full flex justify-center items-center py-3 px-4 rounded-lg shadow-sm text-lg font-medium text-primary bg-accent hover:bg-accent/90 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-accent disabled:opacity-70 transition-all duration-200 font-heading"
          >
            {loading ? (
              <>
                <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-primary" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                </svg>
                Verifying...
              </>
            ) : (
              'Verify & Proceed'
            )}
          </button>
        </form>
      </div>
    </div>
  );
}

export default VerifyPhone;