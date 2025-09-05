import React, { useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { useAuthStore } from '../../Store/auth.store.js';
import { Mail, Lock, Copy, UserPlus, Shield } from 'lucide-react';
import loginImage from './meetImage.avif'; // Make sure to place your image in the assets folder and update the path

function Loginpage() {
  const [formData, setFormData] = useState({
    email: '',
    password: '',
  });

  const { login, loading } = useAuthStore();
  const navigate = useNavigate();

  const handleSubmit = async (e) => {
    e.preventDefault();
    await login({ email: formData.email, password: formData.password });
    if (useAuthStore.getState().user) {
      navigate('/dashboard');
    }
  };

  const handleChange = (e) => {
    setFormData({ ...formData, [e.target.name]: e.target.value });
  };

  const copyToClipboard = (text) => {
    navigator.clipboard.writeText(text);
    // Optional: Add a toast notification here to confirm copy
  };

  return (
    <div className="min-h-screen flex items-center justify-center bg-primary text-text-primary">
      <div className="flex w-full max-w-6xl mx-auto rounded-xl shadow-lg overflow-hidden border border-accent/20 bg-secondary">
        {/* Left side with image */}
        <div className="hidden md:block md:w-1/2 relative">
          <img
            src={loginImage}
            alt="Human DNA helix"
            className="w-full h-full object-cover"
          />
          <div className="absolute inset-0 bg-primary opacity-20"></div>
        </div>

        {/* Right side with Login Form */}
        <div className="w-full md:w-1/2 p-8 flex items-center justify-center">
          <div className="w-full max-w-md">
            <div className="text-center mb-8">
              <h2 className="text-3xl font-bold text-text-primary font-heading mb-2">
                SignIn to Research Portal
              </h2>
              <p className="text-text-secondary font-body">                    Where drug dicovery meet technology
              </p>
            </div>

            <form className="space-y-6" onSubmit={handleSubmit}>
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="email"
                    className="block text-sm font-medium text-text-primary font-label mb-2"
                  >
                    Researcher Email
                  </label>
                  <div className="relative">
                    <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                      <Mail className="h-5 w-5 text-text-secondary" />
                    </div>
                    <input
                      id="email"
                      name="email"
                      type="email"
                      required
                      className="block w-full pl-10 pr-3 py-3 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent font-body"
                      value={formData.email}
                      onChange={handleChange}
                      placeholder="you@research.institution"
                    />
                  </div>
                </div>

                <div>
                  <label
                    htmlFor="password"
                    className="block text-sm font-medium text-text-primary font-label mb-2"
                  >
                    Access Code
                  </label>
                  <div className="relative">
                    <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                      <Lock className="h-5 w-5 text-text-secondary" />
                    </div>
                    <input
                      id="password"
                      name="password"
                      type="password"
                      required
                      className="block w-full pl-10 pr-3 py-3 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent font-body"
                      value={formData.password}
                      onChange={handleChange}
                      placeholder="••••••••"
                    />
                  </div>
                </div>
              </div>

              <button
                type="submit"
                disabled={loading}
                className="w-full flex justify-center items-center py-3 px-4 rounded-lg shadow-sm text-lg font-medium text-primary bg-accent hover:bg-accent/90 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-accent disabled:opacity-70 transition-all duration-200 font-heading"
              >
                {loading ? (
                  <>
                    <svg
                      className="animate-spin -ml-1 mr-3 h-5 w-5 text-white"
                      xmlns="http://www.w3.org/2000/svg"
                      fill="none"
                      viewBox="0 0 24 24"
                    >
                      <circle
                        className="opacity-25"
                        cx="12"
                        cy="12"
                        r="10"
                        stroke="currentColor"
                        strokeWidth="4"
                      ></circle>
                      <path
                        className="opacity-75"
                        fill="currentColor"
                        d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                      ></path>
                    </svg>
                    Authenticating...
                  </>
                ) : (
                  "Access Research Portal"
                )}
              </button>

              <div className="text-center flex justify-center gap-4 text-sm font-body">
                <span className="text-text-secondary">New researcher? </span>
                <Link
                  to="/signup"
                  className="font-medium text-accent hover:text-accent/80 flex items-center justify-center gap-1"
                >
                  <UserPlus className="h-4 w-4" />
                  <span>Sign up</span>
                </Link>
              </div>
            </form>

            {/* Test Credentials Section */}
            <div className="mt-8 border-t border-primary pt-6">
              <h3 className="text-sm font-medium text-text-secondary font-body mb-4">
                Use below test credentials for testing, can't Signup due to use of demo Twilio account.
              </h3>
              <div className="space-y-3">
                <div className="flex items-center justify-between bg-primary/50 p-3 rounded-lg">
                  <div className="flex items-center gap-2">
                    <Mail className="h-4 w-4 text-text-secondary" />
                    <span className="text-sm text-text-primary font-body">
                      test@gmail.com
                    </span>
                  </div>
                  <button
                    onClick={() => copyToClipboard("test@gmail.com")}
                    className="p-1 hover:bg-accent/10 rounded text-text-secondary hover:text-accent"
                    title="Copy email"
                  >
                    <Copy className="h-4 w-4" />
                  </button>
                </div>
                <div className="flex items-center justify-between bg-primary/50 p-3 rounded-lg">
                  <div className="flex items-center gap-2">
                    <Lock className="h-4 w-4 text-text-secondary" />
                    <span className="text-sm text-text-primary font-body">1</span>
                  </div>
                  <button
                    onClick={() => copyToClipboard("1")}
                    className="p-1 hover:bg-accent/10 rounded text-text-secondary hover:text-accent"
                    title="Copy password"
                  >
                    <Copy className="h-4 w-4" />
                  </button>
                </div>
              </div>
            </div>

            <div className="mt-8">
              <div className="relative">
                <div className="absolute inset-0 flex items-center">
                  <div className="w-full border-t border-primary"></div>
                </div>
                <div className="relative flex justify-center">
                  <span className="px-3 bg-secondary text-sm text-text-secondary flex items-center gap-1 font-body">
                    <Shield className="h-4 w-4 text-accent" />
                    Secure AI-powered platform
                  </span>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

export default Loginpage;