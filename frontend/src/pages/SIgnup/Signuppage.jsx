import React, { useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { useAuthStore } from '../../Store/auth.store.js';
import toast from 'react-hot-toast';
import { User, Mail, Phone, Lock, LogIn } from 'lucide-react';
import signupImage from './meetImage.avif';

function Signuppage() {
  const [formData, setFormData] = useState({
    firstName: '',
    lastName: '',
    username: '',
    email: '',
    phoneNumber: '',
    password: '',
    confirmPassword: '',
  });

  const { signup, loading } = useAuthStore();
  const navigate = useNavigate();

  const validateForm = () => {
    // Check for empty fields
    const requiredFields = ['firstName', 'lastName', 'username', 'email', 'phoneNumber', 'password', 'confirmPassword'];
    const emptyFields = requiredFields.filter(field => !formData[field] || formData[field].trim() === '');
    
    if (emptyFields.length > 0) {
      toast.error(`Please fill in all required fields: ${emptyFields.join(', ')}`);
      return false;
    }

    // Username validation (min 3 characters)
    if (formData.username.trim().length < 3) {
      toast.error("Username must be at least 3 characters long");
      return false;
    }

    // Email validation
    const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
    if (!emailRegex.test(formData.email)) {
      toast.error("Please enter a valid email address");
      return false;
    }

    // Phone validation - FIXED to match backend schema
    const phoneRegex = /^\+\d{10,}$/; // Must start with + and have at least 10 digits
    if (!phoneRegex.test(formData.phoneNumber)) {
      toast.error("Please enter a valid phone number starting with + and country code (e.g., +911234567890)");
      return false;
    }

    // Password validation
    if (formData.password.length < 6) {
      toast.error("Password must be at least 6 characters long");
      return false;
    }

    // Password match validation
    if (formData.password !== formData.confirmPassword) {
      toast.error("Passwords do not match");
      return false;
    }

    return true;
  };

  const handleSubmit = async (e) => {
    e.preventDefault();

    console.log('=== FORM SUBMISSION DEBUG ===');
    console.log('Form data:', formData);

    // Validate form
    if (!validateForm()) {
      return;
    }

    try {
      console.log('Attempting signup with data:', formData);
      await signup(formData);
      toast.success("Signup successful! Please verify your OTP");
      navigate('/verifyotp');
    } catch (error) {
      console.error('Signup error in component:', error);
      
      const errorMessage = error.response?.data?.message || "Signup failed. Please try again.";
      
      if (error.response?.status === 400) {
        if (errorMessage.includes("Email already exists")) {
          toast.error("Email already exists. Redirecting to login...");
          setTimeout(() => navigate("/login"), 2000);
        } else if (errorMessage.includes("Username already exists")) {
          toast.error("Username already exists. Please choose a different one.");
        } else {
          toast.error(errorMessage);
        }
      } else {
        toast.error(errorMessage);
      }
    }
  };

  const handleChange = (e) => {
    const { name, value } = e.target;
    
    // Clean up phone number input
    if (name === 'phoneNumber') {
      // Remove any non-digit characters except +
      let cleanValue = value.replace(/[^\d+]/g, '');
      
      // Ensure it starts with + if it has digits
      if (cleanValue.length > 0 && !cleanValue.startsWith('+')) {
        cleanValue = '+' + cleanValue;
      }
      
      setFormData({ ...formData, [name]: cleanValue });
    } else {
      setFormData({ ...formData, [name]: value });
    }
  };

  const renderInput = (id, name, type, placeholder, icon, value) => (
    <div className="relative">
      <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
        {icon}
      </div>
      <input
        id={id}
        name={name}
        type={type}
        required
        className="block w-full pl-10 pr-3 py-3 border border-secondary bg-primary text-text-primary rounded-lg focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent font-body"
        value={value}
        onChange={handleChange}
        placeholder={placeholder}
        // Add specific validations
        minLength={name === 'username' ? 3 : name === 'password' || name === 'confirmPassword' ? 6 : undefined}
        maxLength={name === 'phoneNumber' ? 15 : undefined}
      />
      {/* Show validation hints */}
      {name === 'phoneNumber' && (
        <div className="text-xs text-text-secondary mt-1 ml-2">
          Format: +[country code][number] (e.g., +911234567890)
        </div>
      )}
      {name === 'username' && (
        <div className="text-xs text-text-secondary mt-1 ml-2">
          Minimum 3 characters
        </div>
      )}
      {(name === 'password' || name === 'confirmPassword') && (
        <div className="text-xs text-text-secondary mt-1 ml-2">
          Minimum 6 characters
        </div>
      )}
    </div>
  );

  return (
    <div className="min-h-screen flex items-center justify-center bg-primary text-text-primary">
      <div className="flex w-full max-w-6xl mx-auto rounded-xl shadow-lg overflow-hidden border border-accent/20 bg-secondary">
        {/* Left side with image */}
        <div className="hidden lg:block lg:w-1/2 relative">
          <img
            src={signupImage}
            alt="Molecular structure"
            className="w-full h-full object-cover"
          />
          <div className="absolute inset-0 bg-primary opacity-60"></div>
          <div className="absolute inset-0 flex items-center justify-center">
            <div className="text-center p-8">
              <h1 className="text-4xl font-bold text-white font-heading">Join the Forefront of Discovery</h1>
              <p className="text-text-primary mt-4 text-lg font-body">Accelerate your research with our secure, AI-driven platform. Unlock new possibilities in drug discovery and genetic analysis.</p>
            </div>
          </div>
        </div>

        {/* Right side with Signup Form */}
        <div className="w-full lg:w-1/2 p-8 flex items-center justify-center">
          <div className="w-full max-w-md">
            <div className="text-center mb-8">
              <h2 className="text-3xl font-bold text-text-primary font-heading mb-2">
                Researcher Registration
              </h2>
              <p className="text-text-secondary font-body">Create your account to begin</p>
            </div>

            <form className="space-y-6" onSubmit={handleSubmit}>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {renderInput("firstName", "firstName", "text", "First Name", <User className="h-5 w-5 text-text-secondary" />, formData.firstName)}
                {renderInput("lastName", "lastName", "text", "Last Name", <User className="h-5 w-5 text-text-secondary" />, formData.lastName)}
              </div>
              
              {renderInput("username", "username", "text", "Username", <User className="h-5 w-5 text-text-secondary" />, formData.username)}
              {renderInput("email", "email", "email", "you@research.institution", <Mail className="h-5 w-5 text-text-secondary" />, formData.email)}
              {renderInput("phoneNumber", "phoneNumber", "tel", "+911234567890", <Phone className="h-5 w-5 text-text-secondary" />, formData.phoneNumber)}

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {renderInput("password", "password", "password", "Access Code", <Lock className="h-5 w-5 text-text-secondary" />, formData.password)}
                {renderInput("confirmPassword", "confirmPassword", "password", "Confirm Code", <Lock className="h-5 w-5 text-text-secondary" />, formData.confirmPassword)}
              </div>

              <button
                type="submit"
                disabled={loading}
                className="w-full flex justify-center items-center py-3 px-4 mt-6 rounded-lg shadow-sm text-lg font-medium text-primary bg-accent hover:bg-accent/90 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-accent disabled:opacity-70 transition-all duration-200 font-heading"
              >
                {loading ? (
                  <>
                    <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-primary" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                      <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                      <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                    </svg>
                    Creating Account...
                  </>
                ) : (
                  "Create Account & Send OTP"
                )}
              </button>

              <div className="text-center text-sm font-body pt-4">
                <span className="text-text-secondary">Already have research access? </span>
                <Link
                  to="/login"
                  className="font-medium text-accent hover:text-accent/80 flex items-center justify-center gap-1"
                >
                  <LogIn className="h-4 w-4" />
                  <span>Sign In</span>
                </Link>
              </div>
            </form>
          </div>
        </div>
      </div>
    </div>
  );
}

export default Signuppage;