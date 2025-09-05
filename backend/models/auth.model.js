import mongoose from "mongoose";

const userSchema = new mongoose.Schema(
  {
    firstName: {
      type: String,
      required: [true, "First name is required"],
    },
    lastName: {
      type: String,
      required: [true, "Last name is required"],
    },
    username: {
      type: String,
      required: [true, "Username is required"],
      unique: true,
      trim: true,
      minlength: [3, "Username must be at least 3 characters long"],
    },
    email: {
      type: String,
      required: [true, "Email is required"],
      unique: true,
      trim: true,
      lowercase: true,
      match: [/^\S+@\S+\.\S+$/, "Please use a valid email address"],
    },
    phoneNumber: {
      type: String,
      required: true,
      validate: {
        validator: function (v) {
          return /^\+\d{10,}$/.test(v);
        },
        message: "Please enter a valid phone number with country code (e.g., +911234567890)",
      },
    },
    password: {
      type: String,
      required: [true, "Password is required"],
      minlength: [6, "Password must be at least 6 characters long"],
    },
    confirmPassword: {
      type: String,
      required: [true, "Password is required"],
      minlength: [6, "Password must be at least 6 characters long"],
    },
  
    otp: {
      type: Number,
      required: true,
      validate: {
        validator: function (v) {
          return /^\d{6}$/.test(v.toString());
        },
        message: "OTP must be 6 digits",
      },
    },
    lastLogin: {
      type: Date,
      default: Date.now,
    },
    protienStructures: [
      {
        type: mongoose.Schema.Types.ObjectId,
        ref: "ProteinStructure",
      },
    ],
    NewproteinStructures: [
      { type: mongoose.Schema.Types.ObjectId, ref: "GeneratenewMolecule" },
    ],
    symptoms: [
      {
        type: mongoose.Schema.Types.ObjectId,
        ref: "predictDiseaseSchema",
      },
    ],
    productsmiles:[ 
    {
        type: mongoose.Schema.Types.ObjectId,
        ref: "Productsmiles",
      },
    ]
  },
  { timestamps: true }
);

export const User = mongoose.model("User", userSchema);