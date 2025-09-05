
import { useState } from "react";
import { TextField, InputAdornment } from "@mui/material";
import LocalHospital from "@mui/icons-material/LocalHospital";
import { toast } from "react-hot-toast";

const validateSymptoms = (symptomList) => {
  if (!Array.isArray(symptomList) || symptomList.length === 0) {
    return { valid: false, message: "Please enter at least one valid symptom." };
  }
  if (symptomList.some((s) => typeof s !== "string" || s.trim() === "")) {
    return { valid: false, message: "All symptoms must be non-empty strings." };
  }
  return { valid: true };
};

export default function SymptomInput() {
  const [symptoms, setSymptoms] = useState([]);

  const handleSymptomChange = (e) => {
    let input = e.target.value;
    let symptomArray = input
      ? input
          .split(",")
          .map((s) => s.trim())
          .filter((s) => s)
      : [];

    // Update state to allow typing
    setSymptoms(symptomArray);

    // Validation 1: Check for valid symptoms
    const validation = validateSymptoms(symptomArray);
    if (!validation.valid) {
      if (input) {
        // Only show toast if there's input to avoid spamming on empty input
        toast.error(validation.message);
      }
      return;
    }

    // Validation 2: Check for duplicate symptoms
    const uniqueSymptoms = new Set(symptomArray);
    if (uniqueSymptoms.size !== symptomArray.length) {
      const duplicates = symptomArray.filter(
        (item, index) => symptomArray.indexOf(item) !== index
      );
      toast.error(`Duplicate symptoms found: ${duplicates.join(", ")}`);
      return;
    }

    // Validation 3: Check for maximum 500 characters
    if (input.length > 500) {
      toast.error("Maximum 500 characters allowed");
      return;
    }
  };

  const handleKeyDown = (e) => {
    // Add comma and space after pressing space or enter, but only if the last character isn't already a comma
    if (
      (e.key === " " || e.key === "Enter") &&
      e.target.value &&
      !e.target.value.trim().endsWith(",")
    ) {
      e.preventDefault();
      const newValue = e.target.value.trim() + ", ";
      setSymptoms(newValue.split(",").map((s) => s.trim()).filter((s) => s));
      e.target.value = newValue; // Update input field directly
    }
  };

  return (
    <TextField
      fullWidth
      label="Enter Symptoms"
      value={symptoms.join(", ")}
      onChange={handleSymptomChange}
      onKeyDown={handleKeyDown}
      placeholder="e.g., fever, cough, fatigue"
      variant="outlined"
      sx={{
        "& .MuiOutlinedInput-root": {
          bgcolor: "#0A192F",
          color: "#E0E0E0",
          "& fieldset": { borderColor: "#5E81F4" },
          "&:hover fieldset": { borderColor: "#00F5D4" },
          "&.Mui-focused fieldset": { borderColor: "#00F5D4" },
        },
        "& .MuiInputLabel-root": {
          color: "#A0A0A0",
          fontFamily: "'Lato', sans-serif",
        },
        "& .MuiInputLabel-root.Mui-focused": { color: "#00F5D4" },
      }}
      InputProps={{
        startAdornment: (
          <InputAdornment position="start">
            <LocalHospital sx={{ mr: 1, color: "#A0A0A0" }} />
          </InputAdornment>
        ),
      }}
    />
  );
}
