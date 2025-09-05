  /** @type {import('tailwindcss').Config} */
  export default {
    content: [
      "./index.html",
      "./src/**/*.{js,ts,jsx,tsx}"
    ],
    theme: {
      extend: {
        colors: {
          primary: '#0A192F',           // Main background
          secondary: '#172A45',         // Cards, panels
          accent: '#00F5D4',            // Buttons, highlights
          'accent-secondary': '#5E81F4',// Interactive elements
          error: '#FF4D6D',             // Alerts, warnings
          success: '#70E000',           // Confirmations
          'text-primary': '#E0E0E0',    // Main text
          'text-secondary': '#A0A0A0',  // Subtle text
        },
        fontFamily: {
          heading: ['Inter', 'Barlow', 'sans-serif'],        // For titles, headers
          body: ['Roboto', 'Open Sans', 'sans-serif'],       // For paragraphs, UI
          code: ['Fira Code', 'monospace'],                  // For code/AI output
          label: ['Lato', 'sans-serif'],                     // For scientific labels
        },
     
      },
    },
    plugins: [],
  };
