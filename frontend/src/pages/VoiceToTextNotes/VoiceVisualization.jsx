import React, { useRef, useEffect } from 'react';

const VoiceVisualization = ({ audioLevel = 0, theme }) => {
  const canvasRef = useRef(null);
  
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    
    // Set canvas dimensions
    canvas.width = canvas.offsetWidth * 2; // Higher resolution
    canvas.height = canvas.offsetHeight * 2;
    ctx.scale(2, 2);
    
    // Clear canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    // Create gradient
    const gradient = ctx.createLinearGradient(0, 0, canvas.width / 2, 0);
    gradient.addColorStop(0, '#00F5D4'); // accent
    gradient.addColorStop(0.5, '#5E81F4'); // accent-secondary
    gradient.addColorStop(1, '#00F5D4'); // accent
    
    // Number of bars
    const barCount = 60;
    const barWidth = (canvas.width / 2) / barCount - 3;
    const centerY = (canvas.height / 2) / 2;
    
    // Draw visualization
    for (let i = 0; i < barCount; i++) {
      // Calculate bar height based on audio level and position
      let intensity = Math.sin((i / barCount) * Math.PI * 2) * audioLevel;
      // Add some randomness and wave effect
      intensity *= 0.6 + Math.random() * 0.4;
      intensity += Math.sin(Date.now() * 0.005 + i * 0.2) * 0.1;
      
      // Ensure minimum height
      const barHeight = Math.max(6, (canvas.height / 2) * 0.8 * intensity);
      
      // Position bars symmetrically around center
      const y = centerY - barHeight / 2;
      const x = i * (barWidth + 3);
      
      // Create rounded rectangle path
      ctx.fillStyle = gradient;
      ctx.beginPath();
      ctx.roundRect(x, y, barWidth, barHeight, barWidth / 2);
      ctx.fill();
      
      // Add glow effect
      ctx.shadowColor = '#00F5D4';
      ctx.shadowBlur = 10;
      ctx.fill();
      ctx.shadowBlur = 0;
    }
  }, [audioLevel, theme]);
  
  return (
    <div className="w-full h-16 rounded-2xl overflow-hidden bg-gradient-to-r from-secondary/50 to-secondary/30 border border-secondary/50 backdrop-blur-sm">
      <canvas ref={canvasRef} className="w-full h-full" />
    </div>
  );
};

export default VoiceVisualization;