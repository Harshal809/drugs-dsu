
import { useState, useEffect } from 'react';
import Spline from '@splinetool/react-spline';

function DNA() {
  const [isDesktop, setIsDesktop] = useState(false);
  const [isLoading, setIsLoading] = useState(true);

  useEffect(() => {
    const checkIsDesktop = () => {
      setIsDesktop(window.innerWidth >= 1024);
    };

    checkIsDesktop();
    window.addEventListener('resize', checkIsDesktop);
    
    return () => window.removeEventListener('resize', checkIsDesktop);
  }, []);

  const handleSplineLoad = () => {
    setIsLoading(false);
  };

  // On mobile and tablet, show the animated SVG instead of Spline
  if (!isDesktop) {
    return (
      <div className="w-full h-full flex items-center justify-center min-h-[400px]">
        <div className="relative w-64 h-64 flex items-center justify-center">
          {/* DNA Helix SVG Animation for mobile/tablet */}
          <div className="animate-spin-slow">
            <svg 
              width="200" 
              height="300" 
              viewBox="0 0 200 300" 
              className="text-accent"
              fill="currentColor"
            >
              <defs>
                <linearGradient id="dnaGradient" x1="0%" y1="0%" x2="100%" y2="100%">
                  <stop offset="0%" stopColor="#00F5D4" />
                  <stop offset="50%" stopColor="#5E81F4" />
                  <stop offset="100%" stopColor="#70E000" />
                </linearGradient>
              </defs>
              
              {/* DNA Helix Structure */}
              <path 
                d="M50 50 Q100 100 50 150 Q0 200 50 250" 
                stroke="url(#dnaGradient)" 
                strokeWidth="3" 
                fill="none"
                className="animate-pulse"
              />
              <path 
                d="M150 50 Q100 100 150 150 Q200 200 150 250" 
                stroke="url(#dnaGradient)" 
                strokeWidth="3" 
                fill="none"
                className="animate-pulse"
                style={{ animationDelay: '0.5s' }}
              />
              
              {/* Base pairs */}
              {[...Array(8)].map((_, i) => (
                <line
                  key={i}
                  x1="50"
                  y1={60 + i * 25}
                  x2="150"
                  y2={60 + i * 25}
                  stroke="url(#dnaGradient)"
                  strokeWidth="2"
                  opacity="0.7"
                  className="animate-pulse"
                  style={{ animationDelay: `${i * 0.1}s` }}
                />
              ))}
              
              {/* Molecules */}
              {[...Array(16)].map((_, i) => (
                <circle
                  key={i}
                  cx={i % 2 === 0 ? "50" : "150"}
                  cy={50 + Math.floor(i / 2) * 25}
                  r="4"
                  fill="url(#dnaGradient)"
                  className="animate-pulse"
                  style={{ animationDelay: `${i * 0.05}s` }}
                />
              ))}
            </svg>
          </div>
          
          {/* Floating particles around DNA */}
          <div className="absolute inset-0">
            {[...Array(12)].map((_, i) => (
              <div
                key={i}
                className="absolute w-1 h-1 bg-accent rounded-full animate-float"
                style={{
                  left: `${20 + Math.random() * 60}%`,
                  top: `${20 + Math.random() * 60}%`,
                  animationDelay: `${i * 0.2}s`,
                  animationDuration: `${3 + Math.random() * 2}s`
                }}
              />
            ))}
          </div>
        </div>
      </div>
    );
  }

  // Desktop view - show Spline model
  return (
    <div className="w-full h-full flex items-center justify-center min-h-[500px] lg:min-h-[600px]">
      <div className="spline-container relative w-full h-full">
        {isLoading && (
          <div className="absolute inset-0 flex items-center justify-center bg-primary/50 backdrop-blur-sm rounded-lg z-10">
            <div className="flex flex-col items-center space-y-4">
              <div className="animate-spin w-8 h-8 border-2 border-accent border-t-transparent rounded-full"></div>
              <span className="text-text-secondary text-sm">Loading 3D Model...</span>
            </div>
          </div>
        )}
        <Spline
          scene="https://prod.spline.design/veox8deuEmKQ5EXy/scene.splinecode"
          onLoad={handleSplineLoad}
        />
      </div>
    </div>
  );
}

export default DNA;
