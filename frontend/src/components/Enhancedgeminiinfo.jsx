import React, { useState } from 'react';
import { Icon } from '@iconify/react';

const EnhancedGeminiInfo = ({ information, isLoading, error }) => {
  const [isExpanded, setIsExpanded] = useState(false);

  // Loading state
  if (isLoading) {
    return (
      <div className="animate-pulse bg-secondary rounded-lg border border-secondary p-4">
        <div className="h-5 bg-primary rounded w-40 mb-2 opacity-30"></div>
        <div className="h-3 bg-primary rounded w-60 opacity-20"></div>
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="bg-secondary border border-error rounded-lg p-4 flex items-start">
        <Icon 
          icon="ph:warning-circle-fill" 
          className="w-5 h-5 text-error mr-2 mt-0.5" 
        />
        <div>
          <h4 className="text-error font-semibold text-sm font-heading">Analysis Error</h4>
          <p className="text-text-secondary text-xs font-body">{error}</p>
        </div>
      </div>
    );
  }

  // Empty state
  if (!information || information.includes('Failed to retrieve detailed information')) {
    return (
      <div className="bg-secondary border border-secondary rounded-lg p-4 text-center">
        <Icon icon="ph:atom" className="w-6 h-6 text-accent mx-auto mb-2" />
        <h4 className="text-text-primary font-medium text-sm font-heading">No Analysis Available</h4>
        <p className="text-text-secondary text-xs font-body">Detailed compound information will appear here once generated</p>
      </div>
    );
  }

  // Parse the information into sections and bullet points
  const parseInformation = (text) => {
    const sectionRegex = /(\d+\.\s[A-Za-z\s]+:)/g;
    const sections = text.split(sectionRegex).filter(Boolean);

    const parsedSections = [];
    let currentSection = null;

    for (let i = 0; i < sections.length; i++) {
      const item = sections[i].trim();
      if (item.match(sectionRegex)) {
        if (currentSection) {
          parsedSections.push(currentSection);
        }
        currentSection = { title: item, bulletPoints: [] };
      } else if (currentSection) {
        const bulletPoints = item
          .split('\n')
          .map(line => line.trim())
          .filter(line => line.startsWith('-') || line.length > 0)
          .map(line => line.replace(/^-/, '').trim());
        currentSection.bulletPoints.push(...bulletPoints);
      }
    }

    if (currentSection) {
      parsedSections.push(currentSection);
    }

    return parsedSections.filter(section => section.bulletPoints.length > 0);
  };

  const sections = parseInformation(information);

  return (
    <div className="bg-secondary rounded-lg border border-secondary shadow-lg">
      {/* Dropdown Header */}
      <button
        onClick={() => setIsExpanded(!isExpanded)}
        className={`w-full flex items-center justify-between px-4 py-3 transition-colors ${isExpanded ? 'bg-primary/10' : 'bg-secondary hover:bg-primary/10'}`}
      >
        <div className="flex items-center space-x-2">
          <Icon icon="ph:info-fill" className="w-5 h-5 text-accent" />
          <span className="text-sm font-medium text-text-primary font-heading">
            Compound Analysis <span className='text-xs text-accent-secondary font-label'>(Powered by Gemini)</span>
          </span>
        </div>
        <Icon
          icon={isExpanded ? 'ph:caret-up' : 'ph:caret-down'}
          className="w-5 h-5 text-accent-secondary"
        />
      </button>

      {/* Expanded Content */}
      {isExpanded && (
        <div className="px-4 py-3 border-t border-primary/20 animate-fade-in">
          {sections.length > 0 ? (
            <div className="space-y-4">
              {sections.map((section, index) => (
                <div key={index} className="mb-4 last:mb-0">
                  <h4 className="text-sm font-semibold text-accent mb-2 font-heading tracking-wide">
                    {section.title.replace(':', '')}
                  </h4>
                  <ul className="space-y-2">
                    {section.bulletPoints.map((point, i) => (
                      <li key={i} className="flex items-start text-sm text-text-primary font-body">
                        <span className="text-accent-secondary mr-2 mt-1">•</span>
                        <span className="leading-relaxed">
                          {point || 'No detailed information available for this section.'}
                        </span>
                      </li>
                    ))}
                  </ul>
                </div>
              ))}
            </div>
          ) : (
            <div className="bg-primary/5 p-3 rounded border border-primary/10">
              <p className="text-sm text-text-secondary font-code">
                {information}
              </p>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default EnhancedGeminiInfo;