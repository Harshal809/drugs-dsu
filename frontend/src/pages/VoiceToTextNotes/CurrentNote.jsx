import React, { useState } from 'react';
import { MessageSquare, AlertTriangle, CheckCircle, Info } from 'lucide-react';

const CurrentNote = ({ currentNote, setCurrentNote, sentiment, summary, isEditing, theme }) => {
  const [showMarkdownHelp, setShowMarkdownHelp] = useState(false);
  
  const renderSentimentIcon = () => {
    if (!sentiment) return null;
    
    switch(sentiment) {
      case 'positive':
        return <CheckCircle size={16} className="text-success" />;
      case 'negative':
        return <AlertTriangle size={16} className="text-error" />;
      default:
        return <Info size={16} className="text-accent" />;
    }
  };
  
  const handleNoteChange = (e) => {
    setCurrentNote(e.target.value);
  };
  
  return (
    <div className="p-6 md:p-8 border-b border-secondary/50 bg-gradient-to-r from-secondary/30 to-secondary/20">
      <div className="flex justify-between items-center mb-6">
        <h2 className="flex items-center gap-3 text-xl font-heading font-semibold text-text-primary">
          <div className="p-2 rounded-lg bg-accent/20 border border-accent/30">
            <MessageSquare size={20} className="text-accent" />
          </div>
          {isEditing ? 'Edit Note' : 'Current Note'}
        </h2>
      </div>
      
      <div className="mb-6">
        <textarea
          value={currentNote}
          onChange={handleNoteChange}
          placeholder="Start speaking or type your note here..."
          className="w-full min-h-[200px] p-6 bg-primary/50 border border-secondary/50 text-text-primary placeholder-text-secondary rounded-xl resize-y focus:ring-2 focus:ring-accent/50 focus:border-accent/50 focus:outline-none transition-all duration-300 font-body text-base leading-relaxed backdrop-blur-sm"
        ></textarea>
      </div>
      
      <div className="flex flex-wrap gap-6 text-sm text-text-secondary font-label">
        <div className="flex items-center gap-2 px-3 py-2 bg-secondary/30 rounded-lg border border-secondary/50">
          <span className="font-medium text-text-primary">
            {currentNote.trim() ? currentNote.trim().split(/\s+/).length : 0}
          </span>
          <span>words</span>
        </div>
        
        {sentiment && (
          <div className={`flex items-center gap-2 px-3 py-2 rounded-lg border ${
            sentiment === 'positive' 
              ? 'bg-success/10 border-success/30 text-success' 
              : sentiment === 'negative' 
                ? 'bg-error/10 border-error/30 text-error' 
                : 'bg-accent/10 border-accent/30 text-accent'
          }`}>
            {renderSentimentIcon()}
            <span className="font-medium">Sentiment: {sentiment}</span>
          </div>
        )}
        
        {summary && (
          <div className="flex items-center gap-2 px-3 py-2 bg-accent-secondary/10 border border-accent-secondary/30 rounded-lg text-accent-secondary">
            <span className="font-medium">Summary:</span>
            <span>{summary}</span>
          </div>
        )}
      </div>
    </div>
  );
};

export default CurrentNote;