import { useEffect } from 'react';

export const useSpeechCommands = ({
  isRecording,
  saveNote,
  cancelEditing,
  currentNote,
  setCurrentNote,
}) => {
  useEffect(() => {
    if (!isRecording || !currentNote) return;
    
    // Check for voice commands in the current note
    const lowerNote = currentNote.toLowerCase();
    const words = lowerNote.split(/\s+/);
    const lastWords = words.slice(-3).join(' ');
    
    // Handle commands
    const handleCommands = async () => {
      // Save note command
      if (/(save note|save this|save it)/.test(lastWords)) {
        // Remove the command from the note
        setCurrentNote(prev => {
          const cleanedText = prev.replace(/(save note|save this|save it)\.?$/i, '').trim();
          return cleanedText;
        });
        
        // Execute save after a small delay
        setTimeout(() => {
          saveNote();
        }, 500);
      }
      
      // Cancel command
      else if (/(cancel|cancel edit|go back)/.test(lastWords)) {
        cancelEditing();
      }
      
      // Clear command
      else if (/(clear|clear all|delete all|start over)/.test(lastWords)) {
        setCurrentNote('');
      }
    };
    
    handleCommands();
  }, [currentNote, isRecording]);
  
  return {};
};