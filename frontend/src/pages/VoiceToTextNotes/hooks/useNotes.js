import { useState, useEffect, useRef } from 'react';

export const useNotes = (currentNote, setCurrentNote) => {
  const [notes, setNotes] = useState([]);
  const [sentiment, setSentiment] = useState(null);
  const [summary, setSummary] = useState('');
  const [editingNoteId, setEditingNoteId] = useState(null);
  const [isProcessing, setIsProcessing] = useState(false);
  const [hasUnsavedChanges, setHasUnsavedChanges] = useState(false);
  const [tags, setTags] = useState([]);
  
  const lastSavedContent = useRef('');
  
  // Check for unsaved changes when current note changes
  useEffect(() => {
    if (currentNote !== lastSavedContent.current) {
      setHasUnsavedChanges(true);
      
      // Auto-save draft
      if (currentNote.trim()) {
        localStorage.setItem('voice_notes_draft', currentNote);
      }
    } else {
      setHasUnsavedChanges(false);
      localStorage.removeItem('voice_notes_draft');
    }
  }, [currentNote]);
  
  // Extract all tags from notes
  useEffect(() => {
    const allTags = [];
    
    notes.forEach(note => {
      if (note.tags && Array.isArray(note.tags)) {
        note.tags.forEach(tag => {
          if (!allTags.includes(tag)) {
            allTags.push(tag);
          }
        });
      }
    });
    
    setTags(allTags);
  }, [notes]);
  
  // Mock functions to simulate backend API calls
  const fetchNotes = async (userId) => {
    try {
      // Simulate API call
      const storedNotes = localStorage.getItem('voice_notes');
      return storedNotes ? JSON.parse(storedNotes) : [];
    } catch (error) {
      console.error('Error fetching notes:', error);
      return [];
    }
  };
  
  const saveNoteToAPI = async (noteData) => {
    try {
      // Simulate API call with local storage
      let allNotes = localStorage.getItem('voice_notes');
      allNotes = allNotes ? JSON.parse(allNotes) : [];
      
      if (noteData._id) {
        // Update existing note
        const updatedNotes = allNotes.map(note => 
          note._id === noteData._id ? { ...note, ...noteData } : note
        );
        localStorage.setItem('voice_notes', JSON.stringify(updatedNotes));
        return noteData;
      } else {
        // Create new note
        const newNote = { 
          ...noteData, 
          _id: `note_${Date.now()}`,
          createdAt: new Date().toISOString()
        };
        allNotes.push(newNote);
        localStorage.setItem('voice_notes', JSON.stringify(allNotes));
        return newNote;
      }
    } catch (error) {
      console.error('Error saving note:', error);
      throw new Error('Failed to save note');
    }
  };
  
  const analyzeSentiment = async (text) => {
    // Simplified sentiment analysis
    const lowerText = text.toLowerCase();
    
    const positiveWords = ['good', 'great', 'excellent', 'happy', 'love', 'best', 'awesome'];
    const negativeWords = ['bad', 'terrible', 'awful', 'sad', 'hate', 'worst', 'horrible'];
    
    let positiveCount = 0;
    let negativeCount = 0;
    
    positiveWords.forEach(word => {
      if (lowerText.includes(word)) positiveCount++;
    });
    
    negativeWords.forEach(word => {
      if (lowerText.includes(word)) negativeCount++;
    });
    
    if (positiveCount > negativeCount) return 'positive';
    if (negativeCount > positiveCount) return 'negative';
    return 'neutral';
  };
  
  const summarizeText = async (text) => {
    // Simplified summary - just return first sentence
    const firstSentence = text.split(/[.!?]/).filter(s => s.trim())[0] || '';
    return firstSentence.trim();
  };
  
  // Initialize notes
  useEffect(() => {
    const loadNotes = async () => {
      const loadedNotes = await fetchNotes('demo_user');
      setNotes(loadedNotes);
      
      // Load draft if it exists
      const savedDraft = localStorage.getItem('voice_notes_draft');
      if (savedDraft) {
        setCurrentNote(savedDraft);
        setHasUnsavedChanges(true);
      }
    };
    
    loadNotes();
  }, [setCurrentNote]);
  
  const saveNote = async () => {
    if (!currentNote.trim()) {
      return false;
    }
    
    setIsProcessing(true);
    
    try {
      // Analyze sentiment and generate summary
      const sentimentResult = await analyzeSentiment(currentNote);
      const summaryResult = await summarizeText(currentNote);
      
      // Create note object
      const noteData = {
        content: currentNote,
        sentiment: sentimentResult,
        summary: summaryResult,
        tags: editingNoteId ? (notes.find(note => note._id === editingNoteId)?.tags || []) : []
      };
      
      if (editingNoteId) {
        noteData._id = editingNoteId;
      } else {
        noteData.userId = 'demo_user';
        noteData.createdAt = new Date().toISOString();
      }
      
      // Save note
      const savedNote = await saveNoteToAPI(noteData);
      
      if (editingNoteId) {
        setNotes(notes.map(note => 
          note._id === editingNoteId ? savedNote : note
        ));
        setEditingNoteId(null);
      } else {
        setNotes([...notes, savedNote]);
      }
      
      // Reset current note
      setCurrentNote('');
      setSentiment(null);
      setSummary('');
      lastSavedContent.current = '';
      setHasUnsavedChanges(false);
      
      // Clear draft
      localStorage.removeItem('voice_notes_draft');
      
      return true;
    } catch (error) {
      console.error('Error saving note:', error);
      return false;
    } finally {
      setIsProcessing(false);
    }
  };
  
  const editNote = (noteId) => {
    const note = notes.find((n) => n._id === noteId);
    if (note) {
      setEditingNoteId(noteId);
      setCurrentNote(note.content); // Set currentNote for editing in CurrentNote
      setSentiment(note.sentiment);
      setSummary(note.summary);
      lastSavedContent.current = note.content;
      setHasUnsavedChanges(false);
    }
  };
  
  const updateNote = async (noteId, newContent) => {
    if (!newContent.trim()) {
      return false;
    }
    
    setIsProcessing(true);
    
    try {
      // Analyze sentiment and generate summary
      const sentimentResult = await analyzeSentiment(newContent);
      const summaryResult = await summarizeText(newContent);
      
      // Create updated note object
      const noteData = {
        _id: noteId,
        content: newContent,
        sentiment: sentimentResult,
        summary: summaryResult,
        tags: notes.find(note => note._id === noteId)?.tags || [],
        createdAt: notes.find(note => note._id === noteId)?.createdAt || new Date().toISOString()
      };
      
      // Save updated note
      const savedNote = await saveNoteToAPI(noteData);
      
      setNotes(notes.map(note => 
        note._id === noteId ? savedNote : note
      ));
      setEditingNoteId(null);
      setCurrentNote(''); // Reset currentNote after saving
      setSentiment(null);
      setSummary('');
      setHasUnsavedChanges(false);
      
      return true;
    } catch (error) {
      console.error('Error updating note:', error);
      return false;
    } finally {
      setIsProcessing(false);
    }
  };
  
  const cancelEditing = () => {
    setEditingNoteId(null);
    setCurrentNote(''); // Reset currentNote on cancel
    setSentiment(null);
    setSummary('');
    lastSavedContent.current = '';
    setHasUnsavedChanges(false);
  };
  
  const addTag = (noteId, tag) => {
    // Update note with new tag
    setNotes(notes.map(note => {
      if (note._id === noteId) {
        const updatedTags = note.tags ? [...note.tags] : [];
        if (!updatedTags.includes(tag)) {
          updatedTags.push(tag);
        }
        
        const updatedNote = { ...note, tags: updatedTags };
        saveNoteToAPI(updatedNote);
        return updatedNote;
      }
      return note;
    }));
  };
  
  const exportNotes = (notesToExport = notes) => {
    // Create text content for export
    let content = '# Voice Notes Export\n\n';
    
    notesToExport.forEach(note => {
      content += `## Note from ${new Date(note.createdAt).toLocaleString()}\n\n`;
      content += `${note.content}\n\n`;
      
      if (note.sentiment) {
        content += `Sentiment: ${note.sentiment}\n`;
      }
      
      if (note.summary) {
        content += `Summary: ${note.summary}\n`;
      }
      
      if (note.tags && note.tags.length > 0) {
        content += `Tags: ${note.tags.join(', ')}\n`;
      }
      
      content += '\n---\n\n';
    });
    
    // Create and download file
    const blob = new Blob([content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `voice-notes-export-${new Date().toISOString().slice(0, 10)}.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };
  
  return {
    notes,
    sentiment,
    summary,
    editingNoteId,
    saveNote,
    editNote,
    updateNote,
    cancelEditing,
    isProcessing,
    hasUnsavedChanges,
    addTag,
    exportNotes,
    tags,
  };
};