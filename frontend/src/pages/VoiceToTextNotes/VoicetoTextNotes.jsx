import React, { useState, useEffect } from 'react';
import { Mic, MicOff, AlertCircle, Sun, Moon, Sparkles } from 'lucide-react';
import RecordingControls from './RecordingControls';
import CurrentNote from './CurrentNote';
import NotesList from './NotesList';
import { useVoiceRecognition } from './hooks/useVoiceRecognition';
import { useNotes } from './hooks/useNotes';
import { useSpeechCommands } from './hooks/useSpeechCommands';

function VoiceToTextNotes() {
  const [theme, setTheme] = useState(() => {
    const savedTheme = localStorage.getItem('voice_notes_theme');
    return savedTheme || 'dark';
  });
  const [showAlert, setShowAlert] = useState(false);
  const [alertMessage, setAlertMessage] = useState('');

  const {
    isRecording,
    currentNote,
    setCurrentNote,
    startRecording,
    stopRecording,
    showContinue,
    error,
    audioLevel,
    isEditing,
    setShowContinue,
    setWasLastActionSave,
  } = useVoiceRecognition();

  const {
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
  } = useNotes(currentNote, setCurrentNote);

  // Setup speech commands
  useSpeechCommands({
    isRecording,
    saveNote,
    cancelEditing,
    currentNote,
    setCurrentNote
  });

  // Toggle theme
  const toggleTheme = () => {
    const newTheme = theme === 'light' ? 'dark' : 'light';
    setTheme(newTheme);
    localStorage.setItem('voice_notes_theme', newTheme);
    document.documentElement.classList.toggle('dark', newTheme === 'dark');
  };

  // Set theme class on document
  useEffect(() => {
    document.documentElement.classList.toggle('dark', theme === 'dark');
  }, [theme]);

  // Show alert
  const displayAlert = (message) => {
    setAlertMessage(message);
    setShowAlert(true);
    setTimeout(() => setShowAlert(false), 3000);
  };

  // Handle save note with feedback
  const handleSaveNote = async () => {
    if (isRecording) {
      stopRecording();
    }
    const result = await saveNote();
    if (result) {
      displayAlert('Note saved successfully!');
      localStorage.removeItem('voice_notes_draft');
      setCurrentNote('');
      setShowContinue(false);
      setWasLastActionSave(true);
    } else {
      displayAlert('Failed to save note. Please try again.');
    }
  };

  // Unsaved changes warning
  useEffect(() => {
    const handleBeforeUnload = (e) => {
      if (hasUnsavedChanges && currentNote.trim()) {
        e.preventDefault();
        e.returnValue = '';
        return '';
      }
    };

    window.addEventListener('beforeunload', handleBeforeUnload);
    return () => window.removeEventListener('beforeunload', handleBeforeUnload);
  }, [hasUnsavedChanges, currentNote]);

  return (
    <div className="min-h-screen bg-gradient-to-br transition-all duration-500">
      <div className="max-w-7xl mx-auto p-6 md:p-8">
        <div className="rounded-3xl overflow-hidden shadow-2xl bg-gradient-to-br from-secondary/60 to-secondary/40 backdrop-blur-xl border border-secondary/50">
          <header className="px-8 py-6 border-b border-secondary/50 bg-gradient-to-r from-secondary/40 to-secondary/60 backdrop-blur-sm">
            <div className="flex justify-between items-center">
              <div className="flex items-center gap-4">
                <div className="relative">
                  <div className="p-3 rounded-2xl bg-gradient-to-br from-accent/20 to-accent-secondary/20 border border-accent/30">
                    {isRecording ? (
                      <Mic className="h-8 w-8 text-accent animate-pulse" />
                    ) : (
                      <MicOff className="h-8 w-8 text-text-secondary" />
                    )}
                  </div>
                  {isRecording && (
                    <div className="absolute -top-1 -right-1 h-4 w-4 bg-error rounded-full animate-ping"></div>
                  )}
                </div>
                <div>
                  <h1 className="text-3xl font-heading font-bold bg-gradient-to-r from-accent via-accent-secondary to-accent bg-clip-text text-transparent">
                    Voice-to-Text Notes
                  </h1>
                  <p className="text-text-secondary font-body mt-1">
                    Capture your thoughts with AI-powered voice recognition
                  </p>
                </div>
              </div>
              
              <div className="flex items-center gap-3">
                <div className="flex items-center gap-2 px-4 py-2 bg-secondary/50 rounded-xl border border-secondary/50">
                  <Sparkles size={16} className="text-accent" />
                  <span className="text-sm text-text-secondary font-body">
                    {notes.length} notes saved
                  </span>
                </div>
              </div>
            </div>
          </header>

          {showAlert && (
            <div className={`mx-8 mt-6 px-6 py-4 rounded-2xl flex items-center gap-3 backdrop-blur-sm border transition-all duration-300 ${
              alertMessage.includes('Failed') 
                ? 'bg-error/10 text-error border-error/30' 
                : 'bg-success/10 text-success border-success/30'
            }`}>
              <div className="p-1 rounded-full bg-current/20">
                <AlertCircle size={18} />
              </div>
              <span className="font-medium font-body">{alertMessage}</span>
            </div>
          )}

          <RecordingControls 
            isRecording={isRecording}
            startRecording={() => startRecording(currentNote, false)}
            stopRecording={stopRecording}
            saveNote={handleSaveNote}
            showContinue={showContinue}
            error={error}
            isProcessing={isProcessing}
            hasUnsavedChanges={hasUnsavedChanges}
            currentNote={currentNote}
            audioLevel={audioLevel}
            theme={theme}
            isEditing={isEditing}
          />

          <CurrentNote 
            currentNote={currentNote}
            setCurrentNote={setCurrentNote}
            sentiment={sentiment}
            summary={summary}
            isEditing={!!editingNoteId}
            theme={theme}
          />

          <NotesList 
            notes={notes}
            editNote={(noteId) => {
              if (isRecording) stopRecording();
              editNote(noteId);
            }}
            editingNoteId={editingNoteId}
            updateNote={updateNote}
            cancelEditing={cancelEditing}
            startEditingWithRecording={(noteId) => {
              if (isRecording) stopRecording();
              const note = notes.find(n => n._id === noteId);
              if (note) {
                editNote(noteId);
                startRecording(note.content, true);
              }
            }}
            addTag={addTag}
            exportNotes={exportNotes}
            tags={tags}
            theme={theme}
            currentNote={currentNote}
            setCurrentNote={setCurrentNote}
            stopRecording={stopRecording}
            setShowContinue={setShowContinue}
            setWasLastActionSave={setWasLastActionSave}
          />
        </div>
      </div>
    </div>
  );
}

export default VoiceToTextNotes;