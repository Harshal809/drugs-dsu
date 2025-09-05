import { useState, useRef, useEffect } from 'react';

export const useVoiceRecognition = () => {
  const [isRecording, setIsRecording] = useState(false);
  const [currentNote, setCurrentNote] = useState('');
  const [showContinue, setShowContinue] = useState(false);
  const [error, setError] = useState('');
  const [audioLevel, setAudioLevel] = useState(0);
  const [isEditing, setIsEditing] = useState(false);
  const [wasLastActionSave, setWasLastActionSave] = useState(false); // Track if last action was save

  const recognitionRef = useRef(null);
  const streamRef = useRef(null);
  const audioContextRef = useRef(null);
  const analyserRef = useRef(null);
  const isCleaningUp = useRef(false);
  const lastTranscriptRef = useRef('');
  const unsavedChangesRef = useRef(false);

  // Initialize Web Speech API
  const createRecognition = () => {
    const SpeechRecognition = window.SpeechRecognition || window.webkitSpeechRecognition;
    if (!SpeechRecognition) {
      throw new Error('Speech Recognition API not supported in this browser.');
    }
    const recognition = new SpeechRecognition();
    recognition.lang = 'en-US';
    recognition.interimResults = false;
    recognition.maxAlternatives = 1;
    recognition.continuous = true;
    return recognition;
  };

  // Configure recognition handlers
  const setupRecognition = (recognition) => {
    recognition.onresult = (event) => {
      const speechResult = event.results[event.results.length - 1][0].transcript.trim();
      if (speechResult && speechResult !== lastTranscriptRef.current) {
        console.log('Captured transcript<SUPERSCRIPT>1</SUPERSCRIPT>:', speechResult);
        setCurrentNote((prev) => {
          unsavedChangesRef.current = true;
          return (prev ? prev + ' ' : '') + speechResult;
        });
        lastTranscriptRef.current = speechResult;
      }
    };

    recognition.onend = () => {
      console.log('Recognition ended, isRecording:', isRecording, 'isCleaningUp:', isCleaningUp.current);
      if (isRecording && !isCleaningUp.current) {
        setTimeout(() => {
          try {
            recognition.start();
            console.log('Speech recognition restarted');
          } catch (error) {
            console.error('Error restarting recognition:', error);
            stopRecording();
            setError('Failed to continue recording: ' + error.message);
          }
        }, 300);
      }
    };

    recognition.onerror = (event) => {
      console.error('Speech recognition error:', event.error);
      if (event.error === 'no-speech' || event.error === 'audio-capture') {
        console.log('Ignoring error, continuing recording:', event.error);
        return;
      }

      if (isRecording && !isCleaningUp.current) {
        setTimeout(() => {
          try {
            recognition.start();
            console.log('Retrying speech recognition');
          } catch (error) {
            console.error('Retry failed:', error);
            stopRecording();
            setError('Speech recognition error: ' + error.message);
          }
        }, 300);
      }
    };

    recognition.onstart = () => {
      console.log('Speech recognition started');
      setIsRecording(true);
    };
  };

  // Initialize recognition and handle tab visibility
  useEffect(() => {
    try {
      recognitionRef.current = createRecognition();
      setupRecognition(recognitionRef.current);
    } catch (error) {
      setError(error.message);
    }

    const handleVisibilityChange = () => {
      if (document.hidden && isRecording) {
        console.log('Tab hidden, pausing recording');
        stopRecording();
        if (currentNote.trim() && !wasLastActionSave) {
          // Only set showContinue if the last action wasn't a save
          setShowContinue(true);
          localStorage.setItem('voice_notes_draft', currentNote);
          console.log('Draft saved, showContinue set to true (tab hidden)');
        }
      }
    };

    document.addEventListener('visibilitychange', handleVisibilityChange);

    return () => {
      document.removeEventListener('visibilitychange', handleVisibilityChange);
      isCleaningUp.current = true;
      stopRecording();

      if (recognitionRef.current) {
        recognitionRef.current.stop();
        recognitionRef.current = null;
      }
    };
  }, []);

  // Check for saved draft on component mount
  useEffect(() => {
    const savedDraft = localStorage.getItem('voice_notes_draft');
    if (savedDraft && !wasLastActionSave) {
      setCurrentNote(savedDraft);
      setShowContinue(true);
    }
  }, []);

  // Setup audio visualization
  useEffect(() => {
    if (isRecording && streamRef.current) {
      try {
        if (!audioContextRef.current) {
          audioContextRef.current = new (window.AudioContext || window.webkitAudioContext)();
          analyserRef.current = audioContextRef.current.createAnalyser();
          analyserRef.current.fftSize = 256;

          const source = audioContextRef.current.createMediaStreamSource(streamRef.current);
          source.connect(analyserRef.current);
        }

        const dataArray = new Uint8Array(analyserRef.current.frequencyBinCount);

        const updateAudioLevel = () => {
          if (!isRecording) return;

          analyserRef.current.getByteFrequencyData(dataArray);
          let sum = 0;

          for (let i = 0; i < dataArray.length; i++) {
            sum += dataArray[i];
          }

          const average = sum / dataArray.length / 255; // Normalize to 0-1
          setAudioLevel(average);

          if (isRecording) {
            requestAnimationFrame(updateAudioLevel);
          }
        };

        updateAudioLevel();
      } catch (error) {
        console.error('Error setting up audio visualization:', error);
      }
    }

    return () => {
      if (audioContextRef.current) {
        if (audioContextRef.current.state !== 'closed') {
          try {
            analyserRef.current.disconnect();
          } catch (e) {
            console.log('Error disconnecting audio:', e);
          }
        }
      }
    };
  }, [isRecording]);

  const startRecording = async (existingText = '', isEditingNote = false) => {
    if (isRecording) {
      console.log('Already recording, ignoring start request');
      return;
    }

    try {
      const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
      streamRef.current = stream;

      setIsEditing(isEditingNote);
      if (!isEditingNote) {
        setCurrentNote(existingText);
      } else {
        console.log('Continuing recording for edited note:', currentNote);
      }

      lastTranscriptRef.current = '';
      isCleaningUp.current = false;
      setShowContinue(false);
      setError('');
      setWasLastActionSave(false); // Reset save flag on new recording

      setTimeout(() => {
        try {
          recognitionRef.current.start();
        } catch (error) {
          console.error('Error starting recognition:', error);
          stopRecording();
          setError('Failed to start recording: ' + error.message);
        }
      }, 100);
    } catch (error) {
      console.error('Error starting recording:', error);
      setError('Failed to start recording: Check microphone permissions.');
    }
  };

  const stopRecording = () => {
    if (recognitionRef.current) {
      recognitionRef.current.stop();
    }

    if (streamRef.current) {
      streamRef.current.getTracks().forEach((track) => track.stop());
      streamRef.current = null;
    }

    setIsRecording(false);
    setIsEditing(false);
    isCleaningUp.current = true;
    console.log('Recording stopped');
  };

  return {
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
    setWasLastActionSave, // Expose setWasLastActionSave for use in NotesList
  };
};