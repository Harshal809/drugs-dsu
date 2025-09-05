import React from 'react';
import { Mic, MicOff, Save, Play, Loader, Volume2 } from 'lucide-react';
import VoiceVisualization from './VoiceVisualization';

const RecordingControls = ({
  isRecording,
  startRecording,
  stopRecording,
  saveNote,
  showContinue,
  error,
  isProcessing,
  hasUnsavedChanges,
  currentNote,
  audioLevel,
  theme,
  isEditing,
}) => {
  return (
    <div className="p-6 md:p-8 border-b border-secondary/50 bg-gradient-to-r from-secondary/20 to-secondary/30">
      <div className="flex flex-wrap gap-4 justify-center mb-6">
        {!showContinue && !isRecording && !isEditing && (
          <button
            onClick={startRecording}
            className="flex items-center gap-3 px-8 py-4 rounded-2xl font-medium text-white bg-gradient-to-r from-success to-success/80 hover:from-success/90 hover:to-success/70 transition-all duration-300 shadow-lg disabled:opacity-50 disabled:cursor-not-allowed font-body hover:scale-105 border border-success/30"
            aria-label="Start recording"
            title="Start recording"
            disabled={!!error}
          >
            <div className="p-1 rounded-full bg-white/20">
              <Mic size={20} />
            </div>
            <span className="text-lg">Start Recording</span>
          </button>
        )}

        {isRecording && (
          <button
            onClick={stopRecording}
            className="flex items-center gap-3 px-8 py-4 rounded-2xl font-medium text-white bg-gradient-to-r from-error to-error/80 hover:from-error/90 hover:to-error/70 transition-all duration-300 shadow-lg disabled:opacity-50 disabled:cursor-not-allowed font-body hover:scale-105 border border-error/30"
            aria-label="Stop recording"
            title="Stop recording"
          >
            <div className="p-1 rounded-full bg-white/20">
              <MicOff size={20} />
            </div>
            <span className="text-lg">Stop Recording</span>
          </button>
        )}

        {showContinue && currentNote.trim() && !isRecording && !isEditing && (
          <button
            onClick={startRecording}
            className="flex items-center gap-3 px-8 py-4 rounded-2xl font-medium text-white bg-gradient-to-r from-accent-secondary to-accent-secondary/80 hover:from-accent-secondary/90 hover:to-accent-secondary/70 transition-all duration-300 shadow-lg disabled:opacity-50 disabled:cursor-not-allowed font-body hover:scale-105 border border-accent-secondary/30"
            aria-label="Continue recording"
            title="Continue recording"
            disabled={!!error}
          >
            <div className="p-1 rounded-full bg-white/20">
              <Play size={20} />
            </div>
            <span className="text-lg">Continue Recording</span>
          </button>
        )}

        <button
          onClick={saveNote}
          className={`flex items-center gap-3 px-8 py-4 rounded-2xl font-medium text-white bg-gradient-to-r from-accent to-accent/80 hover:from-accent/90 hover:to-accent/70 transition-all duration-300 shadow-lg disabled:opacity-50 disabled:cursor-not-allowed font-body hover:scale-105 border border-accent/30 ${
            hasUnsavedChanges ? 'animate-pulse' : ''
          }`}
          aria-label="Save note"
          title="Save note"
          disabled={!currentNote.trim() || !!error}
        >
          <div className="p-1 rounded-full bg-white/20">
            {isProcessing ? <Loader className="animate-spin\" size={20} /> : <Save size={20} />}
          </div>
          <span className="text-lg">Save Note</span>
        </button>
      </div>

      {isRecording && (
        <div className="flex flex-col items-center space-y-4">
          <VoiceVisualization audioLevel={audioLevel} theme={theme} />
          <div className="flex items-center text-base gap-3 text-error font-body">
            <div className="relative">
              <span className="h-3 w-3 rounded-full bg-error animate-pulse"></span>
              <span className="absolute h-3 w-3 rounded-full bg-error animate-ping"></span>
            </div>
            <span className="font-medium">Recording in progress...</span>
            <Volume2 size={18} className="text-error animate-pulse" />
          </div>
        </div>
      )}

      {error && (
        <div className="mt-6 p-4 rounded-xl text-sm bg-error/10 text-error border border-error/30 backdrop-blur-sm font-body">
          <div className="flex items-center gap-2">
            <div className="p-1 rounded-full bg-error/20">
              <AlertTriangle size={16} />
            </div>
            <span className="font-medium">{error}</span>
          </div>
        </div>
      )}

      <div className="mt-6">
        <details className="text-sm">
          <summary className="cursor-pointer text-accent hover:text-accent/80 transition-colors font-medium font-body flex items-center gap-2">
            <Volume2 size={16} />
            Voice Commands
          </summary>
          <div className="mt-4 p-4 rounded-xl bg-secondary/30 border border-secondary/50 backdrop-blur-sm">
            <ul className="space-y-2 text-text-secondary font-body">
              <li className="flex items-center gap-3">
                <span className="px-2 py-1 bg-accent/20 text-accent rounded-lg text-xs font-code">"Save note"</span>
                <span>Saves current note</span>
              </li>
              <li className="flex items-center gap-3">
                <span className="px-2 py-1 bg-accent/20 text-accent rounded-lg text-xs font-code">"Cancel"</span>
                <span>Cancels editing</span>
              </li>
              <li className="flex items-center gap-3">
                <span className="px-2 py-1 bg-accent/20 text-accent rounded-lg text-xs font-code">"Clear"</span>
                <span>Clears current text</span>
              </li>
              <li className="flex items-center gap-3">
                <span className="px-2 py-1 bg-accent/20 text-accent rounded-lg text-xs font-code">"Stop recording"</span>
                <span>Stops the recording</span>
              </li>
            </ul>
          </div>
        </details>
      </div>
    </div>
  );
};

export default RecordingControls;