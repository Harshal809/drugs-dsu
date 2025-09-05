import React, { useState } from 'react';
import { Edit, Mic, Save, X, Tag, Download, Search, Calendar,MessagesSquare } from 'lucide-react';

const NotesList = ({
  notes,
  editNote,
  editingNoteId,
  updateNote,
  cancelEditing,
  startEditingWithRecording,
  addTag,
  exportNotes,
  tags,
  theme,
  currentNote,
  setCurrentNote,
  stopRecording,
  setShowContinue,
}) => {
  const [searchTerm, setSearchTerm] = useState('');
  const [selectedTag, setSelectedTag] = useState('');
  const [newTag, setNewTag] = useState('');

  // Filter notes based on search term and selected tag
  const filteredNotes = notes.filter(note => {
    const matchesSearch = note.content.toLowerCase().includes(searchTerm.toLowerCase());
    const matchesTag = selectedTag ? (note.tags && note.tags.includes(selectedTag)) : true;
    return matchesSearch && matchesTag;
  });

  const handleAddTag = (noteId) => {
    if (newTag.trim()) {
      addTag(noteId, newTag.trim());
      setNewTag('');
    }
  };

  const formatDate = (dateString) => {
    if (!dateString || isNaN(new Date(dateString))) {
      return 'Unknown Date';
    }
    const date = new Date(dateString);
    return new Intl.DateTimeFormat('en-US', {
      month: 'short',
      day: 'numeric',
      year: 'numeric',
      hour: '2-digit',
      minute: '2-digit'
    }).format(date);
  };

  return (
    <div className="p-6 md:p-8 bg-gradient-to-b from-primary/50 to-primary/30">
      <div className="flex flex-wrap justify-between items-center mb-8 gap-6">
        <h2 className="text-2xl font-heading font-semibold text-text-primary flex items-center gap-3">
          <div className="p-2 rounded-lg bg-accent-secondary/20 border border-accent-secondary/30">
            <Calendar size={24} className="text-accent-secondary" />
          </div>
          Your Notes
          <span className="text-sm font-body text-text-secondary bg-secondary/50 px-3 py-1 rounded-full border border-secondary/50">
            {filteredNotes.length} notes
          </span>
        </h2>

        <div className="flex flex-wrap items-center gap-4">
          <div className="relative flex items-center rounded-xl bg-secondary/50 border border-secondary/50 backdrop-blur-sm">
            <Search size={18} className="absolute left-4 text-text-secondary" />
            <input
              type="text"
              placeholder="Search notes..."
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              className="py-3 pl-12 pr-4 rounded-xl w-full sm:w-64 bg-transparent text-text-primary placeholder-text-secondary border-none focus:outline-none focus:ring-2 focus:ring-accent/50 font-body"
            />
          </div>

          <select
            value={selectedTag}
            onChange={(e) => setSelectedTag(e.target.value)}
            className="py-3 px-4 rounded-xl bg-secondary/50 border border-secondary/50 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent/50 backdrop-blur-sm"
          >
            <option value="">All Tags</option>
            {tags.map(tag => (
              <option key={tag} value={tag}>{tag}</option>
            ))}
          </select>

          <button
            onClick={() => exportNotes(filteredNotes)}
            className="flex items-center gap-2 py-3 px-6 rounded-xl bg-accent/20 hover:bg-accent/30 text-accent border border-accent/30 transition-all duration-300 font-medium font-body backdrop-blur-sm hover:scale-105"
            title="Export notes"
          >
            <Download size={18} />
            <span className="hidden sm:inline">Export</span>
          </button>
        </div>
      </div>

      {filteredNotes.length === 0 && (
        <div className="text-center py-16 text-text-secondary font-body">
          <div className="p-6 rounded-2xl bg-secondary/30 border border-secondary/50 inline-block backdrop-blur-sm">
            <MessagesSquare size={48} className="mx-auto mb-4 text-text-secondary/50" />
            <p className="text-lg">
              {searchTerm || selectedTag ? 'No matching notes found.' : 'No notes saved yet.'}
            </p>
            <p className="text-sm mt-2">Start recording to create your first note!</p>
          </div>
        </div>
      )}

      <div className="grid gap-6 sm:grid-cols-1 lg:grid-cols-2 xl:grid-cols-3">
        {filteredNotes.map((note) => (
          <div 
            key={note._id} 
            className={`rounded-2xl overflow-hidden border backdrop-blur-sm transition-all duration-300 hover:scale-105 ${
              editingNoteId === note._id 
                ? 'border-accent/50 ring-2 ring-accent/30 bg-secondary/60' 
                : 'border-secondary/50 hover:border-accent/30 bg-secondary/40 hover:bg-secondary/60'
            }`}
          >
            {editingNoteId === note._id ? (
              <div className="p-6">
                <textarea
                  value={currentNote}
                  onChange={(e) => setCurrentNote(e.target.value)}
                  className="w-full min-h-[150px] p-4 mb-4 bg-primary/50 border border-secondary/50 text-text-primary rounded-xl resize-y focus:outline-none focus:ring-2 focus:ring-accent/50 font-body leading-relaxed"
                ></textarea>

                <div className="flex items-center gap-3 mb-4">
                  <input
                    type="text"
                    placeholder="Add tag..."
                    value={newTag}
                    onChange={(e) => setNewTag(e.target.value)}
                    className="flex-1 p-3 text-sm bg-primary/50 border border-secondary/50 text-text-primary rounded-xl focus:outline-none focus:ring-2 focus:ring-accent/50 font-body"
                  />
                  <button
                    onClick={() => handleAddTag(note._id)}
                    className={`p-3 rounded-xl transition-all duration-300 ${
                      !newTag.trim()
                        ? 'opacity-50 cursor-not-allowed bg-secondary/30'
                        : 'bg-accent/20 text-accent hover:bg-accent/30 border border-accent/30'
                    }`}
                    disabled={!newTag.trim()}
                  >
                    <Tag size={16} />
                  </button>
                </div>

                <div className="flex gap-3">
                  <button
                    onClick={async () => {
                      stopRecording();
                      const success = await updateNote(note._id, currentNote);
                      if (success) {
                        cancelEditing();
                        setShowContinue(false);
                        localStorage.removeItem('voice_notes_draft');
                      }
                    }}
                    className="flex items-center gap-2 py-2 px-4 rounded-xl font-medium bg-success/20 text-success hover:bg-success/30 border border-success/30 transition-all duration-300 font-body"
                  >
                    <Save size={16} />
                    <span>Save</span>
                  </button>

                  <button
                    onClick={() => startEditingWithRecording(note._id)}
                    className="flex items-center gap-2 py-2 px-4 rounded-xl font-medium bg-accent-secondary/20 text-accent-secondary hover:bg-accent-secondary/30 border border-accent-secondary/30 transition-all duration-300 font-body"
                  >
                    <Mic size={16} />
                    <span>Record</span>
                  </button>

                  <button
                    onClick={cancelEditing}
                    className="flex items-center gap-2 py-2 px-4 rounded-xl font-medium bg-error/20 text-error hover:bg-error/30 border border-error/30 transition-all duration-300 font-body"
                  >
                    <X size={16} />
                    <span>Cancel</span>
                  </button>
                </div>
              </div>
            ) : (
              <div className="p-6 h-full flex flex-col">
                <div className="mb-4 flex-grow">
                  <p className="text-text-primary font-body leading-relaxed">
                    {note.content.length > 250 
                      ? `${note.content.substring(0, 250)}...` 
                      : note.content}
                  </p>
                </div>

                <div className="flex items-center gap-2 text-xs text-text-secondary mb-4 font-label">
                  <Calendar size={14} />
                  <span>{formatDate(note.createdAt)}</span>
                </div>

                {note.tags && note.tags.length > 0 && (
                  <div className="flex flex-wrap gap-2 mb-4">
                    {note.tags.map(tag => (
                      <span 
                        key={tag} 
                        className="text-xs px-3 py-1 rounded-full bg-accent/20 text-accent border border-accent/30 font-label"
                      >
                        {tag}
                      </span>
                    ))}
                  </div>
                )}

                <div className="flex gap-3 mt-auto pt-4 border-t border-secondary/30">
                  <button
                    onClick={() => editNote(note._id)}
                    className="flex-1 p-3 rounded-xl bg-accent-secondary/20 text-accent-secondary hover:bg-accent-secondary/30 border border-accent-secondary/30 transition-all duration-300 flex items-center justify-center gap-2 font-body"
                    title="Edit note"
                  >
                    <Edit size={16} />
                    <span>Edit</span>
                  </button>

                  <button
                    onClick={() => startEditingWithRecording(note._id)}
                    className="flex-1 p-3 rounded-xl bg-accent/20 text-accent hover:bg-accent/30 border border-accent/30 transition-all duration-300 flex items-center justify-center gap-2 font-body"
                    title="Continue recording"
                  >
                    <Mic size={16} />
                    <span>Record</span>
                  </button>
                </div>
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
};

export default NotesList;