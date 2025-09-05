import React from 'react';
import PropTypes from 'prop-types';

const StructureList = ({ structures, onSelect, selected, loading }) => {
  // Defensive check to ensure structures is an array
  const safeStructures = Array.isArray(structures) ? structures : [];

  // Loading state when no structures are present yet
  if (loading && safeStructures.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center py-8">
        <div className="animate-spin rounded-full h-8 w-8 border-t-2 border-b-2 border-accent"></div>
        <p className="mt-4 text-text-secondary font-body">Loading structures...</p>
      </div>
    );
  }

  // Empty state when no structures exist
  if (safeStructures.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center py-8 text-center">
        <svg
          className="w-12 h-12 text-text-secondary mb-4"
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
          xmlns="http://www.w3.org/2000/svg"
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            strokeWidth={1.5}
            d="M9.172 16.172a4 4 0 015.656 0M9 10h.01M15 10h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"
          />
        </svg>
        <h3 className="text-text-primary font-heading font-medium mb-1">No structures yet</h3>
        <p className="text-text-secondary font-body text-sm">Generate your first structure to get started</p>
      </div>
    );
  }

  return (
    <div className="max-h-[calc(100vh-300px)] overflow-y-auto pr-2">
      <ul className="space-y-2">
        {safeStructures.map((structure) => (
          <li
            key={structure._id || Math.random().toString(36).substr(2)}
            className={`p-3 rounded-lg transition-all duration-200 cursor-pointer ${
              selected?._id === structure._id
                ? 'bg-accent-secondary/20 border-l-4 border-accent-secondary'
                : 'hover:bg-secondary/70 border-l-4 border-transparent'
            }`}
            onClick={() => onSelect(structure)}
            aria-selected={selected?._id === structure._id}
          >
            <div className="flex items-start justify-between">
              <div className="flex-1 min-w-0">
                <p className="text-sm font-heading font-medium text-text-primary truncate">
                  {structure.name || 'Unnamed Structure'}
                </p>
                <p className="text-xs font-label text-text-secondary mt-1">
                  {structure.created
                    ? new Date(structure.created).toLocaleDateString('en-US', {
                        year: 'numeric',
                        month: 'short',
                        day: 'numeric',
                        hour: '2-digit',
                        minute: '2-digit',
                      })
                    : 'Date unknown'}
                </p>
              </div>
              <div className="ml-4 flex-shrink-0">
                <span className="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-heading font-medium bg-primary text-accent">
                  {structure.generatedStructures?.length || 0} variants
                </span>
              </div>
            </div>
            {structure.propertyName && (
              <div className="mt-2 flex items-center">
                <span className="text-xs font-label text-text-secondary mr-2">Property:</span>
                <span className="text-xs font-heading font-medium text-accent">
                  {structure.propertyName}
                  {structure.minimize ? ' (min)' : ' (max)'}
                </span>
              </div>
            )}
          </li>
        ))}
      </ul>
    </div>
  );
};

StructureList.propTypes = {
  structures: PropTypes.array.isRequired,
  onSelect: PropTypes.func.isRequired,
  selected: PropTypes.object,
  loading: PropTypes.bool.isRequired,
};

export default React.memo(StructureList);