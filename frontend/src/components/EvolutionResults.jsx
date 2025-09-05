import React, { useEffect, useRef } from 'react';

const EvolutionResults = ({ results, rdkitLoaded }) => {
  const molRefs = useRef({});

  useEffect(() => {
    if (rdkitLoaded && window.RDKit && results?.length) {
      results.forEach((result, index) => {
        const ref = molRefs.current[index];
        if (ref && result.smiles) {
          const mol = window.RDKit.get_mol(result.smiles);
          if (mol) {
            ref.innerHTML = mol.get_svg(160, 120);
            mol.delete();
          }
        }
      });
    }
  }, [rdkitLoaded, results]);

  if (!results || results.length === 0) {
    return <p className="text-gray-500">No evolution results available.</p>;
  }

  return (
    <div>
      <h3 className="text-lg font-semibold mb-4">Evolution Results</h3>
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
        {results.map((result, index) => (
          <div
            key={index}
            className="border border-gray-200 rounded-lg p-4 hover:shadow-md"
          >
            <div className="h-40 flex items-center justify-center">
              {rdkitLoaded && window.RDKit && result.smiles ? (
                <div ref={(el) => (molRefs.current[index] = el)} />
              ) : (
                <p className="text-gray-500 text-xs">Visualization unavailable</p>
              )}
            </div>
            <div className="mt-2">
              <p className="text-sm"><strong>SMILES:</strong> {result.smiles}</p>
              <p className="text-sm"><strong>Stability:</strong> {result.stability?.toFixed(3) || 'N/A'}</p>
              <p className="text-sm"><strong>Solubility:</strong> {result.solubility?.toFixed(3) || 'N/A'}</p>
              <p className="text-sm"><strong>Efficacy:</strong> {result.efficacy?.toFixed(3) || 'N/A'}</p>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};

export default EvolutionResults;