import React, { useState } from 'react';

const MolecularEvolution = ({ onEvolve, loading, initialSmiles }) => {
  const [evolutionParams, setEvolutionParams] = useState({
    smiles: initialSmiles || '',
    iterations: 5,
    populationSize: 20,
    optimizeFor: ['stability', 'solubility', 'efficacy'],
    mutationRate: 0.1,
  });

  const handleChange = (e) => {
    const { name, value, type, checked } = e.target;
    if (type === 'checkbox') {
      setEvolutionParams((prev) => ({
        ...prev,
        optimizeFor: checked
          ? [...prev.optimizeFor, value]
          : prev.optimizeFor.filter((opt) => opt !== value),
      }));
    } else {
      setEvolutionParams((prev) => ({ ...prev, [name]: value }));
    }
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    onEvolve(evolutionParams);
  };

  return (
    <form onSubmit={handleSubmit} className="space-y-4">
      <div>
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Starting SMILES
        </label>
        <input
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700"
          name="smiles"
          value={evolutionParams.smiles}
          onChange={handleChange}
          placeholder="Enter SMILES string"
          required
        />
      </div>
      <div>
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Iterations
        </label>
        <input
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700"
          name="iterations"
          type="number"
          min="1"
          max="50"
          value={evolutionParams.iterations}
          onChange={handleChange}
        />
      </div>
      <div>
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Population Size
        </label>
        <input
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700"
          name="populationSize"
          type="number"
          min="10"
          max="100"
          value={evolutionParams.populationSize}
          onChange={handleChange}
        />
      </div>
      <div>
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Mutation Rate
        </label>
        <input
          className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700"
          name="mutationRate"
          type="number"
          min="0"
          max="1"
          step="0.01"
          value={evolutionParams.mutationRate}
          onChange={handleChange}
        />
      </div>
      <div>
        <label className="block text-gray-700 text-sm font-bold mb-2">
          Optimize For
        </label>
        {['stability', 'solubility', 'efficacy'].map((opt) => (
          <div key={opt} className="flex items-center">
            <input
              type="checkbox"
              name="optimizeFor"
              value={opt}
              checked={evolutionParams.optimizeFor.includes(opt)}
              onChange={handleChange}
              className="mr-2"
            />
            <span>{opt.charAt(0).toUpperCase() + opt.slice(1)}</span>
          </div>
        ))}
      </div>
      <button
        type="submit"
        className={`bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded w-full ${
          loading ? 'opacity-50 cursor-not-allowed' : ''
        }`}
        disabled={loading}
      >
        {loading ? 'Evolving...' : 'Start Evolution'}
      </button>
    </form>
  );
};

export default MolecularEvolution;