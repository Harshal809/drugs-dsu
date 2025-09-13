import React, { useState } from 'react';
import debounce from 'lodash/debounce';

const StructureForm = ({ onSubmit, loading }) => {
  const [formData, setFormData] = useState({
    name: '',
    smiles: '',
    algorithm: 'CMA-ES',
    numMolecules: 30,
    propertyName: 'QED',
    minimize: false,
    minSimilarity: 0.3,
    particles: 30,
    iterations: 10,
  });

  const [suggestions, setSuggestions] = useState([]);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const [isFetchingSmiles, setIsFetchingSmiles] = useState(false);

  const fetchSuggestions = async (query) => {
    if (!query) {
      setSuggestions([]);
      setShowSuggestions(false);
      return;
    }

    try {
      const response = await fetch(
        `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(query)}/json`
      );
      const data = await response.json();
      if (data.dictionary_terms && data.dictionary_terms.compound) {
        setSuggestions(data.dictionary_terms.compound.slice(0, 5));
        setShowSuggestions(true);
      } else {
        setSuggestions([]);
        setShowSuggestions(false);
      }
    } catch (error) {
      console.error('Error fetching suggestions:', error);
      setSuggestions([]);
      setShowSuggestions(false);
    }
  };

  const debouncedFetchSuggestions = debounce(fetchSuggestions, 300);

  const handleNameChange = (e) => {
    const value = e.target.value;
    setFormData({ ...formData, name: value });
    debouncedFetchSuggestions(value);
  };

  const handleSuggestionSelect = async (selectedName) => {
  console.log('Selected compound:', selectedName);
  setFormData(prev => ({ ...prev, name: selectedName }));
  
  try {
    setIsFetchingSmiles(true);
    const apiUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(selectedName)}/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON`;
    console.log('Fetching SMILES from:', apiUrl);
    
    const response = await fetch(apiUrl);
    
    if (!response.ok) {
      throw new Error(`PubChem API error: ${response.status}`);
    }
    
    const data = await response.json();
    console.log('PubChem response:', data);
    
    if (data.PropertyTable?.Properties?.[0]) {
      const props = data.PropertyTable.Properties[0];
      // Try to get SMILES from different possible fields
      const smiles = props.CanonicalSMILES || props.ConnectivitySMILES || props.IsomericSMILES || '';
      
      if (smiles) {
        console.log('Found SMILES:', smiles);
        setFormData(prev => ({
          ...prev,
          smiles: smiles
        }));
      } else {
        console.warn('No SMILES found in response');
      }
    } else {
      console.warn('No properties found in response');
    }
  } catch (error) {
    console.error('Error fetching SMILES:', error);
  } finally {
    setIsFetchingSmiles(false);
    setShowSuggestions(false);
  }
};

  const handleChange = (e) => {
    const { name, value, type, checked } = e.target;
    setFormData({
      ...formData,
      [name]: type === 'checkbox' ? checked : value,
    });
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    onSubmit(formData);
  };


  return (
    <form onSubmit={handleSubmit} className="space-y-6">
      <div className="relative">
        <label className="block text-text-primary font-label font-semibold mb-2">
          Structure Name
        </label>
        <div className="relative">
          <input
            className="w-full bg-primary border border-accent-secondary rounded-lg py-3 px-4 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent"
            name="name"
            type="text"
            placeholder="e.g., Ibuprofen"
            value={formData.name}
            onChange={handleNameChange}
            onFocus={() => setShowSuggestions(true)}
            onBlur={() => setTimeout(() => setShowSuggestions(false), 200)}
            autoComplete="off"
          />
          {isFetchingSmiles && (
            <div className="absolute right-3 top-3">
              <div className="animate-spin rounded-full h-5 w-5 border-t-2 border-b-2 border-accent"></div>
            </div>
          )}
        </div>
        {showSuggestions && suggestions.length > 0 && (
          <ul className="absolute z-10 mt-1 w-full bg-primary border border-accent-secondary rounded-lg shadow-lg max-h-60 overflow-auto">
            {suggestions.map((suggestion) => (
              <li
                key={suggestion}
                className="px-4 py-3 hover:bg-secondary text-text-primary font-body cursor-pointer transition-colors duration-150"
                onMouseDown={(e) => e.preventDefault()} // Prevent blur before click
                onClick={() => handleSuggestionSelect(suggestion)}
              >
                {suggestion}
              </li>
            ))}
          </ul>
        )}
      </div>

      <div>
        <label className="block text-text-primary font-label font-semibold mb-2">
          SMILES String
        </label>
        <textarea
          className="w-full bg-primary border border-accent-secondary rounded-lg py-3 px-4 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent"
          name="smiles"
          placeholder={isFetchingSmiles ? "Fetching SMILES..." : "SMILES will be auto-filled after selecting a compound"}
          value={formData.smiles}
          onChange={handleChange}
          rows="3"
          required
          disabled={isFetchingSmiles}
        />
        <p className="text-xs text-text-secondary mt-2">
          SMILES will be fetched from PubChem or enter manually
        </p>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <div>
          <label className="block text-text-primary font-label font-semibold mb-2">
            Algorithm
          </label>
          <select
            className="w-full bg-primary border border-accent-secondary rounded-lg py-3 px-4 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent appearance-none"
            name="algorithm"
            value={formData.algorithm}
            onChange={handleChange}
          >
            <option value="CMA-ES">CMA-ES</option>
            <option value="SSD">Sampling Standard Deviation</option>
          </select>
        </div>

        <div>
          <label className="block text-text-primary font-label font-semibold mb-2">
            Property
          </label>
          <select
            className="w-full bg-primary border border-accent-secondary rounded-lg py-3 px-4 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent appearance-none"
            name="propertyName"
            value={formData.propertyName}
            onChange={handleChange}
          >
            <option value="QED">QED</option>
            <option value="logP">logP</option>
          </select>
        </div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <div>
          <label className="block text-text-primary font-label font-semibold mb-2">
            Number of Molecules
          </label>
          <input
            className="w-full bg-primary border border-accent-secondary rounded-lg py-3 px-4 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent"
            name="numMolecules"
            type="number"
            min="1"
            max="100"
            value={formData.numMolecules}
            onChange={handleChange}
          />
        </div>

        <div>
          <label className="block text-text-primary font-label font-semibold mb-2">
            Min Similarity
          </label>
          <div className="relative">
            <input
              className="w-full bg-primary border border-accent-secondary rounded-lg py-3 px-4 text-text-primary font-body focus:outline-none focus:ring-2 focus:ring-accent"
              name="minSimilarity"
              type="number"
              min="0"
              max="1"
              step="0.1"
              value={formData.minSimilarity}
              onChange={handleChange}
            />
            <div className="absolute inset-y-0 right-0 flex items-center pr-3 pointer-events-none text-text-secondary">
              {formData.minSimilarity}
            </div>
          </div>
        </div>
      </div>

      <div className="flex items-center">
        <input
          id="minimize"
          name="minimize"
          type="checkbox"
          className="h-5 w-5 text-accent rounded focus:ring-accent border-accent-secondary bg-primary"
          checked={formData.minimize}
          onChange={handleChange}
        />
        <label htmlFor="minimize" className="ml-3 block text-sm text-text-primary font-body">
          Minimize property (instead of maximize)
        </label>
      </div>

      <button
        className={`w-full bg-accent hover:bg-accent/90 text-primary font-heading font-bold py-3 px-6 rounded-lg transition-all duration-200 ${
          loading ? 'opacity-70 cursor-not-allowed' : 'hover:shadow-lg'
        }`}
        type="submit"
        disabled={loading}
      >
        {loading ? (
          <span className="flex items-center justify-center">
            <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-primary" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
            </svg>
            Generating...
          </span>
        ) : (
          'Generate Structure'
        )}
      </button>
    </form>
  );
};

export default StructureForm;