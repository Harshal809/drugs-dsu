import React, { useState, useEffect, useRef } from 'react';
import { jsPDF } from 'jspdf';
import { Canvg } from 'canvg';
import axios from 'axios';
import { ClipLoader } from 'react-spinners';
import EnhancedGeminiInfo from './Enhancedgeminiinfo.jsx';

const StructureDetails = ({ structure, rdkitLoaded }) => {
  const [activeTab, setActiveTab] = useState('parent');
  const [selectedVariant, setSelectedVariant] = useState(null);
  const [variantInfo, setVariantInfo] = useState(null);
  const [loadingVariantInfo, setLoadingVariantInfo] = useState(false);
  const [variantInfoError, setVariantInfoError] = useState(null);
  const [showParent3D, setShowParent3D] = useState(false);
  const [showVariant3D, setShowVariant3D] = useState(false);
  const [is3DmolLoaded, setIs3DmolLoaded] = useState(false);
  const [parentStyle, setParentStyle] = useState('ballstick');
  const [variantStyle, setVariantStyle] = useState('ballstick');
  const [showLabels, setShowLabels] = useState(false);
  const [aiDescription, setAiDescription] = useState('');
  const [spinAnimation, setSpinAnimation] = useState(false);
  const [highlightedAtom, setHighlightedAtom] = useState(null);
  const parentMolRef = useRef(null);
  const variantMolRef = useRef(null);
  const variantMolRefs = useRef({});
  const parent3DRef = useRef(null);
  const variant3DRef = useRef(null);
  const parentViewerRef = useRef(null);
  const variantViewerRef = useRef(null);
  const animationFrameId = useRef(null);

  useEffect(() => {
    if (showParent3D || showVariant3D) {
      if (!window.$3Dmol && !is3DmolLoaded) {
        const script = document.createElement('script');
        script.src = 'https://3dmol.csb.pitt.edu/build/3Dmol-min.js';
        script.async = true;
        script.onload = () => {
          console.log('3Dmol.js loaded successfully');
          setIs3DmolLoaded(true);
        };
        script.onerror = () => {
          console.error('Failed to load 3Dmol.js');
          setIs3DmolLoaded(false);
        };
        document.head.appendChild(script);
      } else if (window.$3Dmol) {
        setIs3DmolLoaded(true);
      }
    }
  }, [showParent3D, showVariant3D]);

  useEffect(() => {
    // Only reset variant-related state if we're not on the variant tab
    if (activeTab !== 'variant') {
      setSelectedVariant(null);
      setVariantInfo(null);
      setVariantInfoError(null);
    }
    setShowParent3D(false);
    setShowVariant3D(false);
    setAiDescription('');
    setSpinAnimation(false);
    setHighlightedAtom(null);
    if (parent3DRef.current) parent3DRef.current.innerHTML = '';
    if (variant3DRef.current) variant3DRef.current.innerHTML = '';
    if (animationFrameId.current) {
      cancelAnimationFrame(animationFrameId.current);
      animationFrameId.current = null;
    }
  }, [structure, activeTab]);

  useEffect(() => {
    if (activeTab === 'parent' && structure?.smiles) {
      console.log('Fetching AI enhancements for parent:', structure.smiles);
      fetchAiEnhancements(structure.smiles);
    } else if (activeTab === 'variant' && selectedVariant?.smiles) {
      console.log('Fetching AI enhancements for variant:', selectedVariant.smiles);
      fetchAiEnhancements(selectedVariant.smiles);
    }
  }, [activeTab, structure, selectedVariant]);

  useEffect(() => {
    if (rdkitLoaded && window.RDKit && structure?.smiles && activeTab === 'parent') {
      console.log('Rendering parent molecule:', structure.smiles);
      const mol = window.RDKit.get_mol(structure.smiles);
      if (mol && parentMolRef.current) {
        parentMolRef.current.innerHTML = mol.get_svg(600, 400);
        mol.delete();
      } else {
        console.error('Failed to render parent molecule');
      }
    }
  }, [rdkitLoaded, structure, activeTab]);

  useEffect(() => {
    if (rdkitLoaded && window.RDKit && selectedVariant?.smiles && activeTab === 'variant') {
      console.log('Rendering variant molecule:', selectedVariant.smiles);
      const variantMol = window.RDKit.get_mol(selectedVariant.smiles);
      if (variantMol && variantMolRef.current) {
        variantMolRef.current.innerHTML = variantMol.get_svg(600, 400);
        variantMol.delete();
      } else {
        console.error('Failed to render variant molecule');
      }
    }
  }, [rdkitLoaded, selectedVariant, activeTab]);

  useEffect(() => {
    if (rdkitLoaded && window.RDKit && activeTab === 'variants' && structure?.generatedStructures) {
      console.log('Rendering variant thumbnails');
      structure.generatedStructures.forEach((variant, index) => {
        const ref = variantMolRefs.current[index];
        if (ref && variant.smiles) {
          const mol = window.RDKit.get_mol(variant.smiles);
          if (mol) {
            ref.innerHTML = mol.get_svg(160, 120);
            mol.delete();
          }
        }
      });
    }
  }, [rdkitLoaded, structure, activeTab]);

  const fetchAiEnhancements = async (smiles) => {
    const GEMINI_API_KEY = import.meta.env.VITE_GEMINI_API_KEY;
    if (!GEMINI_API_KEY) {
      console.error('Gemini API key not set');
      setAiDescription('Gemini API key not set');
      return;
    }

    try {
      const response = await axios.post(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=${GEMINI_API_KEY}`,
        { contents: [{ parts: [{ text: `Generate a creative, concise description of potential applications for a molecule with SMILES "${smiles}".` }] }] },
        { headers: { 'Content-Type': 'application/json' }, withCredentials: false }
      );
      if (response.data.candidates?.[0]?.content?.parts?.[0]?.text) {
        setAiDescription(response.data.candidates[0].content.parts[0].text);
      } else {
        setAiDescription('No AI description available');
      }
    } catch (error) {
      console.error('Failed to fetch AI enhancements:', error);
      setAiDescription('Failed to fetch AI description');
    }
  };

  const handleAtomClick = (viewer, atom) => {
    if (highlightedAtom && highlightedAtom.index === atom.index) {
      viewer.setStyle({ index: atom.index }, { sphere: { scale: 0.3 }, stick: { radius: 0.2 } });
      setHighlightedAtom(null);
    } else {
      viewer.setStyle({ index: atom.index }, { sphere: { scale: 0.5, color: 'yellow' }, stick: { radius: 0.3, color: 'yellow' } });
      setHighlightedAtom({ element: atom.elem, x: atom.x, y: atom.y, z: atom.z, index: atom.index });
    }
    viewer.render();
  };

  const spinModel = (viewer) => {
    if (!spinAnimation || !viewer) return;
    const animate = () => {
      viewer.rotate(2, 'y');
      viewer.render();
      animationFrameId.current = requestAnimationFrame(animate);
    };
    animationFrameId.current = requestAnimationFrame(animate);
  };

  useEffect(() => {
    if (showParent3D && parent3DRef.current && structure?.smiles && is3DmolLoaded && window.$3Dmol) {
      parent3DRef.current.innerHTML = '';
      const mol = window.RDKit.get_mol(structure.smiles);
      const molData = mol.get_molblock();
      mol.delete();

      const viewer = window.$3Dmol.createViewer(parent3DRef.current, { backgroundColor: 'white' });
      if (!viewer) {
        console.error('Failed to create parent 3D viewer');
        return;
      }
      parentViewerRef.current = viewer;
      viewer.addModel(molData, 'mol');

      const styleConfig = parentStyle === 'ballstick' ? { stick: { radius: 0.2 }, sphere: { scale: 0.3 } } :
        parentStyle === 'stick' ? { stick: { radius: 0.2 } } :
          parentStyle === 'sphere' ? { sphere: { scale: 0.5 } } :
            { surface: { opacity: 0.9, color: 'grey' } };
      viewer.setStyle({}, styleConfig);

      if (showLabels) {
        const model = viewer.getModel(0);
        const atoms = model.selectedAtoms({});
        atoms.forEach((atom, index) => {
          viewer.addLabel(
            `${atom.elem}${index + 1}`,
            { fontSize: 10, backgroundColor: 'black', backgroundOpacity: 0.8 },
            { x: atom.x, y: atom.y, z: atom.z }
          );
        });
      }

      viewer.setClickable({}, true, (atom) => handleAtomClick(viewer, atom));
      viewer.zoomTo();
      viewer.render();

      if (spinAnimation) {
        if (animationFrameId.current) cancelAnimationFrame(animationFrameId.current);
        spinModel(viewer);
      } else if (animationFrameId.current) {
        cancelAnimationFrame(animationFrameId.current);
        animationFrameId.current = null;
      }

      return () => {
        if (animationFrameId.current) {
          cancelAnimationFrame(animationFrameId.current);
          animationFrameId.current = null;
        }
      };
    }
  }, [showParent3D, structure, is3DmolLoaded, parentStyle, showLabels, spinAnimation]);

  useEffect(() => {
    if (showVariant3D && variant3DRef.current && selectedVariant?.smiles && is3DmolLoaded && window.$3Dmol) {
      variant3DRef.current.innerHTML = '';
      const mol = window.RDKit.get_mol(selectedVariant.smiles);
      const molData = mol.get_molblock();
      mol.delete();

      const viewer = window.$3Dmol.createViewer(variant3DRef.current, { backgroundColor: 'white' });
      if (!viewer) {
        console.error('Failed to create variant 3D viewer');
        return;
      }
      variantViewerRef.current = viewer;
      viewer.addModel(molData, 'mol');

      const styleConfig = variantStyle === 'ballstick' ? { stick: { radius: 0.2 }, sphere: { scale: 0.3 } } :
        variantStyle === 'stick' ? { stick: { radius: 0.2 } } :
          variantStyle === 'sphere' ? { sphere: { scale: 0.5 } } :
            { surface: { opacity: 0.9, color: 'grey' } };
      viewer.setStyle({}, styleConfig);

      if (showLabels) {
        const model = viewer.getModel(0);
        const atoms = model.selectedAtoms({});
        atoms.forEach((atom, index) => {
          viewer.addLabel(
            `${atom.elem}${index + 1}`,
            { fontSize: 10, backgroundColor: 'black', backgroundOpacity: 0.8 },
            { x: atom.x, y: atom.y, z: atom.z }
          );
        });
      }

      viewer.setClickable({}, true, (atom) => handleAtomClick(viewer, atom));
      viewer.zoomTo();
      viewer.render();

      if (spinAnimation) {
        if (animationFrameId.current) cancelAnimationFrame(animationFrameId.current);
        spinModel(viewer);
      } else if (animationFrameId.current) {
        cancelAnimationFrame(animationFrameId.current);
        animationFrameId.current = null;
      }

      return () => {
        if (animationFrameId.current) {
          cancelAnimationFrame(animationFrameId.current);
          animationFrameId.current = null;
        }
      };
    }
  }, [showVariant3D, selectedVariant, is3DmolLoaded, variantStyle, showLabels, spinAnimation]);

  const download3DView = (viewerRef, filename) => {
    if (viewerRef.current) {
      viewerRef.current.render(() => {
        const canvas = (parent3DRef.current || variant3DRef.current).querySelector('canvas');
        if (canvas) {
          const pngData = canvas.toDataURL('image/png');
          const link = document.createElement('a');
          link.href = pngData;
          link.download = `${filename}.png`;
          link.click();
        } else {
          console.error('Canvas not found for download');
          alert('Failed to download image: Canvas not found');
        }
      });
    }
  };

  const fetchVariantInfo = async (smiles, retries = 3) => {
    try {
      setLoadingVariantInfo(true);
      setVariantInfoError(null);
      const GEMINI_API_KEY = 'AIzaSyABAk3NKj4eqXhP9SkZIX8lYhfiZINP6ZQ';
      if (!GEMINI_API_KEY) throw new Error('Gemini API key not set');
      const response = await axios.post(
        `https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash-latest:generateContent?key=${GEMINI_API_KEY}`,
        {
          contents: [{
            parts: [{
              text: ` 

You are a pharmaceutical chemistry expert tasked with analyzing a drug molecule represented by the SMILES string "${smiles}". Provide a detailed report with the following sections, ensuring clarity, conciseness, and relevance. Each section should be 3-5 sentences long, unless specified otherwise, and include specific, detailed information, even if it requires predictions or estimations based on computational methods, structural analogies, or literature insights.

1. Structural Details:
- Determine the molecular formula from the SMILES string (e.g., C9H8O4 for aspirin).
- List any known synonyms or alternative names for the compound, citing databases like PubChem with specific identifiers (e.g., PubChem CID: 2244 for aspirin).
- Calculate the exact molecular weight and provide it with units (e.g., 180.16 g/mol).
- Generate the IUPAC name and InChI based on the molecular structure; if exact structure confirmation is unavailable, provide a predicted name and InChI based on the SMILES string.
- Note the provided SMILES string for reference (e.g., SMILES: CC(=O)OC1=CC=CC=C1C(=O)O).

2. Chemical Properties:
- Predict or provide the lipophilicity (logP) value, citing sources if available (e.g., logP: 1.2, predicted via ChemAxon); if unavailable, estimate based on functional groups.
- Provide pKa values for each ionizable group (e.g., -COOH, -NH2) with units (e.g., pKa: 3.5 for carboxylic acid); predict values if experimental data is unavailable.
- Count and list the number of hydrogen bond donors (HBD) and acceptors (HBA) (e.g., 2 HBD, 4 HBA).
- Identify and list significant functional groups present in the molecule (e.g., carboxylic acid, ester, aromatic ring).
- Determine the number of aromatic and non-aromatic rings, specifying each (e.g., 1 aromatic benzene ring, 0 non-aromatic rings).

3. Physical Properties:
- Predict solubility in aqueous and organic solvents, explaining the basis (e.g., moderately soluble in water due to polar groups, highly soluble in ethanol due to logP of 1.2).
- Discuss membrane permeability, relating it to logP and polar surface area (e.g., moderate permeability with logP of 1.2 and PSA of 60 Å²).
- Calculate the polar surface area (PSA) and provide it with units (e.g., PSA: 60 Å², calculated via ChemAxon).
- Count the number of rotatable bonds in the molecular structure (e.g., 3 rotatable bonds).
- Discuss chemical stability, noting any known issues (e.g., stable at neutral pH, hydrolyzes in acidic conditions), and mention crystallinity and polymorphism if data is available or predict based on structure (e.g., likely crystalline due to planar aromatic system).

4. Spectral Information:
- Provide detailed spectral data such as IR, NMR, or mass spectrometry, based on typical expectations for the structure (e.g., IR: 1700 cm⁻¹ for C=O stretch; NMR: 12 ppm for -COOH proton; MS: m/z 181 [M+H]+).
- If experimental data is unavailable, predict values or provide typical ranges based on functional groups, citing prediction methods (e.g., predicted via ChemAxon’s spectral prediction tools).

5. Biological Activity and Pharmacology:
- Identify potential biological targets (e.g., receptors, enzymes like COX-1) and predict binding affinities (e.g., 7.5 kcal/mol), including mechanism of action (e.g., inhibits COX-1 by blocking active site).
- Discuss structure-activity relationship (SAR) data, if known, or predict based on structural features (e.g., the -COOH group enhances COX inhibition).
- Suggest potential therapeutic applications based on structural similarities to known drugs (e.g., similar to aspirin, may treat pain and inflammation).
- Describe pharmacodynamic properties (e.g., inhibits prostaglandin synthesis) and pharmacokinetic properties (e.g., half-life of 4 hours, metabolized by CYP2C9, excreted via kidneys), including interactions and pathways (e.g., interacts with CYP2C9, affects prostaglandin pathway).
- Include any biological test results, such as assay data, or predict likely outcomes (e.g., likely IC50 of 10 µM against COX-1 based on structural analogy to aspirin).

6. Clinical and Therapeutic Information:
- List approved indications, dosage forms, and therapeutic uses, citing clinical sources if applicable (e.g., approved for pain relief, available as 325 mg tablets, per DrugBank DB00945).
- Summarize any clinical trial results, ongoing studies, or comparisons with similar drugs (e.g., reduces pain by 50% in trials, comparable to ibuprofen, per ClinicalTrials.gov NCT123456).
- Identify associated disorders and diseases, explaining the connection (e.g., treats arthritis due to anti-inflammatory effects).
- If clinical data is unavailable, predict potential uses based on biological activity (e.g., likely used for inflammation based on COX inhibition).

7. Safety and Toxicity:
- Discuss handling precautions and safety hazards, citing safety data sheets if possible (e.g., avoid inhalation, may cause skin irritation, per PubChem safety data).
- Provide toxicological information, including LD50 values (e.g., LD50: 200 mg/kg in rats, predicted via QSAR models) and any known or predicted toxicological concerns (e.g., potential gastrointestinal irritation).
- Highlight potential risks, such as hERG inhibition or cytotoxicity (e.g., possible hERG inhibition due to aromatic system, predicted via computational models).

8. Manufacturing and Synthesis:
- Outline a possible synthetic route to obtain the compound, referencing established methods in literature (e.g., acetylate salicylic acid with acetic anhydride, per Organic Syntheses Vol. 1).
- Highlight any challenges or considerations in the synthesis process, including manufacturing details (e.g., requires controlled temperature to avoid side products, scalable for industrial production).

9. References and Patents:
- Cite relevant scientific literature or databases for any data provided, including specific identifiers (e.g., PubChem CID: 2244, DrugBank DB00945, ChEMBL CHEMBL25).
- Note any patents related to the compound, including patent numbers if known (e.g., US Patent 1234567 for synthesis method).

10. Classification and Taxonomy:
- Classify the compound in relevant chemical or biological taxonomies, such as drug classes or chemical categories (e.g., belongs to salicylates, non-steroidal anti-inflammatory drugs).
- Provide standard classifications, noting any regulatory categories if applicable (e.g., FDA-approved analgesic, per PubChem).

**Output Requirements:**
- Use numbered headings for each section (e.g., "1. Structural Details:", "2. Chemical Properties:", etc.).
- Use dashes ("-") for all bullet points, and do not use stars ("*"), bullet symbols ("•"), or other markers.
- Do not use hashtags ("#") or stars ("*") anywhere in the response body.
- Keep each section concise (3-5 sentences, unless specified otherwise) and avoid extraneous information.
- Ensure all numerical values are accompanied by their units (e.g., g/mol for molecular weight, Å² for PSA, cm⁻¹ for IR peaks).
- Cite relevant databases (e.g., PubChem, ChEMBL, DrugBank) or prediction tools (e.g., ChemAxon, QSAR models) for additional context or data sources, if applicable.
- Ensure all text is plain and free of markdown formatting (e.g., no bold, italics, or tables) except for the specified headings and bullet points.
- If exact data is unavailable, provide a prediction or estimation based on computational methods, structural analogies, or literature insights, and note the method used (e.g., "predicted via ChemAxon"); only use "No data available" if the information is completely inapplicable.


` }]
          }]
        },
        { headers: { 'Content-Type': 'application/json' }, withCredentials: false }
      );
      if (response.data.candidates?.[0]?.content?.parts?.[0]?.text) {
        setVariantInfo(response.data.candidates[0].content.parts[0].text);
      } else {
        throw new Error('No content from Gemini API');
      }
    } catch (error) {
      if (retries > 0) {
        console.warn(`Retrying fetchVariantInfo (${retries} attempts left)...`);
        await new Promise((resolve) => setTimeout(resolve, 1000));
        return fetchVariantInfo(smiles, retries - 1);
      }
      setVariantInfoError('Failed to fetch variant info after multiple attempts: ' + error.message);
      setVariantInfo(null);
    } finally {
      setLoadingVariantInfo(false);
    }
  };

  const handleVariantSelect = (variant) => {
    console.log('Selected variant:', variant);
    setSelectedVariant(variant);
    setActiveTab('variant');
    setVariantInfo(null);
    setVariantInfoError(null);
    setShowVariant3D(false);
    if (variant.smiles) {
      fetchVariantInfo(variant.smiles);
      fetchAiEnhancements(variant.smiles);
    }
  };

  const handleTabSwitch = (tab) => {
    console.log('Switching to tab:', tab);
    setActiveTab(tab);
    setShowParent3D(false);
    setShowVariant3D(false);
    setHighlightedAtom(null);
    setSpinAnimation(false);
    if (parent3DRef.current) parent3DRef.current.innerHTML = '';
    if (variant3DRef.current) variant3DRef.current.innerHTML = '';
    if (animationFrameId.current) {
      cancelAnimationFrame(animationFrameId.current);
      animationFrameId.current = null;
    }
  };

  const renderMolecule = (smiles, ref) => {
    if (!rdkitLoaded || !smiles || !window.RDKit) {
      return (
        <div className="border-2 border-dashed border-gray-300 rounded-lg p-8 flex items-center justify-center">
          <p className="text-gray-500">Molecule visualization unavailable</p>
        </div>
      );
    }
    return <div ref={ref} className="flex justify-center" />;
  };



  const exportToPDF = async () => {
    const doc = new jsPDF();
    const pageWidth = doc.internal.pageSize.getWidth();
    const pageHeight = doc.internal.pageSize.getHeight();
    const margin = 10;
    let y = margin;
    const lineHeight = 3;
  
    const checkPageBreak = (additionalHeight) => {
      if (y + additionalHeight > pageHeight - margin) {
        doc.addPage();
        y = margin;
      }
    };
  
    const addText = (text, x, size = 12) => {
      doc.setFontSize(size);
      const lines = doc.splitTextToSize(text, pageWidth - 2 * margin);
      const textHeight = lines.length * size * 0.5;
      checkPageBreak(textHeight);
      doc.text(lines, x, y);
      y += textHeight + lineHeight;
      return y;
    };
  
    const addSvgToPDF = async (svgString, x, yPos, width, height) => {
      if (!svgString) return yPos;
      checkPageBreak(height + 10);
      const canvas = document.createElement('canvas');
      canvas.width = width * 2;
      canvas.height = height * 2;
      const ctx = canvas.getContext('2d');
      const v = Canvg.fromString(ctx, svgString);
      await v.render();
      const imgData = canvas.toDataURL('image/png');
      doc.addImage(imgData, 'PNG', x, yPos, width, height);
      return yPos + height + 5;
    };
  
    // Parse the information into sections
    const parseInformation = (text) => {
      const expectedHeadings = [
        "Structural Details",
        "Chemical Properties",
        "Physical Properties",
        "Spectral Information",
        "Biological Activity and Pharmacology",
        "Clinical and Therapeutic Information",
        "Safety and Toxicity",
        "Manufacturing and Synthesis",
        "References and Patents",
        "Classification and Taxonomy"
      ];
  
      // Regex to match numbered headings like "1. Structural Details:", "2. Chemical Properties:", etc.
      const sectionRegex = /(\d+\.\s[A-Za-z\s]+:)/g;
      const sections = [];
      let lastIndex = 0;
      let match;
  
      // Find all matches of the section headings
      const matches = [...text.matchAll(sectionRegex)];
      matches.forEach((m, index) => {
        const startIndex = m.index;
        const nextMatch = matches[index + 1];
        const endIndex = nextMatch ? nextMatch.index : text.length;
        const content = text.slice(lastIndex, startIndex).trim();
        if (content && lastIndex !== startIndex) {
          sections.push({
            id: sections.length,
            title: expectedHeadings[sections.length] || `Section ${sections.length + 1}`,
            content: content
          });
        }
        lastIndex = startIndex;
      });
  
      // Add the last section
      const lastContent = text.slice(lastIndex).trim();
      if (lastContent) {
        const titleMatch = lastContent.match(/^(\d+\.\s[A-Za-z\s]+:)/);
        let title = titleMatch ? titleMatch[1].replace(/:$/, '') : expectedHeadings[sections.length] || `Section ${sections.length + 1}`;
        title = title.replace(/^\d+\.\s*/, '').trim();
        const content = titleMatch ? lastContent.slice(titleMatch[0].length).trim() : lastContent;
        sections.push({
          id: sections.length,
          title: title,
          content: content
        });
      }
  
      return sections.map((section, index) => ({
        id: index,
        title: expectedHeadings[index] || section.title,
        content: section.content.trim()
      }));
    };
  
    // Add document header
    addText(`${structure.name} Details`, margin, 16);
    addText(`Created: ${new Date(structure.created).toLocaleString()}`, margin, 10);
  
    if (activeTab === 'parent') {
      addText('Parent Structure', margin, 14);
      addText(`SMILES: ${structure.smiles}`, margin);
  
      if (rdkitLoaded && window.RDKit && structure.smiles) {
        const mol = window.RDKit.get_mol(structure.smiles);
        if (mol) {
          const svg = mol.get_svg(300, 200);
          y = await addSvgToPDF(svg, margin, y, 90, 60);
          mol.delete();
        }
      }
  
      if (structure.properties) {
        addText('Properties:', margin, 12);
        addText(`QED: ${structure.properties?.qed?.toFixed(3) || 'N/A'}`, margin + 5, 10);
        addText(`LogP: ${structure.properties?.logp?.toFixed(3) || 'N/A'}`, margin + 5, 10);
      }
  
      // Add detailed information sections
      if (structure.information && structure.information.trim() !== '') {
        const sections = parseInformation(structure.information);
        sections.forEach((section) => {
          addText(section.title, margin, 12);
          const paragraphs = section.content.split('\n\n');
          paragraphs.forEach((paragraph) => {
            const lines = paragraph.split('\n').map(line => {
              // The prompt already uses dashes, so no need to replace bullet points
              return line;
            }).join('\n');
            y = addText(lines, margin + 5, 10);
          });
        });
      } else {
        addText('Detailed Information:', margin, 12);
        y = addText('No detailed information available.', margin + 5, 10);
      }
  
      if (aiDescription && aiDescription.trim() !== '') {
        addText('AI Insight:', margin, 12);
        y = addText(aiDescription, margin + 5, 10);
      }
  
    } else if (activeTab === 'variant' && selectedVariant) {
      addText('Selected Variant', margin, 14);
      addText(`SMILES: ${selectedVariant.smiles}`, margin);
  
      if (rdkitLoaded && window.RDKit && selectedVariant.smiles) {
        const mol = window.RDKit.get_mol(selectedVariant.smiles);
        if (mol) {
          const svg = mol.get_svg(300, 200);
          y = await addSvgToPDF(svg, margin, y, 90, 60);
          mol.delete();
        }
      }
  
      addText('Property Comparison:', margin, 12);
      addText(`Parent QED: ${structure.properties?.qed?.toFixed(3) || 'N/A'}`, margin + 5, 10);
      addText(`Parent LogP: ${structure.properties?.logp?.toFixed(3) || 'N/A'}`, margin + 5, 10);
      addText(`Variant QED: ${selectedVariant.properties?.qed?.toFixed(3) || 'N/A'}`, margin + 5, 10);
      addText(`Variant LogP: ${selectedVariant.properties?.logp?.toFixed(3) || 'N/A'}`, margin + 5, 10);
      addText(`Similarity: ${(selectedVariant.similarity * 100).toFixed(1)}%`, margin + 5, 10);
  
      if (variantInfo && variantInfo.trim() !== '') {
        const sections = parseInformation(variantInfo);
        sections.forEach((section) => {
          addText(section.title, margin, 12);
          const paragraphs = section.content.split('\n\n');
          paragraphs.forEach((paragraph) => {
            const lines = paragraph.split('\n').map(line => {
              return line;
            }).join('\n');
            y = addText(lines, margin + 5, 10);
          });
        });
      } else {
        addText('Detailed Information:', margin, 12);
        y = addText(variantInfoError || 'No detailed information available.', margin + 5, 10);
      }
  
      if (aiDescription && aiDescription.trim() !== '') {
        addText('AI Insight:', margin, 12);
        y = addText(aiDescription, margin + 5, 10);
      }
  
    } else if (activeTab === 'variants') {
      addText('Generated Variants', margin, 14);
      if (structure.generatedStructures) {
        for (let i = 0; i < structure.generatedStructures.length; i++) {
          const variant = structure.generatedStructures[i];
          addText(`Variant ${i + 1}:`, margin, 12);
          addText(`SMILES: ${variant.smiles}`, margin + 5);
          addText(`QED: ${variant.properties?.qed?.toFixed(3) || 'N/A'}`, margin + 5, 10);
          addText(`LogP: ${variant.properties?.logp?.toFixed(3) || 'N/A'}`, margin + 5, 10);
          addText(`Similarity: ${(variant.similarity * 100).toFixed(1)}%`, margin + 5, 10);
  
          if (rdkitLoaded && window.RDKit && variant.smiles) {
            const mol = window.RDKit.get_mol(variant.smiles);
            if (mol) {
              const svg = mol.get_svg(160, 120);
              y = await addSvgToPDF(svg, margin, y, 90, 60);
              mol.delete();
            }
          }
          y += 3;
        }
      }
    }
  
    doc.save(`${structure.name}-${activeTab}.pdf`);
  };
   return (
    <div className="bg-primary text-text-primary p-6 rounded-lg shadow-lg">
      <div className="flex justify-between items-center mb-6">
        <h2 className="text-2xl font-bold font-heading">{structure.name}</h2>
        <span className="text-sm text-text-secondary">
          Created: {new Date(structure.created).toLocaleString()}
        </span>
      </div>

      <div className="border-b border-accent mb-6">
        <nav className="-mb-px flex space-x-8">
          <button
            className={`pb-4 px-1 border-b-2 font-medium text-sm font-body ${
              activeTab === 'parent'
                ? 'border-accent text-accent'
                : 'border-transparent text-text-secondary hover:text-accent-secondary hover:border-accent-secondary'
            }`}
            onClick={() => handleTabSwitch('parent')}
          >
            Parent Structure
          </button>
          <button
            className={`pb-4 px-1 border-b-2 font-medium text-sm font-body ${
              activeTab === 'variants'
                ? 'border-accent text-accent'
                : 'border-transparent text-text-secondary hover:text-accent-secondary hover:border-accent-secondary'
            }`}
            onClick={() => handleTabSwitch('variants')}
          >
            Generated Variants ({structure.generatedStructures?.length || 0})
          </button>
          {selectedVariant && (
            <button
              className={`pb-4 px-1 border-b-2 font-medium text-sm font-body ${
                activeTab === 'variant'
                  ? 'border-accent text-accent'
                  : 'border-transparent text-text-secondary hover:text-accent-secondary hover:border-accent-secondary'
              }`}
              onClick={() => handleTabSwitch('variant')}
            >
              Selected Variant
            </button>
          )}
        </nav>
      </div>

      {activeTab === 'parent' && (
        <div>
          <div className="bg-secondary p-4 rounded-lg mb-6">
            <h3 className="text-lg font-semibold mb-2 font-heading">SMILES</h3>
            <p className="font-code text-sm break-all bg-primary p-2 rounded">{structure.smiles}</p>
          </div>

          {renderMolecule(structure.smiles, parentMolRef)}

          <div className="mt-4 flex justify-center">
            <button
              onClick={() => setShowParent3D(!showParent3D)}
              className="bg-accent hover:bg-accent-secondary text-primary font-bold py-2 px-4 rounded transition-colors duration-300 font-body"
              disabled={!rdkitLoaded}
            >
              {showParent3D ? 'Hide 3D View' : 'Show 3D Visualization'}
            </button>
          </div>

          {showParent3D && (
            <>
              <div
                ref={parent3DRef}
                className="w-full h-96 mt-4 border-2 border-accent rounded-lg relative bg-secondary flex items-center justify-center"
              >
                {!is3DmolLoaded && (
                  <div className="flex flex-col items-center">
                    <ClipLoader color="#00F5D4" size={50} />
                    <p className="mt-2 text-text-secondary">Loading 3D viewer...</p>
                  </div>
                )}
              </div>
              {is3DmolLoaded && (
                <div className="mt-4">
                  <div className="flex justify-center flex-wrap gap-3">
                    <select
                      value={parentStyle}
                      onChange={(e) => setParentStyle(e.target.value)}
                      className="border border-accent rounded p-2 bg-secondary text-text-primary font-body"
                    >
                      <option value="ballstick">Ball & Stick</option>
                      <option value="stick">Stick</option>
                      <option value="sphere">Sphere</option>
                      <option value="surface">Surface</option>
                    </select>
                    <button
                      onClick={() => setShowLabels(!showLabels)}
                      className="bg-accent-secondary hover:bg-accent text-primary py-2 px-4 rounded transition-colors duration-300 font-body"
                    >
                      {showLabels ? 'Hide Labels' : 'Show Labels'}
                    </button>
                    <button
                      onClick={() => setSpinAnimation(!spinAnimation)}
                      className="bg-accent hover:bg-accent-secondary text-primary py-2 px-4 rounded transition-colors duration-300 font-body"
                    >
                      {spinAnimation ? 'Stop Spin' : 'Spin Model'}
                    </button>
                    <button
                      onClick={() => download3DView(parentViewerRef, `${structure.name}-3D`)}
                      className="bg-accent-secondary hover:bg-accent text-primary py-2 px-4 rounded transition-colors duration-300 font-body"
                    >
                      Download PNG
                    </button>
                  </div>
                  {highlightedAtom && (
                    <div className="mt-4 text-center text-sm bg-secondary text-text-primary p-3 rounded-lg border border-accent">
                      <strong className="font-heading">Highlighted Atom:</strong> {highlightedAtom.element} (x: {highlightedAtom.x.toFixed(2)}, y: {highlightedAtom.y.toFixed(2)}, z: {highlightedAtom.z.toFixed(2)})
                    </div>
                  )}
                  {aiDescription && (
                    <div className="mt-4 text-center text-sm bg-secondary text-text-primary p-3 rounded-lg border border-accent">
                      <strong className="font-heading">AI Insight:</strong> {aiDescription}
                    </div>
                  )}
                </div>
              )}
            </>
          )}

          <div className="mt-6">
            <h3 className="text-lg font-semibold mb-2 font-heading">Detailed Information</h3>
            <EnhancedGeminiInfo
              information={structure.information}
              isLoading={false}
              error={null}
              rdkitLoaded={rdkitLoaded}
            />
          </div>

          <div className="mt-6 flex justify-center">
            <button
              className="bg-success hover:bg-green-700 text-primary font-bold py-2 px-6 rounded transition-colors duration-300 font-body"
              onClick={exportToPDF}
            >
              Export to PDF
            </button>
          </div>
        </div>
      )}

      {activeTab === 'variants' && (
        <div>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {structure.generatedStructures?.map((variant, index) => (
              <div
                key={index}
                className="border-2 border-accent-secondary rounded-lg p-4 hover:shadow-lg cursor-pointer transition-all duration-300 bg-secondary hover:border-accent"
                onClick={() => handleVariantSelect(variant)}
              >
                <div className="h-40 flex items-center justify-center bg-primary rounded-lg">
                  {rdkitLoaded && window.RDKit && variant.smiles ? (
                    <div
                      ref={(el) => (variantMolRefs.current[index] = el)}
                      className="w-full h-full"
                    />
                  ) : (
                    <div className="w-full h-full flex items-center justify-center rounded">
                      <p className="text-text-secondary text-xs">Visualization unavailable</p>
                    </div>
                  )}
                </div>
                <div className="mt-4 space-y-2">
                  <div className="flex justify-between text-sm font-label">
                    <span className="font-medium">QED:</span>
                    <span>{variant.properties?.qed?.toFixed(3) || 'N/A'}</span>
                  </div>
                  <div className="flex justify-between text-sm font-label">
                    <span className="font-medium">LogP:</span>
                    <span>{variant.properties?.logp?.toFixed(3) || 'N/A'}</span>
                  </div>
                  <div className="flex justify-between text-sm font-label">
                    <span className="font-medium">Similarity:</span>
                    <span>{(variant.similarity * 100).toFixed(1)}%</span>
                  </div>
                </div>
              </div>
            ))}
          </div>
          <div className="mt-6 flex justify-center">
            <button
              className="bg-success hover:bg-green-700 text-primary font-bold py-2 px-6 rounded transition-colors duration-300 font-body"
              onClick={exportToPDF}
            >
              Export to PDF
            </button>
          </div>
        </div>
      )}

      {activeTab === 'variant' && (
        <div>
          {!selectedVariant ? (
            <div className="text-center text-text-secondary p-6 bg-secondary rounded-lg">
              <p>No variant selected. Please select a variant from the Generated Variants tab.</p>
            </div>
          ) : (
            <>
              <div className="bg-secondary p-4 rounded-lg mb-6">
                <h3 className="text-lg font-semibold mb-2 font-heading">Variant SMILES</h3>
                <p className="font-code text-sm break-all bg-primary p-2 rounded">{selectedVariant.smiles || 'N/A'}</p>
              </div>

              {renderMolecule(selectedVariant.smiles, variantMolRef)}

              <div className="mt-4 flex justify-center">
                <button
                  onClick={() => setShowVariant3D(!showVariant3D)}
                  className="bg-accent hover:bg-accent-secondary text-primary font-bold py-2 px-4 rounded transition-colors duration-300 font-body"
                  disabled={!rdkitLoaded}
                >
                  {showVariant3D ? 'Hide 3D View' : 'Show 3D Visualization'}
                </button>
              </div>

              {showVariant3D && (
                <>
                  <div
                    ref={variant3DRef}
                    className="w-full h-96 mt-4 border-2 border-accent rounded-lg relative bg-secondary flex items-center justify-center"
                  >
                    {!is3DmolLoaded && (
                      <div className="flex flex-col items-center">
                        <ClipLoader color="#00F5D4" size={50} />
                        <p className="mt-2 text-text-secondary">Loading 3D viewer...</p>
                      </div>
                    )}
                  </div>
                  {is3DmolLoaded && (
                    <div className="mt-4">
                      <div className="flex justify-center flex-wrap gap-3">
                        <select
                          value={variantStyle}
                          onChange={(e) => setVariantStyle(e.target.value)}
                          className="border border-accent rounded p-2 bg-secondary text-text-primary font-body"
                        >
                          <option value="ballstick">Ball & Stick</option>
                          <option value="stick">Stick</option>
                          <option value="sphere">Sphere</option>
                          <option value="surface">Surface</option>
                        </select>
                        <button
                          onClick={() => setShowLabels(!showLabels)}
                          className="bg-accent-secondary hover:bg-accent text-primary py-2 px-4 rounded transition-colors duration-300 font-body"
                        >
                          {showLabels ? 'Hide Labels' : 'Show Labels'}
                        </button>
                        <button
                          onClick={() => setSpinAnimation(!spinAnimation)}
                          className="bg-accent hover:bg-accent-secondary text-primary py-2 px-4 rounded transition-colors duration-300 font-body"
                        >
                          {spinAnimation ? 'Stop Spin' : 'Spin Model'}
                        </button>
                        <button
                          onClick={() => download3DView(variantViewerRef, `${structure.name}-variant-3D`)}
                          className="bg-accent-secondary hover:bg-accent text-primary py-2 px-4 rounded transition-colors duration-300 font-body"
                        >
                          Download PNG
                        </button>
                      </div>
                      {highlightedAtom && (
                        <div className="mt-4 text-center text-sm bg-secondary text-text-primary p-3 rounded-lg border border-accent">
                          <strong className="font-heading">Highlighted Atom:</strong> {highlightedAtom.element} (x: {highlightedAtom.x.toFixed(2)}, y: {highlightedAtom.y.toFixed(2)}, z: {highlightedAtom.z.toFixed(2)})
                        </div>
                      )}
                      {aiDescription && (
                        <div className="mt-4 text-center text-sm bg-secondary text-text-primary p-3 rounded-lg border border-accent">
                          <strong className="font-heading">AI Insight:</strong> {aiDescription}
                        </div>
                      )}
                    </div>
                  )}
                </>
              )}

              <div className="mt-6">
                <h3 className="text-lg font-semibold mb-2 font-heading">Property Comparison</h3>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                  <div className="bg-secondary p-4 rounded-lg border border-accent">
                    <h4 className="text-md font-medium font-heading">Parent</h4>
                    <div className="mt-2 space-y-2 font-label">
                      <p><strong>QED:</strong> {structure.properties?.qed?.toFixed(3) || 'N/A'}</p>
                      <p><strong>LogP:</strong> {structure.properties?.logp?.toFixed(3) || 'N/A'}</p>
                    </div>
                  </div>
                  <div className="bg-secondary p-4 rounded-lg border border-accent">
                    <h4 className="text-md font-medium font-heading">Variant</h4>
                    <div className="mt-2 space-y-2 font-label">
                      <p><strong>QED:</strong> {selectedVariant.properties?.qed?.toFixed(3) || 'N/A'}</p>
                      <p><strong>LogP:</strong> {selectedVariant.properties?.logp?.toFixed(3) || 'N/A'}</p>
                      <p><strong>Similarity:</strong> {(selectedVariant.similarity * 100).toFixed(1)}%</p>
                    </div>
                  </div>
                </div>
              </div>

              <div className="mt-6">
                <h3 className="text-lg font-semibold mb-2 font-heading">Detailed Information</h3>
                <EnhancedGeminiInfo
                  information={variantInfo}
                  variantData={{
                    smiles: selectedVariant.smiles,
                    properties: {
                      parent: {
                        qed: structure.properties?.qed?.toFixed(3) || 'N/A',
                        logp: structure.properties?.logp?.toFixed(3) || 'N/A',
                      },
                      variant: {
                        qed: selectedVariant.properties?.qed?.toFixed(3) || 'N/A',
                        logp: selectedVariant.properties?.logp?.toFixed(3) || 'N/A',
                      },
                      similarity: `${(selectedVariant.similarity * 100).toFixed(1)}%`,
                    },
                  }}
                  isLoading={loadingVariantInfo}
                  error={variantInfoError}
                  rdkitLoaded={rdkitLoaded}
                />
              </div>

              <div className="mt-6 flex justify-center">
                <button
                  className="bg-success hover:bg-green-700 text-primary font-bold py-2 px-6 rounded transition-colors duration-300 font-body"
                  onClick={exportToPDF}
                >
                  Export to PDF
                </button>
              </div>
            </>
          )}
        </div>
      )}
    </div>
  );

};

export default StructureDetails;