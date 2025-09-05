// src/components/Summary.jsx
import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { jsPDF } from 'jspdf';
import { toast } from 'react-hot-toast';
import { useAuthStore } from '../../Store/auth.store.js';
import { Download, FileDown, Loader2, AlertCircle, ChevronDown, ChevronUp, Search, X } from 'lucide-react';
import debounce from 'lodash/debounce';
import { v4 as uuidv4 } from 'uuid';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'http://localhost:5000/api';
const axiosInstance = axios.create({
  baseURL: import.meta.mode === 'development' ? API_BASE_URL : '/api',
  withCredentials: true,
});

const Summary = () => {
  const { user } = useAuthStore();
  const [summaryData, setSummaryData] = useState(null);
  const [savedItems, setSavedItems] = useState({});
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState('summary');
  const [expandedSections, setExpandedSections] = useState({});
  const [downloadFormat, setDownloadFormat] = useState('pdf');
  const [searchQuery, setSearchQuery] = useState('');
  const [downloadingItem, setDownloadingItem] = useState(null);
  const [itemDownloadFormats, setItemDownloadFormats] = useState({});
  const [moleculeDownloadFormats, setMoleculeDownloadFormats] = useState({});
  const [currentMoleculeId, setCurrentMoleculeId] = useState(null);
  const [pdfGeneratedItems, setPdfGeneratedItems] = useState({});
  const [generationProgress, setGenerationProgress] = useState(0);

  useEffect(() => {
    fetchSummaryData();
    fetchSavedItems();
  }, []);

  const fetchSummaryData = async () => {
    try {
      const response = await axiosInstance.get('/summary/summary');
      const data = response.data.data;

      // Filter all fields to exclude saved items
      const filterField = (fieldData, type) => {
        if (!fieldData) return [];
        const savedIds = Object.values(savedItems)
          .flatMap(items => items.filter(item => item.type === type).map(item => item.data._id))
          .filter(id => id);
        return fieldData.filter(item => !savedIds.includes(item._id));
      };

      const filteredData = {
        newDrugs: filterField(data.newDrugs, 'newDrug'),
        costEstimations: filterField(data.costEstimations, 'costEstimation'),
        drugNames: filterField(data.drugNames, 'drugName'),
        researchPapers: filterField(data.researchPapers, 'researchPaper'),
        generatedResearchPapers: filterField(data.generatedResearchPapers, 'generatedResearchPaper'),
        predictDiseases: filterField(data.predictDiseases, 'predictDisease'),
        targetProteins: filterField(data.targetProteins, 'targetProtein'),
        // toxicityResults: filterField(data.toxicityResults, 'toxicityResult'),
      };

      setSummaryData(filteredData);
      if (!currentMoleculeId) {
        setCurrentMoleculeId(uuidv4());
      }
    } catch (err) {
      const message = err.response?.data?.message || 'Failed to fetch summary data';
      setError(message);
      toast.error(message);
    } finally {
      setLoading(false);
    }
  };

  const fetchSavedItems = async () => {
    try {
      const response = await axiosInstance.get('/summary/saved');
      setSavedItems(response.data.data);
    } catch (err) {
      toast.error('Failed to fetch saved items: ' + (err.response?.data?.message || err.message));
    }
  };

  // Calculate progress based on non-empty summaryData fields
  useEffect(() => {
    if (summaryData) {
      const fields = [
        summaryData.newDrugs,
        summaryData.costEstimations,
        summaryData.drugNames,
        summaryData.researchPapers,
        summaryData.generatedResearchPapers,
        summaryData.predictDiseases,
        summaryData.targetProteins,
        // summaryData.toxicityResults,
      ];
      const nonEmptyFields = fields.filter(field => field?.length > 0).length;
      const progress = Math.round((nonEmptyFields / fields.length) * 100);
      setGenerationProgress(progress);
    } else {
      setGenerationProgress(0);
    }
  }, [summaryData]);

  useEffect(() => {
    if (summaryData) {
      const filterField = (fieldData, type) => {
        if (!fieldData) return [];
        const savedIds = Object.values(savedItems)
          .flatMap(items => items.filter(item => item.type === type).map(item => item.data._id))
          .filter(id => id);
        return fieldData.filter(item => !savedIds.includes(item._id));
      };

      setSummaryData(prev => ({
        newDrugs: filterField(prev.newDrugs, 'newDrug'),
        costEstimations: filterField(prev.costEstimations, 'costEstimation'),
        drugNames: filterField(prev.drugNames, 'drugName'),
        researchPapers: filterField(prev.researchPapers, 'researchPaper'),
        generatedResearchPapers: filterField(prev.generatedResearchPapers, 'generatedResearchPaper'),
        predictDiseases: filterField(prev.predictDiseases, 'predictDisease'),
        targetProteins: filterField(prev.targetProteins, 'targetProtein'),
        // toxicityResults: filterField(prev.toxicityResults, 'toxicityResult'),
      }));
    }
  }, [savedItems]);

  const saveSummaryData = async () => {
    try {
      const savePromises = [];
      const addItems = (items, type) => {
        items?.forEach(item => {
          if (item._id && !pdfGeneratedItems[`${type}-${item._id}`]) {
            savePromises.push({ item, type, moleculeId: currentMoleculeId });
          }
        });
      };
      addItems(summaryData?.newDrugs, 'newDrug');
      addItems(summaryData?.costEstimations, 'costEstimation');
      addItems(summaryData?.drugNames, 'drugName');
      addItems(summaryData?.researchPapers, 'researchPaper');
      addItems(summaryData?.generatedResearchPapers, 'generatedResearchPaper');
      addItems(summaryData?.predictDiseases, 'predictDisease');
      addItems(summaryData?.targetProteins, 'targetProtein');
      // addItems(summaryData?.toxicityResults, 'toxicityResult');

      if (savePromises.length > 0) {
        await axiosInstance.post('/summary/bulk-save', savePromises);
        toast.success('Summary data saved successfully');
        await fetchSavedItems();
        setCurrentMoleculeId(uuidv4());
        setGenerationProgress(0);
        setSummaryData({
          newDrugs: [],
          costEstimations: [],
          drugNames: [],
          researchPapers: [],
          generatedResearchPapers: [],
          predictDiseases: [],
          targetProteins: [],
          // toxicityResults: [],
        });
      } else {
        toast.error('No data to save');
      }
    } catch (err) {
      toast.error('Failed to save summary data: ' + (err.response?.data?.message || err.message));
    }
  };

  const handleTabChange = (tab) => {
    setActiveTab(tab);
  };

  const flattenDataForCSV = (data) => {
    const rows = [];
    const addSection = (section, items, prefix) => {
      items?.forEach((item, index) => {
        const row = { Section: `${prefix} ${index + 1}` };
        Object.entries(item).forEach(([key, value]) => {
          if (Array.isArray(value)) {
            row[key] = value.join('; ');
          } else if (typeof value === 'object' && value) {
            row[key] = JSON.stringify(value);
          } else {
            row[key] = value || 'N/A';
          }
        });
        rows.push(row);
      });
    };
    addSection('newDrugs', data.newDrugs, 'New Drug');
    addSection('costEstimations', data.costEstimations, 'Cost Estimation');
    addSection('drugNames', data.drugNames, 'Drug Name');
    addSection('researchPapers', data.researchPapers, 'Research Paper');
    addSection('generatedResearchPapers', data.generatedResearchPapers, 'Generated Research Paper');
    addSection('predictDiseases', data.predictDiseases, 'Predict Disease');
    addSection('targetProteins', data.targetProteins, 'Target Protein');
    // addSection('toxicityResults', data.toxicityResults, 'Toxicity Result');
    return rows;
  };

  const generateCSV = (data, filename = 'drug-discovery-summary.csv') => {
    if (!data) {
      toast.error('No data available to generate CSV');
      return;
    }
    try {
      const rows = flattenDataForCSV(data);
      if (rows.length === 0) {
        toast.error('No data to export as CSV');
        return;
      }
      const headers = ['Section', ...new Set(rows.flatMap(row => Object.keys(row)))];
      const csvContent = [
        headers.join(','),
        ...rows.map(row =>
          headers.map(header => `"${(row[header] || '').toString().replace(/"/g, '""')}"`).join(',')
        ),
      ].join('\n');
      const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.click();
      URL.revokeObjectURL(url);
      toast.success('CSV downloaded successfully!');
    } catch (err) {
      toast.error('Failed to generate CSV: ' + err.message);
    }
  };

  const generateTXT = (data, filename = 'drug-discovery-summary.txt') => {
    if (!data) {
      toast.error('No data available to generate TXT');
      return;
    }
    try {
      const content = JSON.stringify(data, null, 2);
      const blob = new Blob([content], { type: 'text/plain;charset=utf-8;' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.click();
      URL.revokeObjectURL(url);
      toast.success('TXT downloaded successfully!');
    } catch (err) {
      toast.error('Failed to generate TXT: ' + err.message);
    }
  };

  const generateJSON = (data, filename = 'drug-discovery-summary.json') => {
    if (!data) {
      toast.error('No data available to generate JSON');
      return;
    }
    try {
      const jsonStr = JSON.stringify(data, null, 2);
      const blob = new Blob([jsonStr], { type: 'application/json' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.click();
      URL.revokeObjectURL(url);
      toast.success('JSON downloaded successfully!');
    } catch (err) {
      toast.error('Failed to generate JSON: ' + err.message);
    }
  };

  const generatePDF = async (data, filename = 'drug-discovery-summary.pdf') => {
    if (!data) {
      throw new Error('No data available to generate PDF');
    }
    try {
      const doc = new jsPDF({
        orientation: 'portrait',
        unit: 'mm',
        format: 'a4',
      });
      const margin = 15;
      const pageWidth = doc.internal.pageSize.width;
      const pageHeight = doc.internal.pageSize.height;
      const contentWidth = pageWidth - 2 * margin;
      const lineHeight = 7;
      const sectionGap = 10;
      let y = margin;

      const styles = {
        title: { size: 24, bold: true, color: [10, 25, 47] }, // primary
        sectionTitle: { size: 18, bold: true, color: [23, 42, 69] }, // secondary
        subsectionTitle: { size: 14, bold: true, color: [55, 65, 81] },
        tableHeader: { size: 10, bold: true, color: [255, 255, 255], fill: [94, 129, 244] }, // accent-secondary
        tableRow: { size: 10, color: [0, 0, 0] },
        footer: { size: 8, color: [160, 160, 160] }, // text-secondary
      };

      const addText = (text, options = {}) => {
        const { size = 12, bold = false, x = margin, color = [0, 0, 0], maxWidth = contentWidth, spacing = 1.2, align = 'left' } = options;
        doc.setFontSize(size);
        doc.setFont(bold ? 'Inter-bold' : 'Inter', 'normal', bold ? 'bold' : 'normal');
        doc.setTextColor(...color);
        const words = text.split(' ');
        let currentLine = '';
        const lines = [];
        for (const word of words) {
          const testLine = currentLine ? `${currentLine} ${word}` : word;
          const testWidth = doc.getStringUnitWidth(testLine) * size / doc.internal.scaleFactor;
          if (testWidth > maxWidth) {
            lines.push(currentLine);
            currentLine = word;
          } else {
            currentLine = testLine;
          }
        }
        if (currentLine) lines.push(currentLine);
        const neededHeight = lines.length * size * 0.352777778 * spacing;
        if (y + neededHeight > pageHeight - margin) {
          doc.addPage();
          y = margin;
        }
        lines.forEach((line, index) => {
          const lineY = y + (index + 1) * (size * 0.352777778 * spacing);
          doc.text(line, x, lineY, { align });
        });
        y += neededHeight + (lines.length > 1 ? size * 0.2 : 0);
        return y;
      };

      const addTable = (headers, rows, columnWidths, options = {}) => {
        const { headerColor = styles.tableHeader.color, headerFill = styles.tableHeader.fill, cellPadding = 3, fontSize = styles.tableRow.size, lineHeight = 1.2, alternateRows = true, startY = y } = options;
        let currentY = startY;
        const headerRowHeight = fontSize * lineHeight + cellPadding * 2;
        const calculateTextHeight = (text, width) => {
          const lines = doc.splitTextToSize(String(text), width - cellPadding * 2);
          return lines.length * fontSize * lineHeight + cellPadding * 2;
        };
        const processedRows = rows.map(row => {
          const heights = row.map((cell, i) => calculateTextHeight(cell, columnWidths[i]));
          return { cells: row, height: Math.max(...heights) };
        });
        const drawHeader = () => {
          let x = margin;
          doc.setFontSize(fontSize);
          doc.setFont('Inter-bold', 'normal', 'bold');
          doc.setTextColor(...headerColor);
          doc.setFillColor(...headerFill);
          headers.forEach((header, i) => {
            doc.rect(x, currentY, columnWidths[i], headerRowHeight, 'F');
            doc.text(header, x + cellPadding, currentY + cellPadding + fontSize * 0.352777778);
            x += columnWidths[i];
          });
          currentY += headerRowHeight;
          doc.setTextColor(...styles.tableRow.color);
          doc.setFont('Inter', 'normal', 'normal');
        };
        drawHeader();
        processedRows.forEach((row, rowIndex) => {
          if (currentY + row.height > pageHeight - margin) {
            doc.addPage();
            currentY = margin;
            drawHeader();
          }
          let x = margin;
          if (alternateRows && rowIndex % 2 === 0) {
            doc.setFillColor(23, 42, 69); // secondary
            doc.rect(margin, currentY, columnWidths.reduce((a, b) => a + b, 0), row.height, 'F');
          }
          row.cells.forEach((cell, i) => {
            doc.rect(x, currentY, columnWidths[i], row.height);
            const text = String(cell || 'N/A');
            const lines = doc.splitTextToSize(text, columnWidths[i] - cellPadding * 2);
            doc.text(lines, x + cellPadding, currentY + cellPadding + fontSize * 0.352777778, { maxWidth: columnWidths[i] - cellPadding * 2, lineHeightFactor: lineHeight });
            x += columnWidths[i];
          });
          currentY += row.height;
        });
        y = currentY + sectionGap;
        return y;
      };

      const safeStringify = (obj) => {
        if (typeof obj === 'string') return obj;
        if (Array.isArray(obj)) return obj.map(item => (typeof item === 'object' ? JSON.stringify(item, null, 2) : item)).join('\n');
        try {
          return typeof obj === 'object' ? JSON.stringify(obj, null, 2) : String(obj);
        } catch {
          return String(obj);
        }
      };

      doc.setFontSize(styles.title.size);
      doc.setFont('Inter-bold', 'normal', 'bold');
      doc.setTextColor(...styles.title.color);
      doc.text('Drug Discovery Report', pageWidth / 2, 60, { align: 'center' });
      doc.setFontSize(14);
      doc.setTextColor(160, 160, 160); // text-secondary
      doc.text(`Generated on: ${new Date().toLocaleString()}`, pageWidth / 2, 90, { align: 'center' });

      if (data.newDrugs?.length > 0) {
        doc.addPage();
        y = margin;
        addText('New Drug Discovery', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        data.newDrugs.forEach((drug, index) => {
          if (y > pageHeight - 100) {
            doc.addPage();
            y = margin;
          }
          addText(`Compound ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
          const rows = [
            ['Title', drug.newmoleculetitle || 'N/A'],
            ['SMILES', drug.newSmiles || 'N/A'],
            ['IUPAC Name', drug.newIupacName || 'N/A'],
            ['Potential Diseases', drug.potentialDiseases || 'N/A'],
            ['Conversion Details', safeStringify(drug.conversionDetails)],
            ['Additional Info', safeStringify(drug.information)],
          ];
          y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
        });
      }

      if (data.costEstimations?.length > 0) {
        doc.addPage();
        y = margin;
        addText('Cost Estimations', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        const rows = data.costEstimations.map(cost => [
          cost.smiles || 'N/A',
          cost.estimatedcost || 'N/A',
          safeStringify(cost.information),
        ]);
        y = addTable(['SMILES', 'Estimated Cost', 'Information'], rows, [60, 40, contentWidth - 100], { startY: y });
      }

      if (data.drugNames?.length > 0) {
        doc.addPage();
        y = margin;
        addText('AI Generated Names', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        data.drugNames.forEach((name, index) => {
          if (y > pageHeight - 100) {
            doc.addPage();
            y = margin;
          }
          addText(`Name Suggestion ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
          const rows = [
            ['Suggested Name', name.suggestedName || 'N/A'],
            ['Molecule Title', name.moleculeTitle || 'N/A'],
            ['SMILES', name.smiles || 'N/A'],
            ['Status', name.status || 'N/A'],
            ['Details', safeStringify(name.namingDetails)],
          ];
          y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
        });
      }

      if (data.researchPapers?.length > 0) {
        doc.addPage();
        y = margin;
        addText('Research Papers', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        data.researchPapers.forEach((paper, index) => {
          if (y > pageHeight - 100) {
            doc.addPage();
            y = margin;
          }
          addText(`Paper ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
          addText(`Molecule: ${paper.molecule?.title || 'N/A'}`);
          addText(`SMILES: ${paper.molecule?.smiles || 'N/A'}`);
          if (paper.papers?.length > 0) {
            paper.papers.forEach((p, pIndex) => {
              if (y > pageHeight - 100) {
                doc.addPage();
                y = margin;
              }
              const rows = [
                ['Title', p.title || 'N/A'],
                ['Authors', p.authors || 'N/A'],
                ['Year', p.year || 'N/A'],
                ['Abstract', p.abstract || 'N/A'],
                ['DOI', p.doi || 'N/A'],
                ['URL', p.url || 'N/A'],
              ];
              y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
            });
          }
        });
      }

      if (data.generatedResearchPapers?.length > 0) {
        doc.addPage();
        y = margin;
        addText('Generated Research Papers', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        data.generatedResearchPapers.forEach((paper, index) => {
          if (y > pageHeight - 100) {
            doc.addPage();
            y = margin;
          }
          addText(`Generated Paper ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
          addText(`Molecule: ${paper.molecule?.title || 'N/A'}`);
          addText(`SMILES: ${paper.molecule?.smiles || 'N/A'}`);
          const rows = [
            ['Title', paper.paper?.title || 'N/A'],
            ['Authors', paper.paper?.authors || 'N/A'],
            ['Abstract', paper.paper?.abstract || 'N/A'],
            ['Keywords', safeStringify(paper.paper?.keywords)],
            ['Introduction', paper.paper?.introduction || 'N/A'],
            ['Methodology', paper.paper?.methodology || 'N/A'],
            ['Results and Discussion', paper.paper?.resultsAndDiscussion || 'N/A'],
            ['Conclusion', paper.paper?.conclusion || 'N/A'],
            ['References', safeStringify(paper.paper?.references)],
          ];
          y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
        });
      }

      if (data.predictDiseases?.length > 0) {
        doc.addPage();
        y = margin;
        addText('Disease Predictions', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        data.predictDiseases.forEach((prediction, index) => {
          if (y > pageHeight - 100) {
            doc.addPage();
            y = margin;
          }
          addText(`Prediction ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
          addText(`SMILES: ${prediction.smiles || 'N/A'}`);
          const rows = [
            ['Disease', prediction.disease || 'N/A'],
            ['Confidence', prediction.confidence || 'N/A'],
            ['Details', safeStringify(prediction.details)],
          ];
          y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
        });
      }

      if (data.targetProteins?.length > 0) {
        doc.addPage();
        y = margin;
        addText('Target Proteins', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
        data.targetProteins.forEach((target, index) => {
          if (y > pageHeight - 100) {
            doc.addPage();
            y = margin;
          }
          addText(`Target ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
          addText(`SMILES: ${target.smiles || 'N/A'}`);
          const rows = [
            ['Protein', target.protein || 'N/A'],
            ['Confidence', target.confidence || 'N/A'],
            ['MOA', target.moa || 'N/A'],
            ['Pathways', safeStringify(target.pathways)],
            ['Diseases', safeStringify(target.diseases)],
            ['Interactions', safeStringify(target.knownInteractions)],
          ];
          y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
        });
      }

      // if (data.toxicityResults?.length > 0) {
      //   doc.addPage();
      //   y = margin;
      //   addText('Toxicity Predictions', { size: styles.sectionTitle.size, bold: true, color: styles.sectionTitle.color });
      //   data.toxicityResults.forEach((toxicity, index) => {
      //     if (y > pageHeight - 100) {
      //       doc.addPage();
      //       y = margin;
      //     }
      //     addText(`Prediction ${index + 1}`, { size: styles.subsectionTitle.size, bold: true });
      //     const rows = [
      //       ['SMILES', toxicity.smiles || 'N/A'],
      //       ['LD50', toxicity.toxicityResult?.acuteToxicity?.LD50 || 'N/A'],
      //       ['Toxicity Class', toxicity.toxicityResult?.acuteToxicity?.toxicityClass || 'N/A'],
      //       ['Hepatotoxicity', toxicity.toxicityResult?.endpoints?.hepatotoxicity || 'N/A'],
      //       ['Carcinogenicity', toxicity.toxicityResult?.endpoints?.carcinogenicity || 'N/A'],
      //       ['Gemini Analysis', safeStringify(toxicity.geminiAnalysis)],
      //     ];
      //     y = addTable(['Field', 'Value'], rows, [50, contentWidth - 50], { startY: y });
      //   });
      // }

      const pageCount = doc.internal.getNumberOfPages();
      for (let i = 1; i <= pageCount; i++) {
        doc.setPage(i);
        doc.setFontSize(styles.footer.size);
        doc.setTextColor(...styles.footer.color);
        doc.text(`Page ${i} of ${pageCount}`, pageWidth - margin, pageHeight - 10, { align: 'right' });
      }

      doc.save(filename);
      return filename;
    } catch (err) {
      console.error('PDF generation failed:', err);
      throw new Error('Failed to generate PDF: ' + err.message);
    }
  };

  const handleDownload = (data, format, filename, key, index) => {
    setDownloadingItem(`${key}-${index}`);
    if (format === 'pdf') {
      setPdfGeneratedItems(prev => ({
        ...prev,
        [`${key}-${index}`]: true,
      }));
    }
    switch (format) {
      case 'pdf':
        generatePDF(data, filename);
        break;
      case 'csv':
        generateCSV(data, filename.replace('.pdf', '.csv'));
        break;
      case 'json':
        generateJSON(data, filename.replace('.pdf', '.json'));
        break;
      case 'txt':
        generateTXT(data, filename.replace('.pdf', '.txt'));
        break;
      default:
        toast.error('Unsupported format');
    }
    setTimeout(() => setDownloadingItem(null), 1000);
  };

  const toggleSection = (section, index) => {
    setExpandedSections(prev => ({
      ...prev,
      [`${section}-${index}`]: !prev[`${section}-${index}`],
    }));
  };

  const renderTableContent = (data, type) => {
    switch (type) {
      case 'newDrug':
        return (
          <table className="w-full border-collapse">
            <thead>
              <tr className="bg-accent-secondary text-white">
                <th className="p-2 border">Field</th>
                <th className="p-2 border">Value</th>
              </tr>
            </thead>
            <tbody>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Title</td><td className="p-2 border text-text-primary">{data.newmoleculetitle || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">SMILES</td><td className="p-2 border text-text-primary">{data.newSmiles || 'N/A'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">IUPAC Name</td><td className="p-2 border text-text-primary">{data.newIupacName || 'Not available'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Potential Diseases</td><td className="p-2 border text-text-primary">{data.potentialDiseases || 'Not specified'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Conversion Details</td><td className="p-2 border text-text-primary">{data.conversionDetails || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Additional Info</td><td className="p-2 border text-text-primary">{data.information || 'N/A'}</td></tr>
            </tbody>
          </table>
        );
      case 'costEstimation':
        return (
          <table className="w-full border-collapse">
            <thead>
              <tr className="bg-accent-secondary text-white">
                <th className="p-2 border">SMILES</th>
                <th className="p-2 border">Estimated Cost</th>
                <th className="p-2 border">Information</th>
              </tr>
            </thead>
            <tbody>
              <tr className="bg-secondary">
                <td className="p-2 border text-text-primary">{data.smiles || 'N/A'}</td>
                <td className="p-2 border text-text-primary">{data.estimatedcost || 'N/A'}</td>
                <td className="p-2 border text-text-primary">{data.information || 'N/A'}</td>
              </tr>
            </tbody>
          </table>
        );
      case 'drugName':
        return (
          <table className="w-full border-collapse">
            <thead>
              <tr className="bg-accent-secondary text-white">
                <th className="p-2 border">Field</th>
                <th className="p-2 border">Value</th>
              </tr>
            </thead>
            <tbody>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Suggested Name</td><td className="p-2 border text-text-primary">{data.suggestedName || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Molecule Title</td><td className="p-2 border text-text-primary">{data.moleculeTitle || 'N/A'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">SMILES</td><td className="p-2 border text-text-primary">{data.smiles || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Status</td><td className="p-2 border text-text-primary">{data.status || 'N/A'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Details</td><td className="p-2 border text-text-primary">{data.namingDetails || 'N/A'}</td></tr>
            </tbody>
          </table>
        );
      case 'researchPaper':
        return (
          <div>
            <h4 className="text-lg font-semibold font-heading text-text-primary">Molecule: {data.molecule?.title || 'N/A'}</h4>
            <p className="font-body text-text-secondary">SMILES: {data.molecule?.smiles || 'N/A'}</p>
            {data.papers?.map((p, i) => (
              <table key={i} className="w-full border-collapse mt-4">
                <thead>
                  <tr className="bg-accent-secondary text-white">
                    <th className="p-2 border">Field</th>
                    <th className="p-2 border">Value</th>
                  </tr>
                </thead>
                <tbody>
                  <tr className="bg-secondary"><td className="p-2 border text-text-primary">Paper Title</td><td className="p-2 border text-text-primary">{p.title || 'N/A'}</td></tr>
                  <tr><td className="p-2 border text-text-primary">Authors</td><td className="p-2 border text-text-primary">{p.authors || 'N/A'}</td></tr>
                  <tr className="bg-secondary"><td className="p-2 border text-text-primary">Year</td><td className="p-2 border text-text-primary">{p.year || 'N/A'}</td></tr>
                  <tr><td className="p-2 border text-text-primary">Abstract</td><td className="p-2 border text-text-primary">{p.abstract || 'N/A'}</td></tr>
                  <tr className="bg-secondary"><td className="p-2 border text-text-primary">DOI</td><td className="p-2 border text-text-primary">{p.doi || 'N/A'}</td></tr>
                  <tr><td className="p-2 border text-text-primary">URL</td><td className="p-2 border text-text-primary">{p.url || 'N/A'}</td></tr>
                </tbody>
              </table>
            ))}
          </div>
        );
      case 'generatedResearchPaper':
        return (
          <div>
            <h4 className="text-lg font-semibold font-heading text-text-primary">Molecule: {data.molecule?.title || 'N/A'}</h4>
            <p className="font-body text-text-secondary">SMILES: {data.molecule?.smiles || 'N/A'}</p>
            <table className="w-full border-collapse mt-4">
              <thead>
                <tr className="bg-accent-secondary text-white">
                  <th className="p-2 border">Field</th>
                  <th className="p-2 border">Value</th>
                </tr>
              </thead>
              <tbody>
                <tr className="bg-secondary"><td className="p-2 border text-text-primary">Paper Title</td><td className="p-2 border text-text-primary">{data.paper?.title || 'N/A'}</td></tr>
                <tr><td className="p-2 border text-text-primary">Authors</td><td className="p-2 border text-text-primary">{data.paper?.authors || 'N/A'}</td></tr>
                <tr className="bg-secondary"><td className="p-2 border text-text-primary">Abstract</td><td className="p-2 border text-text-primary">{data.paper?.abstract || 'N/A'}</td></tr>
                <tr><td className="p-2 border text-text-primary">Keywords</td><td className="p-2 border text-text-primary">{data.paper?.keywords?.join(', ') || 'N/A'}</td></tr>
                <tr className="bg-secondary"><td className="p-2 border text-text-primary">Introduction</td><td className="p-2 border text-text-primary">{data.paper?.introduction || 'N/A'}</td></tr>
                <tr><td className="p-2 border text-text-primary">Methodology</td><td className="p-2 border text-text-primary">{data.paper?.methodology || 'N/A'}</td></tr>
                <tr className="bg-secondary"><td className="p-2 border text-text-primary">Results and Discussion</td><td className="p-2 border text-text-primary">{data.paper?.resultsAndDiscussion || 'N/A'}</td></tr>
                <tr><td className="p-2 border text-text-primary">Conclusion</td><td className="p-2 border text-text-primary">{data.paper?.conclusion || 'N/A'}</td></tr>
                <tr className="bg-secondary"><td className="p-2 border text-text-primary">References</td><td className="p-2 border text-text-primary">{data.paper?.references?.join('\n') || 'N/A'}</td></tr>
              </tbody>
            </table>
          </div>
        );
      case 'predictDisease':
        return (
          <table className="w-full border-collapse">
            <thead>
              <tr className="bg-accent-secondary text-white">
                <th className="p-2 border">Field</th>
                <th className="p-2 border">Value</th>
              </tr>
            </thead>
            <tbody>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">SMILES</td><td className="p-2 border text-text-primary">{data.smiles || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Disease</td><td className="p-2 border text-text-primary">{data.disease || 'N/A'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Confidence</td><td className="p-2 border text-text-primary">{data.confidence || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Details</td><td className="p-2 border text-text-primary">{data.details || 'N/A'}</td></tr>
            </tbody>
          </table>
        );
      case 'targetProtein':
        return (
          <table className="w-full border-collapse">
            <thead>
              <tr className="bg-accent-secondary text-white">
                <th className="p-2 border">Field</th>
                <th className="p-2 border">Value</th>
              </tr>
            </thead>
            <tbody>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">SMILES</td><td className="p-2 border text-text-primary">{data.smiles || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Protein</td><td className="p-2 border text-text-primary">{data.protein || 'N/A'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Confidence</td><td className="p-2 border text-text-primary">{data.confidence || 'N/A'}</td></tr>
              <tr><td className="p-2 border text-text-primary">MOA</td><td className="p-2 border text-text-primary">{data.moa || 'N/A'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Pathways</td><td className="p-2 border text-text-primary">{data.pathways?.join(', ') || 'None'}</td></tr>
              <tr><td className="p-2 border text-text-primary">Diseases</td><td className="p-2 border text-text-primary">{data.diseases?.join(', ') || 'None'}</td></tr>
              <tr className="bg-secondary"><td className="p-2 border text-text-primary">Interactions</td><td className="p-2 border text-text-primary">{data.knownInteractions ? JSON.stringify(data.knownInteractions) : 'N/A'}</td></tr>
            </tbody>
          </table>
        );
      // case 'toxicityResult':
      //   return (
      //     <table className="w-full border-collapse">
      //       <thead>
      //         <tr className="bg-accent-secondary text-white">
      //           <th className="p-2 border">Field</th>
      //           <th className="p-2 border">Value</th>
      //         </tr>
      //       </thead>
      //       <tbody>
      //         <tr className="bg-secondary"><td className="p-2 border text-text-primary">SMILES</td><td className="p-2 border text-text-primary">{data.smiles || 'N/A'}</td></tr>
      //         <tr><td className="p-2 border text-text-primary">LD50</td><td className="p-2 border text-text-primary">{data.toxicityResult?.acuteToxicity?.LD50 || 'N/A'}</td></tr>
      //         <tr className="bg-secondary"><td className="p-2 border text-text-primary">Toxicity Class</td><td className="p-2 border text-text-primary">{data.toxicityResult?.acuteToxicity?.toxicityClass || 'N/A'}</td></tr>
      //         <tr><td className="p-2 border text-text-primary">Hepatotoxicity</td><td className="p-2 border text-text-primary">{data.toxicityResult?.endpoints?.hepatotoxicity || 'N/A'}</td></tr>
      //         <tr className="bg-secondary"><td className="p-2 border text-text-primary">Carcinogenicity</td><td className="p-2 border text-text-primary">{data.toxicityResult?.endpoints?.carcinogenicity || 'N/A'}</td></tr>
      //         <tr><td className="p-2 border text-text-primary">Gemini Analysis</td><td className="p-2 border text-text-primary">{data.geminiAnalysis || 'N/A'}</td></tr>
      //       </tbody>
      //     </table>
      //   );
      default:
        return <p className="text-text-secondary font-body">No content available</p>;
    }
  };

  const renderNewMoleculeSection = () => {
    const canSave = summaryData &&
      summaryData.newDrugs?.length > 0 &&
      summaryData.costEstimations?.length > 0 &&
      summaryData.drugNames?.length > 0 &&
      summaryData.researchPapers?.length > 0 &&
      summaryData.generatedResearchPapers?.length > 0 &&
      summaryData.predictDiseases?.length > 0 &&
      summaryData.targetProteins?.length > 0 
      // summaryData.toxicityResults?.length > 0;

    return (
      <section className="mb-8" aria-label="Summary Data">
        <h2 className="text-2xl font-semibold font-heading text-text-primary mb-4">Summary Data</h2>
        <div className="mb-4">
          <div className="text-sm text-text-secondary font-label mb-2">Summary Data Completeness</div>
          <div className="w-full bg-secondary rounded-full h-2.5">
            <div
              className="bg-accent h-2.5 rounded-full transition-all duration-300"
              style={{ width: `${generationProgress}%` }}
            ></div>
          </div>
          <div className="text-sm text-text-secondary font-label mt-1">{generationProgress}% Complete</div>
        </div>
        <div className="relative group">
          <button
            onClick={saveSummaryData}
            disabled={!canSave}
            className={`mt-4 px-4 py-2 rounded-lg font-body transition-colors ${
              canSave
                ? 'bg-success text-white hover:bg-green-600'
                : 'bg-secondary text-text-secondary cursor-not-allowed'
            }`}
            aria-label={canSave ? 'Save all summary data' : 'All fields must be completed to save'}
            aria-disabled={!canSave}
          >
            Save All
          </button>
          {!canSave && (
            <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
              All fields must be completed to save
            </span>
          )}
        </div>
      </section>
    );
  };

  const renderSection = (title, data, key) => {
    const filteredData = data?.filter((item, index) => !pdfGeneratedItems[`${key}-${item._id}`]);

    return (
      <section className="mb-8" aria-label={title}>
        <h2 className="text-2xl font-semibold font-heading text-text-primary mb-4">{title}</h2>
        {filteredData?.length > 0 ? (
          <div className="space-y-4">
            {filteredData.map((item, index) => (
              <div
                key={index}
                className="bg-secondary p-6 rounded-lg shadow-sm hover:shadow-md transition-shadow"
                aria-label={`${title} item ${index + 1}`}
              >
                <div className="flex justify-between items-center">
                  <div>
                    <h3 className="font-medium font-heading text-text-primary">{item.newmoleculetitle || item.suggestedName || item.smiles || item.molecule?.title || 'Item ' + (index + 1)}</h3>
                    <p className="text-sm font-body text-text-secondary mt-1">SMILES: {item.smiles || item.newSmiles || item.molecule?.smiles || 'N/A'}</p>
                  </div>
                  <div className="flex space-x-2 items-center">
                    <button
                      onClick={() => toggleSection(key, index)}
                      className="p-2 text-accent-secondary hover:text-accent relative group"
                      aria-label={expandedSections[`${key}-${index}`] ? 'Collapse details' : 'Expand details'}
                    >
                      {expandedSections[`${key}-${index}`] ? <ChevronUp size={20} /> : <ChevronDown size={20} />}
                      <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                        {expandedSections[`${key}-${index}`] ? 'Collapse' : 'Expand'}
                      </span>
                    </button>
                    <select
                      value={itemDownloadFormats[`${key}-${index}`] || 'pdf'}
                      onChange={(e) => setItemDownloadFormats(prev => ({
                        ...prev,
                        [`${key}-${index}`]: e.target.value,
                      }))}
                      className="px-2 py-1 border rounded-lg font-body text-sm text-text-primary bg-primary focus:ring-2 focus:ring-accent"
                      aria-label="Select download format for item"
                    >
                      <option value="pdf">PDF</option>
                      <option value="csv">CSV</option>
                      <option value="json">JSON</option>
                      <option value="txt">TXT</option>
                    </select>
                    <button
                      onClick={() => handleDownload(
                        { [key]: [item] },
                        itemDownloadFormats[`${key}-${index}`] || 'pdf',
                        `${key}-${index}.pdf`,
                        key,
                        index
                      )}
                      className="p-2 text-accent-secondary hover:text-accent relative group"
                      disabled={downloadingItem === `${key}-${index}`}
                      aria-label="Download item"
                    >
                      {downloadingItem === `${key}-${index}` ? (
                        <Loader2 size={20} className="animate-spin" />
                      ) : (
                        <FileDown size={20} />
                      )}
                      <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                        Download
                      </span>
                    </button>
                  </div>
                </div>
                {expandedSections[`${key}-${index}`] && (
                  <div className="mt-4 text-sm text-text-secondary font-body bg-primary p-4 rounded-lg">
                    {renderTableContent(item, key)}
                  </div>
                )}
              </div>
            ))}
          </div>
        ) : (
          <p className="text-text-secondary font-body">No {title.toLowerCase()} available</p>
        )}
      </section>
    );
  };

  const renderSavedItems = () => {
    const groupedItems = Object.entries(savedItems).map(([moleculeId, items]) => {
      const moleculeData = items.find(item => item.type === 'newDrug')?.data;
      return {
        moleculeId,
        title: moleculeData?.newmoleculetitle || `Molecule ${moleculeId.slice(0, 8)}`,
        smiles: moleculeData?.newSmiles || 'N/A',
        items,
      };
    });

    return (
      <div>
        <div className="flex justify-between items-center mb-8">
          <h2 className="text-2xl font-bold font-heading text-text-primary">All Saved Items</h2>
          <div className="flex space-x-4 items-center">
            <div className="relative">
              <input
                type="text"
                placeholder="Search saved items..."
                onChange={(e) => debouncedSearch(e.target.value)}
                className="pl-10 pr-4 py-2 border rounded-lg font-body text-sm text-text-primary bg-primary focus:ring-2 focus:ring-accent"
                aria-label="Search saved items"
              />
              <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-5 w-5 text-text-secondary" />
              {searchQuery && (
                <button
                  onClick={() => {
                    setSearchQuery('');
                    debouncedSearch('');
                  }}
                  className="absolute right-3 top-1/2 transform -translate-y-1/2"
                  aria-label="Clear search"
                >
                  <X className="h-5 w-5 text-text-secondary hover:text-text-primary" />
                </button>
              )}
            </div>
            <select
              value={downloadFormat}
              onChange={(e) => setDownloadFormat(e.target.value)}
              className="px-4 py-2 border rounded-lg font-body text-sm text-text-primary bg-primary focus:ring-2 focus:ring-accent"
              aria-label="Select download format for all saved items"
            >
              <option value="pdf">PDF</option>
              <option value="csv">CSV</option>
              <option value="json">JSON</option>
              <option value="txt">TXT</option>
            </select>
            <button
              onClick={() => {
                const data = {
                  newDrugs: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'newDrug').map(item => item.data)),
                  costEstimations: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'costEstimation').map(item => item.data)),
                  drugNames: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'drugName').map(item => item.data)),
                  researchPapers: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'researchPaper').map(item => item.data)),
                  generatedResearchPapers: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'generatedResearchPaper').map(item => item.data)),
                  predictDiseases: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'predictDisease').map(item => item.data)),
                  targetProteins: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'targetProtein').map(item => item.data)),
                  // toxicityResults: Object.values(savedItems).flatMap(items => items.filter(item => item.type === 'toxicityResult').map(item => item.data)),
                };
                handleDownload(
                  data,
                  downloadFormat,
                  `all-saved-items.${downloadFormat}`,
                  'saved',
                  'all'
                );
              }}
              className="flex items-center px-4 py-2 bg-accent-secondary text-white rounded-lg hover:bg-accent transition-colors relative group font-body"
              disabled={downloadingItem === 'saved-all'}
              aria-label={`Download all saved items as ${downloadFormat.toUpperCase()}`}
            >
              {downloadingItem === 'saved-all' ? (
                <Loader2 className="h-5 w-5 mr-2 animate-spin" />
              ) : (
                <Download className="h-5 w-5 mr-2" />
              )}
              Download All
              <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                Download All Saved Items
              </span>
            </button>
          </div>
        </div>
        {groupedItems.length > 0 ? (
          groupedItems.map(({ moleculeId, title, smiles, items }) => (
            <section key={moleculeId} className="mb-8" aria-label={`Molecule ${title}`}>
              <div className="flex justify-between items-center mb-4">
                <div>
                  <h3 className="text-xl font-semibold font-heading text-text-primary">{title}</h3>
                  <p className="text-sm font-body text-text-secondary">SMILES: {smiles}</p>
                </div>
                <div className="flex space-x-2 items-center">
                  <select
                    value={moleculeDownloadFormats[moleculeId] || 'pdf'}
                    onChange={(e) => setMoleculeDownloadFormats(prev => ({
                      ...prev,
                      [moleculeId]: e.target.value,
                    }))}
                    className="px-2 py-1 border rounded-lg font-body text-sm text-text-primary bg-primary focus:ring-2 focus:ring-accent"
                    aria-label="Select download format for molecule report"
                  >
                    <option value="pdf">PDF</option>
                    <option value="csv">CSV</option>
                    <option value="json">JSON</option>
                    <option value="txt">TXT</option>
                  </select>
                  <button
                    onClick={() => {
                      const moleculeData = {
                        newDrugs: items.filter(item => item.type === 'newDrug' && !pdfGeneratedItems[`newDrug-${item.data._id}`]).map(item => item.data),
                        costEstimations: items.filter(item => item.type === 'costEstimation' && !pdfGeneratedItems[`costEstimation-${item.data._id}`]).map(item => item.data),
                        drugNames: items.filter(item => item.type === 'drugName' && !pdfGeneratedItems[`drugName-${item.data._id}`]).map(item => item.data),
                        researchPapers: items.filter(item => item.type === 'researchPaper' && !pdfGeneratedItems[`researchPaper-${item.data._id}`]).map(item => item.data),
                        generatedResearchPapers: items.filter(item => item.type === 'generatedResearchPaper' && !pdfGeneratedItems[`generatedResearchPaper-${item.data._id}`]).map(item => item.data),
                        predictDiseases: items.filter(item => item.type === 'predictDisease' && !pdfGeneratedItems[`predictDisease-${item.data._id}`]).map(item => item.data),
                        targetProteins: items.filter(item => item.type === 'targetProtein' && !pdfGeneratedItems[`targetProtein-${item.data._id}`]).map(item => item.data),
                        // toxicityResults: items.filter(item => item.type === 'toxicityResult' && !pdfGeneratedItems[`toxicityResult-${item.data._id}`]).map(item => item.data),
                      };
                      handleDownload(
                        moleculeData,
                        moleculeDownloadFormats[moleculeId] || 'pdf',
                        `molecule-${title || moleculeId.slice(0, 8)}.${moleculeDownloadFormats[moleculeId] || 'pdf'}`,
                        'molecule',
                        moleculeId
                      );
                    }}
                    className="p-2 text-accent-secondary hover:text-accent relative group"
                    disabled={downloadingItem === `molecule-${moleculeId}`}
                    aria-label={`Download molecule report as ${(moleculeDownloadFormats[moleculeId] || 'pdf').toUpperCase()}`}
                  >
                    {downloadingItem === `molecule-${moleculeId}` ? (
                      <Loader2 size={20} className="animate-spin" />
                    ) : (
                      <FileDown size={20} />
                    )}
                    <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                      Download Molecule Report
                    </span>
                  </button>
                </div>
              </div>
              <div className="space-y-4">
                {items.map((item, index) => (
                  <div
                    key={index}
                    className="bg-secondary p-6 rounded-lg shadow-sm hover:shadow-md transition-shadow"
                    aria-label={`${item.type} item ${index + 1}`}
                  >
                    <div className="flex justify-between items-center">
                      <div>
                        <h4 className="font-medium font-heading text-text-primary">
                          {item.type.charAt(0).toUpperCase() + item.type.slice(1).replace(/([A-Z])/g, ' $1').trim()}
                        </h4>
                        <p className="text-sm font-body text-text-secondary mt-1">
                          Saved on: {new Date(item.createdAt).toLocaleString()}
                        </p>
                      </div>
                      <div className="flex space-x-2 items-center">
                        <button
                          onClick={() => toggleSection(`${moleculeId}-${item.type}`, index)}
                          className="p-2 text-accent-secondary hover:text-accent relative group"
                          aria-label={expandedSections[`${moleculeId}-${item.type}-${index}`] ? 'Collapse details' : 'Expand details'}
                        >
                          {expandedSections[`${moleculeId}-${item.type}-${index}`] ? <ChevronUp size={20} /> : <ChevronDown size={20} />}
                          <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                            {expandedSections[`${moleculeId}-${item.type}-${index}`] ? 'Collapse' : 'Expand'}
                          </span>
                        </button>
                        <select
                          value={itemDownloadFormats[`${moleculeId}-${item.type}-${index}`] || 'pdf'}
                          onChange={(e) =>
                            setItemDownloadFormats((prev) => ({
                              ...prev,
                              [`${moleculeId}-${item.type}-${index}`]: e.target.value,
                            }))
                          }
                          className="px-2 py-1 border rounded-lg font-body text-sm text-text-primary bg-primary focus:ring-2 focus:ring-accent"
                          aria-label="Select download format for saved item"
                        >
                          <option value="pdf">PDF</option>
                          <option value="csv">CSV</option>
                          <option value="json">JSON</option>
                          <option value="txt">TXT</option>
                        </select>
                        <button
                          onClick={() =>
                            handleDownload(
                              { [item.type]: [item.data] },
                              itemDownloadFormats[`${moleculeId}-${item.type}-${index}`] || 'pdf',
                              `saved-${item.type}-${moleculeId}-${index}.pdf`,
                              item.type,
                              index
                            )
                          }
                          className="p-2 text-accent-secondary hover:text-accent relative group"
                          disabled={downloadingItem === `${moleculeId}-${item.type}-${index}`}
                          aria-label="Download saved item"
                        >
                          {downloadingItem === `${moleculeId}-${item.type}-${index}` ? (
                            <Loader2 size={20} className="animate-spin" />
                          ) : (
                            <FileDown size={20} />
                          )}
                          <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                            Download
                          </span>
                        </button>
                      </div>
                    </div>
                    {expandedSections[`${moleculeId}-${item.type}-${index}`] && (
                      <div className="mt-4 text-sm text-text-secondary font-body bg-primary p-4 rounded-lg">
                        {renderTableContent(item.data, item.type)}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            </section>
          ))
        ) : (
          <p className="text-text-secondary font-body">No saved items match your search</p>
        )}
      </div>
    );
  };

  const debouncedSearch = debounce((value) => {
    setSearchQuery(value);
  }, 300);

  const filteredSavedItems = Object.entries(savedItems).reduce((acc, [moleculeId, items]) => {
    const filtered = items.filter(
      (item) =>
        item.type.toLowerCase().includes(searchQuery.toLowerCase()) ||
        JSON.stringify(item.data).toLowerCase().includes(searchQuery.toLowerCase())
    );
    if (filtered.length > 0) {
      acc[moleculeId] = filtered;
    }
    return acc;
  }, {});

  if (loading) {
    return (
      <div
        className="flex items-center justify-center min-h-screen bg-primary"
        aria-live="polite"
      >
        <Loader2 className="h-8 w-8 animate-spin text-accent" />
        <span className="ml-2 text-text-primary font-body">Loading data...</span>
      </div>
    );
  }

  if (error) {
    return (
      <div
        className="flex items-center justify-center min-h-screen bg-primary"
        aria-live="assertive"
      >
        <AlertCircle className="h-8 w-8 text-error" />
        <span className="ml-2 text-error font-body">{error}</span>
      </div>
    );
  }

  if (!user) {
    return (
      <div
        className="flex items-center justify-center min-h-screen bg-primary"
        aria-live="assertive"
      >
        <AlertCircle className="h-8 w-8 text-error" />
        <span className="ml-2 text-error font-body">
          Please log in to view your summary
        </span>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-primary py-12">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="mb-6">
          <div className="border-b border-text-secondary">
            <nav className="-mb-px flex space-x-8" role="tablist">
              <button
                onClick={() => handleTabChange('summary')}
                className={`${
                  activeTab === 'summary'
                    ? 'border-accent text-accent'
                    : 'border-transparent text-text-secondary hover:text-text-primary hover:border-text-secondary'
                } whitespace-nowrap py-4 px-1 border-b-2 font-heading text-sm transition-all duration-200`}
                role="tab"
                aria-selected={activeTab === 'summary'}
                aria-controls="summary-panel"
              >
                Summary
              </button>
              <button
                onClick={() => handleTabChange('saved')}
                className={`${
                  activeTab === 'saved'
                    ? 'border-accent text-accent'
                    : 'border-transparent text-text-secondary hover:text-text-primary hover:border-text-secondary'
                } whitespace-nowrap py-4 px-1 border-b-2 font-heading text-sm transition-all duration-200`}
                role="tab"
                aria-selected={activeTab === 'saved'}
                aria-controls="saved-panel"
              >
                Saved Items
              </button>
            </nav>
          </div>
        </div>

        <div className="bg-secondary rounded-xl shadow-lg p-8">
          {activeTab === 'summary' ? (
            <div id="summary-panel" role="tabpanel">
              <div className="flex justify-between items-center mb-8">
                <h1 className="text-3xl font-bold font-heading text-text-primary">Drug Discovery Summary</h1>
                <div className="flex space-x-4 items-center">
                  {summaryData && (
                    <>
                      <select
                        value={downloadFormat}
                        onChange={(e) => setDownloadFormat(e.target.value)}
                        className="px-4 py-2 border rounded-lg font-body text-sm text-text-primary bg-primary focus:ring-2 focus:ring-accent"
                        aria-label="Select download format"
                      >
                        <option value="pdf">PDF</option>
                        <option value="csv">CSV</option>
                        <option value="json">JSON</option>
                        <option value="txt">TXT</option>
                      </select>
                      <button
                        onClick={() => handleDownload(summaryData, downloadFormat, `drug-discovery-summary.${downloadFormat}`, 'summary', 0)}
                        className="flex items-center px-4 py-2 bg-accent-secondary text-white rounded-lg hover:bg-accent transition-colors relative group font-body"
                        disabled={downloadingItem === 'summary-0'}
                        aria-label={`Download full report as ${downloadFormat.toUpperCase()}`}
                      >
                        {downloadingItem === 'summary-0' ? (
                          <Loader2 className="h-5 w-5 mr-2 animate-spin" />
                        ) : (
                          <Download className="h-5 w-5 mr-2" />
                        )}
                        Download {downloadFormat.toUpperCase()}
                        <span className="absolute hidden group-hover:block bg-primary text-text-primary text-xs font-label rounded py-1 px-2 -top-8 left-1/2 transform -translate-x-1/2">
                          Download Full Report
                        </span>
                      </button>
                    </>
                  )}
                </div>
              </div>
              {summaryData && (!summaryData?.newDrugs?.length ||
                !summaryData?.costEstimations?.length ||
                !summaryData?.drugNames?.length ||
                !summaryData?.researchPapers?.length ||
                !summaryData?.generatedResearchPapers?.length ||
                !summaryData?.predictDiseases?.length ||
                !summaryData?.targetProteins?.length )
                 && (
                <div className="bg-error/10 border-l-4 border-error p-4 mb-6">
                  <div className="flex">
                    <AlertCircle className="h-5 w-5 text-error" />
                    <div className="ml-3">
                      <p className="text-sm text-error font-body">
                        Some data is missing from your summary. Complete the following sections for a comprehensive report:
                      </p>
                      <ul className="mt-2 list-disc list-inside text-sm text-error">
                        {!summaryData?.newDrugs?.length && <li>New Drug Discovery</li>}
                        {!summaryData?.costEstimations?.length && <li>Cost Estimation</li>}
                        {!summaryData?.drugNames?.length && <li>AI Naming Suggestions</li>}
                        {!summaryData?.researchPapers?.length && <li>Research Papers</li>}
                        {!summaryData?.generatedResearchPapers?.length && <li>Generated Research Papers</li>}
                        {!summaryData?.predictDiseases?.length && <li>Disease Predictions</li>}
                        {!summaryData?.targetProteins?.length && <li>Target Proteins</li>}
                        {/* {!summaryData?.toxicityResults?.length && <li>Toxicity Predictions</li>} */}
                      </ul>
                    </div>
                  </div>
                </div>
              )}
              {renderNewMoleculeSection()}
              {summaryData && (
                <>
                  {summaryData.newDrugs?.length > 0 && renderSection('New Drugs', summaryData.newDrugs, 'newDrug')}
                  {summaryData.costEstimations?.length > 0 && renderSection('Cost Estimations', summaryData.costEstimations, 'costEstimation')}
                  {summaryData.drugNames?.length > 0 && renderSection('AI Generated Names', summaryData.drugNames, 'drugName')}
                  {summaryData.researchPapers?.length > 0 && renderSection('Research Papers', summaryData.researchPapers, 'researchPaper')}
                  {summaryData.generatedResearchPapers?.length > 0 && renderSection('Generated Research Papers', summaryData.generatedResearchPapers, 'generatedResearchPaper')}
                  {summaryData.predictDiseases?.length > 0 && renderSection('Disease Predictions', summaryData.predictDiseases, 'predictDisease')}
                  {summaryData.targetProteins?.length > 0 && renderSection('Target Proteins', summaryData.targetProteins, 'targetProtein')}
                  
                </>
              )}
            </div>
          ) : (
            <div id="saved-panel" role="tabpanel">
              {renderSavedItems()}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default Summary;