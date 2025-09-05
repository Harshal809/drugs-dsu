import NewsAPI from 'newsapi';
import News from '../models/news.model.js';

const newsapi = new NewsAPI(process.env.NEWS_API_KEY);

export const getTopHeadlines = async (req, res) => {
  try {
    const response = await newsapi.v2.topHeadlines({
      q: req.query.q || 'drugs OR pharmaceuticals OR medicine OR healthcare', // Focus on drug-related keywords
      category: req.query.category || 'health', // Use health category for drug-related news
      language: 'en',
      country: 'us',
      pageSize: 10, // Ensure at least 10 articles are fetched
    });

    const newsArticles = response.articles.map(article => ({
      title: article.title,
      description: article.description,
      url: article.url,
      source: article.source.name,
      publishedAt: article.publishedAt,
      category: req.query.category || 'health',
    }));

    await News.insertMany(newsArticles);

    res.json(response);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
};

export const getEverything = async (req, res) => {
  try {
    const response = await newsapi.v2.everything({
      q: req.query.q || 'drugs OR pharmaceuticals OR medicine OR healthcare', // Focus on drug-related keywords
      domains: 'bbc.co.uk,techcrunch.com',
      from: req.query.from,
      to: req.query.to,
      language: 'en',
      sortBy: 'relevancy',
      page: req.query.page || 1,
      pageSize: 10, // Ensure at least 10 articles are fetched
    });

    res.json(response);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
};

export const getSources = async (req, res) => {
  try {
    const response = await newsapi.v2.sources({
      category: 'health', // Focus on health sources for drug-related news
      language: 'en',
      country: 'us',
    });
    res.json(response);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
};

export const getSavedNews = async (req, res) => {
  try {
    const news = await News.find({ category: 'health' }) // Filter saved news by health category
      .sort({ publishedAt: -1 })
      .limit(50);
    res.json(news);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
};