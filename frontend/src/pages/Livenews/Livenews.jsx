import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { useAuthStore } from '../../Store/auth.store.js';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000/api";
const axiosInstance = axios.create({
  baseURL: import.meta.env.mode === "development" ? API_BASE_URL : "/api",
  withCredentials: true,
});


function LiveNews() {
  const [liveArticles, setLiveArticles] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  const token = user ? user.token || user._id || JSON.stringify(user) : null;

  useEffect(() => {
    const fetchLiveNews = async () => {
      if (!token) return;
      try {
        setLoading(true);
        const response = await axiosInstance.get(
          '/news/top-headlines?q=drugs&category=health',
          {
            headers: {
              Authorization: `Bearer ${token}`,
            },
          }
        );
        setLiveArticles(response.data.articles || []);
        setLoading(false);
      } catch (err) {
        setError(err.response?.data?.message || 'Failed to fetch live news');
        setLoading(false);
      }
    };

    const initializeLiveNews = async () => {
      await checkAuth();
      if (!useAuthStore.getState().user) {
        setError('Please log in to view live news');
        setLoading(false);
        return;
      }
      fetchLiveNews();
      const interval = setInterval(fetchLiveNews, 300000); // Refetch every 5 minutes
      return () => clearInterval(interval);
    };

    initializeLiveNews();
  }, [checkAuth, token]);

  if (checkingAuth) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-primary">
        <div className="text-center space-y-4">
          <div className="inline-block animate-spin rounded-full h-10 w-10 border-4 border-t-accent border-secondary"></div>
          <p className="text-text-primary text-lg font-body">Verifying authentication...</p>
        </div>
      </div>
    );
  }

  if (!token) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-primary">
        <div className="text-center p-8 bg-secondary rounded-xl shadow-lg max-w-md space-y-6">
          <div className="space-y-4">
            <h2 className="text-3xl font-heading font-bold text-accent">Access Required</h2>
            <p className="text-text-secondary text-lg font-body">Sign in to view real-time health updates</p>
          </div>
          <button
            className="w-full px-6 py-3 bg-accent-secondary text-primary rounded-lg hover:opacity-90 transition-all duration-300 font-semibold font-body shadow-md"
            onClick={() => window.location.href = '/login'}
          >
            Sign In
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="max-w-6xl mx-auto my-12 bg-secondary rounded-xl shadow-lg overflow-hidden border border-accent-secondary">
      <header className="p-6 bg-primary text-white border-b border-accent">
        <div className="flex items-center gap-4">
          <div>
            <h1 className="text-2xl font-heading font-bold tracking-tight text-accent">Drugs News Live</h1>
            <p className="mt-1 text-sm font-body text-text-secondary">Real-time medical breakthroughs and updates</p>
          </div>
        </div>
      </header>

      <main className="relative">
        {loading ? (
          <div className="grid gap-6 p-6">
            {[...Array(5)].map((_, i) => (
              <div key={i} className="animate-pulse space-y-4 bg-primary rounded-lg p-4">
                <div className="flex gap-4">
                  <div className="h-6 w-24 bg-secondary rounded-full"></div>
                  <div className="flex-1 space-y-2">
                    <div className="h-5 bg-secondary rounded w-3/4"></div>
                    <div className="h-4 bg-secondary rounded w-full"></div>
                    <div className="h-4 bg-secondary rounded w-1/2"></div>
                  </div>
                </div>
              </div>
            ))}
          </div>
        ) : error ? (
          <div className="p-6 m-6 bg-error rounded-xl flex items-center gap-4 border border-red-700">
            <svg className="w-8 h-8 text-primary flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2"
                    d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"/>
            </svg>
            <div>
              <h3 className="font-heading text-primary font-medium">Error loading news</h3>
              <p className="text-sm text-primary mt-1 font-body">{error}</p>
            </div>
          </div>
        ) : (
          <div className="divide-y divide-primary">
            {liveArticles.slice(0, 10).map((article, index) => (
              <article
                key={index}
                className="group p-6 hover:bg-primary transition-all duration-200"
              >
                <div className="flex gap-4">
                  <div className="w-28 flex-shrink-0">
                    <time className="text-sm font-label text-text-secondary">
                      {new Date(article.publishedAt).toLocaleTimeString([], {
                        hour: 'numeric',
                        minute: '2-digit',
                        hour12: true,
                      })}
                    </time>
                  </div>
                  <div className="flex-1 space-y-3">
                    <h2 className="text-lg font-heading font-semibold text-text-primary">
                      <a
                        href={article.url}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="hover:text-accent transition-colors duration-200"
                      >
                        {article.title}
                      </a>
                    </h2>
                    {article.urlToImage && (
                      <img
                        src={article.urlToImage}
                        alt={article.title}
                        className="w-full h-48 object-cover rounded-lg mb-3 shadow-sm transition-transform duration-300 group-hover:scale-[1.02] border-2 border-secondary"
                        loading="lazy"
                        onError={(e) => e.target.style.display = 'none'}
                      />
                    )}
                    <p className="text-text-secondary leading-relaxed font-body">
                      {article.description}
                    </p>
                    <div className="flex items-center justify-between text-sm">
                      <span className="font-label text-accent">
                        {article.source?.name}
                      </span>
                      <a
                        href={article.url}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="flex items-center gap-2 text-accent-secondary hover:text-accent transition-colors duration-200 font-body"
                      >
                        Read Full Article
                        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M17 8l4 4m0 0l-4 4m4-4H3"/>
                        </svg>
                      </a>
                    </div>
                  </div>
                </div>
              </article>
            ))}
          </div>
        )}

        <footer className="p-4 text-center text-sm text-text-secondary bg-primary font-body rounded-b-xl">
          <div className="flex items-center justify-center gap-2">
            <span>Limited news access due to API restrictions. Timestamps may be unavailable in the free version.</span>
          </div>
        </footer>
      </main>
    </div>
  );
}

export default LiveNews;