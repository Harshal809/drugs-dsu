import React from 'react';
import { useAuthStore } from '../Store/auth.store.js';
import { Link } from 'react-router-dom';
import {
  Box,
  Typography,
  Grid,
  Button,
  Stack,
  Chip,
} from '@mui/material';
import {
  Home,
  TrendingUp,
  ArrowForward,
  AutoAwesome,
  Psychology,
  Calculate,
  AccountTree,
  ThreeDRotation as ViewIn3D,
  Article,
  Dangerous,
  VoiceChat,
  Summarize,
  Feed,
   
  Science, // More specific for 'discovery'
  Lightbulb, // Better for 'suggestions' or 'ideas'
 
  Hub, // Represents a complex molecular structure well
  
 
  GppBad, // A modern icon for 'bad' or 'unsafe'
  Mic, // More direct for 'audio capture'
 

} from '@mui/icons-material';
import { Toaster } from 'react-hot-toast';

function DashboardHome() {
  const { user } = useAuthStore();

  const quickActions = [
  {
    title: 'New Drug Discovery',
    description: 'Explore New Drugs Discovery',
    icon: <Science sx={{ fontSize: 28 }} />, // Changed from AutoAwesome
    link: '/dashboard/protein-structure-mutation',
    gradient: 'linear-gradient(135deg, #00F5D4 0%, #5E81F4 100%)',
    bgColor: 'rgba(0, 245, 212, 0.1)',
  },
  {
    title: 'AI Naming Suggestions',
    description: 'Let AI suggest the name for your New Drug',
    icon: <Lightbulb sx={{ fontSize: 28 }} />, // Changed from Psychology
    link: '/dashboard/ai-naming',
    gradient: 'linear-gradient(135deg, #5E81F4 0%, #00F5D4 100%)',
    bgColor: 'rgba(94, 129, 244, 0.1)',
  },
  {
    title: 'Cost Estimation',
    description: 'Estimate costs for drug development',
    icon: <Calculate sx={{ fontSize: 28 }} />, // Kept as is, it's perfect
    link: '/dashboard/cost-estimation',
    gradient: 'linear-gradient(135deg, #00F5D4 0%, #70E000 100%)',
    bgColor: 'rgba(0, 245, 212, 0.1)',
  },
  {
    title: 'Protein Structure Generation',
    description: 'Generate and analyze protein structures',
    icon: <Hub sx={{ fontSize: 28 }} />, // Changed from AccountTree
    link: '/dashboard/protein-structure',
    gradient: 'linear-gradient(135deg, #5E81F4 0%, #70E000 100%)',
    bgColor: 'rgba(94, 129, 244, 0.1)',
  },
  {
    title: 'AlphaFold 3D Predictions',
    description: 'Predict and visualize protein structures',
    icon: <ViewIn3D sx={{ fontSize: 28 }} />, // Kept as is, it's perfect
    link: '/dashboard/getalphafoldstrcture',
    gradient: 'linear-gradient(135deg, #70E000 0%, #00F5D4 100%)',
    bgColor: 'rgba(112, 224, 0, 0.1)',
  },
  {
    title: 'Research Paper Generation',
    description: 'Generate research papers using AI',
    icon: <Article sx={{ fontSize: 28 }} />, // Kept as is, it's perfect
    link: '/dashboard/ai-research-paper-generator',
    gradient: 'linear-gradient(135deg, #00F5D4 0%, #5E81F4 100%)',
    bgColor: 'rgba(0, 245, 212, 0.1)',
  },
  {
    title: 'Toxicity & Side Effects',
    description: 'Predict drug toxicity and side effects',
    icon: <GppBad sx={{ fontSize: 28 }} />, // Changed from Dangerous
    link: '/dashboard/toxicityPrediction',
    gradient: 'linear-gradient(135deg, #FF4D6D 0%, #5E81F4 100%)',
    bgColor: 'rgba(255, 77, 109, 0.1)',
  },
  {
    title: 'Audio Note Capture',
    description: 'Voice-to-text note taking system',
    icon: <Mic sx={{ fontSize: 28 }} />, // Changed from VoiceChat
    link: '/dashboard/voice-text-notes',
    gradient: 'linear-gradient(135deg, #5E81F4 0%, #00F5D4 100%)',
    bgColor: 'rgba(94, 129, 244, 0.1)',
  },
  {
    title: 'Summarization',
    description: 'AI-powered content summarization',
    icon: <Summarize sx={{ fontSize: 28 }} />, // Kept as is, it's perfect
    link: '/dashboard/summary',
    gradient: 'linear-gradient(135deg, #00F5D4 0%, #70E000 100%)',
    bgColor: 'rgba(0, 245, 212, 0.1)',
  },
  {
    title: 'NewsFeed',
    description: 'Stay updated with latest research',
    icon: <Feed sx={{ fontSize: 28 }} />, // Kept as is, it's perfect
    link: '/dashboard/live-news',
    gradient: 'linear-gradient(135deg, #70E000 0%, #5E81F4 100%)',
    bgColor: 'rgba(112, 224, 0, 0.1)',
  },
];

  return (
    <Box sx={{ 
      maxWidth: 1400, 
      mx: 'auto', 
      p: { xs: 2, md: 4 },
      minHeight: '100vh',
      background: '#0A192F', // primary color
    }}>
      {/* Header Section */}
      <Box
        sx={{
          mb: 4,
          p: 4,
          background: '#172A45', // secondary color
          border: '1px solid rgba(255, 255, 255, 0.1)',
          borderRadius: 4,
          backdropFilter: 'blur(20px)',
          position: 'relative',
          overflow: 'hidden',
          '&:before': {
            content: '""',
            position: 'absolute',
            top: 0,
            left: 0,
            right: 0,
            height: 4,
            background: 'linear-gradient(90deg, #00F5D4, #5E81F4, #70E000)',
          }
        }}
      >
        <Stack
          direction={{ xs: 'column', md: 'row' }}
          justifyContent="space-between"
          alignItems={{ xs: 'flex-start', md: 'center' }}
          spacing={3}
        >
          <Box>
            <Typography
              variant="h3"
              sx={{
                fontWeight: 800,
                mb: 1,
                background: 'linear-gradient(135deg, #00F5D4, #5E81F4)',
                backgroundClip: 'text',
                WebkitBackgroundClip: 'text',
                WebkitTextFillColor: 'transparent',
                letterSpacing: '-1px',
                fontFamily: 'Inter, Barlow, sans-serif',
              }}
            >
              Welcome back, {user?.firstName || 'Researcher'}!
            </Typography>
            <Typography 
              variant="h6" 
              sx={{
                color: '#A0A0A0', // text-secondary
                fontWeight: 300,
                fontFamily: 'Roboto, Open Sans, sans-serif',
              }}
            >
              Your drug discovery dashboard is ready for today's research.
            </Typography>
          </Box>
          <Button
            component={Link}
            to="/"
            variant="outlined"
            startIcon={<Home />}
            sx={{
              borderRadius: 3,
              px: 4,
              py: 1.5,
              borderColor: 'rgba(0, 245, 212, 0.5)', // accent
              color: '#00F5D4', // accent
              fontWeight: 600,
              transition: 'all 0.3s ease',
              fontFamily: 'Roboto, Open Sans, sans-serif',
              '&:hover': {
                borderColor: '#00F5D4', // accent
                backgroundColor: 'rgba(0, 245, 212, 0.1)',
                transform: 'translateY(-2px)',
                boxShadow: '0 8px 25px rgba(0, 245, 212, 0.3)'
              },
            }}
          >
            Back to Home
          </Button>
        </Stack>
      </Box>

      {/* Quick Actions Section */}
      <Box>
        <Typography
          variant="h4"
          sx={{
            fontWeight: 700,
            color: '#E0E0E0', // text-primary
            mb: 4,
            display: 'flex',
            alignItems: 'center',
            gap: 2,
            letterSpacing: '-0.5px',
            fontFamily: 'Inter, Barlow, sans-serif',
          }}
        >
          <TrendingUp sx={{ 
            color: '#00F5D4', // accent
            fontSize: 36
          }} />
          Quick Actions
        </Typography>
        
        <Grid container spacing={3}>
          {quickActions.map((action, index) => (
            <Grid item xs={12} sm={6} lg={2.4} key={index}>
              <Box
                component={Link}
                to={action.link}
                sx={{
                  display: 'block',
                  height: 220,
                  background: '#172A45', // secondary
                  border: '1px solid rgba(255, 255, 255, 0.1)',
                  borderRadius: 3,
                  p: 3,
                  textDecoration: 'none',
                  position: 'relative',
                  overflow: 'hidden',
                  backdropFilter: 'blur(20px)',
                  transition: 'all 0.4s cubic-bezier(0.175, 0.885, 0.32, 1.275)',
                  cursor: 'pointer',
                  '&:hover': {
                    transform: 'translateY(-12px) scale(1.03)',
                    boxShadow: '0 25px 50px rgba(0, 0, 0, 0.3)',
                    borderColor: 'rgba(255, 255, 255, 0.2)',
                    background: `linear-gradient(135deg, ${action.bgColor}, #172A45)`,
                    '& .action-icon': {
                      transform: 'scale(1.2) rotate(10deg)',
                      background: action.gradient,
                    },
                    '& .action-arrow': {
                      transform: 'translateX(8px) scale(1.1)',
                      opacity: 1,
                    },
                    '& .action-chip': {
                      transform: 'scale(1.05)',
                    }
                  },
                  '&:before': {
                    content: '""',
                    position: 'absolute',
                    top: 0,
                    left: 0,
                    right: 0,
                    height: 4,
                    background: action.gradient,
                    opacity: 0.8,
                  }
                }}
              >
                <Stack spacing={2} sx={{ height: '100%' }}>
                  <Box
                    className="action-icon"
                    sx={{
                      width: 60,
                      height: 60,
                      borderRadius: 2.5,
                      background: 'rgba(255, 255, 255, 0.1)',
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      color: '#E0E0E0', // text-primary
                      transition: 'all 0.4s ease',
                      border: '1px solid rgba(255, 255, 255, 0.1)',
                    }}
                  >
                    {action.icon}
                  </Box>
                  
                  <Box sx={{ flex: 1 }}>
                    <Typography
                      variant="h6"
                      sx={{
                        fontWeight: 700,
                        color: '#E0E0E0', // text-primary
                        mb: 1,
                        lineHeight: 1.3,
                        fontSize: '1.1rem',
                        fontFamily: 'Inter, Barlow, sans-serif',
                      }}
                    >
                      {action.title}
                    </Typography>
                    
                    <Typography
                      variant="body2"
                      sx={{ 
                        color: '#A0A0A0', // text-secondary
                        lineHeight: 1.5,
                        fontSize: '0.9rem',
                        fontFamily: 'Roboto, Open Sans, sans-serif',
                      }}
                    >
                      {action.description}
                    </Typography>
                  </Box>
                  
                  <Stack direction="row" justifyContent="space-between" alignItems="center">
                    <Chip
                      className="action-chip"
                      label="Available"
                      size="small"
                      sx={{
                        backgroundColor: 'rgba(255, 255, 255, 0.1)',
                        color: '#E0E0E0', // text-primary
                        fontWeight: 600,
                        fontSize: '0.75rem',
                        border: '1px solid rgba(255, 255, 255, 0.2)',
                        transition: 'all 0.3s ease',
                        fontFamily: 'Lato, sans-serif',
                      }}
                    />
                    <Box
                      className="action-arrow"
                      sx={{
                        width: 36,
                        height: 36,
                        borderRadius: '50%',
                        backgroundColor: 'rgba(255, 255, 255, 0.1)',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        color: '#00F5D4', // accent
                        transition: 'all 0.3s ease',
                        opacity: 0.7,
                        border: '1px solid rgba(255, 255, 255, 0.2)',
                      }}
                    >
                      <ArrowForward fontSize="small" />
                    </Box>
                  </Stack>
                </Stack>
              </Box>
            </Grid>
          ))}
        </Grid>
      </Box>

      <Toaster 
        position="bottom-right"
        toastOptions={{
          style: {
            background: '#172A45', // secondary
            color: '#E0E0E0', // text-primary
            border: '1px solid rgba(255, 255, 255, 0.2)',
            borderRadius: '16px',
            backdropFilter: 'blur(20px)'
          },
        }}
      />
    </Box>
  );
}

export default DashboardHome;