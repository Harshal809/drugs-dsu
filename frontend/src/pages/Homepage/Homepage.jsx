
import { useState, useEffect, useRef } from "react"
import {
  ArrowRight,
  Atom,
  Beaker,
  Brain,
  Database,
  FileText,
  Microscope,
  Pill,
  Search,
  Shield,
  Mic,
  Users,
  MessageCircle,
  BarChart3,
  Star,
  Mail,
  Phone,
  MapPin,
  Zap,
  Target,
  TrendingUp,
  ChevronDown,
  Play,
  Sparkles,
  Cpu,
  Activity,
  Globe,
  
  DollarSign,
  
 
  BadgeInfo,
 
  BookOpen,
  AlertTriangle,

  ClipboardList,
  Newspaper,
  Bot,
} from "lucide-react"
import DNA from "./dna"
import "./homepage.css"
import Navbar from "../../components/Navbar"
import {Link} from "react-router-dom";

function Homepage() {
  const [activeFeature, setActiveFeature] = useState(null)
  const [mousePosition, setMousePosition] = useState({ x: 0, y: 0 })
  const [counters, setCounters] = useState({ stat1: 0, stat2: 0, stat3: 0, stat4: 0 })
  const [isVisible, setIsVisible] = useState({
    hero: false,
    stats: false,
    features: false,
    testimonials: false,
    cta: false,
  })

  // Refs for scroll animations
  const heroRef = useRef(null)
  const statsRef = useRef(null)
  const featuresRef = useRef(null)
  const testimonialsRef = useRef(null)
  const ctaRef = useRef(null)

  // Intersection Observer for animations
  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const id = entry.target.id
            setIsVisible((prev) => ({ ...prev, [id]: true }))
          }
        })
      },
      { threshold: 0.1 },
    )

    const refs = [heroRef, statsRef, featuresRef, testimonialsRef, ctaRef]
    refs.forEach((ref) => {
      if (ref.current) observer.observe(ref.current)
    })

    return () => observer.disconnect()
  }, [])

  // Mouse tracking for parallax effects (disabled on mobile)
  // useEffect(() => {
  //   const handleMouseMove = (e) => {
  //     if (window.innerWidth >= 768) {
  //       setMousePosition({
  //         x: (e.clientX / window.innerWidth) * 100,
  //         y: (e.clientY / window.innerHeight) * 100,
  //       })
  //     }
  //   }

  //   window.addEventListener("mousemove", handleMouseMove)
  //   return () => window.removeEventListener("mousemove", handleMouseMove)
  // }, [])

  // Counter animations
  useEffect(() => {
    if (isVisible.stats) {
      const animateCounter = (target, key, duration = 2000) => {
        let start = 0
        const increment = target / (duration / 16)
        const timer = setInterval(() => {
          start += increment
          if (start >= target) {
            setCounters((prev) => ({ ...prev, [key]: target }))
            clearInterval(timer)
          } else {
            setCounters((prev) => ({ ...prev, [key]: Math.floor(start) }))
          }
        }, 16)
      }

      animateCounter(57, "stat1")
      animateCounter(60, "stat2")
      animateCounter(3, "stat3")
      setTimeout(() => setCounters((prev) => ({ ...prev, stat4: 24 })), 1000)
    }
  }, [isVisible.stats])

  const features = [
  {
    icon: <Brain className="w-6 h-6" />,
    title: "New Drug Discovery",
    description:
      "Generate new drug candidates with the power of Gemini and RDkit with accuracy range of 82.4 - 86.2%.",
    gradient: "from-accent-secondary to-accent",
  },
  {
    icon: <DollarSign className="w-6 h-6" />,
    title: "Cost Estimation",
    description:
      "Estimate the cost of new drug candidates by complex algorithms and Gemini.",
    gradient: "from-accent to-success",
  },
  {
    icon: <Beaker className="w-6 h-6" />,
    title: "Protein Structure Generation",
    description:
      "Generate drug candidates structures with the power of Gemini and Nvidia.",
    gradient: "from-success to-accent-secondary",
  },
  {
    icon: <Atom className="w-6 h-6" />,
    title: "AlphaFold 3D Predictions",
    description:
      "Get 3D structures of proteins with the power of DeepMind AlphaFold and Gemini.",
    gradient: "from-accent-secondary to-accent",
  },
  {
    icon: <BadgeInfo className="w-6 h-6" />,
    title: "AI Naming Suggestion",
    description:
      "Generate intelligent, systematic names for your novel drug candidates.",
    gradient: "from-accent to-success",
  },
  {
    icon: <FileText className="w-6 h-6" />,
    title: "Research Paper Generation",
    description:
      "Generate research papers in IEEE format for your novel drug candidates with the power of Gemini.",
    gradient: "from-success to-accent-secondary",
  },
 
  {
    icon: <AlertTriangle className="w-6 h-6" />,
    title: "Side Effects Prediction",
    description:
      "Predict potential toxicity and side effects in real-time by Gemini.",
    gradient: "from-accent-secondary to-accent",
  },
  {
    icon: <Mic className="w-6 h-6" />,
    title: "Audio Note Capture",
    description:
      "Capture observations, and AI will auto-transcribe and summarize them for you.",
    gradient: "from-accent to-success",
  },
  {
    icon: <ClipboardList className="w-6 h-6" />,
    title: "Project Summarization",
    description:
      "Summarize your complete process of drug discovery with the power of Gemini.",
    gradient: "from-success to-accent-secondary",
  },
  {
    icon: <Newspaper className="w-6 h-6" />,
    title: "Pharma & Biotech Newsfeed",
    description:
      "Get real-time updates on pharmaceutical and biotech news with the power of Gemini.",
    gradient: "from-accent-secondary to-accent",
  },
  {
    icon: <Bot className="w-6 h-6" />,
    title: "J.A.R.V.I.S.",
    description:
      "Justified AI for Research, Validation & Intelligent Synthesis",
    gradient: "from-accent to-success",
  },
   {
    icon: <Newspaper className="w-6 h-6" />,

    title: "News Feed",
    description:
      "Get real-time updates on pharmaceutical and biotech news",
    gradient: "from-success to-accent-secondary",
  },
  // Optional AI Chatbot block
  // {
  //   icon: <MessageCircle className="w-6 h-6" />,
  //   title: "AI Chatbot",
  //   description: "A Gemini-powered assistant to explain complex molecular structures.",
  //   gradient: "from-success to-accent-secondary",
  // },
];

  return (
    <div className="min-h-screen bg-primary">
      {/* Navigation */}

    <Navbar/>

     
      {/* Floating Particles - Hidden on mobile */}
      <div className="fixed inset-0 pointer-events-none overflow-hidden hidden md:block">
        {[...Array(40)].map((_, i) => (
          <div
            key={i}
            className="particle animate-particle"
            style={{
              left: Math.random() * 100 + "%",
              animationDuration: Math.random() * 15 + 15 + "s",
              animationDelay: Math.random() * 15 + "s",
              opacity: Math.random() * 0.6 + 0.3,
            }}
          />
        ))}
      </div>

      {/* Hero Section */}
      <section
        id="hero"
        ref={heroRef}
        className="relative min-h-screen flex items-center justify-center parallax-bg pt-16"
      >
        {/* Molecular Orbits - Hidden on mobile */}
        <div className="absolute inset-0 overflow-hidden hidden lg:block">
          <div className="molecular-orbit w-96 h-96 top-1/4 left-1/4"></div>
          <div className="molecular-orbit w-80 h-80 top-1/3 right-1/4" style={{ animationDelay: '10s' }}></div>
          <div className="molecular-orbit w-64 h-64 bottom-1/4 left-1/3" style={{ animationDelay: '20s' }}></div>
        </div>

        {/* Floating Background Elements - Hidden on mobile */}
        <div className="absolute inset-0 overflow-hidden hidden lg:block">
          <div
            className="absolute top-20 left-10 floating-icon delay-100"
            style={{
              transform: `translate(${mousePosition.x * 0.02}px, ${mousePosition.y * 0.02}px)`,
            }}
          >
            <Brain className="w-24 h-24 opacity-30 text-accent-secondary" />
          </div>
          <div
            className="absolute top-40 right-20 floating-icon delay-300"
            style={{
              transform: `translate(${mousePosition.x * -0.03}px, ${mousePosition.y * 0.03}px)`,
            }}
          >
            <Beaker className="w-32 h-32 opacity-25 text-accent" />
          </div>
          <div
            className="absolute bottom-40 left-20 floating-icon delay-500"
            style={{
              transform: `translate(${mousePosition.x * 0.025}px, ${mousePosition.y * -0.025}px)`,
            }}
          >
            <Microscope className="w-28 h-28 opacity-30 text-success" />
          </div>
        </div>

        <div className="max-w-7xl mx-auto px-4 mt-5 sm:px-6 lg:px-8 relative">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 lg:gap-12 items-center">
            {/* Left Content */}
            <div className="space-y-6 lg:space-y-8 text-center lg:text-left">
              <div
                className={`inline-flex items-center px-4 lg:px-6 py-2 lg:py-3 rounded-full glass-effect border-glow ${isVisible.hero ? "animate-bounce-in" : "opacity-0"}`}
              >
                <div className="w-3 h-3 rounded-full bg-success animate-pulse-glow mr-3"></div>
                <span className="font-code text-xs lg:text-sm text-text-secondary">
                  Next-Gen AI for Drug Discovery
                </span>
              </div>

              <h1
                className={`hero-title text-3xl sm:text-4xl lg:text-5xl xl:text-7xl leading-tight ${isVisible.hero ? "animate-slide-up delay-200" : "opacity-0"}`}
              >
                <span className="text-text-primary">Revolutionizing</span>
                <br />
                <span className="gradient-text animate-text-glow">Drug Discovery</span>
                <br />
                <span className="text-text-primary">with Gen AI</span>
              </h1>

              <p
                className={`font-body text-lg lg:text-xl leading-relaxed text-text-secondary max-w-lg mx-auto lg:mx-0 ${isVisible.hero ? "animate-slide-up delay-400" : "opacity-0"}`}
              >
                Harness the power of generative AI to accelerate your research, reduce costs, and discover
                breakthrough medications faster than ever before.
              </p>

              <div
                className={`flex flex-col sm:flex-row gap-4 justify-center lg:justify-start ${isVisible.hero ? "animate-slide-up delay-600" : "opacity-0"}`}
              >
                <Link to="/dashboard" className="btn-primary flex items-center justify-center group">
                  <Sparkles className="w-5 h-5 mr-2 animate-wave" />
                  Start Your Discovery Journey
                  <ArrowRight className="w-5 h-5 ml-2 group-hover:translate-x-2 transition-transform duration-300" />
                </Link>
                <button className="btn-secondary flex items-center justify-center group">
                  <Play className="w-5 h-5 mr-2" />
                  Watch Demo
                  <Zap className="w-5 h-5 ml-2 group-hover:rotate-12 transition-transform duration-300" />
                </button>
              </div>

              <div
                className={`flex items-center space-x-4 lg:space-x-6 justify-center lg:justify-start ${isVisible.hero ? "animate-slide-up delay-800" : "opacity-0"}`}
              >
                <div className="flex -space-x-2">
                  {[1, 2, 3, 4].map((i) => (
                    <div
                      key={i}
                      className="w-10 h-10 lg:w-12 lg:h-12 rounded-full border-2 border-accent bg-secondary flex items-center justify-center font-semibold hover-lift cursor-pointer text-accent text-sm lg:text-base"
                      style={{ animationDelay: `${i * 100}ms` }}
                    >
                      {i}
                    </div>
                  ))}
                </div>
                <div className="font-body text-xs lg:text-sm text-text-secondary">
                  <span className="font-semibold gradient-text">500+</span> researchers already using our platform
                </div>
              </div>
            </div>

            {/* Right Content - DNA Component */}
            <div className={`relative ${isVisible.hero ? "animate-slide-in-right delay-300" : "opacity-0"}`}>
              <div className="w-full h-[400px] lg:h-[600px] flex items-center justify-center">
                <DNA />
              </div>
            </div>
          </div>
        </div>

        {/* Scroll Indicator */}
        <div className="scroll-indicator">
          <ChevronDown className="w-6 h-6 text-accent animate-glow-pulse" />
        </div>
      </section>

      {/* Stats Section */}
      {/* <section id="stats" ref={statsRef} className="py-16 lg:py-20 bg-secondary">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-2 lg:grid-cols-4 gap-6 lg:gap-8">
            {[
              {
                key: "stat1",
                value: counters.stat1,
                suffix: "%",
                label: "Average Accuracy",
                icon: <Zap className="w-6 h-6 lg:w-8 lg:h-8" />,
              },
              {
                key: "stat2",
                value: counters.stat2,
                suffix: "%",
                label: "Cost Reduction",
                icon: <TrendingUp className="w-6 h-6 lg:w-8 lg:h-8" />,
              },
              {
                key: "stat3",
                value: counters.stat3,
                suffix: "x",
                label: "More Candidates",
                icon: <Target className="w-6 h-6 lg:w-8 lg:h-8" />,
              },
              {
                key: "stat4",
                value: counters.stat4,
                suffix: "/7",
                label: "AI Assistance",
                icon: <Activity className="w-6 h-6 lg:w-8 lg:h-8" />,
              },
            ].map((stat, index) => (
              <div
                key={index}
                className={`text-center glass-effect rounded-xl p-4 lg:p-6 hover-lift border-glow ${isVisible.stats ? `animate-bounce-in delay-${index * 200}` : "opacity-0"}`}
              >
                <div className="flex justify-center mb-3 lg:mb-4 text-accent animate-glow-pulse">
                  {stat.icon}
                </div>
                <div className="stats-counter text-2xl lg:text-4xl mb-2">
                  {stat.value}
                  {stat.suffix}
                </div>
                <div className="font-body text-sm lg:text-base text-text-secondary">
                  {stat.label}
                </div>
              </div>
            ))}
          </div>
        </div>
      </section> */}

      {/* Features Section */}
      <section id="features" ref={featuresRef} className="py-16 lg:py-20 bg-primary">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className={`text-center mb-12 lg:mb-16 ${isVisible.features ? "animate-slide-up" : "opacity-0"}`}>
            <div className="inline-flex items-center px-4 lg:px-6 py-2 lg:py-3 rounded-full glass-effect border-glow mb-6">
              <Sparkles className="w-4 h-4 mr-2 text-accent animate-wave" />
              <span className="font-code text-xs lg:text-sm text-text-secondary">
                Powerful Capabilities
              </span>
            </div>
            <h2 className="font-heading font-bold text-3xl lg:text-4xl mb-4 gradient-text animate-text-glow">Accelerate Your Research</h2>
            <p className="font-body text-lg lg:text-xl max-w-3xl mx-auto text-text-secondary">
              Our platform combines cutting-edge AI with intuitive tools to supercharge your drug discovery workflow.
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 lg:gap-8">
            {features.map((feature, index) => (
              <div
                key={index}
                className={`feature-card glass-effect rounded-xl p-6 lg:p-8 border-glow cursor-pointer ${isVisible.features ? `animate-slide-up delay-${index * 100}` : "opacity-0"}`}
                onMouseEnter={() => setActiveFeature(index)}
                onMouseLeave={() => setActiveFeature(null)}
              >
                <div
                  className={`inline-flex items-center justify-center w-12 h-12 rounded-lg mb-6 bg-gradient-to-r ${feature.gradient} animate-glow-pulse`}
                >
                  {feature.icon}
                </div>

                <h3
                  className={`font-heading font-semibold text-lg lg:text-xl mb-3 transition-colors duration-300 ${activeFeature === index ? "text-accent" : "text-text-primary"
                    }`}
                >
                  {feature.title}
                </h3>

                <p className="font-body text-sm lg:text-base leading-relaxed mb-4 text-text-secondary">
                  {feature.description}
                </p>

                <div
                  className="flex items-center font-body font-medium text-accent transition-transform duration-300"
                  style={{
                    transform: activeFeature === index ? "translateX(8px)" : "translateX(0)",
                  }}
                >
                  <span className="text-sm">Learn more</span>
                  <ArrowRight className="w-4 h-4 ml-2" />
                </div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* ... keep existing code (testimonials section, CTA section, and footer) the same ... */}
      {/* Testimonials Section */}
      <section id="testimonials" ref={testimonialsRef} className="py-16 lg:py-20 bg-secondary">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className={`text-center mb-12 lg:mb-16 ${isVisible.testimonials ? "animate-slide-up" : "opacity-0"}`}>
            <div className="inline-flex items-center px-4 lg:px-6 py-2 lg:py-3 rounded-full glass-effect border-glow mb-6">
              <Users className="w-4 h-4 mr-2 text-success animate-wave" />
              <span className="font-code text-xs lg:text-sm text-text-secondary">
                Success Stories
              </span>
            </div>
            <h2 className="font-heading font-bold text-3xl lg:text-4xl gradient-text animate-text-glow">What Researchers Are Saying</h2>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 lg:gap-8">
            {[
              {
                quote:
                  "This platform reduced our initial discovery phase from 6 months to just 2 weeks. Unbelievable!",
                author: "Dr. Sarah Chen",
                role: "Lead Researcher",
                affiliation: "Stanford University",
                rating: 5,
              },
              {
                quote: "The AI-generated molecules had better binding affinity than our manually designed compounds.",
                author: "Prof. Raj Patel",
                role: "Principal Investigator",
                affiliation: "MIT Bioengineering",
                rating: 5,
              },
              {
                quote:
                  "Finally, a tool that bridges the gap between computational chemistry and practical drug development.",
                author: "Dr. Elena Rodriguez",
                role: "Senior Scientist",
                affiliation: "Novartis Pharmaceuticals",
                rating: 5,
              },
            ].map((testimonial, index) => (
              <div
                key={index}
                className={`glass-effect rounded-xl p-6 lg:p-8 border-glow hover-lift ${isVisible.testimonials ? `animate-bounce-in delay-${index * 200}` : "opacity-0"}`}
              >
                <div className="flex mb-4">
                  {[...Array(testimonial.rating)].map((_, i) => (
                    <Star
                      key={i}
                      className="w-5 h-5 fill-current text-accent animate-wave"
                      style={{
                        animationDelay: `${i * 100}ms`,
                      }}
                    />
                  ))}
                </div>

                <p className="font-body text-base lg:text-lg mb-6 leading-relaxed text-text-primary">
                  "{testimonial.quote}"
                </p>

                <div className="flex items-center">
                  <div className="w-10 h-10 lg:w-12 lg:h-12 rounded-full flex items-center justify-center font-heading font-semibold mr-4 hover-lift bg-gradient-to-r from-accent to-accent-secondary text-primary">
                    {testimonial.author.charAt(0)}
                  </div>
                  <div>
                    <h4 className="font-heading font-semibold text-sm lg:text-base text-text-primary">
                      {testimonial.author}
                    </h4>
                    <p className="font-body text-xs lg:text-sm text-text-secondary">
                      {testimonial.role}
                    </p>
                    <p className="font-body text-xs lg:text-sm text-text-secondary/70">
                      {testimonial.affiliation}
                    </p>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section
        id="cta"
        ref={ctaRef}
        className="py-16 lg:py-20 relative overflow-hidden bg-gradient-to-r from-accent via-accent-secondary to-success animate-gradient"
      >
        <div className="absolute inset-0 opacity-20 hidden lg:block">
          <div className="absolute top-10 left-10 floating-icon">
            <Atom className="w-16 h-16 lg:w-20 lg:h-20 text-white" />
          </div>
          <div className="absolute bottom-10 right-10 floating-icon delay-300">
            <Brain className="w-20 h-20 lg:w-24 lg:h-24 text-white" />
          </div>
        </div>

        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center relative">
          <div
            className={`inline-flex items-center px-4 lg:px-6 py-2 lg:py-3 rounded-full glass-effect border-glow mb-6 ${isVisible.cta ? "animate-bounce-in" : "opacity-0"}`}
          >
            <div className="w-3 h-3 rounded-full bg-success animate-pulse-glow mr-3"></div>
            <span className="font-code text-xs lg:text-sm text-white">Limited Time Offer</span>
          </div>

          <h2
            className={`font-heading font-bold text-3xl lg:text-4xl text-white mb-6 animate-text-glow ${isVisible.cta ? "animate-slide-up delay-200" : "opacity-0"}`}
          >
            Ready to Transform Drug Discovery?
          </h2>

          <p
            className={`font-body text-lg lg:text-xl text-white/90 mb-8 lg:mb-10 max-w-3xl mx-auto ${isVisible.cta ? "animate-slide-up delay-400" : "opacity-0"}`}
          >
            Join thousands of researchers using our platform to accelerate their drug discovery process with the power
            of generative AI.
          </p>

          <div
            className={`flex flex-col sm:flex-row justify-center gap-4 ${isVisible.cta ? "animate-slide-up delay-600" : "opacity-0"}`}
          >
            <button className="bg-white text-primary px-6 lg:px-8 py-3 lg:py-4 rounded-lg font-heading font-semibold hover-lift flex items-center justify-center group border-2 border-white">
              <Globe className="w-5 h-5 mr-2" />
              Get Started Now
              <ArrowRight className="w-5 h-5 ml-2 group-hover:translate-x-2 transition-transform duration-300" />
            </button>
            <button className="border-2 border-white text-white px-6 lg:px-8 py-3 lg:py-4 rounded-lg font-heading font-semibold hover-lift flex items-center justify-center group hover:bg-white hover:text-primary transition-all duration-300">
              <MessageCircle className="w-5 h-5 mr-2" />
              Contact Our Team
              <Sparkles className="w-5 h-5 ml-2 group-hover:rotate-12 transition-transform duration-300" />
            </button>
          </div>
        </div>
      </section>

      {/* Footer */}
      <footer className="py-12 lg:py-16 bg-primary border-t border-accent/20">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-8 lg:gap-12">
            <div className="space-y-4 col-span-1 sm:col-span-2 lg:col-span-1">
              <div className="flex items-center">
                <div className="animate-rotate-slow">
                  <Atom className="w-8 h-8 mr-3 text-accent animate-glow-pulse" />
                </div>
                <span className="font-heading font-bold text-lg lg:text-xl gradient-text">Drug Discovery AI</span>
              </div>
              <p className="font-body leading-relaxed text-text-secondary text-sm lg:text-base">
                Revolutionizing drug discovery with generative AI to make the process faster, cheaper, and more
                efficient.
              </p>
            </div>

            <div>
              <h3 className="font-heading font-semibold text-base lg:text-lg mb-4 text-text-primary">
                Features
              </h3>
              <ul className="space-y-3">
                {["Molecule Generation", "Visualization Tools", "Research Automation"].map((item, index) => (
                  <li key={index}>
                    <a
                      href="#features"
                      className="font-body hover-lift inline-block transition-all duration-300 text-text-secondary hover:text-accent text-sm lg:text-base"
                    >
                      {item}
                    </a>
                  </li>
                ))}
              </ul>
            </div>

            <div>
              <h3 className="font-heading font-semibold text-base lg:text-lg mb-4 text-text-primary">
                Company
              </h3>
              <ul className="space-y-3">
                {["About Us", "Careers", "Blog"].map((item, index) => (
                  <li key={index}>
                    <a
                      href="#"
                      className="font-body hover-lift inline-block transition-all duration-300 text-text-secondary hover:text-accent text-sm lg:text-base"
                    >
                      {item}
                    </a>
                  </li>
                ))}
              </ul>
            </div>

            <div id="contact">
              <h3 className="font-heading font-semibold text-base lg:text-lg mb-4 text-text-primary">
                Contact
              </h3>
              <div className="space-y-3">
                {[
                  { icon: <Mail className="w-4 h-4" />, text: "alokchaturvedi190@gmail.com" },
                  { icon: <Phone className="w-4 h-4" />, text: "+91 9975175098" },
                  { icon: <Mail className="w-4 h-4" />, text: "tbhangale9@gmail.com" },
                  { icon: <Phone className="w-4 h-4" />, text: "+91 8766816061" },
                  { icon: <MapPin className="w-4 h-4" />, text: "Pune, Maharashtra, India" },
                ].map((contact, index) => (
                  <div
                    key={index}
                    className="flex items-center hover-lift cursor-pointer transition-all duration-300 text-text-secondary hover:text-accent"
                  >
                    <span className="mr-2">{contact.icon}</span>
                    <span className="font-body text-xs lg:text-sm">{contact.text}</span>
                  </div>
                ))}
              </div>
            </div>
          </div>

          <div className="mt-8 lg:mt-12 pt-6 lg:pt-8 border-t border-accent/20 flex flex-col md:flex-row justify-between items-center">
            <p className="font-body text-text-secondary text-sm lg:text-base">
              © {new Date().getFullYear()} Drug Discovery AI. All rights reserved.
            </p>
            <div className="mt-4 md:mt-0">
              <ul className="flex space-x-4 lg:space-x-6 text-sm">
                {["Terms", "Privacy"].map((item, index) => (
                  <li key={index}>
                    <a
                      href="#"
                      className="font-body hover-lift inline-block transition-all duration-300 text-text-secondary hover:text-accent"
                    >
                      {item}
                    </a>
                  </li>
                ))}
              </ul>
            </div>
          </div>
        </div>
      </footer>
    </div>
  )
}

export default Homepage
