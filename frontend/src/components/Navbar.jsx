"use client"

import { useState, useEffect } from "react"
import { Link, useNavigate } from "react-router-dom"
import { useAuthStore } from "../Store/auth.store.js"
import { LogOut, UserPlus, LogIn, Menu, X, Atom, User, ChevronRight, Bell } from "lucide-react"
import { motion, AnimatePresence, useScroll, useMotionValueEvent } from "framer-motion"

function Navbar() {
  const { user, logout } = useAuthStore()
  const navigate = useNavigate()
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)
  const [isScrolled, setIsScrolled] = useState(false)
  const [activeSection, setActiveSection] = useState("")
  const [isLoaded, setIsLoaded] = useState(false)
  const { scrollY } = useScroll()

  // Simplified animation variants - much faster
  const navVariants = {
    hidden: { opacity: 0 },
    visible: { 
      opacity: 1,
      transition: { 
        duration: 0.3, // Reduced from 0.8
        ease: "easeOut"
      }
    }
  }

  const logoVariants = {
    hidden: { scale: 0.8, opacity: 0 },
    visible: { 
      scale: 1, 
      opacity: 1,
      transition: { 
        duration: 0.4, // Reduced from 0.8
        ease: "easeOut"
      }
    }
  }

  const buttonVariants = {
    hover: { 
      scale: 1.05,
      transition: { duration: 0.2 }
    },
    tap: { scale: 0.95 }
  }

  const mobileMenuVariants = {
    hidden: { 
      opacity: 0, 
      y: -20,
      transition: { duration: 0.2 }
    },
    visible: { 
      opacity: 1, 
      y: 0,
      transition: { 
        duration: 0.3,
        ease: "easeOut"
      }
    }
  }

  // Optimized scroll handler with throttling
  useMotionValueEvent(scrollY, "change", (latest) => {
    setIsScrolled(latest > 50)
  })

  // Delayed setup of non-critical features
  useEffect(() => {
    setIsLoaded(true)
    
    // Throttled scroll handler for active section
    let ticking = false
    const handleScroll = () => {
      if (!ticking) {
        requestAnimationFrame(() => {
          const sections = ["features", "testimonials", "contact"]
          const scrollPosition = window.scrollY + 100

          for (const section of sections) {
            const element = document.getElementById(section)
            if (element) {
              const offsetTop = element.offsetTop
              const offsetHeight = element.offsetHeight

              if (scrollPosition >= offsetTop && scrollPosition < offsetTop + offsetHeight) {
                setActiveSection(section)
                break
              }
            }
          }
          ticking = false
        })
        ticking = true
      }
    }

    // Add scroll listener after component is loaded
    const timer = setTimeout(() => {
      window.addEventListener("scroll", handleScroll, { passive: true })
    }, 100)

    return () => {
      clearTimeout(timer)
      window.removeEventListener("scroll", handleScroll)
    }
  }, [])

  const handleLogout = async () => {
    await logout()
    navigate("/")
    setMobileMenuOpen(false)
  }

  const navItems = [
    { href: "#features", label: "Features" },
    { href: "#testimonials", label: "Testimonials" },
    { href: "#contact", label: "Contact" },
  ]

  return (
    <>
      <motion.nav
        className={`fixed w-full z-50 transition-all duration-300 ${
          isScrolled 
            ? "bg-primary/95 backdrop-blur-md shadow-lg border-b border-accent/20" 
            : "bg-primary/90 backdrop-blur-sm shadow-sm"
        }`}
        variants={navVariants}
        initial="hidden"
        animate="visible"
      >
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex justify-between h-16">
            {/* Logo Section */}
            <motion.div
              className="flex items-center"
              variants={logoVariants}
              initial="hidden"
              animate="visible"
            >
              <Link to="/" className="flex items-center gap-3 group">
                <motion.div
                  className="relative flex items-center justify-center w-10 h-10 rounded-xl bg-gradient-to-br from-accent to-accent-secondary text-primary shadow-lg"
                  whileHover={isLoaded ? { 
                    rotate: 360, 
                    scale: 1.1,
                    boxShadow: "0 0 20px rgba(0, 245, 212, 0.5)"
                  } : {}}
                  transition={{ duration: 0.6, ease: "easeInOut" }}
                >
                  <Atom size={22} />
                </motion.div>
                <span className="text-xl md:text-2xl font-bold font-heading bg-gradient-to-r from-accent to-accent-secondary bg-clip-text text-transparent">
                  Drug Discovery AI
                </span>
              </Link>
            </motion.div>

            {/* Desktop Navigation */}
            <div className="hidden md:flex items-center space-x-8">
              {navItems.map((item) => (
                <motion.a
                  key={item.href}
                  href={item.href}
                  className={`relative font-medium font-body transition-all duration-300 group ${
                    activeSection === item.href.slice(1) 
                      ? "bg-gradient-to-r from-accent to-accent-secondary bg-clip-text text-transparent" 
                      : "text-text-primary hover:text-accent"
                  }`}
                  whileHover={isLoaded ? { y: -2 } : {}}
                  transition={{ duration: 0.2 }}
                >
                  {item.label}
                  <motion.span
                    className="absolute bottom-0 left-0 h-0.5 bg-gradient-to-r from-accent to-accent-secondary rounded-full"
                    initial={{ width: 0 }}
                    whileHover={isLoaded ? { width: "100%" } : {}}
                    transition={{ duration: 0.3 }}
                  />
                  {activeSection === item.href.slice(1) && isLoaded && (
                    <motion.span
                      className="absolute bottom-0 left-0 h-0.5 bg-gradient-to-r from-accent to-accent-secondary rounded-full"
                      layoutId="activeSection"
                      initial={{ width: 0 }}
                      animate={{ width: "100%" }}
                      transition={{ duration: 0.3, type: "spring", stiffness: 300 }}
                    />
                  )}
                </motion.a>
              ))}

              {/* User Section */}
              {user ? (
                <div className="flex items-center ml-6 space-x-4">
                  {/* Notifications */}
                  <motion.button
                    className="relative p-2 text-text-secondary hover:text-accent transition-colors duration-200"
                    whileHover={isLoaded ? { scale: 1.1 } : {}}
                    whileTap={{ scale: 0.95 }}
                  >
                    <Bell size={20} />
                    <span className="absolute -top-1 -right-1 w-3 h-3 bg-error rounded-full" />
                  </motion.button>

                  {/* User Profile */}
                  <div className="flex items-center space-x-3">
                    <motion.div
                      className="w-8 h-8 rounded-full bg-gradient-to-br from-accent to-accent-secondary flex items-center justify-center text-primary shadow-md"
                      whileHover={isLoaded ? { 
                        scale: 1.15,
                        boxShadow: "0 0 15px rgba(0, 245, 212, 0.4)"
                      } : {}}
                    >
                      <User size={16} />
                    </motion.div>
                    <div className="hidden lg:block">
                      <div className="text-xs text-text-secondary font-body">Welcome back,</div>
                      <div className="font-semibold text-text-primary font-body">{user.firstName || "User"}</div>
                    </div>
                  </div>

                  {/* Logout Button */}
                  <motion.button
                    className="group relative px-4 py-2 overflow-hidden rounded-lg bg-secondary/80 text-text-primary transition-all duration-300 hover:shadow-md hover:bg-error/20 hover:text-error border border-transparent hover:border-error/30"
                    onClick={handleLogout}
                    variants={buttonVariants}
                    whileHover="hover"
                    whileTap="tap"
                  >
                    <span className="relative font-medium font-body flex items-center">
                      Log Out
                      <LogOut size={16} className="ml-2" />
                    </span>
                  </motion.button>
                </div>
              ) : (
                <div className="flex items-center space-x-4 ml-6">
                  <Link to="/login">
                    <motion.button
                      className="group relative px-6 py-2 overflow-hidden rounded-lg border-2 border-secondary text-text-primary transition-all duration-300 hover:border-accent hover:text-accent"
                      variants={buttonVariants}
                      whileHover="hover"
                      whileTap="tap"
                    >
                      <span className="relative font-medium font-body flex items-center">
                        Login
                        <LogIn size={16} className="ml-2" />
                      </span>
                    </motion.button>
                  </Link>

                  <Link to="/signup">
                    <motion.button
                      className="group relative px-6 py-2 overflow-hidden rounded-lg bg-gradient-to-r from-accent to-accent-secondary text-primary transition-all duration-300 hover:shadow-lg hover:shadow-accent/25"
                      variants={buttonVariants}
                      whileHover={{ 
                        scale: 1.05, 
                        boxShadow: "0 10px 25px -5px rgba(0, 245, 212, 0.4)",
                        y: -2
                      }}
                      whileTap="tap"
                    >
                      <span className="relative font-medium font-body flex items-center">
                        Sign Up
                        <ChevronRight size={16} className="ml-2" />
                      </span>
                    </motion.button>
                  </Link>
                </div>
              )}
            </div>

            {/* Mobile menu button */}
            <div className="md:hidden flex items-center">
              <motion.button
                onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
                className="inline-flex items-center justify-center p-2 rounded-lg text-text-secondary hover:text-accent hover:bg-secondary/50 focus:outline-none transition-all duration-200"
                whileHover={isLoaded ? { scale: 1.1 } : {}}
                whileTap={{ scale: 0.95 }}
              >
                <span className="sr-only">Open main menu</span>
                <AnimatePresence mode="wait">
                  {mobileMenuOpen ? (
                    <motion.div
                      key="close"
                      initial={{ rotate: -90, opacity: 0 }}
                      animate={{ rotate: 0, opacity: 1 }}
                      exit={{ rotate: 90, opacity: 0 }}
                      transition={{ duration: 0.2 }}
                    >
                      <X size={24} />
                    </motion.div>
                  ) : (
                    <motion.div
                      key="menu"
                      initial={{ rotate: 90, opacity: 0 }}
                      animate={{ rotate: 0, opacity: 1 }}
                      exit={{ rotate: -90, opacity: 0 }}
                      transition={{ duration: 0.2 }}
                    >
                      <Menu size={24} />
                    </motion.div>
                  )}
                </AnimatePresence>
              </motion.button>
            </div>
          </div>
        </div>

        {/* Mobile Menu */}
        <AnimatePresence>
          {mobileMenuOpen && (
            <motion.div
              className="md:hidden fixed inset-x-0 top-16 z-50 bg-primary/95 backdrop-blur-md shadow-xl border-b border-accent/20"
              variants={mobileMenuVariants}
              initial="hidden"
              animate="visible"
              exit="hidden"
            >
              <div className="px-4 pt-4 pb-6 space-y-3">
                {/* Navigation Links */}
                {navItems.map((item) => (
                  <motion.a
                    key={item.href}
                    href={item.href}
                    className="block px-4 py-3 rounded-lg text-base font-medium font-body text-text-primary hover:text-accent hover:bg-secondary/50 transition-all duration-200"
                    onClick={() => setMobileMenuOpen(false)}
                    whileHover={isLoaded ? { x: 5, scale: 1.02 } : {}}
                    whileTap={{ scale: 0.98 }}
                  >
                    {item.label}
                  </motion.a>
                ))}

                {/* User Section */}
                {user ? (
                  <div className="pt-4 border-t border-secondary/50">
                    <div className="px-4 py-3 flex items-center space-x-3">
                      <div className="w-10 h-10 rounded-full bg-gradient-to-br from-accent to-accent-secondary flex items-center justify-center text-primary">
                        <User size={18} />
                      </div>
                      <div>
                        <div className="text-xs text-text-secondary font-body">Logged in as</div>
                        <div className="font-semibold text-text-primary font-body">{user.firstName || "User"}</div>
                      </div>
                    </div>

                    <motion.button
                      className="w-full mt-3 flex items-center justify-between px-4 py-3 rounded-lg text-base font-medium font-body text-error bg-error/10 hover:bg-error/20 border border-error/20 hover:border-error/30 transition-colors duration-200"
                      onClick={handleLogout}
                      whileHover={isLoaded ? { scale: 1.02, x: 2 } : {}}
                      whileTap={{ scale: 0.98 }}
                    >
                      <span>Log Out</span>
                      <LogOut size={18} />
                    </motion.button>
                  </div>
                ) : (
                  <div className="space-y-3 pt-4 border-t border-secondary/50">
                    <Link to="/login" onClick={() => setMobileMenuOpen(false)}>
                      <motion.button
                        className="w-full flex items-center justify-between px-4 py-3 rounded-lg text-base font-medium font-body text-text-primary border-2 border-secondary hover:border-accent hover:text-accent hover:bg-secondary/30 transition-all duration-200"
                        whileHover={isLoaded ? { scale: 1.02, x: 2 } : {}}
                        whileTap={{ scale: 0.98 }}
                      >
                        <span>Login</span>
                        <LogIn size={18} />
                      </motion.button>
                    </Link>

                    <Link to="/signup" onClick={() => setMobileMenuOpen(false)}>
                      <motion.button
                        className="w-full flex items-center justify-between px-4 py-3 rounded-lg text-base font-medium font-body text-primary bg-gradient-to-r from-accent to-accent-secondary hover:shadow-lg hover:shadow-accent/25 transition-all duration-200"
                        whileHover={isLoaded ? { 
                          scale: 1.02, 
                          boxShadow: "0 8px 20px -5px rgba(0, 245, 212, 0.4)",
                          x: 2
                        } : {}}
                        whileTap={{ scale: 0.98 }}
                      >
                        <span>Sign Up</span>
                        <UserPlus size={18} />
                      </motion.button>
                    </Link>
                  </div>
                )}
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      </motion.nav>
    </>
  )
}

export default Navbar