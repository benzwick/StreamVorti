# StreamVorti Makefile
# Alternative to CMake build system, compatible with Spack-installed MFEM

# Project configuration
PROJECT_NAME = StreamVorti
VERSION = 1.0.0

# Compiler settings
CXX = g++
CXXSTD = -std=c++17

# Build directories
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
LIB_DIR = $(BUILD_DIR)/lib
BIN_DIR = $(BUILD_DIR)/bin
SRC_DIR = src
INCLUDE_DIR = include

# Spack environment detection
# Users should run: spack load mfem
# This will set up the environment variables
SPACK_ENV := $(shell command -v spack 2> /dev/null)
ifdef SPACK_ENV
    # Try to get MFEM config from spack
    MFEM_PREFIX := $(shell spack location -i mfem 2>/dev/null || echo "")
else
    # Fallback to common Spack installation paths
    MFEM_PREFIX := $(shell find /opt/spack -name "mfem*" -type d 2>/dev/null | head -1)
    ifeq ($(MFEM_PREFIX),)
        MFEM_PREFIX := $(shell find $(HOME)/spack -name "mfem*" -type d 2>/dev/null | head -1)
    endif
endif

# If MFEM_PREFIX is still empty, use environment variable or error
ifeq ($(MFEM_PREFIX),)
    ifdef MFEM_DIR
        MFEM_PREFIX = $(MFEM_DIR)
    else
        $(error MFEM not found. Please run 'spack load mfem' or set MFEM_DIR)
    endif
endif

# MFEM configuration
MFEM_INC = -I$(MFEM_PREFIX)/include
MFEM_LIB = -L$(MFEM_PREFIX)/lib -lmfem

# Try to get MFEM config automatically
MFEM_CONFIG := $(shell find $(MFEM_PREFIX) -name "mfem-config" 2>/dev/null | head -1)
ifneq ($(MFEM_CONFIG),)
    MFEM_CXX_FLAGS := $(shell $(MFEM_CONFIG) --cppflags)
    MFEM_LD_FLAGS := $(shell $(MFEM_CONFIG) --libs)
else
    MFEM_CXX_FLAGS = $(MFEM_INC)
    MFEM_LD_FLAGS = $(MFEM_LIB)
endif

# CGAL configuration (assuming system installation)
CGAL_INC = -I/usr/include
CGAL_LIB = -lCGAL -lCGAL_Core -lgmp -lmpfr

# Eigen3 configuration (assuming system installation)
EIGEN_INC = -I/usr/include/eigen3

# StreamVorti include paths
STREAMVORTI_INC = -I$(INCLUDE_DIR)

# Compiler flags
CXXFLAGS = $(CXXSTD) -O3 -march=native -msse2 -ffast-math -fPIC -Wall
CXXFLAGS += -ansi -frounding-math -fsignaling-nans -mfpmath=sse

# Include paths
INCLUDES = $(STREAMVORTI_INC) $(MFEM_CXX_FLAGS) $(CGAL_INC) $(EIGEN_INC)

# Library flags
LDFLAGS = $(MFEM_LD_FLAGS) $(CGAL_LIB) -lm

# Source files
APPROXIMANTS_SRCS = $(wildcard $(SRC_DIR)/approximants/*.cpp)
SUPPORT_DOMAIN_SRCS = $(wildcard $(SRC_DIR)/support_domain/*.cpp)
LIBRARY_SRCS = $(APPROXIMANTS_SRCS) $(SUPPORT_DOMAIN_SRCS)

# Object files
APPROXIMANTS_OBJS = $(APPROXIMANTS_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
SUPPORT_DOMAIN_OBJS = $(SUPPORT_DOMAIN_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
LIBRARY_OBJS = $(APPROXIMANTS_OBJS) $(SUPPORT_DOMAIN_OBJS)

# Filter out unwanted source files
LIBRARY_OBJS := $(filter-out $(OBJ_DIR)/support_domain/tocompare%,$(LIBRARY_OBJS))

# Target library
STATIC_LIB = $(LIB_DIR)/lib$(PROJECT_NAME)_static.a

# Main executable
MAIN_EXEC = $(BIN_DIR)/MfemRun
MAIN_SRC = mfem_main.cpp

# Default target
.PHONY: all clean install help

all: $(STATIC_LIB) $(MAIN_EXEC)

# Create directories
$(OBJ_DIR) $(LIB_DIR) $(BIN_DIR):
	@mkdir -p $@

$(OBJ_DIR)/approximants $(OBJ_DIR)/support_domain: | $(OBJ_DIR)
	@mkdir -p $@

# Compile approximants objects
$(OBJ_DIR)/approximants/%.o: $(SRC_DIR)/approximants/%.cpp | $(OBJ_DIR)/approximants
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile support_domain objects
$(OBJ_DIR)/support_domain/%.o: $(SRC_DIR)/support_domain/%.cpp | $(OBJ_DIR)/support_domain
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create static library
$(STATIC_LIB): $(LIBRARY_OBJS) | $(LIB_DIR)
	@echo "Creating static library $@..."
	ar rcs $@ $^

# Compile main executable
$(MAIN_EXEC): $(MAIN_SRC) $(STATIC_LIB) | $(BIN_DIR)
	@echo "Compiling main executable $@..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -l$(PROJECT_NAME)_static $(LDFLAGS) -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR)

# Install (optional - copies to system directories)
install: all
	@echo "Installing StreamVorti..."
	@echo "Note: This is a basic install. Adjust paths as needed."
	sudo cp $(STATIC_LIB) /usr/local/lib/
	sudo cp -r $(INCLUDE_DIR)/StreamVorti /usr/local/include/
	sudo cp $(MAIN_EXEC) /usr/local/bin/

# Print configuration info
info:
	@echo "StreamVorti Build Configuration:"
	@echo "  MFEM_PREFIX: $(MFEM_PREFIX)"
	@echo "  MFEM_CONFIG: $(MFEM_CONFIG)"
	@echo "  CXX: $(CXX)"
	@echo "  CXXFLAGS: $(CXXFLAGS)"
	@echo "  INCLUDES: $(INCLUDES)"
	@echo "  LDFLAGS: $(LDFLAGS)"
	@echo "  Source files: $(words $(LIBRARY_SRCS)) files"

# Help target
help:
	@echo "StreamVorti Makefile Usage:"
	@echo ""
	@echo "Prerequisites:"
	@echo "  - Run 'spack load mfem' before building"
	@echo "  - Or set MFEM_DIR environment variable"
	@echo ""
	@echo "Targets:"
	@echo "  all      - Build static library and main executable (default)"
	@echo "  clean    - Remove all build artifacts"
	@echo "  install  - Install library and headers (requires sudo)"
	@echo "  info     - Display build configuration"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Usage examples:"
	@echo "  make                    # Build everything"
	@echo "  make clean all          # Clean and rebuild"
	@echo "  MFEM_DIR=/path make     # Build with specific MFEM path"

# Declare object file dependencies
$(LIBRARY_OBJS): $(wildcard $(INCLUDE_DIR)/StreamVorti/*.hpp) $(wildcard $(INCLUDE_DIR)/StreamVorti/*/*.hpp)