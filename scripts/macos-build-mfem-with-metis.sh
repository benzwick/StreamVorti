#!/bin/bash
#
# Build MFEM 4.8 on macOS with MPI, METIS, HYPRE, and SuiteSparse support
#
# For building MFEM on Will's MacBook. Uses Homebrew to install and
# locate dependencies (metis, hypre, suite-sparse, open-mpi).
#
# This script will:
# 1. Download MFEM 4.8 source (if needed)
# 2. Install dependencies via Homebrew (if needed)
# 3. Configure with all parallel features enabled
# 4. Build and install to ~/local/mfem-4.8
#

set -e  # Exit on error

# Configuration
MFEM_VERSION="4.8"
MFEM_DIR="${HOME}/local/mfem-${MFEM_VERSION}"
BUILD_DIR="${HOME}/mfem-build-${MFEM_VERSION}"
SOURCE_DIR="${BUILD_DIR}/mfem-${MFEM_VERSION}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}MFEM ${MFEM_VERSION} Rebuild Script${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Check if Homebrew dependencies are installed
echo -e "${YELLOW}Checking dependencies...${NC}"
for pkg in metis hypre suite-sparse open-mpi; do
    if brew list --versions $pkg > /dev/null 2>&1; then
        echo -e "  ${GREEN}✓${NC} $pkg: $(brew list --versions $pkg)"
    else
        echo -e "  ${RED}✗${NC} $pkg not installed"
        echo -e "${YELLOW}Installing $pkg via Homebrew...${NC}"
        brew install $pkg
    fi
done
echo ""

# Create build directory
echo -e "${YELLOW}Creating build directory: ${BUILD_DIR}${NC}"
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Download MFEM source if not present
if [ ! -d "${SOURCE_DIR}" ]; then
    echo -e "${YELLOW}Downloading MFEM ${MFEM_VERSION} source...${NC}"
    curl -L "https://github.com/mfem/mfem/archive/refs/tags/v${MFEM_VERSION}.tar.gz" -o "mfem-${MFEM_VERSION}.tar.gz"
    tar -xzf "mfem-${MFEM_VERSION}.tar.gz"
    echo -e "${GREEN}✓ Source downloaded and extracted${NC}"
else
    echo -e "${GREEN}✓ Source already present${NC}"
fi

# Create build subdirectory
mkdir -p "${SOURCE_DIR}/build"
cd "${SOURCE_DIR}/build"

# Find Homebrew package paths
METIS_DIR=$(brew --prefix metis)
HYPRE_DIR=$(brew --prefix hypre)
SUITESPARSE_DIR=$(brew --prefix suite-sparse)
MPI_DIR=$(brew --prefix open-mpi)

echo ""
echo -e "${YELLOW}Dependency paths:${NC}"
echo "  METIS:       ${METIS_DIR}"
echo "  HYPRE:       ${HYPRE_DIR}"
echo "  SuiteSparse: ${SUITESPARSE_DIR}"
echo "  MPI:         ${MPI_DIR}"
echo ""

# Backup existing installation
if [ -d "${MFEM_DIR}" ]; then
    BACKUP_DIR="${MFEM_DIR}.backup.$(date +%Y%m%d_%H%M%S)"
    echo -e "${YELLOW}Backing up existing installation...${NC}"
    echo "  ${MFEM_DIR} -> ${BACKUP_DIR}"
    mv "${MFEM_DIR}" "${BACKUP_DIR}"
fi

# Ensure install directory exists
mkdir -p "$(dirname "${MFEM_DIR}")"

# Configure MFEM with CMake
echo ""
echo -e "${YELLOW}Configuring MFEM with CMake...${NC}"
echo -e "${YELLOW}(This may take a few minutes)${NC}"
echo ""

cmake .. \
    -DCMAKE_INSTALL_PREFIX="${MFEM_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DMFEM_USE_MPI=YES \
    -DMFEM_USE_METIS=YES \
    -DMFEM_USE_METIS_5=YES \
    -DMETIS_DIR="${METIS_DIR}" \
    -DMFEM_USE_HYPRE=YES \
    -DHYPRE_DIR="${HYPRE_DIR}" \
    -DMFEM_USE_SUITESPARSE=YES \
    -DSuiteSparse_DIR="${SUITESPARSE_DIR}" \
    -DMFEM_USE_LAPACK=YES \
    -DCMAKE_CXX_COMPILER="${MPI_DIR}/bin/mpicxx" \
    -DCMAKE_C_COMPILER="${MPI_DIR}/bin/mpicc"

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ CMake configuration failed!${NC}"
    exit 1
fi

echo -e "${GREEN}✓ CMake configuration successful${NC}"
echo ""

# Build MFEM
echo -e "${YELLOW}Building MFEM (using all available cores)...${NC}"
echo -e "${YELLOW}(This will take several minutes)${NC}"
echo ""

make -j$(sysctl -n hw.ncpu)

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Build failed!${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Build successful${NC}"
echo ""

# Install MFEM
echo -e "${YELLOW}Installing to ${MFEM_DIR}...${NC}"
make install

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Installation failed!${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Installation successful${NC}"
echo ""

# Verify installation
echo -e "${YELLOW}Verifying installation...${NC}"
if [ -f "${MFEM_DIR}/share/mfem/config.mk" ]; then
    echo ""
    echo -e "${GREEN}Configuration summary:${NC}"
    grep -E "MFEM_USE_(MPI|METIS|HYPRE|SUITESPARSE)" "${MFEM_DIR}/share/mfem/config.mk"
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}✓ MFEM ${MFEM_VERSION} successfully rebuilt!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Next steps:${NC}"
    echo "  1. Rebuild StreamVorti:"
    echo "     cd /Users/will-li/Coding/Research_project/dev/StreamVorti-1/build"
    echo "     cmake -DMFEM_DIR=${MFEM_DIR} -DCMAKE_BUILD_TYPE=Debug .."
    echo "     make clean && make -j6"
    echo ""
    echo "  2. Test parallel execution:"
    echo "     mpirun -np 4 ./StreamVorti_par -nx 20 -ny 20 -Re 100 -solver cg"
    echo ""
else
    echo -e "${RED}✗ Installation verification failed!${NC}"
    exit 1
fi
