#!/bin/bash
# Build script for DPA Fortran programs

set -e  # Exit on error

echo "=========================================="
echo "DPA Fortran Build Script"
echo "=========================================="
echo ""

# >>> mamba initialize >>>
# !! Contents within this block are managed by 'micromamba shell init' !!
export MAMBA_EXE='/Users/mewall/.local/bin/micromamba';
export MAMBA_ROOT_PREFIX='/Users/mewall/packages/micromamba';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias micromamba="$MAMBA_EXE"  # Fallback on help from micromamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<


# Check if environment exists
if micromamba env list | grep -q "dpa-fortran "; then
    echo "✓ Conda environment 'dpa-fortran' already exists"
    echo "  To recreate it, run: conda env remove -n dpa-fortran"
else
    echo "Creating conda environment 'dpa-fortran'..."
    micromamba create -f environment.yml
    echo "✓ Environment created successfully"
fi

echo ""
echo "Activating conda environment..."
# Note: This needs to be sourced in the parent shell
# For now, we'll provide instructions

echo ""
echo "=========================================="
echo "To compile the programs:"
echo "=========================================="
echo ""
echo "1. Activate the conda environment:"
echo "   conda activate dpa-fortran"
echo ""
echo "2. Compile all programs:"
echo "   make all"
echo ""
echo "Compiled binaries will be in: ./bin/"
echo "=========================================="

micromamba activate dpa-fortran
cd src
make
