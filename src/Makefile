# Makefile for ggsspnma Fortran program
# Compile on macOS

# Compiler
FC = gfortran

# Compiler flags - fixed-form with extended line length for long lines
FFLAGS = -O2 -fallow-argument-mismatch -ffixed-form -ffixed-line-length-132 -fno-range-check

# ggsspnma parallelizes the per-surface-point DPA sweep with OpenMP
GGSSPNMA_FFLAGS = $(FFLAGS) -fopenmp

# Output directory for compiled binaries
BINDIR = ../bin

# Default target
all: ggsspnma extremefit

# Compile ggsspnma
ggsspnma: ggsspnma.f
	$(FC) $(GGSSPNMA_FFLAGS) -o $(BINDIR)/ggsspnma ggsspnma.f
	@echo "✓ Compiled ggsspnma successfully"

extremefit: extremefit.f
	$(FC) $(FFLAGS) -o $(BINDIR)/extremefit dnls.f extremefit.f 
	@echo "✓ Compiled extremefit successfully"

# Clean up
clean:
	rm -f $(BINDIR)/*
	@echo "✓ Cleaned build artifacts"

# Test compilation
test: ggsspnma extremefit
	@echo "✓ Test compilation successful"
	@echo "Executable: ../bin/ggsspnma"

.PHONY: all clean test
