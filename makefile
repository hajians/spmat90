# @file Makefile for GDG90
# @author Soheil Hajian

#
# TODO:
#

FC := gfortran

FCFLAGS := -fbounds-check -fstack-check \
	   -fimplicit-none -fdefault-real-8 \
	   -pg #-Wall

#LIB := -llapack -lblas

SRCDIR := src
BUILDDIR := build
MODDIR := mod
TARGET := bin

SRCEXT = f90

SOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SOURCES:.$(SRCEXT)=.o) )
DEP     := $(SRCDIR)/dependencies.dep

-include $(DEP)

test: $(OBJECTS)
	$(FC) $(FCFLAGS) -o $(TARGET)/$@ $(OBJECTS) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@ -J$(MODDIR) $(LIB)

# $(DEP): $(SRCDIR)/fort_depend.py
# 	@echo "Making dependencies..."
# 	python $(SRCDIR)/fort_depend.py -b $(BUILDDIR) -w -o $(DEP) -f $(SRCDIR)/*.$(SRCEXT)

clean:
	@echo "Cleaning..."
	rm -f $(BUILDDIR)/*.o $(MODDIR)/*.mod $(TARGET)/*


