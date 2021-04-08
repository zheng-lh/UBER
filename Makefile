# UBER library makefile
# Liheng Zheng <zhengliheng@gmail.com>
# William B. Hanson Center for Space Sciences
# The University of Texas at Dallas
# Copyright (C) 2021

include Makefile.in

###################
#   Directories   #
###################

# define directories relative to the project directory (where this Makefile is)
prjdir := .
srcdir := $(prjdir)/src
objdir := $(prjdir)/obj

# all directories that make will search for source files
srcsubdir := params utilities grid domain equation solver interface
VPATH := $(foreach dir, $(srcsubdir), $(addprefix $(srcdir)/, $(dir)))

####################
#   Source files   #
####################

# N.B., the order of these source files defines the dependency relations
src_params    := mod_typedef_params.F90
src_utilities := fornberg.F interpolation.F90 statistics.F90 linear_algebra.F90\
                 functions.F90
src_grid      := mod_grid_typedef.F90 mod_grid_scalar_field.F90\
                 mod_grid_vector_field.F90 mod_grid_tensor_field.F90\
                 mod_grid_io.F90 mod_grid.F90
src_domain    := mod_domain_typedef.F90 mod_domain.F90
src_equation  := mod_equation_typedef.F90 mod_equation.F90
src_solver    := solver_dcmt_omp.c mod_solver_random_number.F90 mod_solver_tree.F90\
                 mod_solver_solution_grid.F90 mod_solver_ito_process.F90\
                 mod_solver_batch.F90 mod_solver_engine.F90 mod_solver_io.F90\
                 mod_solver.F90
src_interface := uber.F90

sources := $(src_params) $(src_utilities) $(src_grid) $(src_domain) $(src_equation)\
           $(src_solver) $(src_interface)

###########################
#   Archive and objects   #
###########################

uber    := uber
libfnm  := lib$(uber).a
objects := $(addprefix $(objdir)/, $(addsuffix .o, $(basename $(sources))))

######################
#   Compiler flags   #
######################

ifeq ($(FC),gfortran)
   OPENMP := -fopenmp
   FFLAGS := -c -O3 -Wall -Wno-unused-dummy-argument -Wno-integer-division\
             -finit-local-zero -J$(objdir) $(OPENMP)
   CFLAGS := -c -O3 -Wall -Wmissing-prototypes -std=c99 $(OPENMP)
else ifeq ($(FC),ifort)
   OPENMP := -qopenmp
   FFLAGS := -c -O3 -warn errors -zero -module $(objdir) $(OPENMP)
   CFLAGS := -c -O3 -Wall -Wmissing-prototypes $(OPENMP)
else
   $(error Unrecognized compiler $(FC): edit Makefile for your compiler and flags)
endif

# in command line, 'make mode=debug' compiles the library using the debug flags
mode =
ifeq ($(mode),debug)
   ifeq ($(FC),gfortran)
      FFLAGS += -Og -fcheck=all -fbacktrace
   else
      FFLAGS += -g -check bounds -traceback -debug all
   endif
endif

#########################
#   Targets and rules   #
#########################

.DELETE_ON_ERROR:
.PHONY: lib oclean clean arc

log := make.log

lib: $(log) $(objdir) $(prjdir)/$(libfnm)

$(objdir):
	@test ! -d $(objdir) && mkdir $(objdir)

$(objdir)/%.o: %.F90
	$(FC) $(FFLAGS) -o $@ $<

$(objdir)/%.o: %.F
	$(FC) $(FFLAGS) -o $@ $<

$(objdir)/%.o: %.c $(incdir)/dc.h
	$(CC) $(CFLAGS) -I$(incdir) -o $@ $<

$(prjdir)/$(libfnm): $(objects)
	ar -crs $@ $^

$(log): FORCE
	@hostname >$@
	@echo -e "\nFC =">>$@
	@$(FC) -v >>$@ 2>&1
	@echo -e "\nCC =">>$@
	@$(CC) -v >>$@ 2>&1
	@echo -e "\nFFLAGS = $(FFLAGS)\n\nCFLAGS = $(CFLAGS)">>$@

oclean:
	-rm -f $(objdir)/*.o

clean:
	-rm -rf $(log) $(objdir) $(prjdir)/$(libfnm)

FORCE: ;


