prjdir := ../..
objdir := $(prjdir)/obj

include $(prjdir)/Makefile.in

main := main.F90
user := user_input.F90
exe  := $(main:%.F90=%.x)

LD := $(FC)

ifeq ($(LD),gfortran)
   OPENMP := -fopenmp
   FFLAGS := -c -O3 -Wall -Wno-unused-dummy-argument -finit-local-zero $(OPENMP)\
             -I$(objdir)
   LFLAGS := -Wall -Wno-unused-dummy-argument $(OPENMP) -I$(objdir)
else ifeq ($(LD),ifort)
   OPENMP := -qopenmp
   FFLAGS := -c -O3 -warn errors -zero $(OPENMP) -I$(objdir)
   LFLAGS := -warn errors $(OPENMP) -I$(objdir)
else
   $(error Unrecognized compiler $(FC): edit Makefile for your compiler and flags)
endif

LIBS := -L$(prjdir) -luber -L$(libdir) -ldcmt -llapack -lrefblas

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ $<

$(exe): $(main:%.F90=%.o) $(user:%.F90=%.o)
	$(LD) $(LFLAGS) -o $@ $^ $(LIBS)

oclean:
	-rm -f *.o

clean: oclean
	-rm -f *.smod *.x

