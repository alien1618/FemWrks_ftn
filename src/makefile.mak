F08 = gfortran
FFLAGS = -O3 -g -Wall -Warray-bounds -ffixed-line-length-none -fbounds-check -fopenmp -O3 -ffast-math -funroll-loops -march=native -mtune=native -J mod
SRCDIR = ./src/femwrks
SLVRDIR = ./src/slvrs
tstDIR = ./src
MODDIR = ./mod
OBJDIR = obj
BINDIR :=  bin

#####################################################################
OBJS =\
          msh.o quad.o eqslvrs.o  krnl.o bc.o prmtrs.o gm.o matfree.o slv_trnsprt.o slv_trnsprt_ebe.o slv_ns.o slv_ns_ebe.o slv_elst.o slv_elst_ebe.o  main.o

##################################################################### 

RUN: $(OBJS)
	$(F08) $(FFLAGS) $(OBJS) -o run

clean:
	rm *.mod *.o *~ run

##################################################################### 
msh.o        	       : $(SRCDIR)/msh.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/msh.f08

quad.o        	       : $(SRCDIR)/quad.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/quad.f08

eqslvrs.o        	       : $(SRCDIR)/eqslvrs.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/eqslvrs.f08
	
krnl.o        	       : $(SRCDIR)/krnl.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl.f08
	
gm.o        	       : $(SRCDIR)/gm.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/gm.f08

matfree.o        	       : $(SRCDIR)/matfree.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/matfree.f08

bc.o        	       : $(SRCDIR)/bc.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/bc.f08

prmtrs.o        	       : $(SRCDIR)/prmtrs.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/prmtrs.f08
##################################################################### 
slv_trnsprt.o        	       : $(SLVRDIR)/slv_trnsprt.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_trnsprt.f08

slv_trnsprt_ebe.o        	       : $(SLVRDIR)/slv_trnsprt_ebe.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_trnsprt_ebe.f08

slv_elst.o        	       : $(SLVRDIR)/slv_elst.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_elst.f08

slv_elst_ebe.o        	       : $(SLVRDIR)/slv_elst_ebe.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_elst_ebe.f08

slv_ns.o        	       : $(SLVRDIR)/slv_ns.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ns.f08

slv_ns_ebe.o        	       : $(SLVRDIR)/slv_ns_ebe.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ns_ebe.f08
##################################################################### 	
main.o        	       : $(tstDIR)/main.f08
	$(F08) $(FFLAGS)  -c $(tstDIR)/main.f08
#####################################################################
