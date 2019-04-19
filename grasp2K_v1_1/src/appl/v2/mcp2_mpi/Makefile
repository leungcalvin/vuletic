EXE = mcp2_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

 
APP_LIBS = ${GRASPBLOCKD}  -L ${GRASPLIB}  -lmpiu ${GRASPLIBS} 

APP_OBJ= \
         genmcpmpi_${FC_MALLOC}.o mcpmpi.o outsdampi.o setdbgmpi.o setmcpmpi.o \
         fndbeg.o getinf.o setdbg.o setmcp.o setsda.o setsum.o settmp.o \
         sort.o strsum.o
 
$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
	-L ${LAPACK_DIR} ${LAPACK_LIB} 

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

