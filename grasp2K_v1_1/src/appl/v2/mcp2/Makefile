EXE = mcp2
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

#FFLAGS = -O3 

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
          genmcp.o  mcp.o outsda.o setmcp2.o \
          fndbeg.o getinf.o setdbg.o setmcp.o setsda.o setsum.o settmp.o \
          sort.o strsum.o
$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

