EXE = mcp1
#BINDIR = ../../../bin
#GRASPLIB = ../../../lib
#GRASPBLOCKD = ../../../lib/graspBD.o
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)


APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
       fndbeg.o genmcp.o getinf.o lodres.o mcp.o setdbg.o setmcp.o \
       setsda.o setsum.o sort.o strsum.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

