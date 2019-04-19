EXE = oscl
#BINDIR = ../../../bin
#GRASPLIB = ../../../lib
#GRASPBLOCKD = ../../../lib/graspBD.o
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
      alclla.o alcnma.o alcnsa.o alcnta.o bessj.o chkplt.o csfm.o \
      getosd.o getrmp.o mctin.o oscl.o   oscl92.o printa.o setdbg.o setsum.o \
      spme.o   strsum.o trsort.o 

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

