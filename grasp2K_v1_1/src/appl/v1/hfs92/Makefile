EXE = rhfs1
#BINDIR = ../../../bin
#GRASPLIB = ../../../lib
#GRASPBLOCKD = ../../../lib/graspBD.o
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)


APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
        chkplt.o gethfd.o hfs.o hfs92.o matelt.o  rinthf.o setdbg.o \
        setsum.o strsum.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS)\
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

