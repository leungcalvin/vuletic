EXE = rci1
#BINDIR = ../../../bin
#GRASPLIB = ../../../lib
#GRASPBLOCKD = $(GRASPLIB)/graspBD.o
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= bessel.o breid.o breit.o brint1.o brint2.o brint3.o brint4.o \
         brint5.o brint6.o brintf.o \
         brra.o chkplt.o cxk.o dalcmc.o dnicmv.o \
         evcout.o fndblk.o funk.o funl.o fzalf.o \
         genintb.o genintrk.o getcid.o hmout.o hovlap.o iabint.o indtpi.o \
         iniestsd.o iniestdm.o keint.o klamaq.o lodmix.o lodres.o maneig.o \
         matrix.o mohr.o ncharg.o qed.o rci92.o \
         rkint.o rkintc.o setdbg.o setham.o setmix.o \
         setres.o setsum.o skint.o snrc.o spicmv.o \
         spodmv.o strsum.o talk.o triangrk.o vac2.o vac4.o \
         vacpol.o vint.o vinti.o vpint.o vpintf.o \
         zkf.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

