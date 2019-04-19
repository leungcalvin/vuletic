EXE = rscf1
BINDIR = $(GRASP)/bin
GRASPLIB = $(GRASP)/lib
GRASPBLOCKD = $(GRASPLIB)/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS} -ldvdson

OBJECTS= alcsca.o chkplt.o consis.o csfwgt.o dacon.o \
         dampck.o dampor.o defcor.o dsubrs.o eigen.o \
         endsum.o estim.o fco.o fndblk.o gco.o \
         getald.o getold.o getscd.o hmout.o improv.o \
         in.o lagcon.o maneig.o matrix.o maxarr.o \
         newco.o newe.o orbout.o out.o \
         outbnd.o prtrsl.o prwf.o rscf92.o scf.o \
         setcof.o setdbg.o setham.o setlag.o setmcp.o \
         setmix.o setrwg.o setsum.o setxuv.o setxv.o \
         setxz.o solve.o spicmv.o strsum.o xpot.o \
         ypot.o orthy.o \
         dsuvec.o fcovec.o gcovec.o

$(BINFILE) : $(OBJECTS)
	$(FC) -o $(BINFILE) $(FC_LD)  $(OBJECTS) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

