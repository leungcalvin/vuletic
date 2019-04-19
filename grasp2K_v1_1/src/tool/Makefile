.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
BIN = $(GRASP)/bin
LIB = $(GRASPLIBS)
BD = $(LIBDIR)/graspBD.o
UTIL = cvtmix  extmix  jsplit mchfmcdf plotmcdf  readrwf \
      rlevels rotate rotate_pair partition

install: EXE 
	cp rsave $(GRASP)/bin

EXE : $(BIN)/cvtmix $(BIN)/extmix $(BIN)/jsplit  \
      $(BIN)/mchfmcdf $(BIN)/plotmcdf $(BIN)/readrwf \
      $(BIN)/rlevels $(BIN)/rotate \
	$(BIN)/rotate_pair $(BIN)/partition
   

OBJ : cvtmix.o  extmix.o jsplit.o  mchfmcdf.o plotmcdf.o  readrwf.o\
      rlevels.o rotate.o format_mix.o rlevelsj.o \
	rotate_pair.o partition.o

$(BIN)/cvtmix  : cvtmix.o
	$(FC) -o $(BIN)/cvtmix  $(FC_LD) -L$(LIBDIR) \
                 cvtmix.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/format_mix  : format_mix.o
	$(FC) -o $(BIN)/format_mix  $(FC_LD) -L$(LIBDIR) \
                 format_mix.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/extmix: extmix.o
	$(FC) -o $(BIN)/extmix $(FC_LD) -L$(LIBDIR)  extmix.o \
                 $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/jsplit: jsplit.o
	$(FC) -o $(BIN)/jsplit $(FC_LD)  -L$(LIBDIR) jsplit.o \
                 $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/mchfmcdf: mchfmcdf.o
	$(FC) -o $(BIN)/mchfmcdf $(FC_LD) -L$(LIBDIR) mchfmcdf.o \
                 $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/plotmcdf: plotmcdf.o 
	$(FC) -o $(BIN)/plotmcdf $(FC_LD) -L$(LIBDIR) \
                 plotmcdf.o $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/readrwf: readrwf.o
	$(FC) -o $(BIN)/readrwf $(FC_LD) -L$(LIBDIR) readrwf.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rlevels: rlevels.o
	$(FC) -o $(BIN)/rlevels $(FC_LD) -L$(LIBDIR) rlevels.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rlevelsj: rlevelsj.o
	$(FC) -o $(BIN)/rlevelsj $(FC_LD) -L$(LIBDIR) rlevelsj.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rotate: rotate.o
	$(FC) -o $(BIN)/rotate $(FC_LD) -L$(LIBDIR) rotate.o $(BD) $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rotate_pair: rotate_pair.o
	$(FC) -o $(BIN)/rotate_pair  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                rotate_pair.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/partition: partition.o
	$(FC) -o $(BIN)/partition partition.o 

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

.f90.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

