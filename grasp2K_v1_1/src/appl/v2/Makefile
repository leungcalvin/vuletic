SUBDIR =  bioscl2  biotra2  cndens2  rhfs2   \
          mcp2  rci2  rscf2   sms2  jj2lsj \
          mcp2_mpi rci2_mpi rscf2_mpi bioscl2_mpi biotra2_mpi

TARGETS	 = install

$(TARGETS):
	@for i in $(SUBDIR) ;\
	do \
		echo "Current directory: $$i" ; \
	                cd $$i ; \
			make; \
echo "....Leaving: $$i" ; echo; echo;\
		cd .. ; \
	done

clean :
	@for i in $(SUBDIR) ;\
	do \
	echo "Current directory: $$i" ; \
           cd $$i ; \
           make clean; \
           cd .. ; \
        done

