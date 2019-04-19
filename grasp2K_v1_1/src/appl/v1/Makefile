SUBDIR = bioscl biotra cndens erwf gencsl geniso  \
         genmcp  hfs92  jjgen jjreduce mrgcsl  oscl92  rci92  \
         rscf92  sms92  xcsl bioscl_mpi biotra_mpi
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

