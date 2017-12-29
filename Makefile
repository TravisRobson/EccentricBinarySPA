CC = gcc-7

LIBS  = gsl gslcblas m 
CCFLAGS += -g -Wall -std=gnu99 -fmax-errors=5 #-Werror 
#CCFLAGS += -g -ffast-math -Wall -O2 -ftree-vectorize -std=gnu99 -fmax-errors=5 #-Werror 

OBJS = Ecc_SPA.o Ecc_Binary.o Ecc_Adiabat_Evol.o Detector.o Ecc_IO.o Ecc_Math.o

all : $(OBJS) mcmc # num_test no_RR_test harm_match spa_harm full_match num_spa_harms fisher 

Ecc_SPA.o : Ecc_SPA.c Ecc_SPA.h Ecc_Binary.h Constants.h
	$(CC) $(CCFLAGS) -c Ecc_SPA.c

Ecc_Binary.o : Ecc_Binary.c Ecc_Binary.h Constants.h
	$(CC) $(CCFLAGS) -c Ecc_Binary.c
	
Ecc_Adiabat_Evol.o : Ecc_Adiabat_Evol.c Ecc_Adiabat_Evol.h Ecc_Binary.h Constants.h
	$(CC) $(CCFLAGS) -c Ecc_Adiabat_Evol.c

Detector.o : Detector.c Detector.h Constants.h
	$(CC) $(CCFLAGS) -c Detector.c
	
Ecc_IO.o : Ecc_IO.c Ecc_IO.h 
	$(CC) $(CCFLAGS) -c Ecc_IO.c
	
Ecc_Math.o : Ecc_Math.c Ecc_Math.h Ecc_IO.h Ecc_Binary.h Detector.h
	$(CC) $(CCFLAGS) -c Ecc_Math.c



# num_test : $(OBJS) Ecc_num_test.c
# 	$(CC) $(CCFLAGS) -o num_test Ecc_num_test.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
# 	
# no_RR_test : $(OBJS) ecc_num_no_RR.c
# 	$(CC) $(CCFLAGS) -o no_RR_test ecc_num_no_RR.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)	
# 
# harm_match : $(OBJS) harmonic_breakout_matches.c
# 	$(CC) $(CCFLAGS) -o harm_match harmonic_breakout_matches.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
# 	
# spa_harm : $(OBJS) harmonic_breakout_num_SPA.c
# 	$(CC) $(CCFLAGS) -o spa_harm harmonic_breakout_num_SPA.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)		
# 
# full_match : $(OBJS) full_ecc_match.c
# 	$(CC) $(CCFLAGS) -o full_match full_ecc_match.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
# 	
# num_spa_harms : $(OBJS) num_SPA_harms.c
# 	$(CC) $(CCFLAGS) -o num_spa_harms num_SPA_harms.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
# 	
# fisher : $(OBJS) fisher_ecc.c
# 	$(CC) $(CCFLAGS) -o fisher fisher_ecc.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
	
mcmc : $(OBJS) mcmc_ecc.c
	$(CC) $(CCFLAGS) -o mcmc mcmc_ecc.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
	
clean: 
	rm *.o mcmc #num_test no_RR_test harm_match spa_harm full_match num_spa_harms fisher 