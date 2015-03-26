
PROGS	=	cnvx_min

HEADS	=	functions.hpp	
CCF	=	

OPTS	=	-Wall -O3  # -m64 -g -fno-inline 
LIBS	=	-lgsl -lgslcblas #-lgomp
LIBSGL	=	#-lGL -lGLU -lSDL

CC	=	g++


cnvx_min: cnvxification_mc.cpp $(HEADS)
	$(CC) $(OPTS) -o $@ cnvxification_mc.cpp $(CCF) $(LIBS)


clean: 
	rm $(PROGS) 
