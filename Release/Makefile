EXEC   = ../../binarios/N-GenIC

OBJS   =  ../main.o ../power.o ../allvars.o ../save.o ../read_param.o  ../read_glass.o ../rw_snapshot.o

INCL   = ../allvars.h ../proto.h    Makefile



#OPT	+=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing

OPT	+=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
OPT	+=  -DCORRECT_CIC  # only switch this on if particles are homogenously distributed over mesh cells (say glass)


OPTIONS =  $(OPT)



SYSTYPE="OPA-Cluster64-Gnu"


CC       =   mpicc       # sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)
MPICHLIB =  -lmpich





OPTIMIZE =  -O3 -g -Wall -m64
GSL_INCL =  -I/usr/include/gsl -I/usr/include/mpich2/ -I/usr/include/libxml2
GSL_LIBS =   -Xlinker -R -Xlinker 
FFTW_INCL= 
FFTW_LIBS=  
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  
OPT      +=  -DNOCALLSOFSYSTEM


FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw -lxml2

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB) $(GSL_LIBS)  -L/usr/lib -lgsl -lgslcblas


CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

all: $(OBJS) $(EXEC)  
$(EXEC) : $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -rf `find .. -name '*.o'` $(EXEC)



