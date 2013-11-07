objects = test_poisson.o lu_solver.o

sources = test_poisson.cpp lu_solver.cpp

    #LINK TO THE COMPILED SUPERLU_4.3 LIBRARY (EARLIER VERSIONS DONT WORK) 
    #MUST BE DOWNLOADED FROM http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
#MYLIBS = /hera/bhs/lib/libsuperlu_4.3.a
    #IF BLAS LIBRARY IS NOT INSTALLED BY DEFAULT
MYLIBS = /hera/bhs/lib/libsuperlu_4.3.a /hera/bhs/lib/blas_LINUX.a
#OTHER LINKS
LIBS =  -lstdc++ -lm -g #-lblas   #UNCOMMENT IF BLAS LIBRARY IS INSTALLED

#LINK TO THE SUPERLU INCLUDE FILES (IN THE INSTALLATION. ZIP). 
INCLS = -I/hera/bhs/lib/include

CC           = g++

LOADER       = $(CC)

CFLAGS       = -DPRNTlevel=0 -O3

FFLAGS       = -O3

FORTRAN	     = g77

#PATH TO THE COMPILED SUPERLU LIBRARY ARCHIVE
SUPERLULIB = /hera/bhs/lib/libsuperlu_4.3.a
#CFLAG SAYING HOW TO USE BLAS FUNCTIONS IN C IMPLEMENTATION
CDEFS        = -DAdd_

TESTPOISSON		= test_poisson.o

all: testpoisson

testpoisson: $(objects)
	$(CC) $(INCLS) $(objects) $(MYLIBS) $(LIBS) -o testpoisson

$(objects): $(sources) 
	$(CC) -c $(INCLS) $(LIBS) $(sources)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(INCLS) $(MYLIBS) $(LIBS) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:
	rm -rf *o testpoisson




