CC = gcc
ROOTPATH = /home/wyao/cosamp/CLAPACK-3.2.1
INCDIRS = -I$(ROOTPATH)/SRC -I$(ROOTPATH) \
          -I$(ROOTPATH)/INCLUDE
F2CDIR  = $(ROOTPATH)/F2CLIBS
LDLIBS  = $(ROOTPATH)/lapack_LINUX.a \
          $(ROOTPATH)/blas_LINUX.a \
	  $(F2CDIR)/libf2c.a -lm

cosamp: cosamp.o
	$(CC) cosamp.o $(LDLIBS) -O3 -o cosamp

cosamp.o: cosamp.c
	$(CC) cosamp.c $(INCDIRS) -c
