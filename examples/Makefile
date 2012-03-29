#       
# Makefile for building a program that calls a routine from
# the CLAPACK Library
#
# Dr. Andrew J. Pounds
# Departments of Chemistry and Computer Science
# Mercer University
#
# November 2004
#

CC = gcc
ROOTPATH = /home/wyao/cosamp/CLAPACK-3.2.1
INCDIRS = -I$(ROOTPATH)/SRC -I$(ROOTPATH) \
          -I$(ROOTPATH)/INCLUDE
F2CDIR  = $(ROOTPATH)/F2CLIBS
LDLIBS  = $(ROOTPATH)/lapack_LINUX.a \
          $(ROOTPATH)/blas_LINUX.a \
	  $(F2CDIR)/libf2c.a -lm

dgesv_example: dgesv_example.o
	$(CC) dgesv_example.o  $(LDLIBS)  -o dgesv_example 

dgesv_example.o: dgesv_example.c
	$(CC) dgesv_example.c  $(INCDIRS) -c 

zgesv_example: zgesv_example.o
	$(CC) zgesv_example.o  $(LDLIBS)  -o zgesv_example 

zgesv_example.o: zgesv_example.c
	$(CC) zgesv_example.c  $(INCDIRS) -c 
