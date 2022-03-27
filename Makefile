#****************************************
#    Makefile for sxrdcalc		*
#   with gcc compiler for linux		*
#****************************************

# Compiler and options 
CC=gcc
#CFLAGS=  -O2 -ffast-math -W -Wall -ansi -pedantic -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fno-common -Wnested-externs  
#CFLAGS=  -O2 -W -Wall -ansi -pedantic -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fno-common -Wnested-externs -g  
CFLAGS=  -O2  -ansi -fshort-enums -fno-common -Wnested-externs -g  
#CFLAGS=  -W -Wall -ansi -pedantic -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fno-common -Wnested-externs  -g

LIBS= -lgsl -lgslcblas -lm
OBJS= sxrdcalc.o sxrdio.o calcs.o
# files
all:	sxrdcalc

sxrdcalc: $(OBJS)
	$(CC) $(CFLAGS) -o sxrdcalc $(OBJS) $(LIBS)

sxrdcalc.o: sxrdcalc.c sxrddefs.h sxrdio.h calcs.h
	$(CC) $(CFLAGS) -c sxrdcalc.c  

sxrdio.o: sxrdio.c sxrddefs.h sxrdio.h 
	$(CC) $(CFLAGS) -c  sxrdio.c

calcs.o: calcs.c sxrddefs.h calcs.h
	$(CC) $(CFLAGS) -c calcs.c

clean: 
	rm -f sxrdcalc sxrdcalc.o sxrdio.o calcs.o
