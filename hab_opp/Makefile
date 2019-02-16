# 
# Sample Makefile for use of libelio.a.
#
# $Revision: 1.1 $
# $Log: main.html,v $
# Revision 1.1  2004/03/10 19:43:21  pturner
# More additions.
#

# -miee for alphas
# CFLAGS = -mieee
# for intel
CFLAGS = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I./
LIB = /usr/local/ace/lib
CC = gcc

hab_opp: hab_opp.o 
	$(CC) $(CFLAGS) hab_opp.c -L$(LIB) -lelio -lm -o hab_opp
