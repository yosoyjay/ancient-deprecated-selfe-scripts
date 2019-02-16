# $Revision$
# 
# Makefile for compute_intrusionlength.c
#
# pturner
#
# $Log$
#
INSTALLDIR = /usr/local/ace

# -miee for alphas
CFLAGS = -mieee -I./
# for intel
CFLAGS = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I/usr/local/ace/include
LIB = /usr/local/ace/lib
LDFLAGS = -L${LIB} -lelio -lm

compute_intrusionlength: compute_intrusionlength.o ${LIB}/libelio.a
	$(CC) $(CFLAGS) -g compute_intrusionlength.o -L${LIB} -lelio -lm -o compute_intrusionlength
