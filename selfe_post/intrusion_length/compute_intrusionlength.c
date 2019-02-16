/*! $Id: compute_intrusionlength.c,v 1.5 2004/11/02 00:43:55 pturner Exp $
 *
 * @file compute_intrusionlength.c
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 * Compute intrusion length from model results.
 * pturner 11-2003
 *
 * $Log: compute_intrusionlength.c,v $
 * Revision 1.5  2004/11/02 00:43:55  pturner
 * Still working on sigma grids and hab opportunity.
 *
 * Revision 1.4  2004/03/03 03:17:13  pturner
 * Worked on documentation.
 *
 * Revision 1.3  2003/11/22 18:30:23  pturner
 * Working version.
 *
 * Revision 1.2  2003/11/21 23:43:09  pturner
 * Generic commit.
 *
 * Revision 1.1  2003/11/20 17:03:05  pturner
 * Port of code to new library.
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "elio.h"

typedef struct _Stations {
    char name[128];
    double x, y;
    double avgdepth;
    double dist;
    double maxsal;
    int elem;
    double h[4];
    int bind[4];
    int sind[4];
    int n;
    double *s;
} Stations;

typedef struct _Var {
    char fname[2048];
    FILE *fp;
    int nsteps;
    float *d;
    ElcircHeader h;
    ElcircTimeStep t;
} Var;

int main(int argc, char **argv)
{
    char stafile[2048], saltfile[2048], outpfile[2048];
    char buf[2048];
    char tz[100]; /* time zone */
    int i, j, k, ilev, nn, nsta, start = 0, step;
    FILE *fp, *fout;	/* used to read in the region and the data set */
    Stations *s;
    Var v;
    int nsteps, err = 0;
    int yy, mo, dd, hh, mm, ss;
    double cstart;
    double salt_intrusion;
    double length_intrusion, lx, ly;

    if (argc != 5) {
	fprintf(stderr, "Usage : %s <station file> <salt file> <output file> <salt intrusion limiting value>\n", argv[0]);
	exit(1);
    }

    sprintf(stafile, "%s", argv[1]);
    sprintf(v.fname, "%s", argv[2]);
    sprintf(outpfile, "%s", argv[3]);
    salt_intrusion = atof(argv[4]);

    if ((fout = fopen(outpfile, "w")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", outpfile);
	exit(1);
    }

    if ((v.fp = fopen(v.fname, "rb")) == NULL) {
	fprintf(stderr, "Could not open file %s \n", v.fname);
	exit(1);
    }
    if (ElioGetFileType(v.fp) == -1) {
	fprintf(stderr, "Incorrect file type for %s\n", v.fname);
	exit(1);
    }
    if (err = ElioGetHeader(v.fname, &v.h)) {
	fprintf(stderr, "Error in GetElcircHeader(): Error # %d\n", err);
	exit(1);
    }
    v.nsteps = ElioGetNStepsInFile(v.fname, &v.h);
    v.t.surfind = (int *) malloc(v.h.np * sizeof(int));
    v.t.e = (float *) malloc(v.h.np * sizeof(float));
    v.t.d = (float *) malloc(v.h.np * v.h.nvrt * v.h.ivs * sizeof(float));
    v.d = (float *) malloc(v.h.nvrt * v.h.ivs * sizeof(float));

    /* get static info, assumed the same for all files (maybe need to check this) */
    nsteps = v.nsteps;
    sscanf(v.h.start_time, "%2d/%2d/%4d %2d:%2d:%2d %s", &mo, &dd, &yy, &hh, &mm, &ss, tz);
    ElioSetCorieTime(1);
    cstart = ElioGetDay(mo, dd, yy, hh, mm, ss);
    //printf("CDay = %lf %d %d %d %d %d %d\n", cstart, mo, dd, yy, hh, mm, ss);

    /* open station file */
    if ((fp = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "Could not open file %s \n", argv[1]);
	exit(1);
    }
    fgets(buf, 255, fp);
    sscanf(buf, "%d", &nsta);
    s = (Stations *) malloc(nsta * sizeof(Stations));
    for (i = 0; i < nsta; i++) {
	fgets(buf, 255, fp);
	sscanf(buf, "%s", s[i].name);
	fgets(buf, 255, fp);
	sscanf(buf, "%lf %lf", &s[i].x, &s[i].y);
	if (i > 0) {		/* compute distances */
	    s[i].dist = s[i - 1].dist + hypot(s[i].x - s[i - 1].x, s[i].y - s[i - 1].y);
	} else {
	    s[i].dist = 0.0;
	}
	s[i].maxsal = 0.0;
	/* find the element that contains the point and compute the shape functions */
	if ((s[i].elem = ElioFindElementXY(&v.h, s[i].x, s[i].y)) < 0) {
	    fprintf(stderr, "elio.ElioFindElementXY(): Unable to locate element at station: %s, x = %lf, y = %lf.\n", s[i].name, s[i].x, s[i].y);
	} else {
	    ElioGetCoefficients(&v.h, s[i].elem, s[i].x, s[i].y, s[i].h);
	    /* compute the average depth for this element */
	    s[i].avgdepth = 0;
	    for (j = 0; j < v.h.etype[s[i].elem]; j++) {
		nn = v.h.icon[j][s[i].elem];
		s[i].avgdepth += v.h.d[nn];
	    }
	    s[i].avgdepth = s[i].avgdepth / v.h.etype[s[i].elem];
	}
	s[i].s = (double *) malloc(v.h.nvrt * sizeof(double));
    }
    fclose(fp);


    for (i = start; i < nsteps; i++) {
	//for (i = 0; i < 1; i++) {
	/* Read a complete time step from each variable's file */
	if (ElioGetTimeStep(v.fp, i, &v.h, &v.t)) {
	    fprintf(stderr, "Error in ElioGetTimeStep()\n");
	    goto out;
	}
	for (j = 0; j < nsta; j++) {
	    if (s[j].elem >= 0) {	/* elem = -1 ==> station not in domain */
		ElioInterpTimeStep(&v.h, s[j].elem, s[j].x, s[j].y, s[j].h, &v.t, s[j].bind, s[j].sind, v.d);
		for (k = imin(v.h.etype[s[j].elem], s[j].bind); k <= imax(v.h.etype[s[j].elem], s[j].sind); k++) {
		    s[j].s[k] = v.d[k];
		}
		for (k = 0; k < v.h.nvrt; k++) {
		    if (imin(v.h.etype[s[j].elem], s[j].bind) <= k 
			&& imax(v.h.etype[s[j].elem], s[j].sind) >= k) {
			if (s[j].maxsal < s[j].s[k]) {
			    s[j].maxsal = s[j].s[k];
			}
		    }
		}
	    }
	}			/* j */

	j = 0;
	while (s[j].maxsal > salt_intrusion && j < nsta) {
	    j++;
	}
	//printf("%d: %lf %lf\n", j, s[j].maxsal, salt_intrusion);
	if (j == 0) {
	    length_intrusion = 0.0;
	} else if (j == nsta) {
	    length_intrusion = s[j - 1].dist;
	} else {
	    length_intrusion = s[j - 1].dist + (s[j].dist - s[j - 1].dist) / (s[j - 1].maxsal - s[j].maxsal) * (s[j - 1].maxsal - salt_intrusion);
	    lx = s[j - 1].x + (s[j].x - s[j - 1].x) / (s[j - 1].maxsal - s[j].maxsal) * (s[j - 1].maxsal - salt_intrusion);
	    ly = s[j - 1].y + (s[j].y - s[j - 1].y) / (s[j - 1].maxsal - s[j].maxsal) * (s[j - 1].maxsal - salt_intrusion);
	}
	/*fprintf(stdout, "time = %f time step = %d intrusion length = %f\n", v.t.t / 86400 + cstart, v.t.it, length_intrusion);*/
	fprintf(fout, "%.12lf %f %f %f\n", (double) ((v.t.t / 86400.0) + cstart), length_intrusion, lx, ly);
	for (j = 0; j < nsta; j++) {
	    s[j].maxsal = 0.0;
	}

    }				/* i */
  out:;
    free(v.t.surfind);
    free(v.t.d);
    free(v.d);
    fclose(v.fp);
    fclose(fout);
    for (i = 0; i < nsta; i++) {
	free(s[i].s);
    }
}

int imin(int n, int *iar)
{
    int i, im = 0;
    im = iar[0];
    for (i = 1; i < n; i++) {
	if (im > iar[i]) {
	    im = iar[i];
	}
    }
    return im;
}

int imax(int n, int *iar)
{
    int i, im = 0;
    im = iar[0];
    for (i = 1; i < n; i++) {
	if (im < iar[i]) {
	    im = iar[i];
	}
    }
    return im;
}
