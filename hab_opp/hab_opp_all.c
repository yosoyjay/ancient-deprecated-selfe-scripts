/*! $Revision: 1.16 $
 *
 * @file hab_opportunity.c
 *
 * Copyright 1990-2004 Oregon Health and Science University
 *
 * Habitat opportunity.
 *
 * @author Original by ??, modified to support hybrid grid by pturner 11-2003
 * 
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "/usr/local/ace/include/elio.h"

#define MIN_ELV   0.1
#define MAX_ELV   2.0
#define MIN_SAL   0.0
#define MAX_SAL   5.0
#define MIN_VEL   0.0
#define MAX_VEL   0.25
#define MIN_TMP   0.0
#define MAX_TMP   19.0

#define min(X, Y) ((X) < (Y) ? (X) : (Y))

typedef struct _Region {
    char id[10];
    int nelems;
    int *ele_map;
    float area;
    float hab_elv;
    float hab_sal2;
    float hab_vel2;
    float hab_tmp2;
    float hab_all;
    float min_val;
    float elv_tmp;
    float sal_tmp;
    float tmp_tmp;
    float vel_tmp;
} Region;

main(int argc, char *argv[])
{
    FILE *felv, *fvel, *fsal, *ftmp, *fhab_elv, *fhab_vel2, *fhab_sal2, *freg, *freg_hab, *fhab_tmp2, *fnod_hab, *fhab_all;
    char fname[200], grtitle[200], inputDIR[200], outputDIR[200];
    double u, v, max_sal, min_sal, max_vel, min_vel, max_tmp, min_tmp, mag, ttime, lyr_dpth, *el_areas;
    double avg_sal, avg_vel, avg_tmp;
    double *hab_elv, *hab_sal2, *hab_vel2, *hab_tmp2, *hab_all;
    double *x, *y, cx, cy, x1, x2, x3, y1, y2, y3, min_val, elv_tmp, sal_tmp, tmp_tmp, vel_tmp;
    int Nvalues, file_count;
    int err, n1, n2, n3, npts;
    int i, j, k, l, time_step, nreg, lyr_count, Nfiles, header_skip, element;
    char buf[2048];
    Region *hab_reg;
    ElcircTimeStep elv, sal, vel, tmp;
    ElcircTimeStep out;
    ElcircHeader helv, hsal, hvel, htmp, oh;

    if (argc != 5) {
	fprintf(stderr, "Usage: %s <input directory> <output directory> <region file>  <number of input files> \n", argv[0]);
	fprintf(stderr, "Example: %s /home/workspace/ccalmr/hindcasts/2004-30-13/run /home/workspace/ccalmr/hindcasts/2004-30-13/process habitat_regions3.dat \
7\n", argv[0]);
	exit(1);
    }
   
    sprintf(inputDIR, "%s", argv[1]);
    sprintf(outputDIR, "%s", argv[2]);
    Nfiles = atoi(argv[4]);

    // read header info  
    sprintf(fname, "%s/%d_elev.61", inputDIR, 1);
    if ((felv = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", fname);
	exit(1);
    }
    if (ElioGetFileType(felv) == -1) {
	fprintf(stderr, "Incorrect file type for %s\n", fname);
	exit(1);
    }
    if (err = ElioGetHeader(fname, &helv)) {
	fprintf(stderr, "Error in GetElcircHeader(): Error # %d\n", err);
	exit(1);
    }

    sprintf(fname, "%s/%d_salt.63", inputDIR, 1);
    if ((fsal = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", fname);
	exit(1);
    }
    if (ElioGetFileType(fsal) == -1) {
	fprintf(stderr, "Incorrect file type for %s\n", fname);
	exit(1);
    }
    if (err = ElioGetHeader(fname, &hsal)) {
	fprintf(stderr, "Error in GetElcircHeader(): Error # %d\n", err);
	exit(1);
    }

    sprintf(fname, "%s/%d_temp.63", inputDIR, 1);
    if ((ftmp = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", fname);
	exit(1);
    }
    if (ElioGetFileType(ftmp) == -1) {
	fprintf(stderr, "Incorrect file type for %s\n", fname);
	exit(1);
    }
    if (err = ElioGetHeader(fname, &htmp)) {
	fprintf(stderr, "Error in GetElcircHeader(): Error # %d\n", err);
	exit(1);
    }

    sprintf(fname, "%s/%d_hvel.64", inputDIR, 1);
    if ((fvel = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", fname);
	exit(1);
    }
    if (ElioGetFileType(fvel) == -1) {
	fprintf(stderr, "Incorrect file type for %s\n", fname);
	exit(1);
    }
    if (err = ElioGetHeader(fname, &hvel)) {
	fprintf(stderr, "Error in GetElcircHeader(): Error # %d\n", err);
	exit(1);
    }

    // read region information  
    //sprintf(fname, "%s/%s", inputDIR, argv[3]);
    sprintf(fname, "%s", argv[3]);
    if ((freg = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "Could not find region file %s\n", fname);
	exit(1);
    }
    fgets(buf, 255, freg);	// one line of alpha  
    fgets(buf, 255, freg);	// number of regions  
    sscanf(buf, "%d", &nreg);
    hab_reg = (Region *) malloc(nreg * sizeof(Region));
    el_areas = (double *) malloc(helv.ne * sizeof(double));
    for (j = 0; j < helv.ne; j++) {
	el_areas[j] = ElioGetElementArea(&helv, j);
    }

    for (j = 0; j < nreg; j++) {
	fgets(buf, 255, freg);
	sscanf(buf, "%d %*d %s", &npts, hab_reg[j].id);	// number of point in polygon  
	hab_reg[j].area = 0.0;
	// Read the polygon for this region  
	x = (double *) malloc(npts * sizeof(double));
	y = (double *) malloc(npts * sizeof(double));
	if (x == NULL || y == NULL) {
	    fprintf(stderr, "%s: Can't allocate %d points for region %s.\n", argv[0], npts, hab_reg[j].id);
	    exit(1);
	}
	for (k = 0; k < npts; k++) {
	    fgets(buf, 255, freg);
	    sscanf(buf, "%lf %lf", &x[k], &y[k]);
	}
	// Read the polygon for this region  
	hab_reg[j].nelems = 0;
	hab_reg[j].ele_map = (int *) malloc(helv.ne * sizeof(int));
	for (k = 0; k < helv.ne; k++) {
	    if (ElioGetElementCenter(&helv, k, &cx, &cy) == ELIO_OK) {
		if (ElioInPolygon(cx, cy, npts, x, y)) {
		    hab_reg[j].ele_map[hab_reg[j].nelems] = k;
		    hab_reg[j].nelems++;
		    hab_reg[j].area += el_areas[k];
		    //printf("In %d %d\n", k, j);
		} else {
		    //printf("%d at (%lf, %lf)\n", k, cx, cy);
		}
	    } else {
	    }
	}
	free(x);
	free(y);
	//printf("%d\n", hab_reg[j].nelems);
    }
    fclose(freg);

    // prepare output files for isoline images  
    oh = helv;
    sprintf(fname, "%s/%s", outputDIR, "habitat_elv.61");
    if ((fhab_elv = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
    }
    if (ElioPutHeader(fhab_elv, &oh)) {
		fprintf(stderr, "Unable to write header to %s\n", fname);
		exit(1);
    }

    sprintf(fname, "%s/%s", outputDIR, "habitat_sal_AVG.61");
    if ((fhab_sal2 = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
    }
    if (ElioPutHeader(fhab_sal2, &oh)) {
		fprintf(stderr, "Unable to write header to %s\n", fname);
		exit(1);
    }

    sprintf(fname, "%s/%s", outputDIR, "habitat_tmp_AVG.61");
    if ((fhab_tmp2 = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
    }
    if (ElioPutHeader(fhab_tmp2, &oh)) {
		fprintf(stderr, "Unable to write header to %s\n", fname);
		exit(1);
    }

    sprintf(fname, "%s/%s", outputDIR, "habitat_vel_AVG.61");
    if ((fhab_vel2 = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
    }
    if (ElioPutHeader(fhab_vel2, &oh)) {
		fprintf(stderr, "Unable to write header to %s\n", fname);
		exit(1);
    }
    sprintf(fname, "%s/%s", outputDIR, "habitat_all_AVG.61");
    if ((fhab_all = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
    }
    if (ElioPutHeader(fhab_all, &oh)) {
		fprintf(stderr, "Unable to write header to %s\n", fname);
		exit(1);
    }

// to compile values for each habitat region  
    hab_elv = (double *) malloc(helv.np * sizeof(double));
    hab_sal2 = (double *) malloc(hsal.np * sizeof(double));
    hab_tmp2 = (double *) malloc(htmp.np * sizeof(double));
    hab_vel2 = (double *) malloc(hvel.np * sizeof(double));
    hab_all = (double *) malloc(hvel.np * sizeof(double));

    ElioAllocateTimeStep(&helv, &elv);
    ElioAllocateTimeStep(&hsal, &sal);
    ElioAllocateTimeStep(&htmp, &tmp);
    ElioAllocateTimeStep(&hvel, &vel);
    ElioAllocateTimeStep(&oh, &out);

// Loop over the files (1 file per day)  
    for (file_count = 1; file_count <= Nfiles; file_count++) {

// One time initialization  
	if (file_count != 1) {
	    sprintf(fname, "%s/%d_elev.61", inputDIR, file_count);
	    if ((felv = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
	    }

	    sprintf(fname, "%s/%d_salt.63", inputDIR, file_count);
	    if ((fsal = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
	    }

	    sprintf(fname, "%s/%d_temp.63", inputDIR, file_count);
	    if ((ftmp = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
	    }

	    sprintf(fname, "%s/%d_hvel.64", inputDIR, file_count);
	    if ((fvel = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "Unable to open file %s\n", fname);
		exit(1);
	    }
	}

// Loop over each time step  

	for (j = 0; j < helv.nsteps; j++) {
	    if (ElioGetTimeStep(felv, j, &helv, &elv)) {
	    }
	    if (ElioGetTimeStep(fsal, j, &hsal, &sal)) {
	    }
	    if (ElioGetTimeStep(ftmp, j, &htmp, &tmp)) {
	    }
	    if (ElioGetTimeStep(fvel, j, &hvel, &vel)) {
	    }
	    //fprintf(stdout, "time= %f step = %d\n", elv.t, elv.it);

// Loop over each node and accumulate hab op quantity in hours (for 7 days, 168 max)  
	    for (k = 0; k < helv.np; k++) {
// elevation  

//                  elv_tmp = 0.0;     // moved down below the wet/dry check under "Accumulate others"
// 		if ((elv.d[k] + helv.d[k]) >= MIN_ELV && (elv.d[k] + helv.d[k]) <= MAX_ELV) {  
// 		    elv_tmp = helv.timestep / 3600.0;  
// 		    hab_elv[k] += elv_tmp;  
// 		    ++dbg1;  
// 		}  
// 		else{  
// 		    elv_tmp = 0.0;  
// 		}  
// Accumulate others  
//		if ((elv.d[k] + helv.d[k]) > 0) {   // elv values can be greater than the depth even when dry.  Thus, check
	                                            // that the element is wet using salinity instead
	      if (sal.d[hsal.no[k]] >= 0) {
                    elv_tmp = 0.0;
	            if ((elv.d[k] + helv.d[k]) >= MIN_ELV && (elv.d[k] + helv.d[k]) <= MAX_ELV) {
		    elv_tmp = helv.timestep / 3600.0;
		        hab_elv[k] += elv_tmp;
		    }
		    else{
		        elv_tmp = 0.0;
		    }

		    max_sal = 0.0;
		    min_sal = 50.0;
		    avg_sal = 0.0;
		    max_tmp = 0.0;
		    min_tmp = 150.0;
		    avg_tmp = 0.0;
		    max_vel = 0.0;
		    min_vel = 10.0;
		    avg_vel = 0.0;
		    lyr_count = 0;
	    //fprintf(stdout, "%d %d\n", sal.surfind[k], hsal.bi[k]);
// Loop over the vertical column  
		    for (i = 0; i < (sal.surfind[k] - 1 - hsal.bi[k]); i++) {
	    //fprintf(stdout, "%f\n", sal.d[hsal.no[k] + i]);
			 if (sal.d[hsal.no[k] + i] >= 0.0) {
			    if (sal.d[hsal.no[k] + i] > max_sal) {
				max_sal = sal.d[hsal.no[k] + i];
			    }
			    if (sal.d[hsal.no[k] + i] < min_sal) {
				min_sal = sal.d[hsal.no[k] + i];
			    }
			    if (tmp.d[htmp.no[k] + i] > max_tmp) {
				max_tmp = tmp.d[htmp.no[k] + i];
			    }
			    if (tmp.d[htmp.no[k] + i] < min_tmp) {
				min_tmp = tmp.d[htmp.no[k] + i];
			    }
			    u = vel.d[hvel.no[k] + i * 2];
			    v = vel.d[hvel.no[k] + i * 2 + 1];
			    mag = sqrt(pow(u, 2.0) + pow(v, 2.0));
			    if (mag > max_vel) {
				max_vel = mag;
			    }
			    if (mag < min_vel) {
				min_vel = mag;
			    }
			    avg_sal += (double) sal.d[hsal.no[k] + i];
			    avg_tmp += (double) tmp.d[htmp.no[k] + i];
			    avg_vel += (double) mag;
			    lyr_count++;
			}
		    }
// Sum up  
		    if (lyr_count > 0) {
			avg_sal /= (double) lyr_count;
			avg_tmp /= (double) lyr_count;
			avg_vel /= (double) lyr_count;
			if (avg_sal <= MAX_SAL && avg_sal >= MIN_SAL) {
			    sal_tmp = hsal.timestep / 3600.0;
			    hab_sal2[k] += sal_tmp;
			}
			else{
			    sal_tmp = 0.0;
			}
			if (avg_tmp <= MAX_TMP && avg_tmp >= MIN_TMP) {
			    tmp_tmp = hsal.timestep / 3600.0;
			    hab_tmp2[k] += tmp_tmp;
			}
			else{
			    tmp_tmp = 0.0;
			}
			if (avg_vel <= MAX_VEL && avg_vel >= MIN_VEL) {
			    vel_tmp = hsal.timestep / 3600.0;
			    hab_vel2[k] += vel_tmp;
			}
			else{
			    vel_tmp = 0.0;
			}
// Do the combined analysis (based on vertically averaged hab params) 
                  	min_val = min(sal_tmp,elv_tmp);
                  	min_val = min(min_val,tmp_tmp);
                  	min_val = min(min_val,vel_tmp);
                  	
                  	//if (k == 3542){
                  	//fprintf(stderr,"hab_elv[k] is %f\n", hab_elv[k]);
                  	//fprintf(stderr,"hab_sal2[k] is %f\n", hab_sal2[k]);
                  	//fprintf(stderr,"hab_tmp2[k] is %f\n", hab_tmp2[k]);
                  	//fprintf(stderr,"hab_vel2[k] is %f\n", hab_vel2[k]);
                  	//fprintf(stderr,"hab_all[k] starts at %f\n", hab_all[k]);
                  	//}
                  	
                  	hab_all[k] += min_val;
                  	
                  	//if (k == 3542){
                  	//fprintf(stderr,"    And now is %f\n", hab_all[k]);
                  	//}
                  
//printf("A %d %d %.2f %.2f %.2lf %.2lf\n", j, k, elv.d[k], helv.d[k], avg_sal, avg_vel);
//printf("B %d %d %.2lf %.2lf %.2lf\n", j, k, hab_elv[k], hab_sal2[k], hab_vel2[k]);

		    }
//printf("%d %d %f %lf\n", j, k, elv.d[k], avg_tmp);
		}
	  }
	fclose(fsal);
	fclose(ftmp);
	fclose(felv);
	fclose(fvel);

    }

    // Write the total pho hours per node file 
    sprintf(fname, "%s/%s", outputDIR, "node_pho.ascii");
    if ((fnod_hab = fopen(fname, "w")) == NULL) {
	fprintf(stderr, "Could not open file %s \n", fname);
	exit(1);
    }
    fprintf(fnod_hab, "%d\n", helv.np);
    fprintf(fnod_hab,"Elevation \t Salinity \t Temperature \t Velocity \t Combined\n");
    for (k = 0; k <= helv.np; k++){
      fprintf(fnod_hab,"%f %f %f %f %f\n",hab_elv[k], hab_sal2[k], hab_tmp2[k], hab_vel2[k], hab_all[k]);
      //fprintf(stderr,"%f %f %f %f %f\n",hab_elv[k], hab_sal2[k], hab_tmp2[k], hab_vel2[k], hab_all[k]);
    }
    fclose(fnod_hab);
   

    // Write out one time step for each file  
    //printf("Writing time steps...\n");
    out.t = 0.0;
    out.it = 1;
    for (k = 0; k < helv.np; k++) {
		out.d[k] = hab_elv[k];
	    ElioPutTimeStep(fhab_elv, &oh, out);
    }
    for (k = 0; k < helv.np; k++) {
		out.d[k] = hab_sal2[k];
	    ElioPutTimeStep(fhab_sal2, &oh, out);
    }
    for (k = 0; k < helv.np; k++) {
		out.d[k] = hab_tmp2[k];
	    ElioPutTimeStep(fhab_tmp2, &oh, out);
    }
    for (k = 0; k < helv.np; k++) {
		out.d[k] = hab_vel2[k];
	    ElioPutTimeStep(fhab_vel2, &oh, out);
    }
    for (k = 0; k < helv.np; k++) {
        out.d[k] = hab_all[k];
	    ElioPutTimeStep(fhab_all, &oh, out);
    }
    fclose(fhab_elv);
    fclose(fhab_sal2);
    fclose(fhab_tmp2);
    fclose(fhab_vel2);
    fclose(fhab_all);

// Write out summary quantities  

    for (j = 0; j < nreg; j++) {
		hab_reg[j].hab_elv = 0.0;
		hab_reg[j].hab_sal2 = 0.0;
		hab_reg[j].hab_tmp2 = 0.0;
		hab_reg[j].hab_vel2 = 0.0;
		hab_reg[j].hab_all = 0.0;
		for (k = 0; k < hab_reg[j].nelems; k++) {
		    int nn;
		    element = hab_reg[j].ele_map[k];
// Sum up for each node in the element 
		    for (nn = 0; nn < helv.etype[element]; nn++) {
			n1 = helv.icon[nn][element];
			hab_reg[j].hab_elv += hab_elv[n1] / helv.etype[element] * el_areas[element];
			hab_reg[j].hab_sal2 += hab_sal2[n1] / helv.etype[element] * el_areas[element];
			hab_reg[j].hab_tmp2 += hab_tmp2[n1] / helv.etype[element] * el_areas[element];
			hab_reg[j].hab_vel2 += hab_vel2[n1] / helv.etype[element] * el_areas[element];
			hab_reg[j].hab_all += hab_all[n1] / helv.etype[element] * el_areas[element];
		    }
//if (j == 0) printf("%d %d %d %.0lf %.0lf %.0lf %.0lf %.0lf\n", j, k, element, hab_reg[j].hab_elv, hab_reg[j].hab_sal2, hab_reg[j].hab_vel2, hab_reg[j].hab_all, hab_reg[j].area);
		}

//printf("%d %d %.0lf %.0lf %.0lf %.0lf %.0lf %.0lf\n", j, hab_reg[j].nelems, hab_reg[j].hab_elv, hab_reg[j].hab_tmp2, hab_reg[j].hab_sal2, hab_reg[j].hab_vel2, hab_reg[j].hab_all, hab_reg[j].area);

// Adjust by the area of the region 
		hab_reg[j].hab_elv  /= hab_reg[j].area;
		hab_reg[j].hab_sal2 /= hab_reg[j].area;
		hab_reg[j].hab_tmp2 /= hab_reg[j].area;
		hab_reg[j].hab_vel2 /= hab_reg[j].area;
		hab_reg[j].hab_all  /= hab_reg[j].area;
    }

    sprintf(fname, "%s/%s", outputDIR, "regional_habitat.ascii");
    if ((freg_hab = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", fname);
		exit(1);
    }
    fprintf(freg_hab, "%d\n", nreg);
    fprintf(freg_hab,
	        "Region \t\
			 Habitat (Elevation criterion) \t\
			 Habitat (Average Salinity Criterion) \t\
			 Habitat (Average Velocity Criterion) \t\
			 Habitat (Average Temperature Criterion) \t\
			 Habitat (Average Combined Criterion) \n");
    for (j = 0; j < nreg; j++) {
		fprintf(freg_hab, "%s %lf %lf %lf %lf %lf\n", 
		        hab_reg[j].id, 
				hab_reg[j].hab_elv, 
				hab_reg[j].hab_sal2, 
				hab_reg[j].hab_vel2, 
				hab_reg[j].hab_tmp2, 
				hab_reg[j].hab_all);
    }
    fclose(freg_hab);

  }
}
