/*  June 17, 2002 Charles Seaton Version 1.0
Uses projects.h geospatial library to convert between grid coordiantes and nad83 y-x coordiantes

Generates the rotation at each node involved in converting from geographic to projected coordinates. Converts between NAD83 y-x and NAD27 ORSPCS-NAD27 meters */

/* cc -c -I/usr/local/include rotate_wind_spcs2ll.c ;
 cc -o rotate_wind_spcs2ll rotate_wind_spcs2ll.o -L/usr/local/lib -lproj -lm
*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <projects.h>

int debug=0;
   static  PJ *spcs, *ll, *WO;
  static  char *parms_spcs[] = {"init=nad27:3601"};
  static char *parms_ll[] = {"proj=latlong","datum=WGS84"};
//  static char *parms_WO[] = {"proj=latlong","nadgrids=WO"};
    static int n_parm_spcs = 1;
    static int n_parm_ll = 2;
    static int n_parm_WO = 2;
/* use proj cartographic library functions to convert high precision nad83 y-x data from Washington, Oregon or Northern California into nad27 Oregon (North) state plane coordinate system. DEG_TO_RAD and RAD_TO_DEG are constants which convert degrees to radians and radians to degrees*/
int convert2ll (double *y, double *x) 
{
if (debug) {	printf("converting %f %f to y-x\n",y[0],x[0]);}
    char teststr[20]; 
    double x_tmp[1],y_tmp[1];
   /*initialize conversion table from nad27 y-x to nad27 OR-N SPCS*/ 
/* convert to y-x from OR-N SPCS */
   x_tmp[0]=x[0];
   y_tmp[0]=y[0];
   /* pj_transform(spcs,WO,1,1,x_tmp,y_tmp,NULL);
	if (debug){
   		 printf ("converted SPCS to high precision %f %f\n",x[0]*RAD_TO_DEG,y[0]*RAD_TO_DEG);
	}
    if (x_tmp[0]!=x[0] & y_tmp[0]!=y[0]){
	x[0]=x_tmp[0];
	y[0]=y_tmp[0];
    }
    else { */
    	pj_transform(spcs,ll,1,1,x,y,NULL);
	if (debug){
   		 printf ("converted SPCS to NAD83 %f %f\n",x[0]*RAD_TO_DEG,y[0]*RAD_TO_DEG);
	}
   // }
 x[0] =	x[0]*RAD_TO_DEG,
   y[0] = y[0]*RAD_TO_DEG;
    return (0);
}

/* use proj cartographic library functions to convert high precision nad83 y-x data from Washington, Oregon or Northern California into nad27 Oregon (North) state plane coordinate system. DEG_TO_RAD and RAD_TO_DEG are constants which convert degrees to radians and radians to degrees*/
int convert2spcs (double *y, double *x) 
{
if (debug){	printf("converting lat %f lon %f to spcs\n",y[0],x[0]);}
    char teststr[200];
    double z[1],x_tmp[1],y_tmp[1],z_tmp[1];
      x[0] = x[0] * DEG_TO_RAD;
      y[0] = y[0] * DEG_TO_RAD;
      z[0] = 0;
    x_tmp[0]=x[0];
    y_tmp[0]=y[0];
    z_tmp[0]=z[0];
if (debug){	printf("convert values to radians %f %f\n",y[0],x[0]);}
    /*Washington and Oregon area GPS lat-long is in WO high precision datum, 
     * convert directly to spcs and test if results are valid, if not , then convert from standard NAD87 to spcs
     */ 

   /* pj_transform(WO,spcs,1,1,x_tmp,y_tmp,z_tmp);
    if (debug){printf ("converted WO region to SPCS %f %f\n",y_tmp[0],x_tmp[0]);}
    if (x_tmp[0]!=x[0] & y_tmp[0]!=y[0]){
	x[0]=x_tmp[0];
	y[0]=y_tmp[0];
	z[0]=z_tmp[0];
    }
    else { */
     	pj_transform(ll,spcs,1,0,x,y,z);
	if (debug){printf ("converted NAD83 to SPCS %f %f\n",y[0],x[0]);}
     // }
    return (0);
} 

void usage(char *progname)
{
fprintf(stderr, "Usage: %s -debug -input [input file] -output [output file] <-ll2spcs or -spcs2ll>\n reads in a grid file and generates the rotation created by converting from ORSPCS-N NAD27 to NAD 83 lat-lon (-spcs2ll) or the other way (-ll2spcs)\nThe input grid should be in the coordinate system that the rotation is from\ne.g %s -input hgrid.ll -output hgrid_rot.gr3 -ll2spcs\nWill return the rotation from lat-lon to ORSPCS-N\n",progname,progname);

    exit(1);
}


int main (int argc, char **argv)
{


	int convert2spcs (double * y, double * x);
	int convert2ll (double * y, double * x);
	void usage(char *progname);
	double y[1], x[1], xshift[1], yshift[1];
	double tmpx, tmpy;
	double *gridy, *gridx;
	float *depth;
	int *ne;
	char buf[256];
	FILE *fgri, *fgrout;

	int convert, input=0, output=0, cnv_arraytest=0;

	int nnodes,nelems;
	int numtmp, i;
	char fname[100], fout[100], grtitle[256];
	/*if (! (WO = pj_init(n_parm_WO,parms_WO))){
		printf("initialization of WO failed because:  %s \n",pj_strerrno(pj_errno));
		return(1);
	} */	
	if (! (ll = pj_init(n_parm_ll,parms_ll))){
		printf("initialization of ll failed because:  %s \n",pj_strerrno(pj_errno));
		return(1);
	}	
	if (! (spcs = pj_init(n_parm_spcs,parms_spcs))){
		printf("initialization of spcs failed because:%s \n",pj_strerrno(pj_errno));
		return(1);
	}	

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if ((strcmp(argv[i], "-debug")) == 0) {
				debug = 1;
			}             
			else if ((strcmp(argv[i],"-ll2spcs"))==0){
				convert=1;
				cnv_arraytest=1;
			}
			else if ((strcmp(argv[i],"-spcs2ll"))==0){
				convert=0;
				cnv_arraytest=1;
			} else if ((strcmp(argv[i], "-input")) == 0) {
				i++;
				strcpy(fname, argv[i]);
				input = 1;
			} else if ((strcmp(argv[i], "-output")) == 0) {
				i++;
				strcpy(fout, argv[i]);
				output = 1;
			} else {
				fprintf(stderr, "Unknown commandline argument: %s\n", argv[i]);
				usage(argv[0]);
			}
		} else {
			fprintf(stderr, "Missing flag or commandline argument: %s\n", argv[i
					]);
			usage(argv[0]);
		}
	}
	if ((output *input * cnv_arraytest)==0)  {
		fprintf(stderr, "Missing flag or commandline argument: %s\n", argv[i
				]);
		usage(argv[0]);
	}  
	if (debug) {printf("spcs parameters are:\n");fflush(stdout);pj_pr_list(spcs);fflush(stdout);}
	if (debug) {printf("ll parameters are:\n");fflush(stdout);pj_pr_list(ll);fflush(stdout);}
	//if (debug) {printf("WO parameters are:\n");fflush(stdout);pj_pr_list(WO);fflush(stdout);}
	if ((fgri = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "Unable to open grid file %s\n", fname);
		exit(1);
	}
	fgets(grtitle,255, fgri);
	fprintf(stderr, "Grid title is %s\n", grtitle);
	fscanf(fgri, "%i %i", &nelems, &nnodes);
	gridx = (double *) calloc(nnodes, sizeof(double));
	gridy = (double *) calloc(nnodes, sizeof(double));
	ne = (int *) calloc(nelems, sizeof(int));
	depth = (float *) calloc(nnodes, sizeof(float));
	fprintf(stderr, "%i nodes, %i elements\n", nnodes, nelems);
	for (i=0; i<nnodes; i++) {
		fscanf(fgri, "%i %lf %lf %f", &numtmp, &gridx[i], &gridy[i], &depth[i]);		
		if(debug){printf("%i %18.8f %18.8f %f\n",numtmp,gridx[i],gridy[i],depth[i]);}
		if (debug) {fprintf(stderr,"reading node %d\n",i);}
		fflush(stderr);
	}
	if (debug) {fprintf(stderr,"finished reading nodes\n");}
	fflush(stderr);
	/*close grid file without reading table of elements*/
	fclose(fgri);
	if ((fgrout = fopen(fout, "w")) == NULL) {
		fprintf(stderr, "Unable to write grid file %s\n", fout);
		exit(1);
	}   
	fprintf(fgrout,"%s %d %d\n",grtitle,nelems,nnodes);
	for (i=0; i<nnodes; i++){
		if (debug){     printf("%d\n",i);}
		if (convert==1){
			x[0] = gridx[i];
			y[0] = gridy[i];
			xshift[0] = gridx[i];
			yshift[0] = gridy[i]+0.00001;

			if(debug){printf("x %f y %f\n",gridx[i],gridy[i]);} 
			if(debug){printf("%f %f\n",x[0],y[0]);} 
			convert2spcs(y,x);
			convert2spcs(yshift,xshift);
			if (debug){	printf("%f %f\n",x[0],y[0]);} 
			if (debug){	printf("%f %f\n",xshift[0],yshift[0]);} 
			fprintf(fgrout,"%i %12.8f %12.8f %f\n",i+1,gridx[i],gridy[i],atan2(-x[0] + xshift[0],-y[0] + yshift[0])*180/M_PI);
			if (debug){     printf("%f %f %f %f\n",-x[0] + xshift[0],-y[0] + yshift[0], (-x[0] + xshift[0])/(-y[0] + yshift[0]), atan2(-x[0] + xshift[0],-y[0] + yshift[0])*180/M_PI);}
		}
		else if (convert == 0) {
			x[0] = gridx[i];
			y[0] = gridy[i];
			if(debug){printf("x %f y %f\n",gridx[i],gridy[i]);} 
			if(debug){printf("spcs %f %f\n",x[0],y[0]); }
			convert2ll(y,x);
			if (debug ) {printf("ll %f %f\n",x[0],y[0]);} 
			xshift[0] = x[0];
			yshift[0] = y[0]+0.00001;
			if (debug ) {printf("ll shifted %f %f\n",xshift[0],yshift[0]);} 
			convert2spcs(yshift,xshift);
			if (debug ) {printf("spcs shifted %f %f\n",xshift[0],yshift[0]);} 
			fprintf(fgrout,"%i %16.6f %16.6f %f\n",i+1,gridx[i],gridy[i],-1*atan2(-gridx[i] + xshift[0],-gridy[i] + yshift[0])*180/M_PI);
			if (debug){     printf("%f %f %f %f\n",-gridx[i] + xshift[0],-gridy[i] + yshift[0], (-gridx[i] + xshift[0])/(-gridy[i] + yshift[0]), atan2(-gridx[i] + xshift[0],-gridy[i] + yshift[0])*180/M_PI);}
		}
	}
	if (debug ) {printf("starting to write table of elements\n");} 
	/* reopen grid file, ignore table of nodes, write table of elemetns to out file*/
	if ((fgri = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "Unable to open grid file %s\n", fname);
		exit(1);
	}
	fgets(buf, 127, fgri);
	fgets(buf,127,fgri);
	for (i=0; i<nnodes; i++) {
		fgets(buf,127,fgri);
	}
	for (i=0; i<nelems; i++) {
		fgets(buf,127,fgri);
		fprintf(fgrout,"%s",buf);
	}
	fclose(fgri);
	fclose(fgrout);
}
