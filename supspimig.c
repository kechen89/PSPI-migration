/* Copyright (c) Signal Analysis and Imaging Group (SAIG), University of Alberta, 2014.*/
/* All rights reserved.*/
/* supspimig:  $Date: 2014/02/24 */

#include "su.h"         /* include file for SU programs */ /*interface*/
#include "segy.h"       /* include file for SEGY traces */
#include "header.h"     /* include file for segy sizes */
#include <signal.h>     /* include library that defines signal handling functions */
#include <fftw3.h>      /* include FFTW library */
                        /* SU has its own complex number library*/
/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                            ",
" SUPSPIMIG - phase-shift plus interpolation migration (depth migration)     ",
"             for zero-offset data                                           ",
"                                                                            ",
" supspimig <infile >outfile vfile= [optional parameters]                    ",
"                                                                            ",
" Required Parameters:							     ",
" nz:		number of depth samples					     ",
" dz:		depth sampling interval					     ",
" vfile:	name of file containing velocities			     ",
"		                                                	     ",
" Optional parameters:                                                       ",
" dt=tr.dt:       time sampling interval (from header)                       ",
" nx=ntr:         number of traces (counted from data)                       ",
" dx:             midpoint sampling interval (from header)                   ",
" tmpdir: 	 if non-empty, use the value as a directory path	     ",
"		 prefix for storing temporary files; else if the	     ",
"	         the CWP_TMPDIR environment variable is set use		     ",
"	         its value for the path; else use tmpfile()		     ",
NULL};

/* 
 * Credits: 
 * Ke Chen
 * Last changes: Feb, 2014
 * Trace header fields accessed: ns,dt,d2
 */
/**************** end self doc ********************************/


/**************************prototypes for functions used internally****************************/
void pspimig(float **data,complex **image,float **v,int nt,int nx,int nz,float dt,float dx, float dz);
static void closefiles(void);


/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	             /*filename for the file of traces*/
char headerfile[BUFSIZ];             /*filename for the file of headers*/
FILE *tracefp;		             /* fp for trace storage file */
FILE *headerfp;		             /* fp for header storage file */

/*segy trace*/
segy tr;                             /*define the type of header variable*/ /*typedef struct*/

int main(int argc, char **argv)      /*argc, argv - the arguments to the main() function*/
{ 
int nt;                              /*number of time samples*/
int nz;			             /*number of migrated depth samples*/
int nx;                              /*number of midpoints (traces)*/
int ix;
int iz;

float dt;                            /*time sampling interval*/                
float dx;                            /*spatial sampling interval*/
float dz;                            /*migrated depth sampling interval*/           
float **data;                        /*input seismic data*/
complex **image;                     /*migrated image*/      
float **rimage;                      /*migrated image*/ 
float **v;                           /*velocity model*/
FILE *vfp;

char *vfile="";                      /*name of velocity file*/
int verbose=1;
char *tmpdir;		             /* directory path for tmp files*/
cwp_Bool istmpdir=cwp_false;         /* true for user-given path*/

/******************************* Intialize *********************************************/
initargs(argc,argv);
requestdoc(1);

/********************************* Get parameters **************************************/
/*get info from first trace*/
if (!gettr(&tr))  err("can't get first trace");  /*fgettr: get a fixed-length segy trace from a file by file pointer*/
nt = tr.ns;                         /*nt*/       /*gettr: macro using fgettr to get a trace from stdin*/

if (!getparfloat("dt", &dt)) {      /*dt*/
if (tr.dt) { 
dt = ((double) tr.dt)/1000000.0;
} 
else {err("dt is not set");}
}

if (!getparfloat("dx", &dx)) {       /*dx*/
if (tr.d2) {
dx = tr.d2;
} 
else {
err("dx is not set");
}
}

/*get optional parameters*/
if (!getparint("nz",&nz)) err("nz must be specified"); 
if (!getparfloat("dz",&dz)) err("dz must be specified");
if (!getparstring("vfile", &vfile)) err("velocity file must be specified");
if (!getparint("verbose", &verbose)) verbose = 0;
/****************************************************************************************/

/* Look for user-supplied tmpdir */
if (!getparstring("tmpdir",&tmpdir) &&
 !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
err("you can't write in %s (or it doesn't exist)", tmpdir);
checkpars();

/**************************** Count trace number nx ******************************/
/* store traces and headers in tempfiles while getting a count */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} 
     else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}

	nx = 0;
	do {
		 ++nx;                                   /*get the number of traces nx*/
		efwrite(&tr,HDRBYTES,1,headerfp);
		efwrite(tr.data, FSIZE, nt, tracefp);
	} while (gettr(&tr));

	erewind(tracefp);                    /*Set position of stream to the beginning*/
	erewind(headerfp);

/******************************************************************************************/

/*allocate memory*/
data = alloc2float(nt,nx);                   /*2D array nx by nt*/
image = alloc2complex(nz,nx);                /*2D array nx by nz*/
rimage = alloc2float(nz,nx);                 /*2D array nx by nz*/
v= alloc2float(nz,nx);                       /*2D array, in Fortran the velocity model is nz by nx 2D array*/ 
                                             /*in binary, it is actually 1D*/

/* load traces into the zero-offset array and close tmpfile */
efread(*data, FSIZE, nt*nx, tracefp);        /*read traces to data*/
efclose(tracefp);                 

/*load velicoty file*/
vfp=efopen(vfile,"r");	
efread(v[0],FSIZE,nz*nx,vfp);                    /*load velocity*/
efclose(vfp);			

/***********************finish reading data*************************************************/
/* call pspi migration function*/
pspimig(data,image,v,nt,nx,nz,dt,dx,dz);

/*get real part of image*/
for (iz=0;iz<nz;iz++){
for (ix=0;ix<nx;ix++){
rimage[ix][iz] = image[ix][iz].r;
}
}

/* restore header fields and write output */
for (ix=0; ix<nx; ix++) {
efread(&tr,HDRBYTES,1,headerfp);
tr.ns = nz;
tr.d1 = dz;
memcpy( (void *) tr.data, (const void *) rimage[ix],nz*FSIZE);
puttr(&tr);
}
	
/* Clean up */
efclose(headerfp);
if (istmpdir) eremove(headerfile);
if (istmpdir) eremove(tracefile);
return(CWP_Exit());	
}

/************************PSPI migration************************/
void pspimig(float **data,complex **image,float **v,int nt,int nx,int nz,float dt,float dx, float dz)
{
int ntfft;                              /*number of samples of the zero padded trace*/
int nxfft;
int nw;                                 /*number of temporal freqs.*/
int it;                                 /*loop index over time sample*/
int ix;                                 /*loop index over midpoint sample*/     
int iw;                                 /*loop index over frequency*/                 
int ik;                                 /*loop index over wavenumber*/
int iz;                                /*loop index over migrated depth samples*/
int iv;                                /*loop index over reference velocities*/
int nvref_max=8;                      /*number of reference velocities in each layer*/
int nvref;
int i1;                               /*nearest reference velocity*/
int i2;
int iw1;
int iw2;


float **vref;                           /*2D reference velocity array*/
float w0;                               /*first frequency sample*/
float w;                                /*frequency*/
float dw;                               /*frequency sampling interval*/                 
float k0;                               /*first wavenumber*/
float k;
float dk;                             /*wave number sampling interval in x*/
float dv;                              /*velocity interval*/
float *in;                             /*input 1D data in FFTW*/
float **cpdata;
float phase;
float wv;
float vmin;
float vmax;
float f1 = 1.0;
float f2 = 35.0;

complex *out;                           /*output of 1D FFT using FFTW*/
complex **datawx;                       /*data in frequency-wavenumber domain*/
complex *in2;
complex *out2;
complex *in3;
complex *out3;
complex **Pkv;                           /*wavefield in k-v*/
complex **Pxv;                           /*wavefield in x-v*/
complex **Pwx;                           /*wavefield in w-x*/

complex cshift;
complex tmp_a;
complex tmp_b;

fftwf_plan p1; 
fftwf_plan p2;
fftwf_plan p3;

/*allocate memory of reference velocities*/
vref = alloc2float(nz,nvref_max);     /*allocate 2D array for reference velocity to avoid changing memory of vector*/
                                      /*nz by nvref_max*/
for (ix=0;ix<nx;ix++){
for (iz=0;iz<nz;iz++){
image[ix][iz] = cmplx(0.0,0.0);       /*nz by nz*/
}
}


/*devide velocity by 2 for downward continuation*/
for (iz=0;iz<nz;iz++){
for (ix=0;ix<nx;ix++){
v[ix][iz]=v[ix][iz]/2.0;
}
}
 /*fprintf(stderr,"nx=%d nz=%d",nx,nz);
 FILE *Fvp = fopen("vel.bin", "wb");
 fwrite(v[0],1,4*nx*nz,Fvp);
 fclose(Fvp);*/


/*zero padding in termporal direction*/
ntfft = 1.0*exp2(ceil(log2(nt)));        /*number of zero padded trace in FFT*/
nw = ntfft/2+1;                          /*number of points of frequency axis after FFTW*/
cpdata = alloc2float(ntfft,nx);          /*data after zero padding*/  

for (ix=0;ix<nx;ix++){
 for (it=0;it<ntfft;it++){
   if(it<=nt) cpdata[ix][it] = data[ix][it];
   else {cpdata[ix][it] = 0.0;
         }
 		          }
	             }


/*allocate memory for w-x domain data*/
datawx = alloc2complex(nw,nx);           /*w-x data nx by nw*/
iw1 = floor(f1*dt*ntfft)+1;
iw2 = floor(f2*dt*ntfft)+1;

/*define plans for FFT using FFTW*/
/*plan 1 from t-x to w-x*/
in = alloc1float(ntfft);
out = alloc1complex(nw);
p1 = fftwf_plan_dft_r2c_1d(ntfft,in,(fftwf_complex*)out,FFTW_ESTIMATE);      /*real to complex*/

/*plan 2 from w-x to w-k*/
nxfft = 1.0*exp2(ceil(log2(nx)));
in2 =  alloc1complex(nxfft);
out2 = alloc1complex(nxfft);
p2 = fftwf_plan_dft_1d(nxfft,(fftwf_complex*)in2,(fftwf_complex*)out2,FFTW_FORWARD,FFTW_ESTIMATE);


/*plan 3 from w-k to w-x*/
in3 =  alloc1complex(nxfft);
out3 = alloc1complex(nxfft);
p3 = fftwf_plan_dft_1d(nxfft,(fftwf_complex*)in3,(fftwf_complex*)out3,FFTW_BACKWARD,FFTW_ESTIMATE);

Pxv = alloc2complex(nvref_max,nxfft);

Pkv = alloc2complex(nvref_max,nxfft);

Pwx = alloc2complex(nw,nxfft); 

/*apply first 1-D Fourier transform on data from t-x to w-x using FFTW package*/
for (ix=0;ix<nx;ix++){
for(it=0;it<ntfft;it++){
in[it] = cpdata[ix][it];                    /*assign one trace to a vector*/
}
fftwf_execute(p1);

for(iw=0;iw<nw;iw++){
datawx[ix][iw] = cdiv(out[iw], cmplx(sqrt(ntfft), 0.0));
}            /*w*/

}            /*x*/

fftwf_destroy_plan(p1);


/*determine frequency and wavenumber axis*/
dw = 2.0*PI/(ntfft*dt);                    /*frequency sampling interval*/
w0 = 0.0;                                  /*first frequency sample*/

dk = 2.0*PI/(nxfft*dx);                   /*wavenumber sampling interval*/
k0 = 0.0;                                 /*first wavenumber sample*/



/*initialization of downward wavefield*/
for (iw=0;iw<nw;iw++){
for (ix=0;ix<nxfft;ix++){
if (ix<nx){Pwx[ix][iw] =  datawx[ix][iw];}     
else{Pwx[ix][iw] = cmplx(0.0,0.0);}
}
}

/*loop over depth z*/
for (iz=0;iz<nz;iz++){                
fprintf(stderr,"depth sample %d\n",iz);

/*calculate reference velocities of each layer*/
vmin = v[0][iz];
vmax = v[0][iz]; 
for (ix=0;ix<nx;ix++){
if(v[ix][iz]>=vmax) vmax=v[ix][iz];       /*get the maximum velocity*/
if(v[ix][iz]<=vmin) vmin=v[ix][iz];       /*get the minimum velocity*/
}                  

dv = (vmax-vmin)/(nvref_max-1);

if(dv/vmax<=0.001){
nvref = 1;
vref[0][iz]=(vmin+vmax)/2;
}
else
{
nvref = nvref_max;
for (iv=0;iv<nvref_max;iv++)
{
vref[iv][iz] = vmin+dv*iv;
}
}


/*loop over frequencies*/
w = w0;
for (iw=iw1;iw<=iw2;iw++){
w = w0 + iw*dw;                               /*frequency axis (important)*/

/*apply phase-shift in w-x (optional)*/

/*datawx*/

/*Apply second FFT to tranform w-x data to w-k domain using FFTW*/
for (ix=0;ix<nxfft;ix++){
in2[ix] = Pwx[ix][iw];
}

fftwf_execute(p2);
for (ik=0;ik<nxfft;ik++){
out2[ik] = cdiv(out2[ik], cmplx(sqrt(nxfft), 0.0));
}

/*loop over wavenumbers*/
k = k0;
for (ik=0;ik<nxfft;ik++){

if (ik<=nxfft/2){
k = ik*dk;                           /*wavenumber axis (important)*/
}
else{
k = (ik-nxfft)*dk;
}
 
/*loop over reference velocities*/
for (iv=0;iv<nvref;iv++){
wv = w/vref[iv][iz];
if(wv>fabs(k)){                     /*note that k can be negative*/
phase = sqrt(wv*wv-k*k)*dz;
cshift = cmplx(cos(phase),sin(phase));
}
else{
cshift = cmplx(0.0,0.0);
}
Pkv[ik][iv] = cmul(out2[ik],cshift); 
}                               /*end for v*/

}                               /*end for k*/

   
/*from w-k go back to w-x domain*/
for (iv=0;iv<nvref;iv++){  /*inverse FFT for each velocity*/

for (ik=0;ik<nxfft;ik++){
in3[ik] = Pkv[ik][iv];
}    /*end for k*/

fftwf_execute(p3);


for (ix=0;ix<nxfft;ix++){
Pxv[ix][iv] = cdiv(out3[ix], cmplx(sqrt(nxfft), 0.0));
}    /*end for x*/

}    /*end for v*/     /*Pxv ix by iv*/


/*interpolation of wavefield in w-x*/
if (nvref==1){
for (ix=0;ix<nx;ix++){
Pwx[ix][iw] = Pxv[ix][0];
}
}
else
{
for (ix=0;ix<nx;ix++){ 
if (v[ix][iz]==vmax){i1=(v[ix][iz]-vmin)/dv-1;}
else
{i1 = (v[ix][iz]-vmin)/dv;}    /*find nearest reference velocity and wavefield*/
i2 = i1+1;
tmp_a = cadd(crmul(Pxv[ix][i1], vref[i2][iz]-v[ix][iz]) , crmul(Pxv[ix][i2], v[ix][iz]-vref[i1][iz]));
tmp_b = cmplx(vref[i2][iz]-vref[i1][iz], 0.0);
Pwx[ix][iw] = cdiv(tmp_a,tmp_b);
}          /*interpolate wavefield*/
}          /*end else*/

/*imaging condition*/
for (ix=0;ix<nx;ix++){
image[ix][iz] = cadd(image[ix][iz],Pwx[ix][iw]);
}

/*zero padding*/
for (ix=nx;ix<nxfft;ix++){
Pwx[ix][iw] = cmplx(0.0,0.0);
}

}       /*w*/
}       /*z*/

fftwf_destroy_plan(p2);
fftwf_destroy_plan(p3);

}   /*end pspimig migration function*/


static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}
