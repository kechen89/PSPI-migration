/* Copyright (c) Signal Analysis and Imaging Group (SAIG), University of Alberta, 2014.*/
/* All rights reserved.*/
/* suprepspimig:  $Date: 2014/04/07 */

#include "su.h"         /* include file for SU programs */ /*interface*/
#include "segy.h"       /* include file for SEGY traces */
#include "header.h"     /* include file for segy sizes */
#include <signal.h>     /* include library that defines signal handling functions */
#include <fftw3.h>      /* include FFTW library */
                        /* SU has its own complex number library*/
/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                            ",
" SUPREPSPIMIG - phase-shift plus interpolation migration (depth migration)  ",
"             for common-shot data                                           ",
"                                                                            ",
" suprepspimig <infile >outfile vfile= [optional parameters]                 ",
"                                                                            ",
" Required Parameters:	                                                     ",
" nshot:	number of shot gathers to be migrated			     ",
" nxi:          number of total horizontal samples of image                  ",
"               i.e. number of horizontal samples of velocity model          ",
" nz:		number of depth samples of image          	             ",
"               i.e. number of depth samples of velocity model	             ",
" dz:		depth sampling interval of image (velocity model)            ",
" vfile:	name of file containing velocities			     ",
"		                                                	     ",
" Optional parameters:                                                       ",
" dx (get from header)           horizontal sampling interval	             ", 
" freq(=25):          the central frequency of Ricker wavelet                ",
" f1(=1);             the frequency band to be migrated                      ",        
" f2(=35);            the frequency band to be migrated                      ",
" stau(=0):           the delay time of source                               ",
NULL};

/* 
 * Credits: 
 * Ke Chen, ke7@ualberta.ca
 * Last changes: April, 2014
 * Trace header fields accessed: ns,dt,d2
 */
/**************** end self doc ********************************/


/**************************prototypes for functions used internally****************************/
void prepspimig(complex **datasw, complex **datarw,float **image,float **v,int ntfft,int iw1,int iw2,int nxi,int nz,float dt,float dx, float dz);
float *ricker(float freq,float dt, int *ntw); /*synthesize ricker wavelet of central frequency fcent*/

/*segy trace*/
segy tr,tro;                             /*define the type of header variable, tr is a variable in segy struct*/ 
                                     /*trace identification header*/
int main(int argc, char **argv)      /*argc, argv - the arguments to the main() function*/
{ 
int nshot;                           /*number of shots to be migrated*/
int nt;                              /*number of time samples*/
int nz;			             /*number of depth samples of image*/
int nxi,nxii;                        /*number of horizontal samples of image*/
int nx;                              /*number of traces in each gather*/
int it;                              /*loop index for time*/
int iw;                              /*loop index for frequency*/
int ix;                              /*loop index for horizontal samples*/
int iz;                              /*loop index for depth samples*/
int is;                              /*loop index for shots*/
int ntfft;                           /*number of sample of FFT*/
int nw,iw1,iw2;                      /*number of frequency samples,starting and ending sample number of frequency*/
int ntw;                             /*number of time samples of ricker wavelet*/                            
int sindex;

float gx;
float sgxmin;
float *r=NULL;                       /*ricker wavelet*/     
float *in;                           /*input of FFTW*/
float dw;                            /*frequency sampling interval*/
float freq;                          /*center frequency of ricker wavelet*/
float f1;                            /*the starting frequency to be migrated*/
float f2;                            /*the ending frequency to be migrated*/               
float stau;                          /*delay time of source*/                      
float dt;                            /*temporal sampling interval*/                
float dx;                            /*spatial sampling interval*/
float dz;                            /*migrated depth sampling interval*/           
float sx;                            /*source location in x*/
float sxold;         
float **dr,**drr;                    /*receiver side wavefield*/
float **image;                      /*migrated image*/      
float **simage;                      /*stacked image*/
float **v;                           /*velocity model*/
float *str;                          /*source wavefield trace with ricker wavelet*/

complex *out;                        /*output of FFTW*/
complex *strw;
complex **dsw;                       /*source side wavefield*/
complex **drw;                       /*receiver side wavefield in w-x*/

char *vfile="";                      /*name of velocity file*/
                                     /*vfile is a pointer that points to a string*/
FILE *vfp;

fftwf_plan p1;


/*hook up getpar to handle the parameters*/
initargs(argc,argv);
requestdoc(1);

/*get required parameters*/
if (!getparint("nxi",&nxi)) err("number of horizontal samples of image nxi must be specified");
if (!getparint("nz",&nz)) err("number of vertical samples of image nz must be specified"); 
if (!getparfloat("dz",&dz)) err("dz must be specified");
if (!getparint("nshot",&nshot)) err("nshot must be specified");
if (!getparstring("vfile", &vfile)) err("velocity file must be specified");

/*get optional parameters*/
if (!getparfloat("freq",&freq)) freq = 25.0;          /*center frequency of ricker wavelet*/
if (!getparfloat("f1",&f1)) f1 = 1.0;                /*center frequency of ricker wavelet*/
if (!getparfloat("f2",&f2)) f2 = 35.0;                /*center frequency of ricker wavelet*/

if (!getparfloat("stau",&stau)) stau = 0.0;     /*delay time of source*/

/*get info from first trace*/
if (!gettr(&tr))  err("can't get first trace");   /*fgettr: get a fixed-length segy trace from a file by file pointer*/
nt = tr.ns;                         /*nt*/        /*gettr: macro using fgettr to get a trace from stdin*/

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

checkpars();

/*allocate memory*/
ntfft = 1.0*exp2(ceil(log2(nt)));          /*number of zero padded trace in FFT*/
nw = ntfft/2+1;                            /*number of points of frequency axis after FFTW*/
iw1 = floor(f1*dt*ntfft)+1;                /*the starting frequency sample to be migrated*/
iw2 = floor(f2*dt*ntfft)+1;                /*the end frequency sample to be migrated*/

image = alloc2float(nz,nxi);               /*image 2D array nxi by nz*/
simage = alloc2float(nz,nxi);              /*stacked image 2D array nxi by nz*/
v = alloc2float(nz,nxi);                   /*2D array, in Fortran the velocity model is nz by nx 2D array*/ 
nxii = nxi+1;                                             /*in binary, it is actually 1D*/ 
drr = alloc2float(ntfft,nxii);             /*auxiliary matrix for stripping receiver side wavefield*/
dr = alloc2float(ntfft,nxi);               /*receiver side wavefield in t-x*/
drw = alloc2complex(nw,nxi);               /*receiver side wavefield in w-x*/
dsw = alloc2complex(nw,nxi);               /*source side wavefield in w-x*/

str = alloc1float(ntfft);                  /*trace that constains ricker wavelet*/                                                       
strw = alloc1complex(nw);

/*load velicoty file*/
vfp=efopen(vfile,"r");	
efread(v[0],FSIZE,nz*nxi,vfp);             /*load velocity*/
efclose(vfp);			


/*make FFTW plan*/
dw = 2.0*PI/(ntfft*dt);                    /*frequency sampling interval*/
in = alloc1float(ntfft);
out = alloc1complex(nw);
p1 = fftwf_plan_dft_r2c_1d(ntfft,in,(fftwf_complex*)out,FFTW_ESTIMATE);      /*real to complex*/

/*************synthesize ricker wavelet**************************/
r=ricker(freq,dt,&ntw);                    /*ricker wavelet centered at zero*/

/*inject the ricker wavelet*/
for (it=0;it<ntfft;it++)
{str[it]=0.0;                              /*initialize trace*/
}

for (it=0;it<ntw;it++)
{str[it]=r[it];                            /*the ricker wavelet starts at 0*/
}

/*shift the ricker wavelet in frequency domain to make it peaks at stau*/
for(it=0;it<ntfft;it++){   
in[it] = str[it];                          /*assign one trace to a vector*/
}
fftwf_execute(p1);                         /*transform to frequency domain*/

for(iw=0;iw<nw;iw++)
strw[iw] = cdiv(out[iw],cmplx(sqrt(ntfft),0.0)); 
for(iw=0;iw<nw;iw++)
strw[iw] = cmul(strw[iw],cmplx(cos(stau*iw*dw),sin(stau*iw*dw)));        /*w(t)->w(t+stau) time shift -> phase shift*/


/*initialize image*/
for (iz=0;iz<nz;iz++){
for (ix=0;ix<nxi;ix++){
simage[ix][iz] = 0.0;
}
}

sx = (float) tr.sx;         /*get source coordinate of first trace*/
gx = (float) tr.gx;
sgxmin = MIN(sx,gx);

/*************************loop over shots*****************************/

for(is=0;is<nshot;is++){    /*loop over shots, there are nshot shots*/
sx = (float) tr.sx;         /*get source coordinate of first trace*/
sindex = (sx - sgxmin)/dx;
fprintf(stderr,"shot %d\n",sindex+1);
nx = 0;                     /*number of traces in each shot gather*/
sxold = sx;

/*synthesize source side wavefiled in w-x*/
for (ix=0;ix<nxi;ix++){
for (iw=0;iw<nw;iw++){
if (ix==sindex) dsw[ix][iw]=strw[iw];
else
dsw[ix][iw]=cmplx(0.0,0.0);
}
}

/*strip receiver side wavefiled in w-x*/

for(it=0;it<ntfft;it++){ 
for(ix=0;ix<nxi+1;ix++){
drr[ix][it]=0.0;
if(ix<nxi){dr[ix][it]=0.0;}
}
}

do {                        /*loop over traces in each shot*/
memcpy( (void *) drr[nx], (const void *) tr.data,nt*FSIZE);    /*get one trace from su data*/
sx = (float) tr.sx;
nx++;                      
}while(gettr(&tr) && sx==sxold);    /*if sx changes, the do loop ends. divide shots*/    


for(it=0;it<ntfft;it++){ 
for(ix=0;ix<nxi;ix++){
dr[ix][it]=drr[ix][it];
}
}

for (ix=0;ix<nxi;ix++){
for (iw=0;iw<nw;iw++){
drw[ix][iw]=cmplx(0.0,0.0);       /*initialize w-x receiver side wavefield*/
}
}

for (ix=0;ix<nxi;ix++){
for(it=0;it<ntfft;it++){
in[it] = dr[ix][it];                    
}
fftwf_execute(p1);                   /*transform to frequency domain*/
for(iw=0;iw<nw;iw++){
drw[ix][iw] = cdiv(out[iw],cmplx(sqrt(ntfft),0.0)); 
}        /*it*/
}        /*ix*/

/*migrate shot using PSPI*/
prepspimig(dsw,drw,image,v,nt,iw1,iw2,nxi,nz,dt,dx,dz);

/*stack the partial images*/
for (iz=0;iz<nz;iz++){
for (ix=0;ix<nxi;ix++){
simage[ix][iz] += image[ix][iz];
}
}

}  /*end the loop for shots*/

/* restore header fields and write output */
for (ix=0; ix<nxi; ix++) {
tro.ns = nz;
tro.d1 = dz;
memcpy( (void *) tro.data, (const void *) simage[ix],nz*FSIZE);
puttr(&tro);
}

fftwf_destroy_plan(p1);
fftwf_free(in); 
fftwf_free(out);

free1float(r);
free2complex(dsw);                       
free2float(dr);
free2float(drr);                   
free2complex(drw);                     
free2float(image);                        
free2float(simage);                     
free2float(v);                          
free1float(str);                       
free1complex(strw);

return(CWP_Exit());	
}

void prepspimig(complex **datasw,complex **datarw,float **image,float **v,int nt,int iw1,int iw2,int nxi,int nz,float dt,float dx, float dz)
{
int ntfft;
int nw;                                 /*number of frequency samples*/
int nxfft;                              /*number of horizontal samples after padding*/
int ix;                                 /*loop index over horizontal sample*/     
int iw;                                 /*loop index over frequency*/                 
int ik;                                 /*loop index over wavenumber*/
int iz;                                 /*loop index over migrated depth samples*/
int iv;                                 /*loop index over reference velocities*/
int nvref_max=2;                        /*number of reference velocities in each layer*/
int nvref;
int i1;                                 /*nearest reference velocity*/
int i2;

float **vref;                           /*2D reference velocity array*/
float w0;                               /*first frequency sample*/
float w;                                /*frequency*/
float dw;                               /*frequency sampling interval*/                 
float k0;                               /*first wavenumber*/
float k;
float dk;                               /*wave number sampling interval in x*/
float dv;                               /*velocity interval*/
float phase;
float wv;
float vmin;
float vmax;

complex tmp;
complex Scshift,Rcshift;                  /*phase shift*/
complex tmp_a;
complex tmp_b;

complex *in2;
complex *out2;
complex *in3;
complex *out3;
complex *Sout2;
complex *Rout2;
complex **SPkv;                           /*source wavefield in k-v*/
complex **SPxv;                           /*source wavefield in x-v*/
complex **SPwx;                           /*source wavefield in w-x*/
complex **RPkv;                           /*receiver wavefield in k-v*/
complex **RPxv;                           /*receiver wavefield in x-v*/
complex **RPwx;                           /*receiver wavefield in w-x*/

fftwf_plan p2;
fftwf_plan p3;

ntfft = 1.0*exp2(ceil(log2(nt)));
nw = ntfft/2+1; 
nxfft = 1.0*exp2(ceil(log2(nxi)));
Sout2 =  alloc1complex(nxfft);
Rout2 = alloc1complex(nxfft);

/*allocate memory of reference velocities*/
vref = alloc2float(nz,nvref_max);     /*allocate 2D array for reference velocity to avoid changing memory of vector*/
                                      /*nvref_max by nz*/

/*initialize image*/
for (ix=0;ix<nxi;ix++){
for (iz=0;iz<nz;iz++){
image[ix][iz] = 0.0;       /*nz by nz*/
}
}

/*plan 2 from w-x to w-k*/
in2 =  alloc1complex(nxfft);
out2 = alloc1complex(nxfft);
p2 = fftwf_plan_dft_1d(nxfft,(fftwf_complex*)in2,(fftwf_complex*)out2,FFTW_FORWARD,FFTW_ESTIMATE);

/*plan 3 from w-k to w-x*/
in3 =  alloc1complex(nxfft);
out3 = alloc1complex(nxfft);
p3 = fftwf_plan_dft_1d(nxfft,(fftwf_complex*)in3,(fftwf_complex*)out3,FFTW_BACKWARD,FFTW_ESTIMATE);

SPxv = alloc2complex(nvref_max,nxfft);       /*extrapolated source wavefield in x-v*/
RPxv = alloc2complex(nvref_max,nxfft);       /*extrapolated receiver wavefield in x-v*/

SPkv = alloc2complex(nvref_max,nxfft);       /*extrapolated source wavefield in k-v*/
RPkv = alloc2complex(nvref_max,nxfft);       /*extrapolated receiver wavefield in k-v*/

SPwx = alloc2complex(nw,nxfft); 
RPwx = alloc2complex(nw,nxfft); 

/*determine frequency and wavenumber axis*/
dw = 2.0*PI/(ntfft*dt);                    /*frequency sampling interval*/
w0 = 0.0;                                  /*first frequency sample*/

dk = 2.0*PI/(nxfft*dx);                    /*wavenumber sampling interval*/
k0 = 0.0;                                  /*first wavenumber sample*/

/*initialization of wavefield*/
for (iw=0;iw<nw;iw++){
for (ix=0;ix<nxfft;ix++){
if (ix<nxi){
SPwx[ix][iw] =  datasw[ix][iw];
RPwx[ix][iw] =  datarw[ix][iw];
}     
else{
SPwx[ix][iw] = cmplx(0.0,0.0);
RPwx[ix][iw] = cmplx(0.0,0.0);
}
}  /*ix*/
}  /*iw*/

/*loop over depth z*/
for (iz=0;iz<nz;iz++){                
/*fprintf(stderr,"depth sample %d\n",iz);*/

/*calculate reference velocities of each layer*/
vmin = v[0][iz];
vmax = v[0][iz]; 
for (ix=0;ix<nxi;ix++){
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

/*Apply second FFT to tranform w-x data to w-k domain using FFTW*/
for (ix=0;ix<nxfft;ix++){
in2[ix] = SPwx[ix][iw];
}
fftwf_execute(p2);
for (ik=0;ik<nxfft;ik++){
Sout2[ik] = cdiv(out2[ik], cmplx(sqrt(nxfft), 0.0));
}

for (ix=0;ix<nxfft;ix++){
in2[ix] = RPwx[ix][iw];
}
fftwf_execute(p2);
for (ik=0;ik<nxfft;ik++){
Rout2[ik] = cdiv(out2[ik], cmplx(sqrt(nxfft), 0.0));
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
Scshift = cmplx(cos(phase),-sin(phase));
Rcshift = cmplx(cos(phase),sin(phase));
}
else{
Scshift = cmplx(0.0,0.0);
Rcshift = cmplx(0.0,0.0);
}
SPkv[ik][iv] = cmul(Sout2[ik],Scshift); 
RPkv[ik][iv] = cmul(Rout2[ik],Rcshift); 
}                               /*end for v*/

}                               /*end for k*/
 
/*source from w-k to w-x domain*/
for (iv=0;iv<nvref;iv++){  /*inverse FFT for each velocity*/
for (ik=0;ik<nxfft;ik++){
in3[ik] = SPkv[ik][iv];
}    /*end for k*/
fftwf_execute(p3);
for (ix=0;ix<nxfft;ix++){
SPxv[ix][iv] = cdiv(out3[ix], cmplx(sqrt(nxfft), 0.0));
}    /*end for x*/
}    /*end for v*/     /*Pxv ix by iv*/

/*receiver from w-k to w-x domain*/
for (iv=0;iv<nvref;iv++){  /*inverse FFT for each velocity*/
for (ik=0;ik<nxfft;ik++){
in3[ik] = RPkv[ik][iv];
}    /*end for k*/
fftwf_execute(p3);
for (ix=0;ix<nxfft;ix++){
RPxv[ix][iv] = cdiv(out3[ix], cmplx(sqrt(nxfft), 0.0));
}    /*end for x*/
}    /*end for v*/     /*Pxv ix by iv*/

/*interpolation of wavefield in w-x*/
if (nvref==1){
for (ix=0;ix<nxi;ix++){
SPwx[ix][iw] = SPxv[ix][0];
RPwx[ix][iw] = RPxv[ix][0];
}
}
else
{
for (ix=0;ix<nxi;ix++){ 
if (v[ix][iz]==vmax){
i1=(v[ix][iz]-vmin)/dv-1;
}
else
{
i1 = (v[ix][iz]-vmin)/dv;
}    /*find nearest reference velocity and wavefield*/
i2 = i1+1;

/*interpolate source wavefield*/
tmp_a = cadd(crmul(SPxv[ix][i1], vref[i2][iz]-v[ix][iz]) , crmul(SPxv[ix][i2], v[ix][iz]-vref[i1][iz]));
tmp_b = cmplx(vref[i2][iz]-vref[i1][iz], 0.0);
SPwx[ix][iw] = cdiv(tmp_a,tmp_b);

/*interpolate receiver wavefield*/
tmp_a = cadd(crmul(RPxv[ix][i1], vref[i2][iz]-v[ix][iz]) , crmul(RPxv[ix][i2], v[ix][iz]-vref[i1][iz]));
tmp_b = cmplx(vref[i2][iz]-vref[i1][iz], 0.0);
RPwx[ix][iw] = cdiv(tmp_a,tmp_b);
}          /*interpolate wavefield*/
}          /*end else*/

/*cross-correlation imaging condition*/
for (ix=0;ix<nxi;ix++){
tmp = cmul(conjg(SPwx[ix][iw]),RPwx[ix][iw]);
image[ix][iz] += tmp.r;
}

/*zero padding*/
for (ix=nxi;ix<nxfft;ix++){
SPwx[ix][iw] = cmplx(0.0,0.0);
RPwx[ix][iw] = cmplx(0.0,0.0);
}

}       /*w*/
}       /*z*/

fftwf_destroy_plan(p2);
fftwf_free(in2);
fftwf_free(out2);
fftwf_destroy_plan(p3);
fftwf_free(in3);
fftwf_free(out3);

free2float(vref);
free1complex(Sout2);
free1complex(Rout2);
free2complex(SPkv);                          
free2complex(SPxv);                           
free2complex(SPwx);                           
free2complex(RPkv);                           
free2complex(RPxv);                          
free2complex(RPwx); 

}   /*end prepspimig migration function*/


/*****************************Synthesize Ricker Wavelet*************************************/
float *ricker(float freq,float dt,int *pntw)  /*synthesize ricker wavelet*/
{
int ncw;                          /*center sample*/
int itw;                          /*time sample index*/
float alpha;
float beta;
float *r;
ncw = floor(1.35*sqrt(6.0)/PI/freq/dt);     /*ncw is k*/
*pntw = 2*ncw+1;                     /*ntw is an odd number 2k+1*/
r = alloc1float(*pntw);              /*ricker wavelet*/

for (itw=0;itw<*pntw;itw++)
{
alpha = (itw-ncw)*freq*dt*PI;      /*-ncw:0:ncw*/
beta = alpha*alpha;
r[itw] = (1-2.0*beta)*exp(-beta);
}                /*end for loop*/

return r;
}                /*end ricker wavelet function*/
