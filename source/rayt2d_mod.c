/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* RAYT2D $Revision: 1.19 $; Date: 94/10/11 $	*/

#include "par.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" RAYT2D_mod - traveltime Tables calculated by 2D paraxial RAY tracing	",
" With shot positions at arbitrary depth	",
"									",
"     rayt2d vfile= tfile= [optional parameters]			",
"									",
" Required parameters:							",
" vfile=stdin		file containning velocity v[nx][nz]		",
" tfile=stdout		file containning traveltime tables		",
"			t[nxs][nxo][nzo]				",
"									",
" Optional parameters							",
" dt=0.008		time sample interval in ray tracing		",
" nt=401		number of time samples in ray tracing		",
"									",
" fz=0			first depth sample in velocity			",
" nz=101		number of depth samples in velocity		",
" dz=100		depth interval in velocity			",
" fx=0			first lateral sample in velocity		",
" nx=101		number of lateral samples in velocity		",
" dx=100		lateral interval in velocity			",
"									",
" fzo=fz		first depth sample in traveltime table		",
" nzo=nz		number of depth samples in traveltime table	",
" dzo=dz		depth interval in traveltime table		",
" fxo=fx		first lateral sample in traveltime table	",
" nxo=nx		number of lateral samples in traveltime table	",
" dxo=dx		lateral interval in traveltime table		",
"									",
" surf=\"0,0;99999,0\"  reflecting surface \"x1,z1;x2,z2;x3,z3;...\"	",
" fxs=fx		x-coordinate of first source			",
" fzs=0			z-coordinate of the source/sources			",
" nxs=1			number of sources				",
" dxs=2*dxo		x-coordinate increment of sources		",
" aperx=0.5*nx*dx  	ray tracing aperature in x-direction		",
"									",
" fa=-60		first take-off angle of rays (degrees)		",
" na=61			number of rays  				",
" da=2			increment of take-off angle  			",
" amin=0		minimum emergence angle 			",
" amax=90		maximum emergence angle 			",
"									",
" fac=0.01		factor to determine radius for extrapolation	",
" ek=1			flag of implementing eikonal in shadow zones 	",
" ms=10			print verbal information at every ms sources	",
" restart=n		job is restarted (=y yes; =n no)		",
" npv=0			flag of computing quantities for velocity analysis",
" if npv>0 specify the following three files				",
" pvfile=pvfile		input file of velocity variation pv[nxo][nzo]	",
" tvfile=tvfile		output file of traveltime variation tables  	",
"			tv[nxs][nxo][nzo]				",
" csfile=csfile		output file of cosine tables cs[nxs][nxo][nzo]	",
"									",
" Notes:								",
" 1. Each traveltime table is calculated by paraxial ray tracing; then 	",
"    traveltimes in shadow zones are compensated by solving eikonal	",
"    equation.								",
" 2. Input velocity is uniformly sampled and smooth one preferred.	",
" 3. Traveltime table and source ranges must be within velocity model.	",
" 4. Ray tracing aperature can be chosen as sum of migration aperature	",
"    plus half of maximum offset.					",
" 5. Memory requirement for this program is about			",
"      [nx*nz+4*mx*nz+3*nxo*nzo+na*(nx*nz+mx*nz+3*nxo*nzo)]*4 bytes	",
"    where mx = min(nx,2*(1+aperx/dx)).					",
"									",
" Note: spatial units of v(z,x) must be the same as those of dx. 	",
" v(z,x) is represented numerically in C-style binary floats v[xn][zn],	",
" where the depth direction is the fast direction in the data. Such     ",
" models can be created with unif2 or makevel.    			",  
"									",
NULL};
/*
 * Author:  Zhenyue Liu, 10/11/94,  Colorado School of Mines 
 *
 *          Trino Salinas, 01/01/96 included the option to handle nonflat
 *          reference surfaces.
 *          Subroutines from Dave Hale's modeling library were adapted in
 *          this code to define topography using cubic splines.
 *
 * References:
 *
 * Beydoun, W. B., and Keho, T. H., 1987, The paraxial ray method:
 *   Geophysics, vol. 52, 1639-1653.
 *
 * Cerveny, V., 1985, The application of ray tracing to the numerical
 *   modeling of seismic wavefields in complex structures, in Dohr, G.,
 *   ED., Seismic shear waves (part A: Theory): Geophysical Press,
 *   Number 15 in Handbook of Geophysical Exploration, 1-124.
 * 
 */
/**************** end self doc ********************************/

/* define initial value for traveltime as an arbitrary large number */
#define INITIAL_T 999999

/* Structures used internally */
/* one  ray */
typedef struct RayStruct {
	float t;		/* time */
	float px,pz;		/* slowness componets */
	float x,z;		/* x,z coordinates */
	float q2,p2;	/* Cerveny's dynamic ray tracing solution */
	float v,dvdx,dvdz;	/* velocity and its derivatives */
	float dtv;	/* traveltime variation */
} Ray;

/* ray tracing geometry parameters */
typedef struct RaytrStruct {
	float xs;		/* source lateral location */	 
	float zs;		/* source depth location */	 
	int nt;			/* number of time samples in ray tracing */
	float dt;		/* time sample interval in ray tracing  */	 
	int na;			/* number of rays */
	float fa;		/* first take-off angle of rays */
	float da;		/* increment of take-off angles */
	float amin;		/* minimum emergence angle  */	 
	float amax;		/* maximum emergence angle  */	 
	float fac;		/* factor to determine extrapolation radius*/ 
} Raytr;

/* 2D section geometry parameters */
typedef struct Geo2DStruct {
	int nx;			/* number of lateral samples  */
	float fx;		/* first lateral sample  */	 
	float dx;		/* lateral interval   */	 
	float odx;		/* 1 over dx   */	 
	int nz;			/* number of depth samples */
	float fz;		/* first depth sample  */	 
	float dz;		/* depth interval  */	
	float odz;		/* 1 over dz   */	 
}Geo2d;

/* Geometry of the recording surface */
typedef struct SurfaceSegmentStruct {
        float x;        /* x coordinate of segment midpoint */
        float z;        /* z coordinate of segment midpoint */
        float s;        /* x component of unit-normal-vector */
        float c;        /* z component of unit-normal-vector */
} SurfaceSegment;

typedef struct SurfaceStruct {
        int ns;               /* number of surface segments */
        float ds;             /* segment length */
        SurfaceSegment *ss;   /* array[ns] of surface segments */
} Surface;

typedef struct ReflectingBoundary {
	float x,z;			/* X,Z for ray to encounter */
	float t;            /* Time it took ray to reach from source */
	float angle;		/* Angle of incidence for ray */
	int flag;			/* Flag,0 initial,1 after hitting boundary,2 at source */       
} Boundary;

/* Prototypes of functions used internally */
void makeRay (Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz,
	Raytr raytr,float a0,int *nrs,Ray *ray,int npv,float *pv);

void raytime2d(Raytr raytr,Geo2d geo2dv,float *vt,Geo2d geo2dt,float *t,
	int npv,int ixs,int ir, int ian,int nr,float *vo,float *pv,float *tv,float *cs, Boundary ***bnd,int dir);

void voint(Geo2d geo2dv,float *v,Geo2d geo2dt,float *ov2,int npv,float *vo);

void dv2(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz);

void eiknl(Geo2d geo2dt,float *t,float *ov2,float tmax);

void trans(int nx,int nz,int nxt,int nx0,float *v,float *vt);

void vel2Interp(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz,
     float x,float z,
     float *u,float *ux,float *uz,float *uxx,float *uxz,float *uzz);

void vintp(Geo2d geo2dv,float *v,float x,float z,float *u);

void decodeSurfaces(int *nrPtr,int **nxzPtr,float ***xPtr,float ***zPtr);

int decodeSurface(char *string,int *nxzPtr,float **xPtr,float **zPtr);

void makesurf(float dsmax,int nr,int *nu,float **xu,float **zu,
             Surface **r);
void zcoorTopog(float fxs,float dxs,int nxs,Surface *srf,float *sz,
                float *nangl);

void decodeBoundaries(int *nrPtr,int **nxzPtr,float ***xPtr,float ***zPtr);

void make_boundaries(int nr,int nxs, int na,float fa, float da,int *nu,float **xu,float **zu,
        Boundary ****r);

float round2grid(float input, float factor);

int
main(int argc, char **argv)
{
	int	na,nat,nt,nxs,nxo,nzo,nx,nz,nxt,nx0,mx,npv,nsrf,*nxzsrf,dir;	
	float	dt,xs,fxs,fzs,dxs,exs,fxo,fzo,dxo,dzo,exo,ezo,fa,ea,amin,	
		amax,da,fat,fac,tmax,aperx,temp,fx,fz,dx,dz,ex,ez,fxt,ext,an;
	float	*v,*vt,*t,*t_down,*ov2,*pv=NULL,*pvt=NULL,*tv=NULL,*cs=NULL,*vo=NULL,
		**xsrf,**zsrf,*szi,*nangl;
	int	ixs,ixs0,nsize,isize,ek,ms,ir,ian,ixo,izo,i,cb=INITIAL_T;
	int ib=INITIAL_T,ia=INITIAL_T;

	Geo2d geo2dv, geo2dt;
	Raytr raytr;
    Surface *srf;
	Boundary ***bnd;
	char *vfile="stdin",*tfile="stdout",*jpfile,*pvfile,*tvfile,*csfile;
	char *restart;
	FILE *vfp, *tfp, *jpfp, *pvfp=NULL, *tvfp=NULL, *csfp=NULL;

	/*hook up getpar to handle the parameters	*/
	initargs(argc,argv);
	requestdoc(1);

	/* get velocity      information */
	if( !getparstring("vfile",&vfile)) {
		vfp = stdin;
	} else  
		vfp = fopen(vfile,"r");
	if(!getparint("nx",&nx)) nx = 101;
	if(!getparint("nz",&nz)) nz = 101;
	if(nx<3 || nz<3) err("nx and nz must be not less than 3!\n");
	if(!getparfloat("fx",&fx)) fx = 0.;
	if(!getparfloat("fz",&fz)) fz = 0.;
	if(!getparfloat("dx",&dx)) dx = 100.;
	if(!getparfloat("dz",&dz)) dz = 100.;
	ex = fx+(nx-1)*dx;
	ez = fz+(nz-1)*dz;
	geo2dv.nx = nx;		geo2dv.fx = fx;	
	geo2dv.dx = dx;		geo2dv.odx = 1.0/dx;
	geo2dv.nz = nz;		geo2dv.fz = fz;	
	geo2dv.dz = dz;		geo2dv.odz = 1.0/dz;

	if(!getparint("nxo",&nxo)) nxo = nx;
	if(!getparint("nzo",&nzo)) nzo = nz;
	if(nxo<3 || nzo<3) err("nxo and nzo must be not less than 3!\n");
	if(!getparfloat("fxo",&fxo)) fxo = fx;
	if(!getparfloat("fzo",&fzo)) fzo = fz;
	if(!getparfloat("dxo",&dxo)) dxo = dx;
	if(!getparfloat("dzo",&dzo)) dzo = dz;
	exo = fxo+(nxo-1)*dxo;
	ezo = fzo+(nzo-1)*dzo;
	geo2dt.nx = nxo;	geo2dt.fx = fxo;	
	geo2dt.dx = dxo;	geo2dt.odx = 1.0/dxo;
	geo2dt.nz = nzo;	geo2dt.fz = fzo;	
	geo2dt.dz = dzo;	geo2dt.odz = 1.0/dzo;

	if(!getparint("nxs",&nxs)) nxs = 1;
	if(!getparfloat("fxs",&fxs)) fxs = fx;
	if(!getparfloat("fzs",&fzs)) fzs = 0.; 		//New
	if(!getparfloat("dxs",&dxs)) dxs = 2*dxo;
	exs = fxs+(nxs-1)*dxs;
	if( !getparfloat("aperx",&aperx)) aperx = 0.5*(ex-fx);
	if(nxs>1) aperx += dxs;
	mx = 2.0*aperx/dx+2.99;
	if(mx>nx) mx = nx;

	if(!getparint("nt",&nt)) nt = 401;
	if(!getparfloat("dt",&dt)) dt = 0.008;
	tmax = (nt-1)*dt;
	if(!getparfloat("fa",&fa)) fa = -60.;
	if(!getparfloat("da",&da)) da = 2.;
	if(!getparint("na",&na)) na = 61;
	if(!getparfloat("amin",&amin)) amin = 0.;
	if(!getparfloat("amax",&amax)) amax = 90.;
	if(amax>180 || amin<0)  
	err("amin and amx must be within 0 to 180 degrees!\n");
	fa *= PI/180;
	da *= PI/180;
	ea = fa+(na-1)*da; //last angle
	amin *= PI/180;
	amax *= PI/180;
	if(!getparfloat("fac",&fac)) fac = 0.01;
	raytr.nt = nt;		raytr.dt = dt;
	raytr.da = da;		raytr.fac = fac;
	raytr.amin = amin;	raytr.amax = amax;

	if(!getparint("ek",&ek)) ek = 1;
	if( !getparstring("restart",&restart)) restart = "n";
	if(!getparint("ms",&ms)) ms = 10;
	if(!getparint("npv",&npv)) npv = 0;

	ixs0 = 0;
	isize = 0;
	nsize = nzo*nxo;
	if( !getparstring("tfile",&tfile)) {
		tfp = stdout;
	} else {
		if((tfp = fopen(tfile,"r"))!=NULL) {
			fclose(tfp);
			tfp = fopen(tfile,"r+");
		} else {
			tfp = fopen(tfile,"w");
		} 
		if(restart[0] =='y') {

			efseeko(tfp,(off_t) 0,SEEK_END);
			isize = (int) (eftello(tfp)/sizeof(float)/nsize);

			ixs0 = isize;
		}else {
			efclose(tfp);
			tfp = efopen(tfile,"w");
		} 
		efseeko(tfp,(off_t) (isize*nsize*sizeof(float)),SEEK_SET);
	}

	if(!getparstring("jpfile",&jpfile)) {
		jpfp = stderr;
	} else {
		jpfp = fopen(jpfile,"w");
	}
	
	if(npv){ 
	    if( !getparstring("pvfile",&pvfile)) 
		pvfile = "pvfile";
	    pvfp = fopen(pvfile,"r");
	    if( !getparstring("tvfile",&tvfile)) 
		tvfile = "tvfile";
	    tvfp = fopen(tvfile,"w");
	    if( !getparstring("csfile",&csfile)) 
		csfile = "csfile";
	    csfp = fopen(csfile,"w");
	}

        checkpars();

	fprintf(jpfp,"\n");
	fprintf(jpfp," RAYT2D parameters\n");
	fprintf(jpfp," ================\n");
	fprintf(jpfp," vfile=%s \n",vfile);
	fprintf(jpfp," tfile=%s \n",tfile);
	fprintf(jpfp," one-way time: dt=%g nt=%d \n",dt,nt);
	fprintf(jpfp," nz=%d fz=%g dz=%g\n",nz,fz,dz);
	fprintf(jpfp," nx=%d fx=%g dx=%g\n",nx,fx,dx);
	fprintf(jpfp," nzo=%d fzo=%g dzo=%g\n",nzo,fzo,dzo);
	fprintf(jpfp," nxo=%d fxo=%g dxo=%g\n",nxo,fxo,dxo);
 	fprintf(jpfp," nxs=%d fxs=%g dxs=%g\n",nxs,fxs,dxs);
	fprintf(jpfp," fzs=%g",fzs);
 	fprintf(jpfp," mx=%d aperx=%g \n",mx,aperx);
 	fprintf(jpfp," na=%d fa=%g da=%g\n",na,fa*180/PI,da*180/PI);
 	fprintf(jpfp," amin=%g amax=%g \n",amin*180/PI,amax*180/PI);
 	fprintf(jpfp," ek=%d ms=%d npv=%d restart=%s\n",ek,ms,npv,restart);
	if(npv)
 	  	fprintf(jpfp," pvfile=%s tvfile=%s csfile=%s\n",
			pvfile,tvfile,csfile);
	fflush(jpfp);

	/* ensure sources and output are in velocity grid	*/
	if(fx>fxs || ex<exs || fz>0) {
            warn("This condition must NOT be satisfied: fx>fxs || ex<exs || fz>0");
            warn("fx=%g fxs=%g ex=%g exs=%g fz=%g",
                    fx,fxo,ex,exo,fz);
		err("source lies outside of specified x grid  \n");
        }
	if(fx>fxo || ex<exo || fz>fzo || ez<ezo) {
            warn("This condition must NOT be satisfied: fx>fxo || ex<exo || fz>fzo || ez<ezo");
            warn("fx=%g fxo=%g ex=%g exo=%g fz=%g fzo=%g,ez=%g ezo=%g",
                    fx,fxo,ex,exo,fz,fzo,ez,ezo);
		err("output lies outside of specified x grid  \n");
        }

	/* allocate space */
	v = alloc1float(nx*nz);
	ov2 = alloc1float(nxo*nzo);
	t = alloc1float(nxo*nzo);
	t_down = alloc1float(nxo*nzo);

	if(npv) {
		pv = alloc1float(nx*nz);
		tv = alloc1float(nxo*nzo);
		vo = alloc1float(nxo*nzo);
		cs = alloc1float(nxo*nzo);
	}
	//printf("number %f\n",*t_up);
	/* read velocity */
	if(fread(v,sizeof(float),nx*nz,vfp)!=nx*nz)
	    err("cannot read %d velocities from file %s",nx*nz,vfp);

	if(npv) 
	    if(fread(pv,sizeof(float),nx*nz,pvfp)!=nx*nz)
	    	err("cannot read %d values from file %s",nx*nz,pvfp);

	fprintf(jpfp," finish velocity input\n");


	/*interpolate velocity and compute slowness squares  */
	voint(geo2dv,v,geo2dt,ov2,npv,vo);

	fprintf(jpfp," begin traveltime calculation\n");
	fflush(jpfp);

	/* Estimation of the reflecting surface 	*/
	decodeBoundaries(&nsrf,&nxzsrf,&xsrf,&zsrf);
	make_boundaries(nsrf,nxs,na,fa,da,nxzsrf,xsrf,zsrf,&bnd);

	/*Setting flag for down- or up-going rays */
	/* 0 for down, 1 for up going */
	dir = 0;

	/* Z coordinate of the source assignment */
        szi = alloc1float(nxs);
		for(ixs=0;ixs<nxs;ixs++)
		{
			szi[ixs] = fzs;
		}


	/* loop over sources */
	for(ixs=ixs0,xs=fxs+ixs0*dxs;ixs<nxs;ixs++,xs+=dxs){
		temp = fxo;
		if(xs<temp) temp = xs;
		if(xs-aperx>temp) temp = xs-aperx;
		nx0 = (temp-fx)/dx;
		fxt = fx+nx0*dx;
		temp = exo;
		if(xs>temp) temp = xs;
		if(xs+aperx<temp) temp = xs+aperx;
		nxt = 1+(int)((temp-fx)/dx+0.99)-nx0;
		ext = fxt+(nxt-1)*dx;

		vt = alloc1float(nxt*nz);
		if(npv) pvt = alloc1float(nxt*nz);

	    /* reducing velocity volume due to aperture	*/
		trans(nx,nz,nxt,nx0,v,vt);
		if(npv) trans(nx,nz,nxt,nx0,pv,pvt);

		/* determine range of take-off angles	*/
		/*Heavily reduced, we will shot in all directions*/
		fat = fa ;
		nat = na ;


		/* update geometry information	*/
		raytr.xs = xs;		raytr.zs = szi[ixs];
		raytr.na = nat;		raytr.fa = fat;		
		geo2dv.nx = nxt;	geo2dv.fx = fxt;
		
	    /* compute traveltime by paraxial ray tracing	*/
		raytime2d(raytr,geo2dv,vt,geo2dt,t,npv,ixs,ib,ia,nsrf,vo,pvt,tv,cs,bnd,dir);

		/* Let's make a copy of a downgoing traveltimes, nice to have around */
		for(ixo=0; ixo<nxo; ++ixo) 
		for(izo=0; izo<nzo; ++izo){
			i = izo+ixo*nzo;
			t_down[i] = t[i]; //t_down is copy
		};

		/* Loop through all the boundary encounters and do raytracing for upgoing waves */
		for (ir=0; ir<nsrf; ++ir){
			//fprintf(jpfp,"Layer number=%d\n",ir+1);
			for(ian=0;ian<na;++ian){

				//fprintf(jpfp,"Ray number=%d\n",ian+1);

				/* Switch geometry for upgoing ray	*/
				/* If there was crossing of ray */
				if (bnd[ir][ixs][ian].flag == 1){

					raytr.xs = bnd[ir][ixs][ian].x;		
					raytr.zs = bnd[ir][ixs][ian].z;
					raytr.na = 1;						
					 if (bnd[ir][ixs][ian].angle > 0) raytr.fa =  (PI / 2) + bnd[ir][ixs][ian].angle;
					 if	(bnd[ir][ixs][ian].angle < 0) raytr.fa = -(PI / 2) + bnd[ir][ixs][ian].angle;

					dir = 1;
					
					//I need to put rayt2d here, but hide a lot behind flag check
					raytime2d(raytr,geo2dv,vt,geo2dt,t,npv,ixs,ir,ian,nsrf,vo,pvt,tv,cs,bnd,dir);

					dir = 0;
				};	
			};	
		};
		

	
		/*make up in shadow zones by eikonal equation	*/
		if(ek) eiknl(geo2dt,t_down,ov2,tmax);
	
		/*write traveltime, downgoing field	*/
		fwrite(t_down,sizeof(float),nxo*nzo,tfp);

		/*write quantities for velocity analysis	*/
		if(npv) {
			fwrite(tv,sizeof(float),nxo*nzo,tvfp);
			fwrite(cs,sizeof(float),nxo*nzo,csfp);
		}

		if(ixs%ms==0) {
		  fprintf(jpfp," traveltime computed at source xs=%g\n",xs);
		  fflush(jpfp);
		}
		free1float(vt);
		if(npv) free1float(pvt);
	}
 	fprintf(jpfp," finish program rayt2d\n\n");
	
	fclose(tfp);
	fclose(vfp);
	fclose(jpfp);
	
	free1float(v);
	free1float(t);
	free1float(t_down);

	if(npv){
		free1float(pv);
		free1float(tv);
		free1float(cs);
	}

	/*Freeing boundary*/
	free1(bnd);

	return(CWP_Exit());
}

void decodeBoundaries(int *nrPtr,int **nxzPtr,float ***xPtr,float ***zPtr)
/*************************************************************************
decodeBoundaries - parse boundary
**************************************************************************
Output:
nrPtr           pointer to nr an int specifying number of surfaces = 2
nxzPtr          pointer to nxz specifying number of (x,z) pairs defining the
                surfaces
xPtr            pointer to array[x][nr] of x values for each surface
zPtr            array[z][nr] of z values for each surface
**************************************************************************/
{
        int nr,*nxz,ir;
        float **x,**z;
		char t[6144],*s;

        /* count surfaces */
        nr = countparname("refl");
        if (nr==0) return;

        /* allocate space */
        nxz = ealloc1(nr,sizeof(int));
        x = ealloc1(nr,sizeof(float*));
        z = ealloc1(nr,sizeof(float*));

        /* get surfaces */
        for (ir=0; ir<nr; ++ir) {
                if (!getnparstring(ir+1,"refl",&s)) s = "0,0;99999,0";
                strcpy(t,s);
                if (!decodeSurface(t,&nxz[ir],&x[ir],&z[ir]))
                        err("Surface number %d specified "
                                "incorrectly!\n",ir+1);
        }

        /* set output parameters before returning */
        *nrPtr = nr;
        *nxzPtr = nxz;
        *xPtr = x;
        *zPtr = z;
}

/*****************************************************************************
Make array of boundaries holding parameters for each boundary,source,angle
******************************************************************************
Input:
nr              number of surfaces = 1
nxs				number of sources = 1
na				number of rays	= 61
nu              array[nr] of numbers of (x,z) pairs; u = 0, 1, ..., nu[ir]
xu              array[nr][nu[ir]] of surface x coordinates x(u)
zu              array[nr][nu[ir]] of surface z coordinates z(u)

Output:
r               array[boundary_no][source_no][ray_no] of boundaries
******************************************************************************/
void make_boundaries(int nr,int nxs,int na,float fa,float da,int *nu,float **xu,float **zu,
        Boundary ****r)
{
	int ir,is,ian;
	float fa0;
	Boundary ***rr;

	*r = rr = (Boundary ***)(ealloc3(na,nxs,nr,sizeof(Boundary)));

	for(ian=0,fa0=fa;ian<na;++ian,fa0+=da){
		for(is=0;is<nxs;++is){
			for (ir=0; ir<nr; ++ir){	
				rr[ir][is][ian].z=zu[ir][0];
				rr[ir][is][ian].x=INITIAL_T;
				rr[ir][is][ian].angle=fa0;
				rr[ir][is][ian].flag=0;
			};
		};
	};

}



int decodeSurface (char *string,int *nxzPtr, float **xPtr, float **zPtr)
/**************************************************************************
decodeSurface - parse one surface specification
***************************************************************************
Input:
string          string representing surface

Output:
nxzPtr          pointer to number of x,z pairs
xPtr            array of x values for one surface
zPtr            array of z values for one surface
**************************************************************************/
{
        int nxz,ixz;
        float *x,*z;
        char *s,*t;

        s = string;
        s = strtok(s,",;\0");
        /* count x and z values, while splitting string into tokens */
        for (t=s,nxz=0; t!=NULL; ++nxz)
                t = strtok(NULL,",;\0");

        /* total number of values must be even */
        if (nxz%2) return 0;

        /* number of (x,z) pairs */
        nxz /= 2;

        /* 2 or more (x,z) pairs are required */
        if (nxz<2) return 0;

        /* allocate space */
        x = ealloc1(nxz,sizeof(float));
        z = ealloc1(nxz,sizeof(float));

        /* convert (x,z) values */
        for (ixz=0; ixz<nxz; ++ixz) {
                x[ixz] = atof(s);
                s += strlen(s)+1;
                z[ixz] = atof(s);
                s += strlen(s)+1;
        }

        /* set output parameters before returning */
        *nxzPtr = nxz;
        *xPtr = x;
        *zPtr = z;
        return 1;

}

/* Prototypes of functions used interally */
 void dfrungk(float v,float vx,float vz,float vxx,float vxz,float vzz,
     	float x,float z,float px,float pz,float p2,float q2,
      	float *dx,float *dz,float *dpx,float *dpz,float *dp2,float *dq2);
 void sum2(float h,float a1,float a2,float a3,float a4,float a5,float a6,
	float b1,float b2,float b3,float b4,float b5,float b6,
	float *c1,float *c2,float *c3,float *c4,float *c5,float *c6);

void makeRay(Geo2d geo2dv,float *vel,float *vxx,float *vxz,float *vzz,
    Raytr raytr,float a0,int *nrs,Ray *ray,int npv,float *pv) 
/*****************************************************************************
 Trace rays in uniformly sampled v(nz,nx)
*****************************************************************************
Input:
geo2dv		grid parameters of velocity 
v2		sampled velocity array
raytr		geometry parameters of ray tracing 
vxx,vxz,vzz	sampled second derivatives of velocity  
a0  		shooting angles at source point
npv		flag of velocity analysis
pv		velocity variation if npv>0

Output:
ray		pointer to ray parameters sampled at discrete ray steps
nrs		number of points on ray 
*****************************************************************************
Note:
The ray ends when it runs out of time or with the first step that is out
of the velocity boudary or out of emergence angle boundary
*****************************************************************************
Author: Zhenyue Liu, CSM, 1995.
****************************************************************************/ 
{
	int nx=geo2dv.nx,nz=geo2dv.nz,nt=raytr.nt;
	float fx=geo2dv.fx,fz=geo2dv.fz;
	float dt=raytr.dt,amin=raytr.amin,amax=raytr.amax;

	int it;
	float tzmin,tzmax,ov,lx,lz,h,h2,h6;

	float v,dvdx,dvdz,uxx,uxz,uzz,tzt,dtv=0.0;
	float x,z,px,pz,p2,q2,
	      dx,dz,dpx,dpz,dp2,dq2;
	float xt,zt,pxt,pzt,p2t,q2t,
	    dxt,dzt,dpxt,dpzt,dp2t,dq2t;


	/* velocity and derivatives at source	*/
	x = raytr.xs;
	z = raytr.zs;
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,x,z,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);
	ov = 1.0/v;

	p2 = 1.0;
	q2 = 0.0;
	*nrs = nt;

/*	compute slowness at source	*/
	px = ov*sin(a0);
	pz = ov*cos(a0);

/*	first ray step output	*/
	tzt = v*pz;
	ray[0].x = x;
	ray[0].z = z;
	ray[0].px = px;
	ray[0].pz = pz;
	ray[0].q2 = q2;
	ray[0].p2 = p2;
	ray[0].v = v;
	ray[0].dvdx = dvdx;
	ray[0].dvdz = dvdz;
	if(npv) ray[0].dtv = 0;

/*	compute minimum and maximum z-component of unit ray vector	*/
	tzmin = cos(amax)-0.01;
	tzmax = cos(amin)+0.01;

/*	determine fraction steps for Rung-Kuta	*/
	h = dt;		h2 = h*0.5;	h6 = h/6.0;

	lx = fx+(nx-1)*geo2dv.dx;
	lz = fz+(nz-1)*geo2dv.dz;

/*	loop over time steps	*/
	for(it=1; it<nt; ++it){
		if(x<fx || x>lx || z<fz || z>lz ||
      		   tzt<tzmin || tzt>tzmax) {
			it -= 1; 
			break;
		}

/*	step 1 of 4th-order Rung-Kuta	*/

/*	  compute k1 = f(y0)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		x,z,px,pz,p2,q2,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
          
/*	  compute y1 = y0+0.5*h*k1	*/
	    sum2(h2,x,z,px,pz,p2,q2,
      		dx,dz,dpx,dpz,dp2,dq2,
      		&xt,&zt,&pxt,&pzt,&p2t,&q2t);

/*	velocity interpolation	*/
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  xt,zt,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	step 2 of 4th-order Rung-Kuta	*/

/*	  compute k2 = f(y1)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		xt,zt,pxt,pzt,p2t,q2t,
      		&dxt,&dzt,&dpxt,&dpzt,&dp2t,&dq2t);

/*	  compute y2 = y0+0.5*h*k2	*/
	    sum2(h2,x,z,px,pz,p2,q2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&xt,&zt,&pxt,&pzt,&p2t,&q2t);

/*	  compute k = k1+2.0*k2	*/
	    sum2(2.0,dx,dz,dpx,dpz,dp2,dq2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
 
/*	velocity interpolation	*/
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  xt,zt,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	step 3 of 4th-order Rung-Kuta	*/

/*	  compute k3 = f(y2)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		xt,zt,pxt,pzt,p2t,q2t,
      		&dxt,&dzt,&dpxt,&dpzt,&dp2t,&dq2t);

/*	  compute y3 = y0+h*k2	*/
 	    sum2(h,x,z,px,pz,p2,q2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&xt,&zt,&pxt,&pzt,&p2t,&q2t);

/*	  compute k = k1+2.0*k2+2.0*k3	*/
	    sum2(2.0,dx,dz,dpx,dpz,dp2,dq2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
 
/*	velocity interpolation	*/
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  xt,zt,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	step 4 of 4th-order Rung-Kuta	*/

/*	  compute k4 = f(y3)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		xt,zt,pxt,pzt,p2t,q2t,
      		&dxt,&dzt,&dpxt,&dpzt,&dp2t,&dq2t);

/*	  compute k = k1+2.0*k2+2.0*k3+k4	*/
	    sum2(1.0,dx,dz,dpx,dpz,dp2,dq2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
 
/*	  compute y4 = y0+h*k/6	*/
 	    sum2(h6,x,z,px,pz,p2,q2,
      		dx,dz,dpx,dpz,dp2,dq2,
      		&x,&z,&px,&pz,&p2,&q2);
         
/*	velocity interpolation	*/
	    vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  	x,z,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	velocity itself interpolation	*/
	if(npv)	vintp(geo2dv,pv,x,z,&dtv);

		tzt = v*pz;

/*	save ray parameters */
		ray[it].x = x;
		ray[it].z = z;
		ray[it].px = px;
		ray[it].pz = pz;
		ray[it].p2 = p2;
		ray[it].q2 = q2;
		ray[it].v = v;
		ray[it].dvdx = dvdx;
		ray[it].dvdz = dvdz;
		if(npv)  ray[it].dtv = dtv/v;
	}
	*nrs = it; 

}

/***********************************************************************/
  void dfrungk(float v,float vx,float vz,float vxx,float vxz,float vzz,
     	float x,float z,float px,float pz,float p2,float q2,
      	float *dx,float *dz,float *dpx,float *dpz,float *dp2,float *dq2)
{

	float ov,vv,c,s,v11;

		ov = 0.0*z+ 0.0*x + 1.0/v;
		vv = v*v;
		*dx = vv*px;
		*dz = vv*pz;
		*dpx = -ov*vx;
		*dpz = -ov*vz;

		c = v*pz;
		s = v*px;
		v11 = vxx*c*c-2.0*vxz*c*s+vzz*s*s;
		*dp2 = -ov*v11*q2;
		*dq2 = vv*p2;

}


/************************ c = a+h*b ***********************************/
 void sum2(float h,float a1,float a2,float a3,float a4,float a5,float a6,
	float b1,float b2,float b3,float b4,float b5,float b6,
	float *c1,float *c2,float *c3,float *c4,float *c5,float *c6)
{

	*c1 = a1+h*b1;
	*c2 = a2+h*b2;
	*c3 = a3+h*b3;
	*c4 = a4+h*b4;
	*c5 = a5+h*b5;
	*c6 = a6+h*b6;

}

/* Prototype of function used internally */
void ddt(float p,float q,float c,float s, float *dv,float v,float *d2t,
	float *cuv);

void raytime2d(Raytr raytr,Geo2d geo2dv,float *v,Geo2d geo2dt,float *t,
	int npv,int ixs,int ir, int ian,int nr,float *vo,float *pv,float *tv,float *cs, Boundary ***bnd,int dir) 
/****************************************************************************
 raytime2d - calculate traveltimes by paraxial ray tracing
*****************************************************************************
Input:
geo2dt		grid parameters of traveltime 
geo2dv		grid parameters of velocity 
v2		sampled velocity array
raytr		geometry parameters of ray traycing 

Output:
t		traveltime 
*****************************************************************************
Note: when traveltime has multiple values, the one has the shortest 
 ray path is chosen.
*****************************************************************************
Author: Zhenyue Liu, CSM 1995.
****************************************************************************/ 
{
	int na=raytr.na,nt=raytr.nt,nxo=geo2dt.nx,nzo=geo2dt.nz,
		nx=geo2dv.nx,nz=geo2dv.nz;
	float fxo=geo2dt.fx,fzo=geo2dt.fz,dxo=geo2dt.dx,dzo=geo2dt.dz,
		odxo=geo2dt.odx,odzo=geo2dt.odz;
	float fac=raytr.fac,dt=raytr.dt,fa=raytr.fa,da=raytr.da,xs=raytr.xs;
	float *vxx,*vxz,*vzz,*s,zs=raytr.zs;
	float closest_x =INITIAL_T,closest_z=INITIAL_T,closest_t=INITIAL_T;
	int nrs;
	Ray *ray;

 
 	int i,ia,ixo,izo,ib;
	float a,xo,zo,exo,ezo;

/*	variables used for extrapolation	*/
	int it,jt,nxf,nxe,nzf,nze;
	float tc,xoc,zoc,xc,zc,r1,v0=0.0,norm2,terr,vd1,cuv,
      		sd,sind,cosd,vd,pxd,pzd,r2,t1,t2,r1x,r2x,t1x,t2x,t2xz,
      		p2d,q2d,gradv[2],d2t[3],dtvd,odt=1.0/dt,xcosd=0.0,*tvd;



	vxx = alloc1float(nx*nz);
	vxz = alloc1float(nx*nz);
	vzz = alloc1float(nx*nz);
	tvd = alloc1float(nt);
	s = alloc1float(nxo*nzo);
	ray = (Ray*)alloc1(nt,sizeof(Ray)); 

	/* compute second derivatives of velocity	*/
 	dv2(geo2dv,v,vxx,vxz,vzz); 


/*	maximum range of output points	*/
	exo = fxo+(nxo-1)*dxo;
	ezo = fzo+(nzo-1)*dzo;

/*	initialize traveltime and raypath length	*/
	for(ixo=0; ixo<nxo; ++ixo) 
		for(izo=0; izo<nzo; ++izo){
			i = izo+ixo*nzo;
			t[i] = INITIAL_T;
			s[i] = INITIAL_T;
		     if(npv) tv[i] = cs[i] = 0.0;
	}


/* 	loop over shooting angles at source point 	*/
	for(ia=0,a=fa; ia<na; ++ia,a+=da){  
	
	/* trace rays	*/
		makeRay(geo2dv,v,vxx,vxz,vzz,raytr,a,&nrs,ray,npv,pv);

/*		extropolate to grids near central rays	*/
		    v0 = ray[0].v;
		    r2 = 0.0;
		    xc = raytr.xs;
		    zc = raytr.zs;
		    sd = 0;
		    vd1 = v0;
		    if(npv) {
		      dtvd = 0.0;
		      for(it=1; it<nrs; ++it) {
			dtvd += 0.5*dt*(ray[it].dtv+ray[it-1].dtv);
			tvd[it] = ray[it].v*dtvd;
		      }
		    }
		    for (it=1; it<nrs; ++it) {
			xo = ray[it].x;
			zo = ray[it].z;
			vd = ray[it].v;
			sd = sd+0.5*dt*(vd+vd1);
			vd1 = vd;

/*		    	if the point is within the previous circle	*/
			if(r2 > (xc-xo)*(xc-xo)+(zc-zo)*(zc-zo)
      			  && it != nrs-1) continue;

/*		    	if the point is out of output range	*/
			if(xo<fxo || xo>exo || zo<fzo || zo>ezo) continue;
								
			xc = xo;
			zc = zo;
			q2d = ray[it].q2; 

/*			if caustics	*/
			if(q2d == 0.0) {
				r2 = 0.0;
				continue;
			}

			vd = ray[it].v;
			pxd = ray[it].px;
			pzd = ray[it].pz;

			p2d = ray[it].p2;
			sind = vd*pxd;
			cosd = vd*pzd;
			gradv[0] = ray[it].dvdx;
			gradv[1] = ray[it].dvdz;

/*			calculate second derivatives of traveltime*/
			ddt(p2d,q2d,cosd,sind,gradv,vd,d2t,&cuv);

/*			determine radius for extrapolation	*/
			tc = it*dt;
			terr = tc*fac;
			norm2 = sqrt(d2t[0]*d2t[0]+d2t[2]*d2t[2]+
      				2.0*d2t[1]*d2t[1]);

			r2 = terr/norm2;
			r1 = sqrt(r2);
/*			keep ray-centered coordinate system regular */
			if(r1*cuv > 0.1) r1 = 0.1/cuv;

/*			radius cannot be too large	*/
			if(r1 > 0.1*sd) r1 = 0.1*sd;

			r2 = r1*r1;

			nxf = (xc-r1-fxo)*odxo+0.9;
			if(nxf < 0) nxf = 0;
			nxe = (xc+r1-fxo)*odxo+0.1;
			if(nxe >= nxo) nxe = nxo-1;
/*       	fprintf(stderr,"x,z,t,r=%g %g %g %g\n",
			xc,zc,tc,r1);*/	 

			for(ixo=nxf; ixo<=nxe; ++ixo){
				xoc = fxo-xc+ixo*dxo;
				r2x = r2-xoc*xoc;
				if(r2x < 0) continue;

				r1x = sqrt(r2x);
				t1x = tc+xoc*pxd;
				t2x = tc*xoc*xoc*d2t[0];
				t2xz = 2.0*xoc*d2t[1];
				if(npv) xcosd = pzd+xoc*d2t[1];

				nzf = (zc-r1x-fzo)*odzo+0.9;
				if(nzf < 0) nzf = 0;
				nze = (zc+r1x-fzo)*odzo+0.1;
				if(nze>=nzo) nze = nzo-1;

				i = ixo*nzo;
				for(izo=nzf; izo<=nze; ++izo){ 
/*				    if ray path is shorter	*/
				    if(sd < s[i+izo] ) {
					zoc = fzo-zc+izo*dzo;
					t1 = t1x+zoc*pzd;
					t2 = t2x+zoc*(zoc*d2t[2]+t2xz);
					s[i+izo] = sd;
					t[i+izo] = t1*t1+tc*t2;
				      if(npv) {
					jt = (t1+0.5*t2)*odt+0.5;
					if(jt<0) jt = 0;
					if(jt>nrs-1) jt = nrs-1;
					tv[i+izo] = tvd[jt];
					cs[i+izo] = xcosd+d2t[2]*zoc;
				      }
				    }
 				}
 			}
 		}
		

		/* Reflection magic */

		switch (dir){

			/*Downgoing rays */
			case 0:

				/*Checking if boundary is crossed by ray */
				for(ib=0; ib<nr; ++ib){
					for(it=0; it<nrs; ++it){

						/* Rounding to nearest grid point */
						if(round2grid(ray[it].z,dzo) == round2grid(bnd[ib][ixs][ia].z,dzo)){

							/* If first loop, first values are closest to z */
							if(bnd[ib][ixs][ia].x == INITIAL_T){

									bnd[ib][ixs][ia].x = ray[it].x;
									bnd[ib][ixs][ia].z = ray[it].z;
							}

							/* Checks if the current values are closer to boundary z than the ones from previous loops */
							if(abs(ray[it].z - round2grid(ray[it].z,dzo)) < abs(bnd[ib][ixs][ia].z - round2grid(ray[it].z,dzo))){

									bnd[ib][ixs][ia].x = ray[it].x;
									bnd[ib][ixs][ia].z = ray[it].z;
	
							}
						}		
					}

					/*Rounding closest values, changing flag for encountering reflector */
					/* But only if boundary was crossed */
					if(bnd[ib][ixs][ia].x != INITIAL_T){
						bnd[ib][ixs][ia].x = round2grid(bnd[ib][ixs][ia].x,dxo);
						bnd[ib][ixs][ia].z = round2grid(bnd[ib][ixs][ia].z,dzo);
						bnd[ib][ixs][ia].flag = 1;
					}

					/*Loop for finding time at the x,z of reflection */
					for (ixo=0; ixo<nxo; ++ixo){
						xo = fxo+ixo*dxo;
						i = ixo*nzo;
						for (izo=0; izo<nzo; ++izo){
							zo = fzo+izo*dzo;
							if(bnd[ib][ixs][ia].x == xo && bnd[ib][ixs][ia].z == zo){

								/*Unfortunately with repeating part that happens later for entire timefield */ 
								bnd[ib][ixs][ia].t = t[i+izo];
			  					bnd[ib][ixs][ia].t = MAX(0.0,bnd[ib][ixs][ia].t);
								if(bnd[ib][ixs][ia].t < 999) bnd[ib][ixs][ia].t = sqrt(bnd[ib][ixs][ia].t);
								
							}
 						}
					}
				}


				break;

			/* Upgoing rays */
			case 1:

				for(it=0; it<nrs; ++it){
					if(round2grid(ray[it].z,dzo) == fzo){
						
						/* If this is 1st time it reached surf */
						if (closest_z == INITIAL_T){

							/* Save all parameters */
							closest_z = ray[it].z;
							closest_x = ray[it].x;

						}
						else {
							/* If there is already time stored for surface */
							/* Check if this one is not closer to it and if so, replace */
							if(abs(ray[it].z - fzo) < abs(closest_z - fzo))
								closest_z = ray[it].z;
								closest_x = ray[it].x;
						}
					}

				
				}

				/* Overwrite downgoing times with sum of down- and up-going */
				/*If we already found surface hit */
				if (closest_z != INITIAL_T){
					/*X and Z of the ray hit */
					bnd[ir][ixs][ian].x =  round2grid(closest_x,dxo);
					bnd[ir][ixs][ian].z =  round2grid(closest_z,dzo);
					/*Loop through all time field to get t value for x,z */
					for (ixo=0; ixo<nxo; ++ixo){
						xo = fxo+ixo*dxo;
						i = ixo*nzo;
						for (izo=0; izo<nzo; ++izo){
							zo = fzo+izo*dzo;
							if(bnd[ir][ixs][ian].x == xo && bnd[ir][ixs][ian].z == zo){
								/* Restrictions from orignal code, repeated i am afraid */
								if(t[i+izo] > 0.0 && t[i+izo] < 999){

									//printf("Time for down %f\n",bnd[ir][ixs][ian].t);
									//printf("Time for up %f\n",sqrt(t[i+izo]));
									/*Overwrite of down time with 2WT */
									bnd[ir][ixs][ian].t += sqrt(t[i+izo]); //!
								}
							}
						}
					}

					//printf("Final Z %f\n",bnd[ir][ixs][ian].z);
					//printf("Final X %f\n",bnd[ir][ixs][ian].x);
					//printf("Final Time %f\n",bnd[ir][ixs][ian].t);	

				}
				break;
			
			default:
				printf("Uh oh!\n");
		}
	}

/*	square root of traveltimes */
	for (ixo=0; ixo<nxo; ++ixo){ 
		i = ixo*nzo;
		for (izo=0; izo<nzo; ++izo){  
			t[i+izo] = MAX(0.0,t[i+izo]);
			if(t[i+izo] < 999) t[i+izo] = sqrt(t[i+izo]);
  			if(npv) cs[i+izo] *= vo[i+izo];
		}
	}

/*	compute traveltime near source	*/
	nxf = (xs-fxo)/dxo-1.5;
	if(nxf<0) nxf = 0;
	nxe = (xs-fxo)/dxo+2.5;
	if(nxe>nxo-1) nxe = nxo-1;
	nzf = (zs-fzo)/dzo-0.5;
	if(nzf<0) nzf = 0;
	nze = (zs-fzo)/dzo+1.5;
	if(nze>=nzo) nze = nzo-1;
	for (ixo=nxf; ixo<=nxe; ++ixo){
		xo = fxo+ixo*dxo;
		i = ixo*nzo;
		for (izo=nzf; izo<=nze; ++izo){
			zo = fzo+izo*dzo;
			if(t[i+izo] > 999) 
			  t[i+izo] = sqrt((xo-xs)*(xo-xs)+(zo-zs)*(zo-zs))/v0;
 		}
	}



	free1float(vxx);
	free1float(vxz);
	free1float(vzz);
	free1float(s);
	free1((void*)ray);

}

void ddt(float p,float q,float c,float s, float *dv,float v,float *d2t,
	float *cuv)
/************************************************************************
ddt - calculate second derivatives of traveltime with respect to x,y
************************************************************************/
{
 	float ov2,m11,m12,m22;

	ov2 = 1.0/(v*v);
	m11 = p/q;

/*	calculate 2nd column of m	*/
	m12 = -(dv[0]*c-dv[1]*s)*ov2;
	m22 = -(dv[0]*s+dv[1]*c)*ov2;

/*	compute h*m*h^T	*/
	d2t[0] = m11*c*c+2.0*m12*c*s+m22*s*s;
	d2t[1] = (m22-m11)*c*s+m12*(c*c-s*s);
	d2t[2] = m11*s*s-2.0*m12*c*s+m22*c*c;

/*	compute the curvature of raypath	*/
	*cuv = fabs(m12)*v;

}


void trans(int nx,int nz,int nxt,int nx0,float *v,float *vt) 
{
	int ix,iz,i,i0;

	for(ix=0*nx; ix<nxt; ++ix){
		i = ix*nz;
		i0 = (ix+nx0)*nz;
		for(iz=0; iz<nz; ++iz)
			vt[i+iz] = v[i0+iz];
	}
}

void voint(Geo2d geo2dv,float *v,Geo2d geo2dt,float *ov2,int npv,float *vo)
{

	int nx=geo2dv.nx,nz=geo2dv.nz,nxo=geo2dt.nx,nzo=geo2dt.nz;
	int i,io,jx,jz,ixo,izo;
	float fx=geo2dv.fx,fz=geo2dv.fz,odx=geo2dv.odx,odz=geo2dv.odz;
	float fxo=geo2dt.fx,fzo=geo2dt.fz,dxo=geo2dt.dx,dzo=geo2dt.dz;
	float x,z,ax,az,sx,sz,temp;


	for(ixo=0,x=fxo; ixo<nxo; ++ixo,x+=dxo){
		ax = (x-fx)*odx;
		jx = ax;
		sx = ax-jx;
		if(jx < 0) jx = 0;
		if(jx > nx-2) jx = nx-2;
		for(izo=0,z=fzo; izo<nzo; ++izo,z+=dzo){
			az = (z-fz)*odz;
			jz = az;
			sz = az-jz;
			if(jz < 0) jz = 0;
			if(jz > nz-2) jz = nz-2;
			io = ixo*nzo+izo;
			i = jx*nz+jz;
		    	temp = (1.0-sx)*((1.0-sz)*v[i]+sz*v[i+1])
      			    	+sx*((1.0-sz)*v[i+nz]+sz*v[i+nz+1]);
			if(npv) vo[io] = temp;
			ov2[io] = 1.0/(temp*temp);
		}
	}
}

 
void dv2(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz) 
/*calculate second DeriVatives from a 2D VELocity grid by finite difference*/
/*  	input: 	velocity v 
  	output: vxx,vxz and vzz  */	 
{ 
  	int  nx=geo2dv.nx, nz=geo2dv.nz, ix, iz, i; 	 
 	float   odx=geo2dv.odx,odz=geo2dv.odz,odxx,odzz,odxz; 
   	
	odxx = odx*odx; 
	odzz = odz*odz; 
	odxz = 0.25*odx*odz; 

/*	initialize	*/
	for(ix=0; ix<nx; ++ix)    
 	    for(iz=0; iz<nz; ++iz){
 		i = ix*nz+iz;
		vxx[i] = 0.; 
		vxz[i] = 0.; 
		vzz[i] = 0.; 
	}
	
/*	finite difference 	*/
	for(ix=1; ix<nx-1; ++ix)    
 	    for(iz=1; iz<nz-1; ++iz){
 		i = ix*nz+iz;
 		vxx[i] = odxx*(v[i+nz]-2.0*v[i]+v[i-nz]); 
 		vxz[i] = odxz*(v[i+nz+1]-v[i+nz-1]-v[i-nz+1]+v[i-nz-1]); 
      		vzz[i] = odzz*(v[i+1]-2.0*v[i]+v[i-1]);
	}

}

void eiknl(Geo2d geo2dt,float *t,float *ov2,float tmax) 
/* 	compute traveltime in shadow zones by eikonal equation
  input
   t 		traveltimes form ray tracing
   ov2 		slowness squares
   tmax	 	upper limit of ordinary traveltime value
  output:
   t	traveltime (unchanged if t<=tmax)
*/
{
	int nx=geo2dt.nx,nz=geo2dt.nz,ix,iz;
	float dx=geo2dt.dx,dz=geo2dt.dz;
	float tx2,tz2,t0,tl,tr,temp,odx2,*tt1,*tt2;

	tt1 = alloc1float(nx);
	tt2 = alloc1float(nx);

	odx2 = 1.0/(dx*dx);

	for(ix=0; ix<nx; ++ix)
		tt1[ix] = t[ix*nz];

	for(iz=1; iz<nz; ++iz){ 
	    for(ix=0; ix<nx; ++ix) 
		tt2[ix] = t[ix*nz+iz];

	    for(ix=1; ix<nx; ++ix){
/*	if traveltime is abnormal and the upper is normal	*/
		t0 = tt1[ix];
		if(tt2[ix] > tmax && t0 <= tmax) {
			tl = tr = 0.0; 
			if(ix > 0) tl = tt1[ix-1];
			if(ix < nx-1) tr = tt1[ix+1];
                        tx2 = (tl-t0)*(tl-t0);
                        temp = (tr-t0)*(tr-t0);
                        if(tx2>temp) tx2 = temp;
                        temp = 0.25*(tl-tr)*(tl-tr);
                        if(tx2>temp) tx2 = temp;
                        tx2 *= odx2;

			tz2 = 0.5*(ov2[ix*nz+iz-1]+ov2[ix*nz+iz])-tx2;
			if(tz2 >= 0) tt2[ix] = t0+dz*sqrt(tz2);
		}
	    }
	    for(ix=0; ix<nx; ++ix){
		t[ix*nz+iz] = tt2[ix];
		tt1[ix] = tt2[ix];
	    }
	}

	free1float(tt1);
	free1float(tt2);

}

void vel2Interp(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz,
     float x,float z,
     float *u,float *ux,float *uz,float *uxx,float *uxz,float *uzz) 
/*************************************************************************
vel2Interp - Function to support interpolation of velocity and its
		derivatives	 
*************************************************************************
Input:
x,z		2D coordinate at which to interpolate 
v		velocity array(nz,nx) on unoform grids 
vxx,vxz,vzz	second derivaitve arrays(nz,nx) on uniform grids 

Output:  
u 	v(x,z)
ux  	dv/dx
uz  	dv/dz
uxx	ddv/dxdx
uxz 	ddv/dzdx
uzz	ddv/dzdz
*************************************************************************
Author: Zhenyue Liu, CSM 1995
*************************************************************************/
{
	int k0,k1,k2,k3,jx,jz,nx=geo2dv.nx,nz=geo2dv.nz;
	float fx=geo2dv.fx,dx=geo2dv.dx,odx=geo2dv.odx;
	float fz=geo2dv.fz,dz=geo2dv.dz,odz=geo2dv.odz;
	float ax,az,sx,sz,sxx,szz,a0,a1,a2,a3;
	float g0,g1,g2,g3,gx0,gx1,gx2,gx3,gz0,gz1,gz2,gz3;

/*	determine interpolate coefficients	*/
	    ax = (x-fx)*odx;
	    jx = ax;
	    if(jx<0) jx = 0;
	    if(jx>nx-2) jx = nx-2;
	    sx = ax-jx;
	    az = (z-fz)*odz;
	    jz = az;
	    if(jz<0) jz = 0;
	    if(jz>nz-2) jz = nz-2;
	    sz = az-jz;

	    sxx = 0.5*sx*(1.0-sx)*dx*dx;
	    szz = 0.5*sz*(1.0-sz)*dz*dz;

	    a0 = (1.0-sx)*(1.0-sz);
	    a1 = (1.0-sx)*sz;
	    a2 = sx*(1.0-sz);
	    a3 = sx*sz;

/*	    set the table of indices for interpolation	*/
	    k0 = nz*jx+jz;
	    k1 = k0+1;
 	    k2 = k0+nz;
	    k3 = k2+1;

	    g0 = v[k0];
	    g1 = v[k1];
	    g2 = v[k2];
	    g3 = v[k3];
	    gx0 = vxx[k0];
	    gx1 = vxx[k1];
	    gx2 = vxx[k2];
	    gx3 = vxx[k3];
	    gz0 = vzz[k0];
	    gz1 = vzz[k1];
	    gz2 = vzz[k2];
	    gz3 = vzz[k3];

/*	interpolation	*/
	    *uxx = a0*gx0+a1*gx1+a2*gx2+a3*gx3;
	    *uxz = a0*vxz[k0]+a1*vxz[k1]+a2*vxz[k2]+a3*vxz[k3];
	    *uzz = a0*gz0+a1*gz1+a2*gz2+a3*gz3;

	    *u = a0*g0+a1*g1+a2*g2+a3*g3-(sxx*(*uxx)+szz*(*uzz));
	    *ux = ((1.0-sz)*(g2-g0-sxx*(gx2-gx0)-szz*(gz2-gz0))
      		    +sz*(g3-g1-sxx*(gx3-gx1)-szz*(gz3-gz1)))*odx
      		    +(sx-0.5)*dx*(*uxx);
	    *uz = ((1.0-sx)*(g1-g0-sxx*(gx1-gx0)-szz*(gz1-gz0))
      		    +sx*(g3-g2-sxx*(gx3-gx2)-szz*(gz3-gz2)))*odz
      		    +(sz-0.5)*dz*(*uzz);
}

void vintp(Geo2d geo2dv,float *v,float x,float z,float *u) 
/*************************************************************************
 Function to support interpolation of a single fuction 
*************************************************************************
Input:
x,z		2D coordinate at which to interpolate
v		array(nz,nx) on unoform grids 
Output: 
u		v(x,z)
*************************************************************************
Author: Zhenyue Liu, CSM 1995
*************************************************************************/
{

	int k0,k1,k2,k3,jx,jz,nx=geo2dv.nx,nz=geo2dv.nz;
	float fx=geo2dv.fx,odx=geo2dv.odx;
	float fz=geo2dv.fz,odz=geo2dv.odz;
	float ax,az,sx,sz;

/*	determine interpolate coefficients	*/
	    ax = (x-fx)*odx;
	    jx = ax;
	    if(jx<0) jx = 0;
	    if(jx>nx-2) jx = nx-2;
	    sx = ax-jx;
	    az = (z-fz)*odz;
	    jz = az;
	    if(jz<0) jz = 0;
	    if(jz>nz-2) jz = nz-2;
	    sz = az-jz;

/*	    set the table of indices for interpolation	*/
	    k0 = nz*jx+jz;
	    k1 = k0+1;
 	    k2 = k0+nz;
	    k3 = k2+1;

/*	interpolation	*/
	    *u = (1.0-sx)*((1.0-sz)*v[k0]+sz*v[k1])
      			+sx*((1.0-sz)*v[k2]+sz*v[k3]);

}

/*************************************************************************
Function to round x or z value of ray to the nearest x,z value of grid
*************************************************************************/
float round2grid(float input, float factor)
{
	return round(input/factor) * factor;
}
