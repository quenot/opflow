/*------------------------------------------------------------------------------------*/
/*  Copyright (C) 1995, 1996 and 1997  Georges QUENOT  LIMSI-CNRS                     */
/*  Copyright (C) 1998                 Georges QUENOT  CLIPS-IMAG                     */
/*  Version 2.00 Last revision: may 25, 1998                                          */
/*  This software comes with ABSOLUTELY NO WARRANTY                                   */
/*                                                                                    */
/*  Free use for educational and research purpose only of this software is possible   */
/*  under following conditions:                                                       */
/*  - No redistribution or inclusion in a product with or without modifications,      */
/*  - Maintainance of the copyright notice in (modified or not) source files,         */
/*  - Modifications and enhancements possible for local use only,                     */
/*  - If publications are made about a system using this software or part of it,      */
/*    the source software and author must be referenced. Citing any one paper of:     */
/*      ftp://ftp.limsi.fr/pub/quenot/articles/icassp92.ps                            */
/*      ftp://ftp.limsi.fr/pub/quenot/articles/imacom96.ps                            */
/*      ftp://ftp.limsi.fr/pub/quenot/articles/mva96.ps                               */
/*    is sufficient.                                                                  */
/*                                                                                    */
/*  Use within a commercial product or within the execution of a contract may         */
/*  require a licence agreement or a cooperation convention with LIMSI-CNRS.          */
/*------------------------------------------------------------------------------------*/

#include <strings.h>
#include "glob.h"
#include "tens.h"
#include "imio.h"
#include "inpo.h"
#include "opfl.h"
#include "velo.h"

/*------------------------------------------------------------------------------------*/

#define PROGRAM "opflow"
#define CM1     0.03
#define CM2     1.00
#define EPS     4.0
#define VS      2

/*------------------------------------------------------------------------------------*/
/*--------------------- ARGUMENTS EXTRACTION AND USAGE PROCEDURES --------------------*/
/*------------------------------------------------------------------------------------*/

void usage(FILE *fp) {
    fprintf(fp,"Usage: %s input-path output-path seq-name [-options]\n",PROGRAM);
    fprintf(fp,"\nGeneral options:\n");
    fprintf(fp,"  [-QUIET]              : silent mode, default is verbose\n");
    fprintf(fp,"\nInput sequence specification parameters:\n");
    fprintf(fp,"  [-P n]                : first frame offset relative to -M value,\n");
    fprintf(fp,"                          default value is 1\n");
    fprintf(fp,"  [-N n]                : last frame offset relative to -M value,\n");
    fprintf(fp,"                          default value is 1\n");
    fprintf(fp,"  [-F n]                : first frame number, overrides -P option,\n");
    fprintf(fp,"                          default is -P option.\n");
    fprintf(fp,"  [-L n]                : last frame number, overrides -N option,\n");
    fprintf(fp,"                          default is -N option.\n");
    fprintf(fp,"  [-PGM]                : pgm grayscale file format\n");
    fprintf(fp,"  [-PPM]                : ppm color file format\n");
    fprintf(fp,"  [-B cols rows [bands]]: raw image size, bands defaults to 1\n");
    fprintf(fp,"                          (the (cols x rows x bands) last bytes\n");
    fprintf(fp,"                          of the image file are read; therefore\n");
    fprintf(fp,"                          most image file formats are supported)\n");
    fprintf(fp,"                          sun raster file format is the default\n");
    fprintf(fp,"\nOutput result specification parameters:\n");
    fprintf(fp,"  [-M f]                : frame for which flow is computed\n");
    fprintf(fp,"                          may be float, default value is 0.5\n");
    fprintf(fp,"  [-S n]                : result border shrink, default value is 0\n");
    fprintf(fp,"  [-II]                 : computes an interpolated image\n");
    fprintf(fp,"                          (useful if -M value is not integer)\n");
    fprintf(fp,"  [-EX]                 : extrapolate for interpolated image\n");
    fprintf(fp,"  [-C name]             : correct result filename\n");
    fprintf(fp,"  [-CSCALE f]           : correct result scale\n");
    fprintf(fp,"  [-CM1 f]              : 1D confidence measure threshold (float)\n");
    fprintf(fp,"  [-CM2 f]              : 2D confidence measure threshold (float)\n");
    dpusage(fp,1,(fdata) EPS);
}

/*------------------------------------------------------------------------------------*/

int getsnumbers(int *pac, char *av[], fdata *pm, int *pnf, int *pnl)
{
    /*--------------------------*/
    /* get source image numbers */
    /*--------------------------*/
    int f,l,n,p;
    fdata m;
    if (getiarg(pac,av,"-F",-1073741824,&f) == FAILURE) return(FAILURE);
    if (getiarg(pac,av,"-L",-1073741824,&l) == FAILURE) return(FAILURE);
    if (getfarg(pac,av,"-M",(fdata) 0.5,&m) == FAILURE) return(FAILURE);
    if (getiarg(pac,av,"-P",1,&p) == FAILURE) return(FAILURE);
    if (getiarg(pac,av,"-N",1,&n) == FAILURE) return(FAILURE);
    if (n < 0) {
        fprintf(stderr,"ERROR: invalid next frame offset: %d\n",n);
        fprintf(stderr,"-N value must be non negative\n");
        return(FAILURE);
    }
    if (p < 0) {
        fprintf(stderr,"ERROR: invalid previous frame offset: %d\n",p);
        fprintf(stderr,"-P value must be non negative\n");
        return(FAILURE);
    }
    if ((n + p) == 0) {
        fprintf(stderr,"ERROR: incompatible frame offsets: %d %d\n",n,p);
        fprintf(stderr,"-P and -N values must not be both null\n");
        return(FAILURE);
    }
    if ((m != floor(m)) && ((n + p) == 1)) {
        fprintf(stderr,"ERROR: incompatible frame offsets: %d %d\n",n,p);
        fprintf(stderr,"-P value + -N value must be greater ");
        fprintf(stderr,"than 1 if -M value is not integer\n");
        return(FAILURE);
    }
    if (f != -1073741824) p = f; else p = (int) ceil(m-(fdata)p);
    if (l != -1073741824) n = l; else n = (int) floor(m+(fdata)n);
    if ((n-p) <= 0) {
        fprintf(stderr,"ERROR: invalid extreme frames specification: %d %d\n",n,p);
        fprintf(stderr,"-F or -L values may be incorrect\n");
        return(FAILURE);
    }
    *pm = m;
    *pnf = p;
    *pnl = n;
    return(SUCCESS);
}

int getfnames(int *pac, char *av[], fdata m, int p, int n, char *iname, char *oname,
char *inamei)
{
    /*---------------*/
    /* get filenames */
    /*---------------*/
    int i,ii;
    char str[32];
    ii = getbarg(pac,av,"-II");
    if (*pac < 4) {
        return(FAILURE);
    }
    sprintf(iname,"%s/%s",av[1],av[3]);
    prflp(m,str);
    sprintf(oname,"%s/%s.%s%s",av[2],PROGRAM,av[3],str);
    if (ii) {
        if (m == floor(m)) {
            sprintf(inamei,"%s/%s%d",av[2],av[3],(int)m);
        } else {
            sprintf(inamei,"%s/%s%1.3f",av[2],av[3],m);
        }
    } else {
        inamei[0] = 0;
    }
    for (i = 0; i < *pac-3; i++) av[i] = av[i+3];
    *pac -= 3;
    return(SUCCESS);
}

int getshrink(int *pac, char *av[], int *pns)
{
    if (getiarg(pac,av,"-S",0,pns) == FAILURE) return(FAILURE);
    if (*pns < 0) {
        fprintf(stderr,"ERROR: invalid shrink value: %d\n",*pns);
        fprintf(stderr,"-S value must be non negative\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int getcorrect(int *pac, char *av[], char *cname, fdata *cm1, fdata *cm2, fdata *cscale)
{
    /*-----------------------------*/
    /* get correct result filename */
    /*-----------------------------*/
    if (getsarg(pac,av,"-C",cname) == FAILURE) return(FAILURE);
    if (getfarg(pac,av,"-CM1",(fdata) CM1,cm1) == FAILURE) return(FAILURE);
    if (*cm1 < 0) {
        fprintf(stderr,"ERROR: invalid threshold: %f\n",(float) *cm1);
        fprintf(stderr,"-CM1 value must be non negative\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-CM2",(fdata) CM2,cm2) == FAILURE) return(FAILURE);
    if (*cm2 < 0) {
        fprintf(stderr,"ERROR: invalid threshold: %f\n",(float) *cm2);
        fprintf(stderr,"-CM2 value must be non negative\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-CSCALE",1.0,cscale) == FAILURE) return(FAILURE);
    if (*cscale <= 0) {
        fprintf(stderr,"ERROR: invalid scale: %f\n",(float) *cscale);
        fprintf(stderr,"-CSCALE value must be positive\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int getargs(int ac, char *av[], int *pnx, int *pny, int *pnb, int *pfo, int *pnf,
int *pnl, int *pns, DPCONTROL *dpcontrol, char *iname, char *oname, char *inamei,
char *cname, fdata *cm1, fdata *cm2, fdata *cscale, int *pex)
{
    /*-----------------------*/
    /* get program arguments */
    /*-----------------------*/
    fdata m;
    char inamep[256];
    *pex = getbarg(&ac,av,"-EX");
    if (getsnumbers(&ac,av,&m,pnf,pnl) == FAILURE) return(FAILURE);
    if (getfnames(&ac,av,m,*pnf,*pnl,iname,oname,inamei) == FAILURE) return(FAILURE);
    sprintf(inamep,"%s%d",iname,*pnf);
    if (getimsize(&ac,av,inamep,pnx,pny,pnb,pfo) == FAILURE) return(FAILURE);
    if (getshrink(&ac,av,pns) == FAILURE) return(FAILURE);
    if (getcontrol(&ac,av,dpcontrol,(fdata) EPS,VS,1) == FAILURE) return(FAILURE);
    dpcontrol->m = m;
    if (getcorrect(&ac,av,cname,cm1,cm2,cscale) == FAILURE) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*  Backward reconstruction: build im(t) from im(t-dt) and v(t).                      */
/*  Forward reconstruction: build im(t+dt) from im(t) and v(t).                       */
/*------------------------------------------------------------------------------------*/

/*  If t is integer:                                                                  */
/*    1) Backward reconstruction with all integer values of dt and whole image        */
/*    2) Backward reconstruction with all integer values of dt and partial image      */
/*    3) Forward reconstruction with all integer values of dt and whole image         */
/*    4) Forward reconstruction with all integer values of dt and partial image       */
/*    5)-8) Like 1)-4) but using the correct flow field                               */
/*  If t is not integer:                                                              */
/*    1) Backward reconstruction using the two nearest images and whole image         */
/*    2) Backward reconstruction using the two nearest images and partial image       */

/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/

int confmeas(CIMAGE *cim, CIMAGE *com, fdata thr1, fdata thr2)
{
    int k,nx,ny,mx,my,nb,ii,jj,jm,jp,im,ip;
    uchar *omd,*imd;
    fdata *smd,dx,dy,dxx,dxy,dyy,dd,dax,day,dmin,dmax;
    thr1 = thr1*thr1;
    thr2 = thr2*thr2;
    com->nx = nx = cim->nx;
    com->ny = ny = cim->ny;
    nb = cim->nb;
    imd = cim->im;
    if ((malloc1((void **) &omd,nx*ny,sizeof(uchar)) == FAILURE)
     || (malloc1((void **) &smd,nb*nx*ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: confmeas: fail to allocate buffers\n");
        return(FAILURE);
    }
    com->nb = 1;
    com->im = omd;
    mx = nx-1;
    my = ny-1;
    for (jj = 0; jj < nx*ny; jj += nx) {
        for (ii = 0; ii < nx; ii++) {
            omd[ii+jj] = 0;
        }
    }
    for (k = 0; k < nb; k++) {
        for (jj = 0; jj < nx*ny; jj += nx) {
            jm = imax(0,jj-nx);
            jp = imin(my*nx,jj+nx);
            for (ii = 0; ii < nx; ii++) {
                im = imax(0,ii-1);
                ip = imin(nx-1,ii+1);
                smd[ii+jj] = imd[ii+jj];
/*              smd[ii+jj] = 0.04*imd[im+jm]+0.08*imd[ii+jm]+0.04*imd[ip+jm]+
                             0.08*imd[im+jj]+0.52*imd[ii+jj]+0.08*imd[ip+jj]+
                             0.04*imd[im+jp]+0.08*imd[ii+jp]+0.04*imd[ip+jp]; */
/*              smd[ii+jj] = 0.06*imd[im+jm]+0.13*imd[ii+jm]+0.06*imd[ip+jm]+
                             0.13*imd[im+jj]+0.24*imd[ii+jj]+0.13*imd[ip+jj]+
                             0.06*imd[im+jp]+0.13*imd[ii+jp]+0.06*imd[ip+jp]; */
/*              smd[ii+jj] = 0.11*imd[im+jm]+0.11*imd[ii+jm]+0.11*imd[ip+jm]+
                             0.11*imd[im+jj]+0.12*imd[ii+jj]+0.11*imd[ip+jj]+
                             0.11*imd[im+jp]+0.11*imd[ii+jp]+0.11*imd[ip+jp]; */
            }
        }
        imd += nx*ny;
        smd += nx*ny;
    }
    smd -= nb*nx*ny;
    for (k = 0; k < nb; k++) {
        dmin = 255.0;
        dmax = 0.0;
        for (jj = 0; jj < nx*ny; jj += nx) {
            for (ii = 0; ii < nx; ii++) {
                if (dmin > smd[ii+jj]) dmin = smd[ii+jj];
                if (dmax < smd[ii+jj]) dmax = smd[ii+jj];
            }
        }
        thr1 *= (dmax-dmin)*(dmax-dmin);
        for (jj = 0; jj < nx*ny; jj += nx) {
            jm = imax(0,jj-nx);
            jp = imin(my*nx,jj+nx);
            for (ii = 0; ii < nx; ii++) {
                im = imax(0,ii-1);
                ip = imin(nx-1,ii+1);
                dx = 0.5*(smd[ip+jj]-smd[im+jj]);
                dy = 0.5*(smd[ii+jp]-smd[ii+jm]);
                if (((dd = dx*dx+dy*dy) > thr1) && (omd[ii+jj] < 1)) {
                    omd[ii+jj] = 1;
                    dxx = smd[ip+jj]+smd[im+jj]-2*smd[ii+jj];
                    dyy = smd[ii+jp]+smd[ii+jm]-2*smd[ii+jj];
                    dxy = 0.25*((smd[ip+jp]+smd[im+jm])-(smd[ip+jm]+smd[im+jp]));
                    dax = dy*dxx-dx*dxy;
                    day = dx*dyy-dy*dxy;
                    dd = (dax*dax+day*day)/(dd*dd);
                    if ((dd > thr2) && (omd[ii+jj] < 2)) omd[ii+jj] = 2;     
                }
            }
        }
        smd += nx*ny;
    }
    smd -= nb*nx*ny;
    if (free1((void **) &smd,nb*nx*ny,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: confmeas: fail to free buffers\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int confidence(VIMAGE *vim, CIMAGE *cim, fdata thr1, fdata thr2, int cml)
{
    int i,j,nx,ny;
    uchar *imd;
    fdata **dx,**dy;
    CIMAGE ccm;
    if (confmeas(cim,&ccm,thr1,thr2) == FAILURE) return(FAILURE);
    nx = ccm.nx;
    ny = ccm.ny;
    imd = ccm.im;
    dx = vim->dx;
    dy = vim->dy;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            if (((int) imd[i+nx*j]) < cml) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
    }
    if (free1((void **) &imd,nx*ny,sizeof(uchar)) == FAILURE) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

void errorstat(VIMAGE *vim, VIMAGE *cvim)
{
    int i,j,nx,ny,nx0,ny0,nx1,ny1;
    fdata **dx,**dy,**cdx,**cdy;
    double vx,vy,cvx,cvy,cv,dvx,dvy,dv,c,m,s,m2,s2,m3,s3;
    double ang,scal,min,max,min2,max2,min3,max3,degrad=90.0/asin(1.0);
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    nx0 = vim->nx0;
    ny0 = vim->ny0;
    nx1 = vim->nx1;
    ny1 = vim->ny1;
    cdx = cvim->dx;
    cdy = cvim->dy;
    c = m = s = m2 = s2 = m3 = s3 = 0.0;
    min = min2 = min3 = 1000000000.0;
    max = max2 = max3 = 0.0;
    for (i = nx0; i < nx-nx1; i++) {
        for (j = ny0; j < ny-ny1; j++) {
            vx = dx[i][j];
            vy = -dy[i][j];
            cvx = cdx[i][j];
            cvy = -cdy[i][j];
            if (((vx != 100.0) || (vy != 100.0)) && ((cvx != 100.0) || (cvy != 100.0))) {
                dvx = vx-cvx;
                dvy = vy-cvy;
                dv = sqrt(dvx*dvx+dvy*dvy);
                cv = sqrt(cvx*cvx+cvy*cvy);
                scal = (1+vx*cvx+vy*cvy)/sqrt((1+vx*vx+vy*vy)*(1+cvx*cvx+cvy*cvy));
                ang = degrad*acos(scal);
                c += 1.0;
                m += ang;
                s += ang*ang;
                min = (ang < min) ? ang : min;
                max = (ang > max) ? ang : max;
                m2 += dv;
                s2 += dv*dv;
                min2 = (dv < min2) ? dv : min2;
                max2 = (dv > max2) ? dv : max2;
                m3 += cv;
                s3 += cv*cv;
                min3 = (cv < min3) ? cv : min3;
                max3 = (cv > max3) ? cv : max3;
            }
        }
    }
    m = m/c;
    s = sqrt(fabs(s/c-m*m));
    m2 = m2/c;
    s2 = sqrt(fabs(s2/c-m2*m2));
    m3 = m3/c;
    s3 = sqrt(fabs(s3/c-m3*m3));
    printf("Count: %d  Density: %.1f\n",(int) c,(float) ((100.0*c)/(nx*ny)));
    printf("Velocity Angle Error: %.2f+%.2f  ",(float) m,(float) s);
    printf("Min: %.2f  Max: %.2f\n",(float) min,(float) max);
    printf("Velocity Value Error: %.2f+%.2f  ",(float) m2,(float) s2);
    printf("Min: %.2f  Max: %.2f\n",(float) min2,(float) max2);
    printf("Velocity Value:       %.2f+%.2f  ",(float) m3,(float) s3);
    printf("Min: %.2f  Max: %.2f\n",(float) min3,(float) max3);
}

/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/

void exitfail()
{
    usage(stderr);
    fprintf(stderr,"ERROR: %s: execution failed\n",PROGRAM);
    exit(1);
}

/*------------------------------------------------------------------------------------*/

void dispinfo(DPCONTROL dpcontrol, int nx, int ny, int nb, int fo, int ns,
char *inamef, char *inamel, char *oname, char *inamei, char *cname, fdata cscale)
{
    char ext[6];
    printf("%s: command line OK\n",PROGRAM);
    strcpy(ext,"");
    if (fo == PGM) strcpy(ext,".pgm");
    if (fo == PPM) strcpy(ext,".ppm");
    if (fo == RAW) printf("images format:         raw\n");
    if (fo == PGM) printf("images format:         pgm\n");
    if (fo == PPM) printf("images format:         ppm\n");
    if (fo == RAS) printf("images format:         sun raster\n");
    if (fo == CCF) printf("images format:         calfonc\n");
    printf("images width:          %d\n",nx);
    printf("images height:         %d\n",ny);
    printf("number of bands:       %d\n",nb);
    printf("result shrink:         %d\n",ns);
    printf("first input image:     %s%s\n",inamef,ext);
    printf("last input image:      %s%s\n",inamel,ext);
    printf("velocity file name:    %s\n",oname);
    if (cname[0] != 0) printf("correct flow field:    %s\n",cname);
    if (cname[0] != 0) printf("correct flow scale:    %f\n",(float) cscale);
    if (inamei[0] != 0) printf("interpolated image:    %s%s\n",inamei,ext);
    dpinfo(dpcontrol,1);
}

/*------------------------------------------------------------------------------------*/
/*--------------------------------------- MAIN ---------------------------------------*/
/*------------------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    int i,j,nx,ny,nb,fo,ns,nf,nl,ex;
    fdata cm1,cm2,cscale,**cdx,**cdy,**dx,**dy;
    char iname[256],inamef[256],inamel[256],oname[256],inamei[256],cname[256];
    CIMAGE com;
    SIMAGE sim;
    VIMAGE vim,cvim;
    DPCONTROL dpcontrol;
    /*-----------------------*/
    /* get program arguments */
    /*-----------------------*/
    if (getargs(argc,argv,&nx,&ny,&nb,&fo,&nf,&nl,&ns,&dpcontrol,iname,oname,inamei,
    cname,&cm1,&cm2,&cscale,&ex) == FAILURE) exitfail();
    sprintf(inamef,"%s%d",iname,nf);
    sprintf(inamel,"%s%d",iname,nl);
    /*--------------*/
    /* display info */
    /*--------------*/
    if (dpcontrol.vb)
        dispinfo(dpcontrol,nx,ny,nb,fo,ns,inamef,inamel,oname,inamei,cname,cscale);
    /*--------------------------*/
    /* get input image sequence */
    /*--------------------------*/
    sim.nx = nx;
    sim.ny = ny;
    sim.nb = nb;
    sim.nf = nf;
    sim.nl = nl;
    if (rdsimage(iname,&sim,fo) == FAILURE) exitfail();
    /*---------------------------*/
    /* optical flow search here! */
    /*---------------------------*/
    vim.d = NULL;
    vim.nx0 = vim.ny0 = ns;
    vim.nx1 = vim.ny1 = ns;
    if (opflow(&sim,&vim,&dpcontrol) == FAILURE) exitfail();
    /* if (invdep(&vim,1,&vim,1) == FAILURE) exitfail();
    if (invdep(&vim,1,&vim,1) == FAILURE) exitfail(); */
    /*----------------------------*/
    /* compute interpolated image */
    /*----------------------------*/
    if (iminterp(&sim,&com,&vim,0,dpcontrol.m,ex) == FAILURE) exitfail();
    /*-------------------*/
    /* read correct flow */
    /*-------------------*/
    if ((cname[0] != 0) && (rdvimage(cname,&cvim,UWO,1)) == FAILURE) exitfail();
    if (cname[0] != 0) {
        cdx = cvim.dx;
        cdy = cvim.dy;
        dx = vim.dx;
        dy = vim.dy;
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                if ((cdx[i][j] != 100.0) || (cdy[i][j] != -100.0)) {
                    cdx[i][j] *= cscale;
                    cdy[i][j] *= cscale;
                } else {
                    dx[i][j] = 100.0;
                    dy[i][j] = -100.0;
                }
            }
        }
    }
    /*----------------*/
    /* put velocities */
    /*----------------*/
    strcat(oname,".cm0");
    if ((cname[0] != 0) && (dpcontrol.vb)) {
        printf("\nFull velocity error calculation\n");
        printf("No confidence measure\n");
        errorstat(&vim,&cvim);
    } else {
        errorstat(&vim,&vim);
    }
    if (wrvimage(oname,&vim,UWO,0) == FAILURE) exitfail();
    if (confidence(&vim,&com,cm1,cm2,1) == FAILURE) exitfail();
    oname[strlen(oname)-1] = '1';
    if ((cname[0] != 0) && (dpcontrol.vb)) {
        printf("\nFull velocity error calculation\n");
        printf("1D confidence measure\n");
        errorstat(&vim,&cvim);
    }
    if (wrvimage(oname,&vim,UWO,0) == FAILURE) exitfail();
    if (confidence(&vim,&com,cm1,cm2,2) == FAILURE) exitfail();
    oname[strlen(oname)-1] = '2';
    if ((cname[0] != 0) && (dpcontrol.vb)) {
        printf("\nFull velocity error calculation\n");
        printf("2D confidence measure\n");
        errorstat(&vim,&cvim);
        printf("\n");
    }
    if (wrvimage(oname,&vim,UWO,1) == FAILURE) exitfail();
    /*------------------------*/
    /* put interpolated image */
    /*------------------------*/
    if ((inamei[0] != 0) && (wrcimagec(inamef,inamei,&com,fo,0) == FAILURE))
        exitfail();
    /*--------------*/
    /* free buffers */
    /*--------------*/
    if ((free1((void **) &com.im,nb*nx*ny,sizeof(uchar)) == FAILURE)
     || (frsimage(&sim) == FAILURE)) exitfail();
    if ((cname[0] != 0) && (frvimage(&cvim) == FAILURE)) exitfail();
    /*--------------*/
    /* display info */
    /*--------------*/
    if (dpcontrol.vb) {
        printf("Maximum allocated memory: %8d bytes\n",MaxUsedMemory());
        printf("Still allocated memory:   %8d bytes\n",UsedMemory());
        printf("%s: execution completed\n",PROGRAM);
    }
    exit(0);
}

/*------------------------------------------------------------------------------------*/
