/*------------------------------------------------------------------------------------*/
/*  Copyright (C) 1995, 1996 and 1997  Georges QUENOT  LIMSI-CNRS                     */
/*  Copyright (C) 1998 and 1999        Georges QUENOT  CLIPS-IMAG                     */
/*  Version 2.00 Last revision: november 16, 1999                                     */
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

#include "glob.h"
#include "tens.h"
#include "imio.h"
#include "inpo.h"
#include "velo.h"
#include "opfl.h"

/*------------------------------------------------------------------------------------*/
/*--------------------- ARGUMENTS EXTRACTION AND USAGE PROCEDURES --------------------*/
/*------------------------------------------------------------------------------------*/

void mfusage(FILE *fp)
{
    fprintf(fp,"  [-MF n]               : DP multi-frame method, numer of frames\n");
    fprintf(fp,"                          default is two frames\n");
    fprintf(fp,"  [-SUM]                : DP multi-frame distance SUM option\n");
    fprintf(fp,"                          default is MAX option\n");
    fprintf(fp,"  [-VS n]               : DP variable scale strategy\n");
    fprintf(fp,"                          0: fixed scale,\n");
    fprintf(fp,"                          1: full range variation,\n");
    fprintf(fp,"                          2: adaptive variation,\n");
    fprintf(fp,"                          default is 2\n");
    fprintf(fp,"  [-MD n]               : DP maximun relative displacement\n");
    fprintf(fp,"                          during adaptive variation,\n");
    fprintf(fp,"                          default is %f\n",(float) MD);
}

void dpusage(FILE *fp, int mf, fdata eps)
{
    fprintf(fp,"\nOrthogonal Dynamic Programming control parameters:\n");
    fprintf(fp,"  [-GOA]                : DP enable gain-offset correction\n");
    fprintf(fp,"                          default is no correction\n");
    fprintf(fp,"  [-HX]                 : horizontal shrink of images\n");
    fprintf(fp,"                          (useful for interlaced video fields)\n");
    fprintf(fp,"  [-DNI]                : DP do not interpolate distances\n");
    fprintf(fp,"                          default is interpolate\n");
    fprintf(fp,"  [-BSI]                : DP use bicubic spline interpolation\n");
    fprintf(fp,"                          default is bilinear interpolation\n");
    fprintf(fp,"  [-MB n]               : multiply number of bands with smoothing,\n");
    fprintf(fp,"                          default value is 1.\n");
    fprintf(fp,"  [-EPS f]              : DP parameter epsilon (float)\n");
    fprintf(fp,"                          default value is %f\n",(float) eps);
    fprintf(fp,"  [-DL f]               : DP Ln distance specification (float)\n");
    fprintf(fp,"                          default value is %f\n",(float) DLF);
    fprintf(fp,"  [-STX]                : DP horizontal stereovision option\n");
    fprintf(fp,"  [-STY]                : DP vertical stereovision option\n");
    fprintf(fp,"  [-ID h|v]             : DP initial slicing direction (char)\n");
    fprintf(fp,"                          default value is h (horizontal)\n");
    fprintf(fp,"  [-NP n]               : DP parameter number of passes (int)\n");
    fprintf(fp,"                          default value is %d\n",NP);
    fprintf(fp,"  [-TAU f]              : DP decreasing control parameter (float)\n");
    fprintf(fp,"                          default value is %f\n",(float) TAU);
    fprintf(fp,"  [-SPW f]              : DP strip spacing and width factor (float)\n");
    fprintf(fp,"                          default value is %f\n",(float) SPW);
    fprintf(fp,"  [-WIN f]              : DP search window width factor (float)\n");
    fprintf(fp,"                          default value is %f\n",(float) WIN);
    fprintf(fp,"  [-SMO f]              : DP smoothing factor (float)\n");
    fprintf(fp,"                          default value is %f\n",(float) SMO);
    fprintf(fp,"  [-SMF n]              : DP smoothing function (int)\n");
    fprintf(fp,"                          default value is %d\n", SMF);
    fprintf(fp,"  [-ENW n]              : DP parameter extra width (int)\n");
    fprintf(fp,"                          default value is %d\n",ENW);
    fprintf(fp,"  [-EIT n]              : DP number of extra iterations (int)\n");
    fprintf(fp,"                          default value is %d\n",EIT);
    fprintf(fp,"  [-SC n]               : DP slope constraint\n");
    fprintf(fp,"                          may be 1, 2 or 3, default is %d\n",SC);
    if (mf) mfusage(fp);
}

int getcontrol(int *pac, char *av[], DPCONTROL *dpcontrol, fdata eps, int vs, int mf)
{
    /*------------------------------------------------------------------*/
    /* extracts dynamic programming control parameters from command line*/
    /*------------------------------------------------------------------*/
    char dir;
    dpcontrol->vb = !getbarg(pac,av,"-QUIET");
    if (getbarg(pac,av,"-STY")) {
        dpcontrol->dlin[0] = DXL;
        dpcontrol->dsig[0] = DXS;
        dpcontrol->id = 1;
    } else {
        dpcontrol->dlin[0] = 1.0;
        dpcontrol->dsig[0] = 1.0;
    }
    if (getbarg(pac,av,"-STX")) {
        dpcontrol->dlin[1] = DYL;
        dpcontrol->dsig[1] = DYS;
        dpcontrol->id = 0;
    } else {
        dpcontrol->dlin[1] = 1.0;
        dpcontrol->dsig[1] = 1.0;
    }
    dpcontrol->goa = getbarg(pac,av,"-GOA");
    dpcontrol->hx  = getbarg(pac,av,"-HX");
    dpcontrol->dni = getbarg(pac,av,"-DNI");
    dpcontrol->bsi = getbarg(pac,av,"-BSI");
    if (getiarg(pac,av,"-MB",MB,&dpcontrol->mb) == FAILURE) return(FAILURE);
    if (dpcontrol->mb < 1) {
        fprintf(stderr,"ERROR: invalid multiply bands value: %d\n",dpcontrol->mb);
        fprintf(stderr,"-MB value must be positive\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-EPS",eps,&dpcontrol->eps) == FAILURE) return(FAILURE);
    if (dpcontrol->eps < 0.0) {
        fprintf(stderr,"ERROR: invalid epsilon value: %f\n",(float) dpcontrol->eps);
        fprintf(stderr,"-EPS value must be non negative\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-DL",(fdata) DLF,&dpcontrol->dlf) == FAILURE) return(FAILURE);
    if (dpcontrol->dlf <= 0.0) {
        fprintf(stderr,"ERROR: invalid Ln distance specification: %f\n",
        (float) dpcontrol->dlf);
        fprintf(stderr,"-DL value must be positive\n");
        return(FAILURE);
    }
    if (getcarg(pac,av,"-ID",ID,&dir) == FAILURE) return(FAILURE);
    dpcontrol->id = 0;
    if (dir == 'v') dpcontrol->id = 1;
    else if (dir != 'h') {
        fprintf(stderr,"-ID option must be only h or v\n");
        return(FAILURE);
    }
    if (getiarg(pac,av,"-NP",NP,&dpcontrol->np) == FAILURE) return(FAILURE);
    if (dpcontrol->np <= 0) {
        fprintf(stderr,"ERROR: invalid number of passes: %d\n",dpcontrol->np);
        fprintf(stderr,"-NP value must be positive\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-TAU",(fdata) TAU,&dpcontrol->tau) == FAILURE)
        return(FAILURE);
    if (dpcontrol->tau <= 0.0) {
        fprintf(stderr,"ERROR: invalid tau value: %f\n",(float) dpcontrol->tau);
        fprintf(stderr,"-TAU value must be positive\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-SPW",(fdata) SPW,&dpcontrol->spwfac) == FAILURE)
        return(FAILURE);
    if ((dpcontrol->spwfac <= 0.0) || (dpcontrol->spwfac > 1.0)) {
        fprintf(stderr,"ERROR: invalid spwfac value: %f\n",(float) dpcontrol->spwfac);
        fprintf(stderr,"-SPW value must be within the [0.0,1.0] range\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-WIN",(fdata) WIN,&dpcontrol->winfac) == FAILURE)
        return(FAILURE);
    if (dpcontrol->winfac <= 0.0) {
        fprintf(stderr,"ERROR: invalid winfac value: %f\n",(float) dpcontrol->winfac);
        fprintf(stderr,"-WIN value must be positive\n");
        return(FAILURE);
    }
    if (getfarg(pac,av,"-SMO",(fdata) SMO,&dpcontrol->smofac) == FAILURE)
        return(FAILURE);
    if (dpcontrol->smofac <= 0.0) {
        fprintf(stderr,"ERROR: invalid smofac value: %f\n",(float) dpcontrol->smofac);
        fprintf(stderr,"-SMO value must be positive\n");
        return(FAILURE);
    }
    if (getiarg(pac,av,"-SMF",SMF,&dpcontrol->smf) == FAILURE) return(FAILURE);
    if (getiarg(pac,av,"-ENW",ENW,&dpcontrol->enw) == FAILURE) return(FAILURE);
    if (dpcontrol->enw < 0) {
        fprintf(stderr,"ERROR: invalid extra width value: %d\n",dpcontrol->enw);
        fprintf(stderr,"-ENW value must be non negative\n");
        return(FAILURE);
    }
    if (getiarg(pac,av,"-EIT",EIT,&dpcontrol->eit) == FAILURE) return(FAILURE);
    if (dpcontrol->eit < 0) {
        fprintf(stderr,"ERROR: invalid number of extra iterations: %d\n",
        dpcontrol->eit);
        fprintf(stderr,"-EIT value must be non negative\n");
        return(FAILURE);
    }
    if (getiarg(pac,av,"-SC",SC,&dpcontrol->sc) == FAILURE) return(FAILURE);
    if ((dpcontrol->sc <= 0) || (dpcontrol->sc > 3)) {
        fprintf(stderr,"ERROR: invalid DP slope constraint: %d\n",dpcontrol->sc);
        fprintf(stderr,"-SC value must be 1, 2 or 3\n");
        return(FAILURE);
    }
    if (mf) {
        if (getiarg(pac,av,"-MF",MF,&dpcontrol->mf) == FAILURE) return(FAILURE);
        if (dpcontrol->mf < 2) {
            fprintf(stderr,"ERROR: invalid number of frames: %d\n",dpcontrol->mf);
            fprintf(stderr,"-MF value must be at least 2\n");
            return(FAILURE);
        }
        dpcontrol->sum = getbarg(pac,av,"-SUM");
        if (getiarg(pac,av,"-VS",vs,&dpcontrol->vs) == FAILURE) return(FAILURE);
        if ((dpcontrol->vs < 0) || (dpcontrol->vs > 2)) {
            fprintf(stderr,"ERROR: invalid variable scale selection: %d\n",dpcontrol->vs);
            fprintf(stderr,"-VS value must be 0, 1 or 2\n");
            return(FAILURE);
        }
        if (getfarg(pac,av,"-MD",(fdata) MD,&dpcontrol->md) == FAILURE) return(FAILURE);
        if (dpcontrol->md <= 0) {
            fprintf(stderr,"ERROR: invalid maximum displacement: %f\n",
            (float) dpcontrol->md);
            fprintf(stderr,"-MD value must be positive\n");
            return(FAILURE);
        }
    } else {
        dpcontrol->vs = vs;
        dpcontrol->mf = MF;
        dpcontrol->sum = 0;
        dpcontrol->md = MD;
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

void dpinfo(DPCONTROL dpcontrol, int mf)
{
    char str[32];
    printf("DP first slicing:      %s\n",dpcontrol.id ? "vertical" : "horizontal");
    prflp(dpcontrol.m,str);
    printf("DP frame number:       %s\n",str);
    prflp(dpcontrol.eps,str);
    printf("DP epsilon:            %s\n",str);
    printf("DP extra width:        %d\n",dpcontrol.enw);
    printf("DP gain-offset:        %s\n",dpcontrol.goa ? "ON" : "OFF");
    printf("DP number of passes:   %d\n",dpcontrol.np);
    if (dpcontrol.mf > 2) {
        printf("DP multi-frame option: %s\n",dpcontrol.sum ? "sum" : "max");
    }
    printf("DP slope constraint:   %d\n",dpcontrol.sc);
    printf("DP extra iterations:   %d\n",dpcontrol.eit);
    prflp(dpcontrol.dlf,str);
    printf("DP global distance:    L%s\n",str);
}

/*------------------------------------------------------------------------------------*/
/*------------------------- IMAGE TRANSFORMATION PROCEDURES --------------------------*/
/*------------------------------------------------------------------------------------*/

int imin(int i, int j) {return((i < j) ? i : j);}
int imax(int i, int j) {return((i < j) ? j : i);}

/*------------------------------------------------------------------------------------*/

int xyswap(VIMAGE *vim)
{
    /*--------------------------------------------------------------*/
    /* performs a horizontal / vertical swap of velocity components */
    /*--------------------------------------------------------------*/
    fdata ***d,**dx,**dy,***nd,**ndx,**ndy;
    int nx,ny,i,j;
    nx = vim->nx;
    ny = vim->ny;
    d = vim->d;
    dx = d[0];
    dy = d[1];
    if (malloc3((void ****) &nd,2,ny,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    ndx = nd[0];
    ndy = nd[1];
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            ndx[j][i] = dy[i][j];
            ndy[j][i] = dx[i][j];
        }
    }
    if (free3((void ****) &d,2,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    vim->d = nd;
    vim->dx = nd[0];
    vim->dy = nd[1];
    vim->nx = ny;
    vim->ny = nx;
    return(SUCCESS);
}

int getwwww(SIMAGE *sim, VIMAGE *vim, fdata m, int *pwl, int *pwr, int *pwd, int *pwu)
{
    /*---------------------------------------------------------------*/
    /* Gets the parameters (border shrinks) of a window of a maximum */
    /* size so that the velocity field applied to it keeps it inside */
    /* the image size. The applied velocity field is considered when */
    /* bringing images from planes nf and nl into a plane m.         */
    /*---------------------------------------------------------------*/
    int i,j,wl,wr,wu,wd,nx,ny,mx,my,nf,nl;
    fdata **dx,**dy,el,er,ed,eu,rx,ry,rmx,rmy,l0,l1;
    nf = sim->nf;
    nl = sim->nl;
    l0 = ((fdata) nf)-m;
    l1 = ((fdata) nl)-m;
    nx = vim->nx;
    ny = vim->ny;
    mx = nx-1;
    my = ny-1;
    rmx = (fdata) mx;
    rmy = (fdata) my;
    dx = vim->dx;
    dy = vim->dy;
    wl = -1;
    el = -INFINI;
    while (el < 0) {
        wl++;
        if (wl > mx) {
            fprintf(stderr,"ERROR: getwwww: found velocities are out of bounds1\n");
            return(FAILURE);
        }
        el = 0;
        for (j = 0; j < ny; j++) {
            if ((rx = wl+l0*dx[wl][j]) < el) el = rx;
            if ((rx = wl+l1*dx[wl][j]) < el) el = rx;
        }
    }
    wr = nx;
    er = INFINI;
    while (er > rmx) {
        wr--;
        if (wr < 0) {
            fprintf(stderr,"ERROR: getwwww: found velocities are out of bounds2\n");
            return(FAILURE);
        }
        er = rmx;
        for (j = 0;j < ny; j++) {
            if ((rx = wr+l0*dx[wr][j]) > er) er = rx;
            if ((rx = wr+l1*dx[wr][j]) > er) er = rx;
        }
    }
    wd = -1;
    ed = -INFINI;
    while (ed < 0) {
        wd++;
        if (wd > my) {
            fprintf(stderr,"ERROR: getwwww: found velocities are out of bound3\n");
            return(FAILURE);
        }
        ed = 0;
        for (i = 0; i < nx; i++) {
            if ((ry = wd+l0*dy[i][wd]) < ed) ed = ry;
            if ((ry = wd+l1*dy[i][wd]) < ed) ed = ry;
        }
    }
    wu = ny;
    eu = INFINI;
    while (eu > rmy) {
        wu--;
        if (wu < 0) {
            fprintf(stderr,"ERROR: getwwww: found velocities are out of bound4\n");
            return(FAILURE);
        }
        eu = rmy;
        for (i = 0; i < nx; i++) {
            if ((ry = wu+l0*dy[i][wu]) > eu) eu = ry;
            if ((ry = wu+l1*dy[i][wu]) > eu) eu = ry;
        }
    }
    if (wl > wr) {
        fprintf(stderr,"ERROR: getwwww: found velocities are out of bound5\n");
        return(FAILURE);
    }
    if (wd > wu) {
        fprintf(stderr,"ERROR: getwwww: found velocities are out of bound6\n");
        return(FAILURE);
    }
    *pwl=wl;
    *pwr=wr;
    *pwd=wd;
    *pwu=wu;
    return(SUCCESS);
}

int getwwww2(SIMAGE *sim, VIMAGE *vim, fdata m, int *pwl, int *pwr, int *pwd, int *pwu)
{
    /*---------------------------------------------------------------*/
    /*                          NOT OPTIMUM!                         */
    /*---------------------------------------------------------------*/
    /* Gets the parameters (border shrinks) of a window of a maximum */
    /* size so that the velocity field applied to it keeps it inside */
    /* the image size for at least two consecutive valid frames.     */
    /* The applied velocity field is considered when bringing        */
    /* images from planes from nf to nl into a plane m.              */
    /*---------------------------------------------------------------*/
    int i,j,f,wl,wr,wu,wd,nx,ny,mx,my,nf,nl,*vf,vv,co;
    fdata **dx,**dy,rmx,rmy;
    nf = sim->nf;
    nl = sim->nl;
    vf = sim->vf;
    nx = vim->nx;
    ny = vim->ny;
    mx = nx-1;
    my = ny-1;
    rmx = (fdata) mx;
    rmy = (fdata) my;
    dx = vim->dx;
    dy = vim->dy;
    wl = -1;
    vv = 1;
    while (vv) {
        wl++;
        if (wl > mx) {
            fprintf(stderr,"ERROR: getwwww2: found velocities are out of bounds1\n");
            return(FAILURE);
        }
        vv = 0;
        for (j = 0; j < ny; j++) {
            co = 0;
            for (f = nf; f <= nl; f++) {
                co += (vf[f] && ((wl+(f-m)*dx[wl][j]) >= 0.0));
            }
            if (co < 2) vv = 1;
        }
    }
    wr = nx;
    vv = 1;
    while (vv) {
        wr--;
        if (wr < 0) {
            fprintf(stderr,"ERROR: getwwww2: found velocities are out of bounds2\n");
            return(FAILURE);
        }
        vv = 0;
        for (j = 0; j < ny; j++) {
            co = 0;
            for (f = nf; f <= nl; f++) {
                co += (vf[f] && ((wr+(f-m)*dx[wr][j]) <= rmx));
            }
            if (co < 2) vv = 1;
        }
    }
    wd = -1;
    vv = 1;
    while (vv) {
        wd++;
        if (wd > my) {
            fprintf(stderr,"ERROR: getwwww2: found velocities are out of bound3\n");
            return(FAILURE);
        }
        vv = 0;
        for (i = 0; i < nx; i++) {
            co = 0;
            for (f = nf; f <= nl; f++) {
                co += (vf[f] && ((wd+(f-m)*dy[i][wd]) >= 0.0));
            }
            if (co < 2) vv = 1;
        }
    }
    wu = ny;
    vv = 1;
    while (vv) {
        wu--;
        if (wu < 0) {
            fprintf(stderr,"ERROR: getwwww2: found velocities are out of bound4\n");
            return(FAILURE);
        }
        vv = 0;
        for (i = 0; i < nx; i++) {
            co = 0;
            for (f = nf; f <= nl; f++) {
                co += (vf[f] && ((wu+(f-m)*dy[i][wu]) <= rmy));
            }
            if (co < 2) vv = 1;
        }
    }
    if (wl > wr) {
        fprintf(stderr,"ERROR: getwwww2: found velocities are out of bound5\n");
        return(FAILURE);
    }
    if (wd > wu) {
        fprintf(stderr,"ERROR: getwwww2: found velocities are out of bound6\n");
        return(FAILURE);
    }
    *pwl=wl;
    *pwr=wr;
    *pwd=wd;
    *pwu=wu;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*------------------------------ SMOOTHING PROCEDURES --------------------------------*/
/*------------------------------------------------------------------------------------*/

int setfilter(fdata **v, int n)
{
    /*--------------------------------------------------------------*/
    /* defines a convolution kernel for vector and matrix filtering */
    /*--------------------------------------------------------------*/
    int k;
    fdata *vec,pi=2*asin(1.0);
    if (altens1((void **) v,-n,n,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: setfilter: fail to allocate buffer\n");
        return(FAILURE);
    }
    vec = *v;
    vec[0] = 1.0;
    for (k = 1; k < n; k++) vec[-k] = vec[k] = (1+cos((pi*k)/n))/2;
    vec[-n] = vec[n] = 0.0;
    return(SUCCESS);
}

int smovec(fdata **pv, int smo, int nx, fdata *fil)
{
    /*-----------------------------------------------------*/
    /* performs a low pass filter                          */
    /* vo[i] = ((sum/k)(vi[i+k]*fil[k]))/((sum/k)(fil[k])) */
    /*-----------------------------------------------------*/
    int i,k,kmin,kmax;
    fdata *vi,*vo,s1,s0,filk,*xi;
    vi = *pv;
    if (malloc1((void **) &vo,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    for (i = 0; i < nx; i++) {
        s0 = s1 = 0;
        xi = vi+i;
        kmin = imax(1-smo,-i);
        kmax = imin(smo-1,(nx-1)-i);
        for (k = kmin; k <= kmax; k++) {
            filk = fil[k];
            s0 += filk;
            s1 += xi[k]*filk;
        }
        vo[i] = s1/s0;
    }
    if (free1((void **) &vi,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pv = vo;
    return(SUCCESS);
}

int smomat(fdata ***pm, int smo, int nx, int ny, fdata *fil)
{
    /*-----------------------------------------------------------*/
    /* performs a low pass horizontal filter                     */
    /* mo[i][j] = ((sum/k)(mi[i+k][j]*fil[k]))/((sum/k)(fil[k])) */
    /*-----------------------------------------------------------*/
    int i,j,k,kmin,kmax;
    fdata **mi,**mo,*moi,*mik,filk,mc;
    mi = *pm;
    if (malloc2((void ***) &mo,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    for (i = 0; i < nx; i++) {
        mc = 0;
        moi = mo[i];
        for (j = 0; j < ny; j++) moi[j] = 0;
        kmin = imax(1-smo,-i);
        kmax = imin(smo-1,(nx-1)-i);
        for (k = kmin; k <= kmax; k++) {
            filk = fil[k];
            mc += filk;
            mik = mi[i+k];
            for (j = 0;j < ny; j++) moi[j] += mik[j]*filk;
        }
        for (j = 0; j < ny; j++) moi[j] /= mc;
    }
    if (free2((void ***) &mi,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pm = mo;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int smovec1(fdata **pv, int smo, int nx)
{
    /*----------------------------*/
    /* performs a low pass filter */
    /*----------------------------*/
    int i,w;
    fdata *vi,*vo,*s1,*s0;
    w = smo/2;
    vi = *pv;
    if (altens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s1,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &vo,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    s0[-w] = 0.0;
    s1[-w] = 0.0;
    for (i = -w; i < 0; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i];
    }
    for (i = 0; i < w; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
    }
    for (i = w; i < nx; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
        vo[i-w] = (s1[i+1]-s1[i-2*w])/(s0[i+1]-s0[i-2*w]);
    }
    for (i = nx; i < nx+w; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i];
        vo[i-w] = (s1[i+1]-s1[i-2*w])/(s0[i+1]-s0[i-2*w]);
    }
    if (frtens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &s1,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free1((void **) &vi,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pv = vo;
    return(SUCCESS);
}

int smomat1(fdata ***pm, int smo, int nx, int ny)
{
    /*----------------------------*/
    /* performs a low pass filter */
    /*----------------------------*/
    int i,j,w;
    fdata **mi,**mo,**s1,**ss1,*s0;
    w = (smo < nx) ? smo/2 : nx/2;
    mi = *pm;
    if (malloc2((void ***) &mo,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s1,-w-1,nx+w,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (altens2((void ***) &ss1,-w-1,w,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    s0[-w] = 0.0;
    s1[-w-1] = ss1[-w-1];
    s1[-w] = ss1[-w];
    for (j = 0; j < ny; j++) {
        s1[-w][j] = 0.0;
    }
    for (i = -w; i < 0; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = ss1[i+1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j];
        }
    }
    for (i = 0; i < w; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = ss1[i+1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j]+mi[i][j];
        }
    }
    for (i = w; i < nx; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i-2*w-1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j]+mi[i][j];
            mo[i-w][j] = (s1[i+1][j]-s1[i-2*w][j])/(s0[i+1]-s0[i-2*w]);
        }
    }
    for (i = nx; i < nx+w; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i-2*w-1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j];
            mo[i-w][j] = (s1[i+1][j]-s1[i-2*w][j])/(s0[i+1]-s0[i-2*w]);
        }
    }
    if (frtens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &s1,-w-1,nx+w,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (frtens2((void ***) &ss1,-w-1,w,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free2((void ***) &mi,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pm = mo;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int smovec2(fdata **pv, int smo, int nx)
{
    /*----------------------------*/
    /* performs a low pass filter */
    /*----------------------------*/
    int i,w0,w1;
    fdata *vi,*vo,*s1,*s0,*sl1,*sl0,*sr1,*sr0;
    w1 = smo/2;
    w0 = smo/2;
    vi = *pv;
    if (altens1((void **) &s0,-w0,nx+w0,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sl0,-w1,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sr0,0,nx+w1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s1,-w0,nx+w0,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sl1,-w1,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sr1,0,nx+w1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &vo,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    s0[-w0] = 0.0;
    s1[-w0] = 0.0;
    for (i = -w0; i < 0; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i];
    }
    for (i = 0; i < w0; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
    }
    for (i = w0; i < nx; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
        vo[i-w0] = (s1[i+1]-s1[i-2*w0])/(s0[i+1]-s0[i-2*w0]);
    }
    for (i = nx; i < nx+w0; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i];
        vo[i-w0] = (s1[i+1]-s1[i-2*w0])/(s0[i+1]-s0[i-2*w0]);
    }
    sr0[0] = sl0[-w1] = 0.0;
    sr1[0] = sl1[-w1] = 0.0;
    for (i = -w1; i < 0; i++) {
        sl0[i+1] = sl0[i];
        sl1[i+1] = sl1[i];
    }
    for (i = 0; i < w1; i++) {
        sl0[i+1] = sl0[i]+1.0;
        sr0[i+1] = sr0[i]+sl0[i+1]-sl0[i-w1];
        sl1[i+1] = sl1[i]+vo[i];
        sr1[i+1] = sr1[i]+sl1[i+1]-sl1[i-w1];
    }
    for (i = w1; i < nx; i++) {
        sl0[i+1] = sl0[i]+1.0;
        sr0[i+1] = sr0[i]+sl0[i+1]-sl0[i-w1];
        sl1[i+1] = sl1[i]+vo[i];
        sr1[i+1] = sr1[i]+sl1[i+1]-sl1[i-w1];
        vi[i-w1] = (sr1[i+1]-sr1[i-w1])/(sr0[i+1]-sr0[i-w1]);
    }
    for (i = nx; i < nx+w1; i++) {
        sr0[i+1] = sr0[i];
        sr1[i+1] = sr1[i];
        vi[i-w1] = (sr1[i+1]-sr1[i-w1])/(sr0[i+1]-sr0[i-w1]);
    }
    if (frtens1((void **) &s0,-w0,nx+w0,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sl0,-w1,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sr0,0,nx+w1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &s1,-w0,nx+w0,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sl1,-w1,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sr1,0,nx+w1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free1((void **) &vo,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pv = vi;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int smovec3(fdata **pv, int smo, int nx)
{
    /*----------------------------*/
    /* performs a low pass filter */
    /*----------------------------*/
    int i,w;
    fdata *vi,*vo,*s1,*s0,*sl1,*sl0,*sr1,*sr0;
    w = (smo < nx) ? smo/2 : nx/2;
    vi = *pv;
    if (altens1((void **) &sl0,-w,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sl1,-w,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sr0,0,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sr1,0,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s1,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &vo,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    sr0[0] = sl0[-w] = s0[-w] = 0.0;
    sr1[0] = sl1[-w] = s1[-w] = 0.0;
    for (i = -w; i < 0; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i];
    }
    for (i = 0; i < w; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
        sl0[i-w+1] = sl0[i-w];
        sl1[i-w+1] = sl1[i-w];
    }
    for (i = w; i < 2*w; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
        sl0[i-w+1] = sl0[i-w]+s0[i+1]-s0[i-2*w];
        sl1[i-w+1] = sl1[i-w]+s1[i+1]-s1[i-2*w];
        sr0[i-w+1] = sr0[i-w]+sl0[i-w+1]-sl0[i-2*w];
        sr1[i-w+1] = sr1[i-w]+sl1[i-w+1]-sl1[i-2*w];
    }
    for (i = 2*w; i < nx; i++) {
        s0[i+1] = s0[i]+1.0;
        s1[i+1] = s1[i]+vi[i];
        sl0[i-w+1] = sl0[i-w]+s0[i+1]-s0[i-2*w];
        sl1[i-w+1] = sl1[i-w]+s1[i+1]-s1[i-2*w];
        sr0[i-w+1] = sr0[i-w]+sl0[i-w+1]-sl0[i-2*w];
        sr1[i-w+1] = sr1[i-w]+sl1[i-w+1]-sl1[i-2*w];
        vo[i-2*w] = (sr1[i-w+1]-sr1[i-2*w])/(sr0[i-w+1]-sr0[i-2*w]);
    }
    for (i = nx; i < nx+w; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = s1[i];
        sl0[i-w+1] = sl0[i-w]+s0[i+1]-s0[i-2*w];
        sl1[i-w+1] = sl1[i-w]+s1[i+1]-s1[i-2*w];
        sr0[i-w+1] = sr0[i-w]+sl0[i-w+1]-sl0[i-2*w];
        sr1[i-w+1] = sr1[i-w]+sl1[i-w+1]-sl1[i-2*w];
        vo[i-2*w] = (sr1[i-w+1]-sr1[i-2*w])/(sr0[i-w+1]-sr0[i-2*w]);
    }
    for (i = nx+w; i < nx+2*w; i++) {
        sr0[i-w+1] = sr0[i-w];
        sr1[i-w+1] = sr1[i-w];
        vo[i-2*w] = (sr1[i-w+1]-sr1[i-2*w])/(sr0[i-w+1]-sr0[i-2*w]);
    }
    if (frtens1((void **) &sl0,-w,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sl1,-w,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sr0,0,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sr1,0,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &s1,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free1((void **) &vi,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pv = vo;
    return(SUCCESS);
}

int smomat3(fdata ***pm, int smo, int nx, int ny)
{
    /*----------------------------*/
    /* performs a low pass filter */
    /*----------------------------*/
    int i,j,w;
    fdata **mi,**mo,**s1,**ss1,*s0,**sl1,**ssl1,*sl0,**sr1,**ssr1,*sr0;
    w = (smo < nx) ? smo/2 : nx/2;
    mi = *pm;
    if (malloc2((void ***) &mo,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sl0,-w,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sr0,0,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens1((void **) &s1,-w-1,nx+w,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sl1,-w-1,nx,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (altens1((void **) &sr1,-1,nx+w,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (altens2((void ***) &ss1,-w-1,w,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens2((void ***) &ssl1,-w-1,0,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (altens2((void ***) &ssr1,-1,w,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    sr0[0] = sl0[-w] = s0[-w] = 0.0;
    s1[-w-1] = ss1[-w-1];
    sl1[-w-1] = ssl1[-w-1];
    sr1[-1] = ssr1[-1];
    s1[-w] = ss1[-w];
    sl1[-w] = ssl1[-w];
    sr1[0] = ssr1[0];
    for (j = 0; j < ny; j++) {
        sr1[0][j] = sl1[-w][j] = s1[-w][j] = 0.0;
    }
    for (i = -w; i < 0; i++) {
        s0[i+1] = s0[i];
        s1[i+1] = ss1[i+1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j];
        }
    }
    for (i = 0; i < w; i++) {
        s0[i+1] = s0[i]+1.0;
        sl0[i-w+1] = sl0[i-w];
        s1[i+1] = ss1[i+1];
        sl1[i-w+1] = ssl1[i-w+1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j]+mi[i][j];
            sl1[i-w+1][j] = sl1[i-w][j];
        }
    }
    for (i = w; i < 2*w; i++) {
        s0[i+1] = s0[i]+1.0;
        sl0[i-w+1] = sl0[i-w]+s0[i+1]-s0[i-2*w];
        sr0[i-w+1] = sr0[i-w]+sl0[i-w+1]-sl0[i-2*w];
        s1[i+1] = s1[i-2*w-1];
        sl1[i-w+1] = sl1[i-2*w-1];
        sr1[i-w+1] = ssr1[i-w+1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j]+mi[i][j];
            sl1[i-w+1][j] = sl1[i-w][j]+s1[i+1][j]-s1[i-2*w][j];
            sr1[i-w+1][j] = sr1[i-w][j]+sl1[i-w+1][j]-sl1[i-2*w][j];
        }
    }
    for (i = 2*w; i < nx; i++) {
        s0[i+1] = s0[i]+1.0;
        sl0[i-w+1] = sl0[i-w]+s0[i+1]-s0[i-2*w];
        sr0[i-w+1] = sr0[i-w]+sl0[i-w+1]-sl0[i-2*w];
        s1[i+1] = s1[i-2*w-1];
        sl1[i-w+1] = sl1[i-2*w-1];
        sr1[i-w+1] = sr1[i-2*w-1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j]+mi[i][j];
            sl1[i-w+1][j] = sl1[i-w][j]+s1[i+1][j]-s1[i-2*w][j];
            sr1[i-w+1][j] = sr1[i-w][j]+sl1[i-w+1][j]-sl1[i-2*w][j];
            mo[i-2*w][j] = (sr1[i-w+1][j]-sr1[i-2*w][j])/(sr0[i-w+1]-sr0[i-2*w]);
        }
    }
    for (i = nx; i < nx+w; i++) {
        s0[i+1] = s0[i];
        sl0[i-w+1] = sl0[i-w]+s0[i+1]-s0[i-2*w];
        sr0[i-w+1] = sr0[i-w]+sl0[i-w+1]-sl0[i-2*w];
        s1[i+1] = s1[i-2*w-1];
        sl1[i-w+1] = sl1[i-2*w-1];
        sr1[i-w+1] = sr1[i-2*w-1];
        for (j = 0; j < ny; j++) {
            s1[i+1][j] = s1[i][j];
            sl1[i-w+1][j] = sl1[i-w][j]+s1[i+1][j]-s1[i-2*w][j];
            sr1[i-w+1][j] = sr1[i-w][j]+sl1[i-w+1][j]-sl1[i-2*w][j];
            mo[i-2*w][j] = (sr1[i-w+1][j]-sr1[i-2*w][j])/(sr0[i-w+1]-sr0[i-2*w]);
        }
    }
    for (i = nx+w; i < nx+2*w; i++) {
        sr0[i-w+1] = sr0[i-w];
        sr1[i-w+1] = sr1[i-2*w-1];
        for (j = 0; j < ny; j++) {
            sr1[i-w+1][j] = sr1[i-w][j];
            mo[i-2*w][j] = (sr1[i-w+1][j]-sr1[i-2*w][j])/(sr0[i-w+1]-sr0[i-2*w]);
        }
    }
    if (frtens1((void **) &s0,-w,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sl0,-w,nx,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sr0,0,nx+w,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &s1,-w-1,nx+w,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sl1,-w-1,nx,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (frtens1((void **) &sr1,-1,nx+w,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (frtens2((void ***) &ss1,-w-1,w,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens2((void ***) &ssl1,-w-1,0,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (frtens2((void ***) &ssr1,-1,w,0,ny-1,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free2((void ***) &mi,nx,ny,sizeof(fdata)) == FAILURE) return(FAILURE);
    *pm = mo;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*------------------------------- DYNAMIC PROGRAMMING --------------------------------*/
/*------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/*                          Ln distance matrix calculation                            */
/*------------------------------------------------------------------------------------*/

int xcaldist(fdata **dd, DIMAGE *dim0, DIMAGE *dim1, int na, int km, int kp, int kz,
int wm, int wp, fdata *dfilter, fdata dlf, int dni)
{   
    int nb,l,sl,dl,l0,l1,ky,slmin,slmax,dlmin,dlmax,dli,dni2;
    fdata ***xm0,***xm1,**xml0,**xml1,dist0,dist1,dist2,dist3,filk;
    uchar **in0,**in1,*inl0,*inl1;
    if (dlf == 0.5) dli = 0;
    else if (dlf == 1.0) dli = 1;
    else if (dlf == 2.0) dli = 2;
    else dli = 3;
    dni2 = 2-dni;
    slmin = 2*wm;
    slmax = 2*wp;
    dlmin = -na-1;
    dlmax = na+1;
    for (sl = slmin; sl <= slmax; sl++) {
        for (dl = dlmin; dl <= dlmax; dl++) {
            dd[sl][dl] = INFINI;
        }
    }
    nb = dim0->nb;
    xm0 = dim0->xm;
    in0 = dim0->in;
    xm1 = dim1->xm;
    in1 = dim1->in;
    if (malloc1((void **) &xml0,nb,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &xml1,nb,sizeof(fdata *)) == FAILURE) return(FAILURE);
    for (sl = slmin; sl <= slmax; sl++) {
        for (dl = ((sl+na)&(!dni))-na; dl <= na; dl += dni2) {
            l0 = sl-dl;
            l1 = dl+sl;
            if (!dni) {
                l0 /= 2;
                l1 /= 2;
            }
            inl0 = in0[l0];
            inl1 = in1[l1];
            for (l = 0; l < nb; l++) {
                xml0[l] = xm0[l][l0];
                xml1[l] = xm1[l][l1];
            }
            dist0 = dist1 = 0.0;
            for (ky = km; ky <= kp; ky++) {
                filk = dfilter[ky-kz];
                if ((inl0[ky]) && (inl1[ky])) {
                    dist2 = 0.0;
                    for (l = 0; l < nb; l++) {
                        dist3 = xml0[l][ky]-xml1[l][ky];
                        switch(dli) {
                            case 0: {
                                dist2 += sqrt((dist3 > 0) ? dist3 : -dist3);
                                break;
                            }
                            case 1: {
                                dist2 += (dist3 > 0) ? dist3 : -dist3;
                                break;
                            }
                            case 2: {
                                dist2 += dist3*dist3;
                                break;
                            }
                            case 3: {
                                dist2 += pow(((dist3 > 0) ? dist3 : -dist3),dlf);
                                break;
                            }
                            default: {
                            }
                        }
                    }
                    dist0 += filk;
                    dist1 += filk*dist2;
                }
            }
            dd[sl][dl] = 1.0e-20+dist1/dist0;
        }
    }
    if (free1((void **) &xml0,nb,sizeof(fdata *)) == FAILURE) return(FAILURE);
    if (free1((void **) &xml1,nb,sizeof(fdata *)) == FAILURE) return(FAILURE);
    return(SUCCESS);
}

int xcaldista(fdata **dd, VIMAGE *vim, SIMAGE *sim, int dir, int sum, fdata m,
fdata res, int na, int km, int kp, int kz, int wm, int wp, fdata *dfilter, fdata dlf,
int dni, int bsi)
{   
    int nx,ny,nb,nf,nl,f,l,sl,dl,ky,*vf,slm,slp,slmin,slmax,wx,wy;
    int prev,curr,dli,dni2;
    fdata dist0,dist1,dist2,dist3,filk,rx,ry,rm,rp,rr,rl,rdl;
    fdata *rim,*pim,*immin,*immax,**dx,**dy,m2;
    uchar **im;
    if (dlf == 0.5) dli = 0;
    else if (dlf == 1.0) dli = 1;
    else if (dlf == 2.0) dli = 2;
    else dli = 3;
    dni2 = 2-dni;
    slmin = 2*wm;
    slmax = 2*wp;
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    wx = sim->nx;
    wy = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    im = sim->im;
    vf = sim->vf;
    m2 = 0.5*(nf+nl);
    if (malloc1((void **) &rim,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &pim,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &immin,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (malloc1((void **) &immax,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    for (sl = slmin; sl <= slmax; sl++) {
        dd[sl][-na-1] = INFINI;
        rl = 0.5*sl;
        for (dl = ((sl+na)&(!dni))-na; dl <= na; dl += dni2) {
            dist0 = dist1 = 0.0;
            rdl = res*dl;
            for (ky = km; ky <= kp; ky++) {
                filk = dfilter[ky-kz];
                prev = 0;
                for (l = 0; l < nb; l++) {
                    immin[l] = INFINI;
                    immax[l] = -INFINI;
                }
                for (f = nf; f <= nl; f++) {
                    if (vf[f]) {
                        slm = rr = rl+(f-m2)*rdl;
                        if (slm == rr) {
                            rx = (f-m)*dx[slm][ky]+rr;
                            ry = (f-m)*dy[slm][ky]+ky;
                        } else {
                            slp = slm+(slm < wp);
                            rp = rr-slm;
                            rm = 1.0-rp;
                            rx = (f-m)*(rm*dx[slm][ky]+rp*dx[slp][ky])+rr;
                            ry = (f-m)*(rm*dy[slm][ky]+rp*dy[slp][ky])+ky;
                        }
                        curr = dir ? bsi ? inter2(ry,rx,wx,wy,nb,im[f],rim)
                                         : inter1(ry,rx,wx,wy,nb,im[f],rim)
                                   : bsi ? inter2(rx,ry,wx,wy,nb,im[f],rim)
                                         : inter1(rx,ry,wx,wy,nb,im[f],rim);
                        if (curr) {
                            if (sum) {
                                if (prev) {
                                    dist2 = 0;
                                    for (l = 0; l < nb; l++) {
                                        dist3 = pim[l]-rim[l];
                                        switch(dli) {
                                            case 0: {
                                                dist2 += sqrt((dist3 > 0) ? dist3 : -dist3);
                                                break;
                                            }
                                            case 1: {
                                                dist2 += (dist3 > 0) ? dist3 : -dist3;
                                                break;
                                            }
                                            case 2: {
                                                dist2 += dist3*dist3;
                                                break;
                                            }
                                            case 3: {
                                                dist2 += pow(((dist3 > 0) ? dist3 : -dist3),dlf);
                                                break;
                                            }
                                            default: {
                                            }
                                        }
                                    }
                                    dist0 += filk;
                                    dist1 += filk*dist2;
                                }
                                for (l = 0; l < nb; l++) {
                                    pim[l] = rim[l];
                                }
                            } else {
                                for (l = 0; l < nb; l++) {
                                    if (rim[l] < immin[l]) immin[l] = rim[l];
                                    if (rim[l] > immax[l]) immax[l] = rim[l];
                                }
                            }
                            prev += 1;
                        }
                    }
                }
                if (!sum) {
                    dist2 = 0;
                    for (l = 0; l < nb; l++) {
                        dist3 = immax[l]-immin[l];
                        switch(dli) {
                            case 0: {
                                dist2 += sqrt((dist3 > 0) ? dist3 : -dist3);
                                break;
                            }
                            case 1: {
                                dist2 += (dist3 > 0) ? dist3 : -dist3;
                                break;
                            }
                            case 2: {
                                dist2 += dist3*dist3;
                                break;
                            }
                            case 3: {
                                dist2 += pow(((dist3 > 0) ? dist3 : -dist3),dlf);
                                break;
                            }
                            default: {
                            }
                        }
                    }
                    if (prev >= 2) {
                        dist0 += filk;
                        dist1 += filk*dist2;
                    }
                }
            }
            dd[sl][dl] = 1.0e-20+dist1/dist0;
        }
        dd[sl][na+1] = INFINI;
    }
    if (free1((void **) &rim,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free1((void **) &pim,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free1((void **) &immin,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    if (free1((void **) &immax,nb,sizeof(fdata)) == FAILURE) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*       Dynamic progamming inside the search window around the matrix diagonal       */
/*------------------------------------------------------------------------------------*/

void monodp(int na, int wm, int wp, fdata eps, fdata **dd, fdata **gg, int **bb,
SIMAGE *sim, int sc, float dlf, int dni)
{
    int b2,sl,dl,slmin,slmax,dlmin,dlmax,*bbm0,n,nb;
    fdata g1,g2,g3,cc,*ggm3,*ggm2,*ggm1,*ggm0,*ddm3,*ddm2,*ddm1,*ddm0,ddmax,dist0,dist1;
    fdata *val;
    slmin = 2*wm;
    slmax = 2*wp;
    if (!dni) {
        for (sl = slmin; sl <= slmax; sl++) {
            for (dl = ((sl+na+1)&1)-na; dl <= na; dl += 2) {
                dist1 = dist0 = 0.0;
                if ((sl-1) >= slmin) {
                    dist1 += dd[sl-1][dl];
                    dist0 += 1.0;
                }
                if ((sl+1) <= slmax) {
                    dist1 += dd[sl+1][dl];
                    dist0 += 1.0;
                }
                if ((dl-1) >= -na) {
                    dist1 += dd[sl][dl-1];
                    dist0 += 1.0;
                }
                if ((dl+1) <= na) {
                    dist1 += dd[sl][dl+1];
                    dist0 += 1.0;
                }
                dd[sl][dl] = dist1/dist0;
            }
        }
    }
    nb = sim->nb;
    val = sim->val;
    ddmax = 0.0;
    for (n = 0; n < nb; n++) ddmax += pow((double) (0.4*val[n]),dlf);
    cc = eps*ddmax;
    dlmin = -na-1;
    dlmax = na+1;
    for (sl = slmin; sl <= slmax; sl++) {
        for (dl = dlmin; dl <= dlmax; dl++) {
            gg[sl][dl] = INFINI;
            bb[sl][dl] = -1;
        }
    }
    for (dl = -na; dl <= na; dl++) {
        gg[slmin][dl] = 0;
        bb[slmin][dl] = 0;
        if (sc > 1) {
            gg[slmin+1][dl] = gg[slmin][dl]+dd[slmin][dl]+dd[slmin+1][dl];
            bb[slmin+1][dl] = 2;
        }
        if (sc > 2) {
            gg[slmin+2][dl] = gg[slmin+1][dl]+dd[slmin+1][dl]+dd[slmin+2][dl];
            bb[slmin+2][dl] = 2;
        }
    }
    for (sl = slmin+sc; sl <= slmax; sl++) {
        bbm0 = bb[sl];
        ggm3 = gg[sl-3];
        ggm2 = gg[sl-2];
        ggm1 = gg[sl-1];
        ggm0 = gg[sl];
        ddm3 = dd[sl-3];
        ddm2 = dd[sl-2];
        ddm1 = dd[sl-1];
        ddm0 = dd[sl];
        for (dl = -na; dl <= na; dl++) {
            b2 = 2;
            g2 = ggm1[dl]+ddm1[dl]+ddm0[dl];
            switch(sc) {
                case 1: {
                    g1 = ddm1[dl-1]+ddm0[dl]+cc+ggm1[dl-1];
                    g3 = ddm1[dl+1]+ddm0[dl]+cc+ggm1[dl+1];
                    break;
                }
                case 2: {
                    g1 = ddm2[dl-1]+ddm1[dl-1]+ddm1[dl]+ddm0[dl]+cc+ggm2[dl-1];
                    g3 = ddm2[dl+1]+ddm1[dl+1]+ddm1[dl]+ddm0[dl]+cc+ggm2[dl+1];
                    break;
                }
                case 3: {
                    g1 = ddm3[dl-1]+1.333333*ddm2[dl-1]+0.666667*ddm2[dl]+
                         0.666667*ddm1[dl-1]+1.333333*ddm1[dl]+ddm0[dl]+cc+ggm3[dl-1];
                    g3 = ddm3[dl+1]+1.333333*ddm2[dl+1]+0.666667*ddm2[dl]+
                         0.666667*ddm1[dl+1]+1.333333*ddm1[dl]+ddm0[dl]+cc+ggm3[dl+1];
                    break;
                }
                default: {
                    g1 = g3 = INFINI;
                }
            }
            if (g1 < g2) {
                g2 = g1;
                b2 = 1;
            }
            if (g3 < g2) {
                g2 = g3;
                b2 = 3;
            }
            ggm0[dl] = g2;
            bbm0[dl] = b2;
        }
    }
}

/*------------------------------------------------------------------------------------*/
/*               Backtrack of the optimal path inside the search window               */
/*------------------------------------------------------------------------------------*/

int backtrack(fdata *adjk, int na, int nx, int wp, fdata lambda, fdata res, fdata **gg,
int **bb, int sc)
{
    int l,sl,dl,dlmin;
    fdata rmin,rval,l1,l0,dl1,dl0;
    sl = 2*wp;
    rmin = gg[sl][0];
    dlmin = 0;
    for (dl = 1; dl <= na; dl++) {
        if ((rval = gg[sl][dl]) < rmin) {
            rmin = rval;
            dlmin = dl;
        }
        if ((rval = gg[sl][-dl]) < rmin) {
            rmin = rval;
            dlmin = -dl;
        }
    }
    dl = dlmin;
    l0 = 0.5*sl+(lambda-0.5)*dl;
    dl0 = dl;
    for(l = nx-1; l >= l0; l--) adjk[l] = res*dl0;
    while (bb[sl][dl]) {
        switch(bb[sl][dl]) {
            case 1: {
                sl -= sc;
                dl -= 1;
                break;
            }
            case 2: {
                sl -= 1;
                break;
            }
            case 3: {
                sl -= sc;
                dl += 1;
                break;
            }
            default: {
                fprintf(stderr,"\nError in switch backtrack\n");
                return(FAILURE);
            }
        }
        l1 = l0;
        dl1 = dl0;
        l0 = 0.5*sl+(lambda-0.5)*dl;
        if (l0 < 0.0) l0 = 0.0;
        dl0 = dl;
        for(; l >= l0; l--) adjk[l] = res*(dl0+(l-l0)*(dl1-dl0)/(l1-l0));
    }
    for(; l >= 0; l--) adjk[l] = res*dl0;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int insw(SIMAGE *sim, fdata m, fdata dx, fdata ml)
{
    int f,nf,nl,*vf,co;
    if (ml < 0.0) return(0);
    nf = sim->nf;
    nl = sim->nl;
    vf = sim->vf;
    co = 0;
    for (f = nf; f <= nl; f++) {
        co += (vf[f] && ((f-m)*dx <= ml));
    }
    return(co >= 2);
}

int xopflow(SIMAGE *sim, VIMAGE *vim, DPCONTROL *dpcontrol, int n)
{
    int i,j,k,km,kp,ky,l,nj,na,nl,smo,nw,nx,ny,nb,wl,wr,wd,wu,wm,wp,ne,ql,**bb,dir;
    int ff,lf,mx,my,mw,sum,wx,sc,ar,dni,bsi,smf;
    fdata *sfilter,*dfilter,**dd,**gg,**adj,**dx,**dy,r0,r1,eps,m,res,rna,dlf;
    fdata dlin,dsig,xavr,xval,hih,hil,lambda,scale,tau,spwfac,winfac,smofac;
    DIMAGE dim0,dim1;
    nx = vim->nx;
    ny = vim->ny;
    nb = sim->nb;
    ff = sim->nf;
    lf = sim->nl;
    bsi = dpcontrol->bsi;
    dni = dpcontrol->dni;
    eps = dpcontrol->eps;
    tau = dpcontrol->tau;
    spwfac = dpcontrol->spwfac;
    winfac = dpcontrol->winfac;
    smofac = dpcontrol->smofac;
    scale = exp(tau*n);
    nl = floor(1.0+0.5*spwfac*scale);
    na = floor(1.5+0.5*winfac*spwfac*scale);
    smo = floor(1.5+0.5*smofac*spwfac*scale);
    smf = dpcontrol->smf;
    res = (n >= 0) ? 1.0 : scale;
    rna = 0.5*res*na;
    eps *= sqrt(res);
    res /= lf-ff;
    nw = nl+dpcontrol->enw;
    dlf = dpcontrol->dlf;
    m  = dpcontrol->m;
    dir = dpcontrol->id;
    dlin = nx*dpcontrol->dlin[dir];
    dsig = nx*dpcontrol->dsig[dir];
    sc = dpcontrol->sc;
    ar =  (dpcontrol->mf > 2) || (n < 0);
    sum = dpcontrol->sum;
    lambda = (m-(fdata)ff)/((fdata)(lf-ff));
    mx = nx-1;
    my = ny-1;
    mw = nw-1;
    wx = 2*mx;
    if (setfilter(&sfilter,smo) == FAILURE) return(FAILURE);
    if (setfilter(&dfilter,nw) == FAILURE) return(FAILURE);
    if (getwwww(sim,vim,m,&wl,&wr,&wd,&wu) == FAILURE) {
        fprintf(stderr,"ERROR: xopflow: fail to get window parameters\n");
        return(FAILURE);
    }
    if (dpcontrol->vb) {
        if (dir) {
            printf("   vertical: %2d %2d %2d %2d",wd,(ny-1)-wu,wl,(nx-1)-wr);
            fflush(stdout);
        } else {
            printf("   horizontal: %2d %2d %2d %2d",wl,(nx-1)-wr,wd,(ny-1)-wu);
            fflush(stdout);
        }
    }
    if (getwwww2(sim,vim,m,&wl,&wr,&wd,&wu) == FAILURE) {
        fprintf(stderr,"ERROR: xopflow: fail to get window parameters\n");
        return(FAILURE);
    }
    dx = vim->dx;
    dy = vim->dy;
    if (dlin < nx) {
        xavr = 0;
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                xavr += dx[i][j];
            }
        }
        xavr /= ny*nx;
        hil = xavr-dlin;
        hih = xavr+dlin;
    } else {
        hil = hih = 0;
    }
    ne = (wu-wd)/nl;
    ql = ((wu-wd)-nl*ne)/2;
    ne += 1;
    if ((!ar && (getimage(dir,sim,vim,m,&dim0,&dim1,dpcontrol->goa,dni,bsi) == FAILURE))
     || (altens2((void ***) &dd,0,wx,-na-1,na+1,sizeof(fdata)) == FAILURE)
     || (altens2((void ***) &gg,0,wx,-na-1,na+1,sizeof(fdata)) == FAILURE)
     || (altens2((void ***) &bb,0,wx,-na-1,na+1,sizeof(int)) == FAILURE)
     || (malloc2((void ***) &adj,ne,nx,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: xopflow: fail to allocate buffers\n");
        return(FAILURE);
    }
    for (k = 0; k < ne; k++) {
        ky = k*nl+wd+ql;
        wm = l = (wl+wr)/2;
        l = (wl+wr)/2;
        while ((l >= 0) && (insw(sim,m,-dx[l][ky],l-rna))) wm = l--;
        wp = l = (wl+wr)/2;
        l = (wl+wr)/2;
        while ((l < nx) && (insw(sim,m,dx[l][ky],mx-l-rna))) wp = l++;
        km = imax(0,ky-mw);
        kp = imin(my,ky+mw);
        if ((ar ? 
            xcaldista(dd,vim,sim,dir,sum,m,res,na,km,kp,ky,wm,wp,dfilter,dlf,dni,bsi) :
            xcaldist(dd,&dim0,&dim1,na,km,kp,ky,wm,wp,dfilter,dlf,dni)) == FAILURE) {
            fprintf(stderr,"ERROR: xopflow: fail to compute distance matrix\n");
            return(FAILURE);
        }
        monodp(na,wm,wp,eps,dd,gg,bb,sim,sc,dlf,dni);
        if ((backtrack(adj[k],na,nx,wp,lambda,res,gg,bb,sc) == FAILURE)) {
            fprintf(stderr,"ERROR: xopflow: fail to backtrack\n");
            return(FAILURE);
        }
        if (((smf == 3) && (smovec3(adj+k,smo,nx) != FAILURE))
         || ((smf == 2) && (smovec2(adj+k,smo,nx) != FAILURE))
         || ((smf == 1) && (smovec2(adj+k,smo,nx) != FAILURE))
         || (smovec(adj+k,smo,nx,sfilter) != FAILURE)) {
        } else {
            fprintf(stderr,"ERROR: xopflow: fail to smooth vector\n");
            return(FAILURE);
        }
    }
    for (j = 0; j < ny; j++) {
        nj = j;
        if (j < wd) nj = wd;
        if (j > wu) nj = wu;
        k = ((nj-wd)+(nl-ql))/nl-1;
        if (k < 0) k = 0;
        if (k > (ne-2)) k = ne-2;
        l = (nj-wd)-(ql+k*nl);
        r1 = l;
        r1 /= nl;
        r0 = 1-r1;
        for (i = 0; i < nx; i++) dx[i][j] += r0*adj[k][i]+r1*adj[k+1][i];
    }
    if ((!ar && ((free2((void ***) &dim0.in,dim0.nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim0.xm,nb,dim0.nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim1.in,dim1.nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim1.xm,nb,dim1.nx,ny,sizeof(fdata)) == FAILURE)))
     || (frtens2((void ***) &dd,0,wx,-na-1,na+1,sizeof(fdata)) == FAILURE)
     || (frtens2((void ***) &gg,0,wx,-na-1,na+1,sizeof(fdata)) == FAILURE)
     || (frtens2((void ***) &bb,0,wx,-na-1,na+1,sizeof(int)) == FAILURE)
     || (free2((void ***) &adj,ne,nx,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: xopflow: fail to free buffers\n");
        return(FAILURE);
    }
    if (dlin < nx) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
            if (dx[i][j] < hil) {
                    xval = dx[i][j]-hil;
                    dx[i][j] = hil+(xval*dsig)/(dsig-xval);
                }
            if (dx[i][j] > hih) {
                    xval = dx[i][j]-hih;
                    dx[i][j] = hih+(xval*dsig)/(dsig+xval);
                }
            }
        }
    }
    xyswap(vim);
    dpcontrol->id = 1-dir;
    if (((smf == 3) && (smomat3(&vim->dy,smo,ny,nx) != FAILURE))
     || ((smf == 1) && (smomat1(&vim->dy,smo,ny,nx) != FAILURE))
     || (smomat(&vim->dy,smo,ny,nx,sfilter) != FAILURE)) {
    } else {
        fprintf(stderr,"ERROR: xopflow: fail to smooth matrix\n");
        return(FAILURE);
    }
    vim->d[1] = vim->dy;
    if ((frtens1((void **) &dfilter,-nw,nw,sizeof(fdata)) == FAILURE)
     || (frtens1((void **) &sfilter,-smo,smo,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: xopflow: fail to free buffers\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

fdata maxvr(VIMAGE *vim)
{
    fdata x,y,v,vmax,**dx,**dy,mx,my,mxx,myy,m,mi,mj,mii,mjj,ax,bx,ay,by,xi,xj;
    int i,j,nx,ny;
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    mx = my = mxx = myy = m = mi = mj = mii = mjj = 0.0;
    for (j = 0; j < ny; j++) {
        xj = (fdata) j;
        for (i = 0; i < nx; i++) {
            xi = (fdata) i;
            m += 1;
            mi += xi;
            mj += xj;
            mii += xi*xi;
            mjj += xj*xj;
            mx += dx[i][j];
            my += dy[i][j];
            mxx += i*dx[i][j];
            myy += j*dy[i][j];
        }
    }
    ax = (mx*mi-mxx*m)/(mi*mi-mii*m);
    ay = (my*mj-myy*m)/(mj*mj-mjj*m);
    bx = (mi*mxx-mii*mx)/(mi*mi-mii*m);
    by = (mj*myy-mjj*my)/(mj*mj-mjj*m);
    vmax = 0;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            x = dx[i][j]-(ax*i+bx);
            y = dy[i][j]-(ay*j+by);
            v = x*x+y*y;
            if (vmax < v) vmax = v;
        }
    }
    return(sqrt(vmax/((fdata) (nx*nx+ny*ny))));
}

fdata maxva(VIMAGE *vim)
{
    int i,j,nx,ny;
    fdata x,y,v,vmax,**dx,**dy;
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    vmax = 0;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            x = dx[i][j];
            y = dy[i][j];
            v = x*x+y*y;
            if (vmax < v) vmax = v;
        }
    }
    return(sqrt(vmax/(nx*nx+ny*ny)));
}

int oiter(SIMAGE *sim, VIMAGE *vim, DPCONTROL *dpcontrol, int i, int itc, int nf0,
int nl0, int *pa)
{
    int nf,nl,f,mf;
    float ff,md;
    nf = sim->nf;
    nl = sim->nl;
    mf = dpcontrol->mf;
    md = dpcontrol->md;
    if (mf > 2) md *= 2.0;
    for (f = nf; f <= nl; f++) sim->vf[f] = 0;
    for (ff = 0.0; ff < mf; ff += 1.0) {
        sim->vf[(int) (floor(0.5+nf0+(ff/(mf-1))*(nl0-nf0)))] = 1;
    }
    sim->nf = nf0;
    sim->nl = nl0;
    if (dpcontrol->vb) {
        printf("Orthogonal iteration %2d:",itc);
        fflush(stdout);
    }
    if (xopflow(sim,vim,dpcontrol,i) == FAILURE) return(FAILURE);
    if (xopflow(sim,vim,dpcontrol,i) == FAILURE) return(FAILURE);
    if (dpcontrol->vb) printf("   %d\n",nl0-nf0);
    for (f = nf; f <= nl; f++) sim->vf[f] = 1;
    sim->nf = nf;
    sim->nl = nl;
    if (dpcontrol->vs == 2) *pa = ((((fdata) (nl0+1-nf0))*maxvr(vim)) < md);
    if (dpcontrol->vs == 3) *pa = 0;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int absimage(SIMAGE *sim, int mb)
{
    /*---------------------------------------*/
    /* add smoothed bands to image sequences */
    /*---------------------------------------*/
    int i,l,f,nx,ny,nb,nf,nl;
    uchar **im,*imf1,*imfn;
    CIMAGE cim0,cim1;
    cim1.nx = cim0.nx = nx = sim->nx;
    cim1.ny = cim0.ny = ny = sim->ny;
    cim1.nb = cim0.nb = nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    im = sim->im;
    for (f = nf; f <= nl; f++) {
        if (malloc1((void **) (&imfn),nb*mb*nx*ny,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"ERROR: absimage: fail to allocate buffer\n");
            return(FAILURE);
        }
        imf1 = im[f];
        for (i = 0; i < nx*ny*nb; i++) imfn[i] = imf1[i];
        cim1.im = imfn;
        for (l = 1; l < mb; l++) {
            cim0.im = cim1.im;
            cim1.im += nx*ny*nb;
            if (smocim(&cim0,&cim1) == FAILURE) {
                fprintf(stderr,"ERROR: absimage: fail to smooth image\n");
                return(FAILURE);
            }
        }
        if (free1((void **) (&imf1),nb*nx*ny,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"ERROR: absimage: fail to allocate buffer\n");
            return(FAILURE);
        }
        im[f] = imfn;
    }
    if (free1((void **) (&sim->val),nb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: absimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    if (malloc1((void **) (&sim->val),nb*mb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: absimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    sim->nb = nb*mb;
    return(SUCCESS);
}

int rbsimage(SIMAGE *sim, int mb)
{
    /*--------------------------------------------*/
    /* remove smoothed bands from image sequences */
    /*--------------------------------------------*/
    int i,f,nx,ny,nb,nf,nl;
    uchar **im,*imf1,*imfn;
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    im = sim->im;
    for (f = nf; f <= nl; f++) {
        if (malloc1((void **) (&imf1),nb*nx*ny/mb,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"ERROR: rbsimage: fail to allocate buffer\n");
            return(FAILURE);
        }
        imfn = im[f];
        for (i = 0; i < nx*ny*nb/mb; i++) imf1[i] = imfn[i];
        if (free1((void **) (&imfn),nb*nx*ny,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"ERROR: rbsimage: fail to allocate buffer\n");
            return(FAILURE);
        }
        im[f] = imf1;
    }
    if (free1((void **) (&sim->val),nb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: rbsimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    if (malloc1((void **) (&sim->val),nb/mb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: rbsimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    sim->nb = nb/mb;
    return(SUCCESS);
}

void mkvalsimage(SIMAGE *sim)
{
    int i,f,n,nx,ny,nb,nf,nl;
    uchar *imf,**im;
    fdata s0,s1,s2,imfi;
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    im = sim->im;
    for (n = 0; n < nb; n++) {
        s0 = s1 = s2 = 0.0;
        for (f = nf; f <= nl; f++) {
            imf = im[f]+n*nx*ny;
            for (i = 0; i < nx*ny; i++) {
                s0 += 1.0;
                s1 += imfi = imf[i];
                s2 += imfi*imfi;
            }
        }
        s1 /= s0;
        s2 /= s0;
        sim->val[n] = sqrt(s2-s1*s1);
    }
}

int opflow(SIMAGE *sim, VIMAGE *vim, DPCONTROL *dpcontrol)
{
    /*----------------------------------*/
    /* optical flow detection procedure */
    /*----------------------------------*/
    int i,j,k,n,nx,ny,nb,mb,nr,ns,nf,nl,eit,nf0,nl0,itc,a,wl,wr,wd,wu;
    fdata ***d,**dx,**dy,m,tau,spwfac;
    uchar **wim;
    eit = dpcontrol->eit;
    m = dpcontrol->m;
    tau = dpcontrol->tau;
    spwfac = dpcontrol->spwfac;
    if ((dpcontrol->hx) && (simxshrink(sim,&wim) == FAILURE)) {
        fprintf(stderr,"ERROR: opflow: fail to shrink images\n");
        return(FAILURE);
    }
    if (((mb = dpcontrol->mb) > 1) && (absimage(sim,mb) == FAILURE)) {
        fprintf(stderr,"ERROR: opflow: fail to add bands\n");
        return(FAILURE);
    }
    mkvalsimage(sim);
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    vim->nx = nx;
    vim->ny = ny;
    vim->nb = 2;
    if (vim->d == NULL) {
        if (alvimage(vim) == FAILURE) {
            fprintf(stderr,"ERROR: imalign: fail to allocate buffer for velocities\n");
            return(FAILURE);
        }
        d = vim->d;
        for (k = 0; k < 2; k++)
            for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) d[k][i][j] = 0.0;
    }
    dx = vim->dx;
    dy = vim->dy;
    if (getwwww(sim,vim,m,&wl,&wr,&wd,&wu) == FAILURE) {
        fprintf(stderr,"ERROR: xopflow: fail to get window parameters\n");
        return(FAILURE);
    }
    ns = 1+imax(wu-wd,wr-wl);
    nr = 1+imin(wu-wd,wr-wl);
    if (nr < 32) {
        fprintf(stderr,"ERROR: opflow: images too small\n");
        return(FAILURE);
    }
    if (ns > 4096) {
        fprintf(stderr,"ERROR: opflow: images too big\n");
        return(FAILURE);
    }
    n = log(nr*spwfac)/tau;
    if (dpcontrol->id) xyswap(vim);
    for (j = 0; j < dpcontrol->np; j++) {
        if (dpcontrol->vb) printf("Pass %d:\n",j+1);
        if (j != 0) autonorm(vim);
        itc = 0;
        if (dpcontrol->vs) {
            nf0 = (int)floor(m);
            nl0 = (int)ceil(m);
            if (nf0 < nf) nf0 = nf;
            if (nl0 > nl) nl0 = nl;
            if (nf0 == nl0) {
                if (nf0 > nf) nf0 -=1; else nl0 += 1;
            }
        } else {
            nf0 = nf;
            nl0 = nl;
        }
        for (i = n+1; i >= -eit; i--) {
            a = 1;
            if (oiter(sim,vim,dpcontrol,i,++itc,nf0,nl0,&a) == FAILURE) return(FAILURE);
            if (a && ((nf0 != nf) || (nl0 != nl))) {
                if ((nl-nl0) > (nf0-nf)) nl0 += 1; else nf0 -= 1;
                if (oiter(sim,vim,dpcontrol,i,++itc,nf0,nl0,&a) == FAILURE) return(FAILURE);
            }
        }
    }
    if (dpcontrol->id) xyswap(vim);
    if ((mb > 1) && (rbsimage(sim,mb) == FAILURE)) {
        fprintf(stderr,"ERROR: opflow: fail to remove bands\n");
        return(FAILURE);
    }
    if ((dpcontrol->hx) && (simxexpand(sim,wim) == FAILURE)) {
        fprintf(stderr,"ERROR: opflow: fail to expand images\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int opflow3(SIMAGE *sim, VIMAGE *vim, DPCONTROL *dpcontrol)
{
    uchar **im,*im0;
    int nf;
    im = sim->im;
    nf = sim->nf;
    dpcontrol->m = nf+0.5;
    im0      = im[nf];
    im[nf]   = im[nf+1];
    im[nf+1] = im[nf+2];
    im[nf+2] = im0;
    if (dpcontrol->vb) printf("Matching images %d and %d: ",nf+1,nf+2);
    if (opflow(sim,vim+0,dpcontrol) == FAILURE) return(FAILURE);
    im0      = im[nf];
    im[nf]   = im[nf+1];
    im[nf+1] = im[nf+2];
    im[nf+2] = im0;
    if (dpcontrol->vb) printf("Matching images %d and %d: ",nf+2,nf);
    if (opflow(sim,vim+1,dpcontrol) == FAILURE) return(FAILURE);
    im0      = im[nf];
    im[nf]   = im[nf+1];
    im[nf+1] = im[nf+2];
    im[nf+2] = im0;
    if (dpcontrol->vb) printf("Matching images %d and %d: ",nf,nf+1);
    if (opflow(sim,vim+2,dpcontrol) == FAILURE) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
