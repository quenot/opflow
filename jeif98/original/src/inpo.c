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

#include "glob.h"
#include "tens.h"
#include "imio.h"
#include "inpo.h"
#include "velo.h"

/*------------------------------------------------------------------------------------*/

int inter1(fdata rx, fdata ry, int nx, int ny, int nb, uchar *im, fdata *rim)
{
    int l,mx,my,nxy,lx0,lx1,ly0,ly1,p00,p01,p10,p11,ins=1;
    fdata px0,py0,px1,py1;
    mx = nx-1;
    my = ny-1;
    nxy = nx*ny;
    if (rx < 0) {
        ins = 0;
        lx0 = lx1 = 0;
        px0 = px1 = 0.5;
    } else if (rx < mx) {
        lx0 = (int)rx;
        lx1 = lx0+1;
        px1 = rx-(fdata)lx0;
        px0 = 1.0-px1;
    } else {
        if (rx > mx) ins = 0;
        lx0 = lx1 = mx;
        px0 = px1 = 0.5;
    }
    if (ry < 0) {
        ins = 0;
        ly0 = ly1 = 0;
        py0 = py1 = 0.5;
    } else if (ry < my) {
        ly0 = (int)ry;
        ly1 = ly0+1;
        py1 = ry-(fdata)ly0;
        py0 = 1.0-py1;
    } else {
        if (ry > my) ins = 0;
        ly0 = ly1 = my;
        py0 = py1 = 0.5;
    }
    ly0 *= nx;
    ly1 *= nx;
    p00 = lx0+ly0;
    p10 = lx1+ly0;
    p01 = lx0+ly1;
    p11 = lx1+ly1;
    for (l = 0; l < nb; l++) {
        rim[l] = py0*(px0*((fdata)(im[p00]))+px1*((fdata)(im[p10])))
               + py1*(px0*((fdata)(im[p01]))+px1*((fdata)(im[p11])));
        p00 += nxy;
        p10 += nxy;
        p01 += nxy;
        p11 += nxy;
    }
    return(ins);
}

int inter2(fdata rx, fdata ry, int nx, int ny, int nb, uchar *im, fdata *rim)
{
    int l,mx,my,nxy,lx0,lx1,lx2,lx3,ly0,ly1,ly2,ly3,ins=1;
    int p00,p01,p02,p03,p10,p11,p12,p13,p20,p21,p22,p23,p30,p31,p32,p33;
    fdata px0,px1,px2,px3,py0,py1,py2,py3,c6=1.0/6.0;
    mx = nx-1;
    my = ny-1;
    nxy = nx*ny;
    if (rx < 0) {
        ins = 0;
        lx1 = lx2 = 0;
        px1 = px2 = 0.5;
    } else if (rx < mx) {
        lx1 = (int)rx;
        lx2 = lx1+1;
        px2 = rx-(fdata)lx1;
        px1 = 1.0-px2;
    } else {
        if (rx > mx) ins = 0;
        lx1 = lx2 = mx;
        px1 = px2 = 0.5;
    }
    if (ry < 0) {
        ins = 0;
        ly1 = ly2 = 0;
        py1 = py2 = 0.5;
    } else if (ry < my) {
        ly1 = (int)ry;
        ly2 = ly1+1;
        py2 = ry-(fdata)ly1;
        py1 = 1.0-py2;
    } else {
        if (ry > my) ins = 0;
        ly1 = ly2 = my;
        py1 = py2 = 0.5;
    }
    lx0 = (lx1 == 0) ? 0 : lx1-1;
    lx3 = (lx2 == mx) ? mx : lx2+1;
    ly0 = (ly1 == 0) ? 0 : ly1-1;
    ly3 = (ly2 == my) ? my : ly2+1;
    px0 = c6*px1*(px1*px1-1.0);
    px3 = c6*px2*(px2*px2-1.0);
    px1 += px3-2*px0;
    px2 += px0-2*px3;
    py0 = c6*py1*(py1*py1-1.0);
    py3 = c6*py2*(py2*py2-1.0);
    py1 += py3-2*py0;
    py2 += py0-2*py3;
    ly0 *= nx;
    ly1 *= nx;
    ly2 *= nx;
    ly3 *= nx;
    p00 = lx0+ly0;
    p10 = lx1+ly0;
    p20 = lx2+ly0;
    p30 = lx3+ly0;
    p01 = lx0+ly1;
    p11 = lx1+ly1;
    p21 = lx2+ly1;
    p31 = lx3+ly1;
    p02 = lx0+ly2;
    p12 = lx1+ly2;
    p22 = lx2+ly2;
    p32 = lx3+ly2;
    p03 = lx0+ly3;
    p13 = lx1+ly3;
    p23 = lx2+ly3;
    p33 = lx3+ly3;
    for (l = 0; l < nb; l++) {
        rim[l] = py0*(px0*((fdata)(im[p00]))+px1*((fdata)(im[p10]))+
                      px2*((fdata)(im[p20]))+px3*((fdata)(im[p30])))
               + py1*(px0*((fdata)(im[p01]))+px1*((fdata)(im[p11]))+
                      px2*((fdata)(im[p21]))+px3*((fdata)(im[p31])))
               + py2*(px0*((fdata)(im[p02]))+px1*((fdata)(im[p12]))+
                      px2*((fdata)(im[p22]))+px3*((fdata)(im[p32])))
               + py3*(px0*((fdata)(im[p03]))+px1*((fdata)(im[p13]))+
                      px2*((fdata)(im[p23]))+px3*((fdata)(im[p33])));
        p00 += nxy;
        p10 += nxy;
        p20 += nxy;
        p30 += nxy;
        p01 += nxy;
        p11 += nxy;
        p21 += nxy;
        p31 += nxy;
        p02 += nxy;
        p12 += nxy;
        p22 += nxy;
        p32 += nxy;
        p03 += nxy;
        p13 += nxy;
        p23 += nxy;
        p33 += nxy;
    }
    return(ins);
}

int apply(VIMAGE *vim, SIMAGE *sim, int f, DIMAGE *dim, fdata m, int dir, int dni,
int bsi)
{
    /*---------------------------------------------------------------*/
    /* Builds the transform of an image using a given velocity field */
    /* Builds also a boolean image defining where the velocity field */
    /* points inside the image                                       */
    /*---------------------------------------------------------------*/
    int i,j,l,nx,ny,nb,wx,wy,im,ip;
    fdata **dx,**dy,*dxi,*dyi,*dxm,*dxp,*dym,*dyp,***xm,**xmi;
    fdata rx,ry,*rim,s;
    uchar *imf,**in,*ini;
    s = ((fdata) f)-m;
    dx = vim->dx;
    dy = vim->dy;
    nx = vim->nx;
    ny = vim->ny;
    if (dni) nx = 2*nx-1;
    dim->nx = nx;
    dim->ny = ny;
    wx = sim->nx;
    wy = sim->ny;
    nb = sim->nb;
    dim->nb = nb;
    imf = sim->im[f];
    if ((malloc1((void **) &rim,nb,sizeof(fdata)) == FAILURE)
     || (malloc1((void **) &xmi,nb,sizeof(fdata *)) == FAILURE)
     || (malloc2((void ***) &in,nx,ny,sizeof(uchar)) == FAILURE)
     || (malloc3((void ****) &xm,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: apply: fail to allocate buffers\n");
        return(FAILURE);
    }
    for (i = 0; i < nx; i++) {
        ini = in[i];
        im = i >> 1;
        ip = (i+1) >> 1;
        dxm = dx[im];
        dym = dy[im];
        dxp = dx[ip];
        dyp = dy[ip];
        dxi = dx[i];
        dyi = dy[i];
        for (l = 0; l < nb; l++) xmi[l] = xm[l][i];
        for (j = 0; j < ny; j++) {
            if (dni) {
                rx = 0.5*(i+s*(dxm[j]+dxp[j]));
                ry = j+0.5*s*(dym[j]+dyp[j]);
            } else {
                rx = i+s*dxi[j];
                ry = j+s*dyi[j];
            }
            ini[j] = dir ? bsi ? inter2(ry,rx,wx,wy,nb,imf,rim)
                               : inter1(ry,rx,wx,wy,nb,imf,rim)
                         : bsi ? inter2(rx,ry,wx,wy,nb,imf,rim)
                               : inter1(rx,ry,wx,wy,nb,imf,rim);
            for (l = 0; l < nb; l++) xmi[l][j] = rim[l];
        }
    }
    if ((free1((void **) &rim,nb,sizeof(fdata)) == FAILURE)
     || (free1((void **) &xmi,nb,sizeof(fdata *)) == FAILURE)) {
        fprintf(stderr,"ERROR: apply: fail to free buffer\n");
        return(FAILURE);
    }
    dim->xm = xm;
    dim->in = in;
    return(SUCCESS);
}

void goadj(DIMAGE *dim0, DIMAGE *dim1, fdata lambda)
{
    /*----------------------------------------------------------------*/
    /* Apply a gain/offset correction between two images so that they */
    /* have the same average and dispersion for each of their bands   */
    /*----------------------------------------------------------------*/
    int i,j,l,nx,ny,nb;
    fdata ***xm0,**xm0l,*xm0li,***xm1,**xm1l,*xm1li;
    fdata aa,m0,m1,s0,s1,x0,x1,mm,ss,l0,l1;
    uchar **in0,*in0i,**in1,*in1i;
    l1 = lambda;
    if (l1<0.0) l1 = 0.0;
    if (l1>1.0) l1 = 1.0;
    l0 = 1-l1;
    nx = dim0->nx;
    ny = dim0->ny;
    nb = dim0->nb;
    xm0 = dim0->xm;
    in0 = dim0->in;
    xm1 = dim1->xm;
    in1 = dim1->in;
    for (l = 0; l < nb; l++) {
        aa = m0 = m1 = s0 = s1 = 0;
        xm0l = xm0[l];
        xm1l = xm1[l];
        for (i = 0; i < nx; i++) {
            in0i = in0[i];
            in1i = in1[i];
            xm0li = xm0l[i];
            xm1li = xm1l[i];
            for (j = 0; j < ny; j++) {
                if ((in0i[j]) && (in1i[j])) {
                    aa += 1;
                    m0 += x0 = xm0li[j];
                    m1 += x1 = xm1li[j];
                    s0 += x0*x0;
                    s1 += x1*x1;
                }
            }
        }
        if (aa != 0) {        
            m0 /= aa;
            m1 /= aa;
            s0 /= aa;
            s1 /= aa;
            s0 = sqrt(s0-m0*m0);
            s1 = sqrt(s1-m1*m1);
            ss = l0*s0+l1*s1;
            mm = l0*m0+l1*m1;
            if (ss != 0) {        
                s0 /= ss;
                s1 /= ss;
            } else {
                s0 = s1 = 1;
            }
            m0 -= mm*s0;
            m1 -= mm*s1;
        } else {
            s0 = s1 = 1;
            m0 = m1 = 0;
        }
        for (i = 0; i < nx; i++) {
            xm0li = xm0l[i];
            xm1li = xm1l[i];
            for (j = 0; j < ny; j++) {
                xm0li[j] = (xm0li[j]-m0)/s0;
                xm1li[j] = (xm1li[j]-m1)/s1;
            }
        }
    }
}

int getimage(int dir, SIMAGE *sim, VIMAGE *vim, fdata m, DIMAGE *dim0, DIMAGE *dim1,
int goa, int dni, int bsi)
{
    /*---------------------------------------------------------------------*/
    /* Transforms two images im(nf) and im(nl) defined in planes nl and nf */
    /* using a velocity field into an intermediate plane defines by a      */
    /* parameter m (where matching will occur).                            */
    /* A gain/offset correction is performed optionally.                   */
    /*---------------------------------------------------------------------*/
    if (apply(vim,sim,sim->nf,dim0,m,dir,dni,bsi) == FAILURE) return(FAILURE);
    if (apply(vim,sim,sim->nl,dim1,m,dir,dni,bsi) == FAILURE) return(FAILURE);
    if (goa) goadj(dim0,dim1,m-((fdata)sim->nf)/((fdata)(sim->nl-sim->nf)));
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int rerrorf(SIMAGE *sim, VIMAGE *vim, fdata m, fdata *prerrorf, fdata *prerrorp)
{
    int i,j,l,nx,ny,nb,nf,nl,nf0,nf1;
    fdata ***xm0,***xm1,rerrorf,cerrorf,rerrorp,cerrorp,errornb;
    uchar **in0,**in1;
    DIMAGE dim0,dim1;
    VIMAGE wim[1];
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    nf0 = (int)floor(m);
    nf1 = (int)ceil(m);
    if (nf0 < nf) nf0 = nf;
    if (nf1 > nl) nf1 = nl;
    if (nf0 == nf1) {
            if (nf0 > nf) nf0 -=1; else nf1 += 1;
    }
    if (nx == (2*vim->nx)) {
        if (vimxexpand(vim,0,wim,1) == FAILURE) {
            fprintf(stderr,"ERROR: opflow: fail to expand images\n");
            return(FAILURE);
        }
    } else {
        wim->dx  = vim->dx;
        wim->dy  = vim->dy;
        wim->nx  = vim->nx;
        wim->ny  = vim->ny;
        wim->nx0 = vim->nx0;
        wim->ny0 = vim->ny0;
        wim->nx1 = vim->nx1;
        wim->ny1 = vim->ny1;
    }
    if ((apply(wim,sim,nf0,&dim0,m,0,0,1) == FAILURE)
     || (apply(wim,sim,nf1,&dim1,m,0,0,1) == FAILURE)) {
        fprintf(stderr,"ERROR: rerrorf: fail to transform images\n");
        return(FAILURE);
    }
    xm0 = dim0.xm;
    xm1 = dim1.xm;
    in0 = dim0.in;
    in1 = dim1.in;
    rerrorf = cerrorf = rerrorp = cerrorp = 0.0;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            errornb = 0.0;
            for (l = 0; l < nb; l++) {
                errornb += fabs(xm0[l][i][j]-xm1[l][i][j]);
            }
            cerrorf += 255.0;
            rerrorf += errornb;
            if ((in0[i][j]) && (in1[i][j])) {
                cerrorp += 255.0;
                rerrorp += errornb;
            }
        }
    }
    *prerrorf = rerrorf /= cerrorf*nb;
    *prerrorp = rerrorp /= cerrorp*nb;
    if ((free2((void ***) &dim0.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim0.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim1.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim1.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: rerrorf: fail to free buffers\n");
        return(FAILURE);
    }
    if ((nx == (2*vim->nx)) && ((frvimage(wim)) == FAILURE)) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int iminterp(SIMAGE *sim, CIMAGE *com, VIMAGE *vim, int fr, fdata m, int ex)
{
    int i,j,l,nx,ny,nb,nf,nl,nf0,nf1,omlij;
    fdata l0,l1,***xm0,***xm1,lambda;
    uchar *omd,**in0,**in1;
    DIMAGE dim0,dim1;
    VIMAGE wim[1];
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    if ((m > nf) && (m < nl)) {
        if (m == ((int)m)) {
            nf0 = m-1;
            nf1 = m+1;
        } else {
            nf0 = (int)floor(m);
            nf1 = (int)ceil(m);
        } 
    } else if (m <= nf) {
        nf0 = nf;
        nf1 = nf+1;
    } else {
        nf0 = nl;
        nf1 = nl-1;
    }
    if (nx == (2*vim->nx)) {
        if (vimxexpand(vim,0,wim,1) == FAILURE) {
            fprintf(stderr,"ERROR: opflow: fail to expand images\n");
            return(FAILURE);
        }
    } else {
        wim->d  = vim->d;
        wim->dx  = vim->dx;
        wim->dy  = vim->dy;
        wim->nx  = vim->nx;
        wim->ny  = vim->ny;
        wim->nx0 = vim->nx0;
        wim->ny0 = vim->ny0;
        wim->nx1 = vim->nx1;
        wim->ny1 = vim->ny1;
    }
    if ((apply(wim,sim,nf0,&dim0,m,0,0,1) == FAILURE)
     || (apply(wim,sim,nf1,&dim1,m,0,0,1) == FAILURE)) {
        fprintf(stderr,"ERROR: iminterp: fail to transform images\n");
        return(FAILURE);
    }
    xm0 = dim0.xm;
    xm1 = dim1.xm;
    in0 = dim0.in;
    in1 = dim1.in;
    if (malloc1((void **) &omd,nb*nx*ny,sizeof(uchar)) == FAILURE) {
        fprintf(stderr,"ERROR: iminterp: fail to allocate output buffer\n");
        return(FAILURE);
    }
    com->nb = nb;
    com->ny = ny;
    com->nx = nx;
    com->im = omd; 
    l1 = lambda = (m-(fdata)nf0)/((fdata)(nf1-nf0));
    if (l1 < 0.0) l1 = 0.0;
    if (l1 > 1.0) l1 = 1.0;
    l0 = 1-l1;
    if ((lambda > 0) && (lambda < 1)) {
        for (l = 0; l < nb; l++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    if (in0[i][j] && !in1[i][j] && !ex) {
                        omlij = (int) rint(xm0[l][i][j]);
                    } else if (!in0[i][j] && in1[i][j] && !ex) {
                        omlij = (int) rint(xm1[l][i][j]);
                    } else {
                        omlij = (int) rint(l0*xm0[l][i][j]+l1*xm1[l][i][j]);
                    }
                    if (omlij < 0)   omlij = 0;
                    if (omlij > 255) omlij = 255;
                    *omd++ = (uchar) omlij;
                }
            }
        }
    }
    else if (lambda <= 0) {
        for (l = 0; l < nb; l++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    omlij = (int) rint(xm0[l][i][j]);
                    if (omlij < 0)   omlij = 0;
                    if (omlij > 255) omlij = 255;
                    *omd++ = (uchar) omlij;
                }
            }
        }
    }
    else if (lambda >= 1) {
        for (l = 0; l < nb; l++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    omlij = (int) rint(xm1[l][i][j]);
                    if (omlij < 0)   omlij = 0;
                    if (omlij > 255) omlij = 255;
                    *omd++ = (uchar) omlij;
                }
            }
        }
    }
    if ((free2((void ***) &dim0.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim0.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim1.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim1.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: iminterp: fail to free buffers\n");
        return(FAILURE);
    }
    if ((nx == (2*vim->nx)) && (frvimage(wim) == FAILURE)) return(FAILURE);
    if (fr && (frvimage(vim) == FAILURE)) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int triinterp(SIMAGE *sim, CIMAGE *com, VIMAGE *vim, fdata *l, int ex)
{
    int i,j,k,f,nx,ny,nb,nf,omkij;
    fdata ***xm0,***xm1,***xm2,l0,l1,l2;
    uchar *omd,**in0,**in1,**in2;
    DIMAGE dim0,dim1,dim2;
    VIMAGE wim[3];
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    printf("lambda%d = %2.3f  ",nf+0,(float) l[0]);
    printf("lambda%d = %2.3f  ",nf+1,(float) l[1]);
    printf("lambda%d = %2.3f\n",nf+2,(float) l[2]);
    vcolin(vim+5,vim+7,vim+0,l[1],l[2]);
    vcolin(vim+3,vim+8,vim+1,l[2],l[0]);
    vcolin(vim+4,vim+6,vim+2,l[0],l[1]);
    if (invdep(vim+0,1,vim+0,1) == FAILURE) return(FAILURE);
    if (invdep(vim+1,1,vim+1,1) == FAILURE) return(FAILURE);
    if (invdep(vim+2,1,vim+2,1) == FAILURE) return(FAILURE);
    if (nx == (2*vim[0].nx)) {
        if ((vimxexpand(vim+0,0,wim+0,1) == FAILURE)
         || (vimxexpand(vim+1,0,wim+1,1) == FAILURE)
         || (vimxexpand(vim+2,0,wim+2,1) == FAILURE)) {
            fprintf(stderr,"ERROR: triinterp: fail to expand images\n");
            return(FAILURE);
        }
    } else {
        for (f = 0; f < 3; f++) {
            wim[f].dx  = vim[f].dx;
            wim[f].dy  = vim[f].dy;
            wim[f].nx  = vim[f].nx;
            wim[f].ny  = vim[f].ny;
            wim[f].nx0 = vim[f].nx0;
            wim[f].ny0 = vim[f].ny0;
            wim[f].nx1 = vim[f].nx1;
            wim[f].ny1 = vim[f].ny1;
        }
    }
    if ((apply(wim+0,sim,nf+0,&dim0,nf-1.0,0,0,1) == FAILURE)
     || (apply(wim+1,sim,nf+1,&dim1,nf+0.0,0,0,1) == FAILURE)
     || (apply(wim+2,sim,nf+2,&dim2,nf+1.0,0,0,1) == FAILURE)) {
        fprintf(stderr,"ERROR: triinterp: fail to transform images\n");
        return(FAILURE);
    }
    xm0 = dim0.xm;
    xm1 = dim1.xm;
    xm2 = dim2.xm;
    in0 = dim0.in;
    in1 = dim1.in;
    in2 = dim2.in;
    if (malloc1((void **) &omd,nb*nx*ny,sizeof(uchar)) == FAILURE) {
        fprintf(stderr,"ERROR: triinterp: fail to allocate output buffer\n");
        return(FAILURE);
    }
    com->nb = nb;
    com->ny = ny;
    com->nx = nx;
    com->im = omd;
    if ((l[1] < 0.0) && (l[2] < 0.0)) {
        l0 = 1.0;
        l1 = 0.0;
        l2 = 0.0;
    } else if ((l[2] < 0.0) && (l[0] < 0.0)) {
        l0 = 0.0;
        l1 = 1.0;
        l2 = 0.0;
    } else if ((l[0] < 0.0) && (l[1] < 0.0)) {
        l0 = 0.0;
        l1 = 0.0;
        l2 = 1.0;
    } else if (l[0] < 0.0) {
        l0 = 0.0;
        l1 = l[1]/(l[1]+l[2]);
        l2 = l[2]/(l[1]+l[2]);
    } else if (l[1] < 0.0) {
        l0 = l[0]/(l[2]+l[0]);
        l1 = 0.0;
        l2 = l[2]/(l[2]+l[0]);
    } else if (l[2] < 0.0) {
        l0 = l[0]/(l[0]+l[1]);
        l1 = l[1]/(l[0]+l[1]);
        l2 = 0.0;
    } else {
        l0 = l[0];
        l1 = l[1];
        l2 = l[2];
    }
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                omkij = (int) rint(l0*xm0[k][i][j]
                                  +l1*xm1[k][i][j]
                                  +l2*xm2[k][i][j]);
                if (omkij < 0)   omkij = 0;
                if (omkij > 255) omkij = 255;
                *omd++ = (uchar) omkij;
            }
        }
    }
    if ((free2((void ***) &dim0.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim0.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim1.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim1.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim2.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim2.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: triinterp: fail to free buffers\n");
        return(FAILURE);
    }
    if ((nx == (2*vim[0].nx))
      && ((frvimage(wim+0) == FAILURE)
       || (frvimage(wim+1) == FAILURE)
       || (frvimage(wim+2) == FAILURE))) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int steinterp(SIMAGE *sim, CIMAGE *com, VIMAGE *vim, fdata x, fdata y, fdata z, fdata r,
fdata dxm, fdata dxs, int ex)
{
    int i,j,k,f,nvx,nvy,nx,ny,nb,nf,omkij;
    fdata ***xm0,***xm1,l0,l1,x0,x1,di0,dj0,dk0,di1,dj1,dk1,ic,jc,s;
    fdata **dx0,**dy0,**dx1,**dy1,**dx00,**dy00,**dx11,**dy11;
    uchar *omd,**in0,**in1;
    DIMAGE dim0,dim1;
    VIMAGE wim[2];
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    printf("x = %2.3f  ",(float) x);
    printf("y = %2.3f  ",(float) y);
    printf("z = %2.3f\n",(float) z);
    nvx = vim[0].nx;
    nvy = vim[0].ny;
    dx00 = vim[0].dx;
    dy00 = vim[0].dy;
    dx11 = vim[1].dx;
    dy11 = vim[1].dy;
    dx0 = vim[2].dx;
    dy0 = vim[2].dy;
    dx1 = vim[3].dx;
    dy1 = vim[3].dy;
    x0 = x+0.5;
    x1 = x-0.5;
    ic = 0.5*(nvx-1.0);
    jc = 0.5*(nvy-1.0);
    s = sqrt(3.0/(ic*ic+jc*jc));
    for (i = 0; i < nvx; i++) {
        for (j = 0; j < nvy; j++) {
            di0 = x0*dx0[i][j];
            dj0 = x0*dy0[i][j]-y*dx0[i][j];
            di1 = x1*dx1[i][j];
            dj1 = x1*dy1[i][j]-y*dx1[i][j];
            dk0 = z*(r*(dxm-dx0[i][j])+dxs);
            dk1 = z*(r*(dxm-dx1[i][j])+dxs);
            dx00[i][j] = di0+s*(i+di0-ic)*dk0;
            dy00[i][j] = dj0+s*(j+dj0-jc)*dk0;
            dx11[i][j] = di1+s*(i+di1-ic)*dk1;
            dy11[i][j] = dj1+s*(j+dj1-jc)*dk1;
        }
    }
    if (invdep(vim+0,1,vim+0,1) == FAILURE) return(FAILURE);
    if (invdep(vim+1,1,vim+1,1) == FAILURE) return(FAILURE);
    if (nx == 2*nvx) {
        if ((vimxexpand(vim+0,0,wim+0,1) == FAILURE)
         || (vimxexpand(vim+1,0,wim+1,1) == FAILURE)) {
            fprintf(stderr,"ERROR: steinterp: fail to expand images\n");
            return(FAILURE);
        }
    } else {
        for (f = 0; f < 2; f++) {
            wim[f].dx  = vim[f].dx;
            wim[f].dy  = vim[f].dy;
            wim[f].nx  = vim[f].nx;
            wim[f].ny  = vim[f].ny;
            wim[f].nx0 = vim[f].nx0;
            wim[f].ny0 = vim[f].ny0;
            wim[f].nx1 = vim[f].nx1;
            wim[f].ny1 = vim[f].ny1;
        }
    }
    if ((apply(wim+0,sim,nf+0,&dim0,nf-1.0,0,0,1) == FAILURE)
     || (apply(wim+1,sim,nf+1,&dim1,nf+0.0,0,0,1) == FAILURE)) {
        fprintf(stderr,"ERROR: steinterp: fail to transform images\n");
        return(FAILURE);
    }
    xm0 = dim0.xm;
    xm1 = dim1.xm;
    in0 = dim0.in;
    in1 = dim1.in;
    if (malloc1((void **) &omd,nb*nx*ny,sizeof(uchar)) == FAILURE) {
        fprintf(stderr,"ERROR: steinterp: fail to allocate output buffer\n");
        return(FAILURE);
    }
    com->nb = nb;
    com->ny = ny;
    com->nx = nx;
    com->im = omd;
    l1 = x+0.5;
    if (l1 < 0.0) l1 = 0.0;
    if (l1 > 1.0) l1 = 1.0;
    l0 = 1-l1;
    if ((x > -0.5) && (x < 0.5)) {
        for (k = 0; k < nb; k++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    if (in0[i][j] && !in1[i][j] && !ex) {
                        omkij = (int) rint(xm0[k][i][j]);
                    } else if (!in0[i][j] && in1[i][j] && !ex) {
                        omkij = (int) rint(xm1[k][i][j]);
                    } else {
                        omkij = (int) rint(l0*xm0[k][i][j]+l1*xm1[k][i][j]);
                    }
                    if (omkij < 0)   omkij = 0;
                    if (omkij > 255) omkij = 255;
                    *omd++ = (uchar) omkij;
                }
            }
        }
    }
    else if (x <= -0.5) {
        for (k = 0; k < nb; k++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    omkij = (int) rint(xm0[k][i][j]);
                    if (omkij < 0)   omkij = 0;
                    if (omkij > 255) omkij = 255;
                    *omd++ = (uchar) omkij;
                }
            }
        }
    }
    else if (x >= 0.5) {
        for (k = 0; k < nb; k++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    omkij = (int) rint(xm1[k][i][j]);
                    if (omkij < 0)   omkij = 0;
                    if (omkij > 255) omkij = 255;
                    *omd++ = (uchar) omkij;
                }
            }
        }
    }
    if ((free2((void ***) &dim0.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim0.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim1.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim1.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: steinterp: fail to free buffers\n");
        return(FAILURE);
    }
    if ((nx == 2*nvx)
      && ((frvimage(wim+0) == FAILURE)
       || (frvimage(wim+1) == FAILURE))) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
