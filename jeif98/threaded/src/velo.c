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

/*------------------------------------------------------------------------------------*/
/*  intersect routine: compute (*xi,*yi) intersection of the lines:                   */
/*  (1) passing by (x1,y1) and (x2,y2) and                                            */
/*  (2) passing by (x,y) and directed by (vx,vy).                                     */
/*------------------------------------------------------------------------------------*/

void intersect(double x, double y, double x1, double y1, double x2, double y2,
double vx, double vy, double *xi, double *yi)
{
    double xm1,ym1,x12,y12,lambda,r;
    xm1 = x-x1;
    ym1 = y-y1;
    x12 = x2-x1;
    y12 = y2-y1;
    if (fabs(r = vx*y12-vy*x12) < 1.0e-30) r = 1.0e+30;
    lambda = (xm1*y12-ym1*x12)/r;
    *xi = x-lambda*vx;
    *yi = y-lambda*vy;
}

/*------------------------------------------------------------------------------------*/
/*  get_bar function: computes barycentric coordinates (*r00,*r01,*r10,*r11) of a     */
/*  (x,y) point in a ((x00,y00),(x01,y01),(x10,y10),(x11,y11)) quadrilateral and      */
/*  and returns 1 if the point is inside with an epsr accuracy (0 else).              */
/*------------------------------------------------------------------------------------*/

int get_bar(double x, double y, double x00, double y00, double x01, double y01,
double x10, double y10, double x11, double y11, double *r00, double *r01,
double *r10, double *r11, double epsr)
{
    double xx0,yx0,xx1,yx1,x0y,y0y,x1y,y1y;
    double xvx0,yvx0,xvx1,yvx1,xv0y,yv0y,xv1y,yv1y;
    double xv00,yv00,xv11,yv11,rd00m,rd11m,rd00p,rd11p;
    double xvxx,yvxx,xvyy,yvyy,cxx,cx0,cyy,c0y,rx0,rx1,r0y,r1y;
    double ss00,ss01,ss10,ss11;
    int inside,sec00,sec01,sec10,sec11,singular;
    xvx0 = x10-x00; /* vx0 = M10 - M00 */
    yvx0 = y10-y00;
    xvx1 = x11-x01; /* vx1 = M11 - M01 */
    yvx1 = y11-y01;
    xv0y = x01-x00; /* v0y = M01 - M00 */
    yv0y = y01-y00;
    xv1y = x11-x10; /* v1y = M11 - M10 */
    yv1y = y11-y10;
    xv00 = x-x00;   /* v00 = M - M00 */
    yv00 = y-y00;
    xv11 = x-x11;   /* v11 = M - M11 */
    yv11 = y-y11;
    rd00m = yv0y*xv00-xv0y*yv00; /* (v00 ^ +v0y) */
    rd00p = xvx0*yv00-yvx0*xv00; /* (+vx0 ^ v00) */
    rd11m = xv1y*yv11-yv1y*xv11; /* (v11 ^ -v1y) */
    rd11p = yvx1*xv11-xvx1*yv11; /* (-vx1 ^ v11) */
    ss00 = xvx0*yv0y-yvx0*xv0y;  /* ss00 = (vx0 ^ v0y) */
    ss01 = xvx1*yv0y-yvx1*xv0y;  /* ss01 = (vx1 ^ v0y) */
    ss10 = xvx0*yv1y-yvx0*xv1y;  /* ss10 = (vx0 ^ v1y) */
    ss11 = xvx1*yv1y-yvx1*xv1y;  /* ss11 = (vx1 ^ v1y) */
    sec00 = ((ss00>0)&&(epsr*ss00+rd00m>0)&&(epsr*ss00+rd00p>0));
    sec01 = ((ss01>0)&&(epsr*ss01+rd00m>0)&&(epsr*ss01+rd11p>0));
    sec10 = ((ss10>0)&&(epsr*ss10+rd11m>0)&&(epsr*ss10+rd00p>0));
    sec11 = ((ss11>0)&&(epsr*ss11+rd11m>0)&&(epsr*ss11+rd11p>0));
    if ((ss00 > 0) && (ss01 > 0) && (ss10 > 0) && (ss11 > 0)) {
        inside = (sec00 && sec11);
        singular = 0;
    } else if ((ss00 > 0) && (ss01 > 0) && (ss10 > 0) && (ss11 <= 0)) {
        inside = (sec00 && (sec10 || sec01));
        singular = 1;
    } else if ((ss00 > 0) && (ss01 > 0) && (ss10 <= 0) && (ss11 > 0)) {
        inside = (sec01 && (sec00 || sec11));
        singular = 1;
    } else if ((ss00 > 0) && (ss01 <= 0) && (ss10 > 0) && (ss11 > 0)) {
        inside = (sec10 && (sec10 || sec01));
        singular = 1;
    } else if ((ss00 <= 0) && (ss01 > 0) && (ss10 > 0) && (ss11 > 0)) {
        inside = (sec11 && (sec00 || sec11));
        singular = 1;
    } else if ((ss00 > 0) && (ss01 > 0) && (ss11 <= 0) && (ss10 <= 0)) {
        inside = (sec00 && sec01);
        singular = 2;
    } else if ((ss00 > 0) && (ss01 <= 0) && (ss11 <= 0) && (ss10 > 0)) {
        inside = (sec00 && sec10);
        singular = 2;
    } else if ((ss00 <= 0) && (ss01 <= 0) && (ss11 > 0) && (ss10 > 0)) {
        inside = (sec11 && sec10);
        singular = 2;
    } else if ((ss00 <= 0) && (ss01 > 0) && (ss11 > 0) && (ss10 <= 0)) {
        inside = (sec11 && sec01);
        singular = 2;
    } else if ((ss00 <= 0) && (ss01 <= 0) && (ss10 <= 0) && (ss11 <= 0)) {
        inside = 0;
        singular = 4;
    } else {
        inside = 0;
        singular = 3;
    }
    if (inside) {
        cxx = xvx0*yvx1-yvx0*xvx1; /* (vx0 ^ vx1) */
        cx0 = xv0y*yvx1-yv0y*xvx1; /* (v0y ^ vx1) */
        xvyy = cxx*xv00-cx0*xvx0;  /* direction of M1y - M0y */
        yvyy = cxx*yv00-cx0*yvx0;
        cyy = xv0y*yv1y-yv0y*xv1y; /* (v0y ^ v1y) */
        c0y = xvx0*yv1y-yvx0*xv1y; /* (vx0 ^ v1y) */
        xvxx = cyy*xv00-c0y*xv0y;  /* direction of Mx1 - Mx0 */
        yvxx = cyy*yv00-c0y*yv0y;
        intersect(x,y,x00,y00,x10,y10,xvxx,yvxx,&xx0,&yx0);
        intersect(x,y,x01,y01,x11,y11,xvxx,yvxx,&xx1,&yx1);
        intersect(x,y,x00,y00,x01,y01,xvyy,yvyy,&x0y,&y0y);
        intersect(x,y,x10,y10,x11,y11,xvyy,yvyy,&x1y,&y1y);
        xvxx = xx1-xx0;
        yvxx = yx1-yx0;
        xvyy = x1y-x0y;
        yvyy = y1y-y0y;
        rx0 = (xvxx*(x-xx0)+yvxx*(y-yx0))/(xvxx*xvxx+yvxx*yvxx);
        if ((rx0 < (0.0-2*epsr)) || (rx0 > (1.0+2*epsr))) inside = 0;
        rx1 = 1-rx0;
        r0y = (xvyy*(x-x0y)+yvyy*(y-y0y))/(xvyy*xvyy+yvyy*yvyy);
        if ((r0y < (0.0-2*epsr)) || (r0y > (1.0+2*epsr))) inside = 0;
        r1y = 1-r0y;
        *r00 = rx1*r1y;
        *r01 = rx0*r1y;
        *r10 = rx1*r0y;
        *r11 = rx0*r0y;
    }
    return(inside);
}

/*------------------------------------------------------------------------------------*/
/*  (rdk[][],rdl[][]) defines the displacement for the points with real               */
/*  coordinates (rk[][],rl[][]).						      */
/*------------------------------------------------------------------------------------*/
/*  fillquad computes the interpolation (rdi[][],rdj[][]) for all points with         */
/*  integer coordinates inside the ((rk[i][j],rl[i][j]),(rk[i][j+1],rl[i][j+1]),      */
/*  (rk[i+1][j],rl[i+1][j]),(rk[i+1][j+1],rl[i+1][j+1])) quadrilateral.		      */
/*------------------------------------------------------------------------------------*/
/*  A boolean matrix donex is used to remember if a point has already been computed   */
/*  and save execution time							      */
/*------------------------------------------------------------------------------------*/

void fillquad(int i0, int j0, fdata **rk, fdata **rl, fdata **rdk, fdata **rdl,
fdata **rdi, fdata **rdj, int nx, int ny, int **donex, fdata epsr, int ep)
{
    double x,y,x00,y00,x01,y01,x10,y10,x11,y11,r00,r01,r10,r11,xmin,xmax,ymin,ymax;
    fdata *rdi0,*rdj0,*rk0,*rk1,*rl0,*rl1,*rdk0,*rdk1,*rdl0,*rdl1;
    int *donex0,k,l,mx,my,kmin,kmax,lmin,lmax,first,i1,j1;
    mx = nx-1;
    my = ny-1;
    i1 = i0+1;
    j1 = j0+1;
    rk0 = rk[i0];
    xmax = xmin = x00 = rk0[j0];
    if ((x01 = rk0[j1]) < xmin) xmin = x01;
    if (x01 > xmax) xmax = x01;
    rk1 = rk[i1];
    if ((x10 = rk1[j0]) < xmin) xmin = x10;
    if (x10 > xmax) xmax = x10;
    if ((x11 = rk1[j1]) < xmin) xmin = x11;
    if (x11 > xmax) xmax = x11;
    kmin = ceil(xmin-epsr);
    if (ep) kmin -= epsr/EPS_BAR;
    if (kmin < 0)  kmin = 0;
    kmax = floor(xmax+epsr);
    if (ep) kmax += epsr/EPS_BAR;
    if (kmax > mx) kmax = mx;
    rl0 = rl[i0];
    ymax = ymin = y00 = rl0[j0];
    if ((y01 = rl0[j1]) < ymin) ymin = y01;
    if (y01 > ymax) ymax = y01;
    rl1 = rl[i1];
    if ((y10 = rl1[j0]) < ymin) ymin = y10;
    if (y10 > ymax) ymax = y10;
    if ((y11 = rl1[j1]) < ymin) ymin = y11;
    if (y11 > ymax) ymax = y11;
    lmin = ceil(ymin-epsr);
    if (ep) lmin -= epsr/EPS_BAR;
    if (lmin < 0)  lmin = 0;
    lmax = floor(ymax+epsr);
    if (ep) lmax += epsr/EPS_BAR;
    if (lmax > my) lmax = my;
    first = (epsr == 0.0);
    rdk0 = rdk[i0];
    rdk1 = rdk[i1];
    rdl0 = rdl[i0];
    rdl1 = rdl[i1];
    for (k = kmin; k <= kmax; k++) {
        donex0 = donex[k];
        x = k;
        rdi0 = rdi[k];
        rdj0 = rdj[k];
	for (l = lmin; l <= lmax; l++) {
            if ((first) || (!donex0[l])) {
                y = l;
		if (get_bar(x,y,x00,y00,x01,y01,x10,y10,x11,y11,
	        &r00,&r01,&r10,&r11,(double) epsr)) {
                    donex0[l] += 1;
		    rdi0[l] += r00*rdk0[j0]+r01*rdk0[j1]+r10*rdk1[j0]+r11*rdk1[j1];
		    rdj0[l] += r00*rdl0[j0]+r01*rdl0[j1]+r10*rdl1[j0]+r11*rdl1[j1];
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------------*/
/* vshift computes vim1 which is the velocity field defined in frame m1               */
/* from vim0 which is a velocity field defined in frame m0 and scales it by sc.       */
/*------------------------------------------------------------------------------------*/
	
int vshift(VIMAGE *vim0, int fr, VIMAGE *vim1, int al, fdata m0, fdata m1, fdata sc)
{
    fdata dxmin,dxmax,dymin,dymax,ldx,ldy,leps,epsr,normij;
    fdata **rk,**rl,**rdk,**rdl,*rki,*rli,*rdki,*rdli;
    fdata ***d0,**dx0,*dx0i,**dy0,*dy0i,***d1,**dx1,*dx1i,**dy1,*dy1i,lambda;
    int i,j,nx,ny,nb,mx,my,ei,ej,il,ih,jl,jh,**donex,*donexi,ndone;
    lambda = m1-m0;
    nx = vim0->nx;
    ny = vim0->ny;
    nb = vim0->nb;
    d0 = vim0->d;
    dx0 = d0[0];
    dy0 = d0[1];
    mx = nx-1;
    my = ny-1;
    if (al) {
        if (malloc3((void ****) &d1,nb,nx,ny,sizeof(fdata)) == FAILURE) {
            fprintf(stderr,"ERROR: vshift: failed to allocate buffers\n");
            return(FAILURE);
        }
    } else {
        d1 = vim1->d;
    }
    dx1 = d1[0];
    dy1 = d1[1];
    if (m0 == m1) {
        for (i = 0; i < nx; i++) {
            donexi = donex[i];
            dx0i = dx0[i];
            dy0i = dy0[i];
            dx1i = dx1[i];
            dy1i = dy1[i];
            for (j = 0; j < ny; j++) {
                normij = sc/((fdata) donexi[j]);
                dx1i[j] = sc*dx0i[j];
                dy1i[j] = sc*dy0i[j];
            }
        }
    } else {
        if (malloc2((void ***) &donex,nx,ny,sizeof(int)) == FAILURE) {
            fprintf(stderr,"ERROR: vshift: failed to allocate buffer\n");
            return(FAILURE);        
        }
        dxmin = dxmax = dymin = dymax = 0.0;
        leps = 1+EPS_BAR;
        for (i = 0; i < nx; i++) {
            donexi = donex[i];
            dx0i = dx0[i];
            dy0i = dy0[i];
            for (j = 0; j < ny; j++) {
                donexi[j] = 0;
	        ldx = -lambda*dx0i[j];
	        if (ldx < dxmin) dxmin = ldx;
                if (ldx > dxmax) dxmax = ldx;
                ldy = -lambda*dy0i[j];
	        if (ldy < dymin) dymin = ldy;
                if (ldy > dymax) dymax = ldy;
            }
        }
        il = (int)(floor(dxmin-leps));
        ih = mx+(int)(ceil(dxmax+leps));
        jl = (int)(floor(dymin-leps));
        jh = my+(int)(ceil(dymax+leps));
        if ((altens2((void ***) &rk,il,ih,jl,jh,sizeof(fdata)) == FAILURE)
         || (altens2((void ***) &rl,il,ih,jl,jh,sizeof(fdata)) == FAILURE)
         || (altens2((void ***) &rdk,il,ih,jl,jh,sizeof(fdata)) == FAILURE)
         || (altens2((void ***) &rdl,il,ih,jl,jh,sizeof(fdata)) == FAILURE)) {
            fprintf(stderr,"ERROR: vshift: failed to allocate buffers\n");
            return(FAILURE);        
        }
        for (i = il; i <= ih; i++) {
            ei = i;
            if (ei < 0)  ei = 0;
            if (ei > mx) ei = mx;
            dx0i = dx0[ei];
            dy0i = dy0[ei];
            rki = rk[i];
            rli = rl[i];
            rdki = rdk[i];
            rdli = rdl[i];
            for (j = jl; j <= jh; j++) {
	        ej = j;
                if (ej < 0)  ej = 0;
                if (ej > my) ej = my;
	        rki[j] = i+lambda*dx0i[ej];
	        rli[j] = j+lambda*dy0i[ej];
                rdki[j] = dx0i[ej];
                rdli[j] = dy0i[ej];
            }
        }
        epsr = 0.0;
        ndone = 1;
        for (i = 0; i < nx; i++) {
            dx1i = dx1[i];
            dy1i = dy1[i];
            for (j = 0; j < ny; j++) {
                dx1i[j] = dy1i[j] = 0.0;
            }
        }
        while (ndone) {
	    for (i = il; i < ih; i++) {
	        for (j = jl; j < jh; j++) {
            	    fillquad(i,j,rk,rl,rdk,rdl,dx1,dy1,nx,ny,donex,epsr,0);
                }
            }
            ndone = nx*ny;
            for (i = 0; i < nx; i++) {
                donexi = donex[i];
                for (j = 0; j < ny; j++) {
                    if (donexi[j]) ndone--;
                }
            }
	    if (ndone) {
                if (epsr == 0.0) epsr = EPS_BAR;
                else epsr *= 2;
            }
        }
        for (i = 0; i < nx; i++) {
            donexi = donex[i];
            dx1i = dx1[i];
            dy1i = dy1[i];
            for (j = 0; j < ny; j++) {
                normij = sc/((fdata) donexi[j]);
                dx1i[j] *= normij;
                dy1i[j] *= normij;
            }
        }
        if ((free2((void ***) &donex,nx,ny,sizeof(int)) == FAILURE)
         || (frtens2((void ***) &rk,il,ih,jl,jh,sizeof(fdata)) == FAILURE)
         || (frtens2((void ***) &rl,il,ih,jl,jh,sizeof(fdata)) == FAILURE)
         || (frtens2((void ***) &rdk,il,ih,jl,jh,sizeof(fdata)) == FAILURE)
         || (frtens2((void ***) &rdl,il,ih,jl,jh,sizeof(fdata)) == FAILURE)) {
            fprintf(stderr,"ERROR: vshift: failed to free buffers\n");
            return(FAILURE);        
        }
    }
    if (fr) {
        if (free3((void ****) &d0,nb,nx,ny,sizeof(fdata)) == FAILURE) {
            fprintf(stderr,"ERROR: vshift: failed to free buffers\n");
            return(FAILURE);
        }
    }
    vim1->nx = nx;
    vim1->ny = ny;
    vim1->nb = nb;
    vim1->d = d1;
    vim1->dx = dx1;
    vim1->dy = dy1;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*  invdep routine: invert a displacement field.                           	      */
/*------------------------------------------------------------------------------------*/

int invdep(VIMAGE *vim, int fr, VIMAGE *ivim, int al)
{
    if (vshift(vim,fr,ivim,al,0.0,1.0,-1.0) == FAILURE) {
        fprintf(stderr,"ERROR: invdep: vshift failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*  compdep routine: compose two displacment fields.                           	      */
/*------------------------------------------------------------------------------------*/

float rim(float rx, float ry, float **im)
{
    int nx,ny;
    float px,py,qx,qy;
    nx = (int)floor(rx);
    ny = (int)floor(ry);
    qx = rx-(float)nx;
    qy = ry-(float)ny;
    px = 1.0-qx;
    py = 1.0-qy;
    return(py*(px*im[nx][ny]+qx*im[nx+1][ny])+qy*(px*im[nx][ny+1]+qx*im[nx+1][ny+1]));
}

float rimc(float rx, float ry, float **im, int wx, int wy)
{
    int nx,ny,mx,my;
    float px,py,qx,qy;
    nx = floor(rx);
    ny = floor(ry);
    mx = ceil(rx);
    my = ceil(ry);
    qx = rx-(float)nx;
    qy = ry-(float)ny;
    px = 1.0-qx;
    py = 1.0-qy;
    if (nx < 0) {nx = 0; px = 1.0; qx = 0.0; if (mx < 0) mx = 0;}
    if (mx > (wx-1)) {mx = wx-1; px = 0.0; qx = 1.0; if (nx > (wx-1)) nx = wx-1;}
    if (ny < 0) {ny = 0; py = 1.0; qy = 0.0; if (my < 0) my = 0;}
    if (my > (wy-1)) {my = wy-1; py = 0.0; qy = 1.0; if (ny > (wy-1)) ny = wy-1;}
    return(py*(px*im[nx][ny]+qx*im[mx][ny])+qy*(px*im[nx][my]+qx*im[mx][my]));
}

int compdep(VIMAGE *vim0, int fr0, VIMAGE *vim1, int fr1, VIMAGE *vim, int al)
{
    fdata ***d0,**dx0,**dy0,***d1,**dx1,**dy1,***d,**dx,**dy,x,y;
    int i,j,nx,ny,nb;
    nx = vim0->nx;
    ny = vim0->ny;
    nb = vim0->nb;
    d0 = vim0->d;
    dx0 = d0[0];
    dy0 = d0[1];
    d1 = vim1->d;
    dx1 = d1[0];
    dy1 = d1[1];
    if (al) {
        if (malloc3((void ****) &d,nb,nx,ny,sizeof(fdata)) == FAILURE) {
            fprintf(stderr,"ERROR: vshift: failed to allocate buffers\n");
            return(FAILURE);
        }
    } else {
        d = vim->d;
    }
    dx = d[0];
    dy = d[1];
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            x = i+dx0[i][j];
            y = j+dy0[i][j];
            dx[i][j] = dx0[i][j]+rimc(x,y,dx1,nx,ny);
            dy[i][j] = dy0[i][j]+rimc(x,y,dy1,nx,ny);
        }
    }
    if (fr0) {
        if (free3((void ****) &d0,nb,nx,ny,sizeof(fdata)) == FAILURE) {
            fprintf(stderr,"ERROR: vshift: failed to free buffers\n");
            return(FAILURE);
        }
    }
    if (fr1) {
        if (free3((void ****) &d1,nb,nx,ny,sizeof(fdata)) == FAILURE) {
            fprintf(stderr,"ERROR: vshift: failed to free buffers\n");
            return(FAILURE);
        }
    }
    vim->nb = nb;
    vim->nx = nx;
    vim->ny = ny;
    vim->nx0 = 0;
    vim->ny0 = 0;
    vim->nx1 = 0;
    vim->ny1 = 0;
    vim->d = d;
    vim->dx = d[0];
    vim->dy = d[1];
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

void vscale(VIMAGE *vim, fdata scale)
{
    fdata **dx,*dxi,**dy,*dyi;
    int i,j,nx,ny;
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    for (i = 0; i < nx; i++) {
        dxi = dx[i];
        dyi = dy[i];
        for (j = 0; j < ny; j++) {
            dxi[j] *= scale;
            dyi[j] *= scale;
        }
    }
}

/*------------------------------------------------------------------------------------*/

void vcolin(VIMAGE *vim0, VIMAGE *vim1, VIMAGE *vim, fdata l0, fdata l1)
{
    fdata **dx0,*dx0i,**dy0,*dy0i,**dx1,*dx1i,**dy1,*dy1i,**dx,*dxi,**dy,*dyi;
    int i,j,nx,ny;
    dx0 = vim0->dx;
    dy0 = vim0->dy;
    dx1 = vim1->dx;
    dy1 = vim1->dy;
    dx  = vim->dx;
    dy  = vim->dy;
    nx  = vim->nx;
    ny  = vim->ny;
    for (i = 0; i < nx; i++) {
        dx0i = dx0[i];
        dy0i = dy0[i];
        dx1i = dx1[i];
        dy1i = dy1[i];
        dxi  = dx[i];
        dyi  = dy[i];
        for (j = 0; j < ny; j++) {
            if ((dx0i[j] != 100.0) && (dy0i[j] != -100.0)
             && (dx1i[j] != 100.0) && (dy1i[j] != -100.0)) {
                dxi[j] = l0*dx0i[j]+l1*dx1i[j];
                dyi[j] = l0*dy0i[j]+l1*dy1i[j];
            } else {
                dxi[j] = 100.0;
                dyi[j] = -100.0;
            }
        }
    }
}

/*------------------------------------------------------------------------------------*/

int vsplit(VIMAGE *vim, VIMAGE *vim0, VIMAGE *vim1)
{
    if ((vshift(vim,0,vim0,1,0.5,0.0,1.0) == FAILURE)
     || (vshift(vim,0,vim1,1,0.5,1.0,-1.0) == FAILURE)) {
        fprintf(stderr,"ERROR: vsplit: vshift failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int vsplit2(VIMAGE *vim)
{
    /*   vim01 = vim+2;   vim10 = vim+3;   */
    if ((vshift(vim,0,vim+2,1,0.5,0.0,1.0) == FAILURE)
     || (vshift(vim,0,vim+3,1,0.5,1.0,1.0) == FAILURE)) {
        fprintf(stderr,"ERROR: vsplit2: vshift failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int vsplit3(VIMAGE *vim)
{
    /*   vim12 = vim+3;   vim20 = vim+4;  vim01 = vim+5;   */
    /*   vim21 = vim+6;   vim02 = vim+7;  vim10 = vim+8;   */
    if (vsplit(vim+0,vim+3,vim+6) == FAILURE) return(FAILURE);
    if (vsplit(vim+1,vim+4,vim+7) == FAILURE) return(FAILURE);
    if (vsplit(vim+2,vim+5,vim+8) == FAILURE) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/

int apply0(VIMAGE *vim, CIMAGE *iim, int f, DIMAGE *dim, fdata m, int bsi)
{
    /*---------------------------------------------------------------*/
    /* Builds the transform of an image using a given velocity field */
    /* Builds also a boolean image defining where the velocity field */
    /* points inside the image                                       */
    /*---------------------------------------------------------------*/
    int i,j,l,nx,ny,nvx,nvy,nb;
    fdata **dx,**dy,*dxi,*dyi,***xm,**xmi,rx,ry,*rim,s;
    uchar *imf,**in,*ini;
    s = m-((fdata) f);
    dx = vim->dx;
    dy = vim->dy;
    dim->nx = nvx = vim->nx;
    dim->ny = nvy = vim->ny;
    nx = iim->nx;
    ny = iim->ny;
    nb = iim->nb;
    dim->nb = nb;
    imf = iim->im;
    if ((malloc1((void **) &rim,nb,sizeof(fdata)) == FAILURE)
     || (malloc1((void **) &xmi,nb,sizeof(fdata *)) == FAILURE)
     || (malloc2((void ***) &in,nvx,nvy,sizeof(uchar)) == FAILURE)
     || (malloc3((void ****) &xm,nb,nvx,nvy,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: apply0: fail to allocate buffers\n");
        return(FAILURE);
    }
    for (i = 0; i < nvx; i++) {
        ini = in[i];
        dxi = dx[i];
        dyi = dy[i];
        for (l = 0; l < nb; l++) xmi[l] = xm[l][i];
        for (j = 0; j < nvy; j++) {
            rx = 0.5*(nx-nvx)+i+s*dxi[j];
            ry = 0.5*(ny-nvy)+j+s*dyi[j];
            ini[j] = bsi ? inter2(rx,ry,nx,ny,nb,imf,rim)
                         : inter1(rx,ry,nx,ny,nb,imf,rim);
            for (l = 0; l < nb; l++) xmi[l][j] = rim[l];
        }
    }
    if ((free1((void **) &rim,nb,sizeof(fdata)) == FAILURE)
     || (free1((void **) &xmi,nb,sizeof(fdata *)) == FAILURE)) {
        fprintf(stderr,"ERROR: apply0: fail to free buffer\n");
        return(FAILURE);
    }
    dim->xm = xm;
    dim->in = in;
    return(SUCCESS);
}

int ishift(CIMAGE *iim, CIMAGE *oim, VIMAGE *vim, int f, int fr, fdata m, fdata l1,
fdata noise)
{
    int i,j,l,nx,ny,nb,omlij;
    fdata l0,***xm0,***xm1;
    uchar *omd;
    DIMAGE dim0,dim1;
    l0 = 1.0-l1;
    nx = vim->nx;
    ny = vim->ny;
    nb = iim->nb;
    if (apply0(vim,iim+0,f,&dim0,m,0) == FAILURE) {
        fprintf(stderr,"ERROR: ishift: fail to transform images\n");
        return(FAILURE);
    }
    if (apply0(vim,iim+1,f,&dim1,m,0) == FAILURE) {
        fprintf(stderr,"ERROR: ishift: fail to transform images\n");
        return(FAILURE);
    }
    xm0 = dim0.xm;
    xm1 = dim1.xm;
    if (malloc1((void **) &omd,nb*nx*ny,sizeof(uchar)) == FAILURE) {
        fprintf(stderr,"ERROR: ishift: fail to allocate output buffer\n");
        return(FAILURE);
    }
    oim->nb = nb;
    oim->ny = ny;
    oim->nx = nx;
    oim->im = omd;
    noise *= 255.0/2147483648.0;
    for (l = 0; l < nb; l++) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                omlij = (int) rint(l0*xm0[l][i][j]+l1*xm1[l][i][j]+noise*random());
                if (omlij < 0)   omlij = 0;
                if (omlij > 255) omlij = 255;
                *omd++ = (uchar) omlij;
            }
        }
    }
    if ((free2((void ***) &dim0.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim0.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)
     || (free2((void ***) &dim1.in,nx,ny,sizeof(uchar)) == FAILURE)
     || (free3((void ****) &dim1.xm,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: ishift: fail to free buffers\n");
        return(FAILURE);
    }
    if (fr && ((free2((void ***) &vim->dx,nx,ny,sizeof(float)) == FAILURE)
     || (free2((void ***) &vim->dy,nx,ny,sizeof(float)) == FAILURE))) return(FAILURE);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/

void lininfo(fdata ax, fdata bx, fdata cx, int nx,
             fdata ay, fdata by, fdata cy, int ny,
             fdata *z, fdata *r)
{
    double mx,my,a,b,c,d,u,v,t,au,av,zm,zp,zs,tm,tp,pi=4.0*atan(1.0);
    mx = nx-1;
    my = ny-1;
    a = (bx*ax+by*ay);
    b = (ax*ax+ay*ay)-(bx*bx+by*by);
    c = -(ax*bx+ay*by);
    d = sqrt(b*b-4.0*a*c);
    if ((fabs(a) < 1.0e-6) || (fabs(c) < 1.0e-6)) {
        t = 0.0;
    } else if (a > c) {
        t = atan2(d-b,2.0*a);
    } else {
        t = atan2(2.0*c,d-b);
    }
    if (t < 0.0) t += pi;
    if (t > pi/2) t -= pi/2;
    u = cos(t);
    v = sin(t);
    printf("Check: avv+buv+cuu = %f\n",(float) (a*v*v+b*u*v+c*u*u));
    au = ax*u+bx*v;
    av = ay*u+by*v;
    *r = atan2(av,au)-t;
    zp = hypot(au,av);
    u = sin(t);
    v = -cos(t);
    au = ax*u+bx*v;
    av = ay*u+by*v;
    zm = hypot(au,av);
    if (zm > zp) {
        zs = zm;
        zm = zp;
        zp = zs;
        tm = t;
        tp = t+pi/2;
    } else {
        tm = t+pi/2;
        tp = t;
    }
    *z = (zm+zp)/2.0;
    printf("Center horizontal translation: %3.4f\n",(float) (cx+(ax-1.0)*mx/2+bx*my/2));
    printf("Center vertical translation:   %3.4f\n",(float) (cy+ay*mx/2+(by-1.0)*my/2));
    printf("Rotation angle:                %3.4f\n",(float) (*r));
    if (fabs(zm-zp) < 1.0e-4) {
        printf("Zoom factor:                   %3.4f\n",(float) (*z));
    } else {
        printf("Zoom factor 1:                 %3.4f\n",(float) (zp));
        printf("Zoom axis angle 1:             %3.4f\n",(float) (tp));
        printf("Zoom factor 2:                 %3.4f\n",(float) (zm));
        printf("Zoom axis angle 2:             %3.4f\n",(float) (tm));
    }
}

/*------------------------------------------------------------------------------------*/

void autonorm(VIMAGE *vim)
{
    /*------------------------------------------------------*/
    /* Reduces the velocity field to a linear approximation */
    /*------------------------------------------------------*/
    fdata **dx,**dy,ax,bx,cx,ay,by,cy,mx,my,z,r;
    fdata sdd,sii,sdx,sdy,sdi,sij,six,siy,sdj,sjj,sjx,sjy;
    int i,j,nx,ny;
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    mx = (fdata) (nx-1);
    my = (fdata) (ny-1);
    sdd = sdi = sdj = sii = sij = sjj = sdx = sdy = six = sjx = siy = sjy = 0.0;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            sdd += 1;
            sdi += i;
            sdj += j;
            sii += i*i;
            sij += i*j;
            sjj += j*j;
            sdx += dx[i][j];
            sdy += dy[i][j];
            six += i*dx[i][j];
            siy += i*dy[i][j];
            sjx += j*dx[i][j];
            sjy += j*dy[i][j];
        }
    }
    eqlin3(sdi,sdj,sdd,sdx,&ax,sii,sij,sdi,six,&bx,sij,sjj,sdj,sjx,&cx);
    eqlin3(sdi,sdj,sdd,sdy,&ay,sii,sij,sdi,siy,&by,sij,sjj,sdj,sjy,&cy);
    lininfo(ax+1.0,bx,cx,nx,ay,by+1.0,cy,ny,&z,&r);
    for (i = 0; i< nx; i++) {
        for (j = 0; j < ny; j++) {
            dx[i][j] = ax*i+bx*j+cx;
            dy[i][j] = ay*i+by*j+cy;
        }
    }
    printf("%f %f %f\n",(float) (1+ax),(float) (bx),(float) (cx+ax*mx/2+bx*my/2));
    printf("%f %f %f\n",(float) (ay),(float) (1+by),(float) (cy+ay*mx/2+by*my/2));
}

/*------------------------------------------------------------------------------------*/
