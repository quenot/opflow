/*------------------------------------------------------------------------------------*/
/*  Copyright 1995  Georges QUENOT  LIMSI-CNRS 					      */
/*  Version 2.00 Last revision: june 26 1996					      */
/*  Absolutely NO WARRANTY  							      */
/*------------------------------------------------------------------------------------*/

#include <stdio.h>
#include "glob.h"
#include "tens.h"
#include "imio.h"

#define INFINI 10000000

/*------------------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{   
    int i,j,nx,ny;
    fdata **dx,**dy,gamma,dvx,dvy,dvv,dv;
    VIMAGE vim;
    if (getfarg(&argc,argv,"-GAMMA",1.0,&gamma) == FAILURE) {
	fprintf(stderr,"Error uwogamma: failed to read scale\n");
        exit(1);
    }
    if (rdvimage(argv[1],&vim,UWO,1) == FAILURE) {
	fprintf(stderr,"Error uwogamma: failed to read file: %s\n",argv[1]);
        exit(1);
    }
    nx = vim.nx;
    ny = vim.ny;
    dx = vim.dx;
    dy = vim.dy;
    dv = 0.0;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            dvx = dx[i][j];
            dvy = dy[i][j];
            if ((dvx != 100.0) || (dvy != -100.0)) {
                if ((dvv = dvx*dvx+dvy*dvy) > dv) dv = dvv;
            }
        }
    }
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            dvx = dx[i][j];
            dvy = dy[i][j];
            if ((dvx != 100.0) || (dvy != -100.0)) {
                if ((dvv = dvx*dvx+dvy*dvy) > 0.0)
                    dvv = pow(sqrt((dvx*dvx+dvy*dvy)/dv),gamma-1);
                dx[i][j] = dvx*dvv;
                dy[i][j] = dvy*dvv;
            }
        }
    }
    if (wrvimage(argv[2],&vim,UWO,1) == FAILURE) {
	fprintf(stderr,"Error uwogamma: failed to write file: %s\n",argv[2]);
        exit(1);
    }
    exit(0);
}
