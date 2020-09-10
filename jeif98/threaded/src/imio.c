/*------------------------------------------------------------------------------------*/
/*  Copyright (C) 1995, 1996 and 1997  Georges QUENOT  LIMSI-CNRS                     */
/*  Copyright (C) 1998 and 1999        Georges QUENOT  CLIPS-IMAG                     */
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

/*------------------------------------------------------------------------------------*/

int filelen(char *fname, int *len)
{
    struct stat buf;
    if (stat(fname,&buf)) {
	fprintf(stderr,"Error filelen: can't stat file: %s\n",fname);
	return(FAILURE);
    }
    *len = (int) (buf.st_size);
    return(SUCCESS);
}

void bswap(void *p, size_t l)
{
    int i;
    unsigned *q;
    q = (unsigned *) p;
    for (i = 0; i < l; i++) q[i] = htonl(q[i]);
}

/*------------------------------------------------------------------------------------*/

void prflp(fdata f, char *str)
{
    int i;
    sprintf(str,"%f",(float) f);
    i = 0;
    while (str[i] != 0) i++;
    while ((--i >= 0) && (str[i] == '0')) str[i] = 0;
    if (str[i] == '.') str[i] = 0;
}

/*------------------------------------------------------------------------------------*/
/*                            PROGRAM ARGUMENTS EXTRACTION                            */
/*------------------------------------------------------------------------------------*/

int testbarg(int *pac, char *av[], char *name)
{
    /*----------------------------------------------------*/
    /* Extracts a bolean parameter from the argument list */
    /*----------------------------------------------------*/
    int i,barg;
    barg = 0;
    for (i = 1; i < *pac; i++) {
        if (strcmp(av[i],name) == 0) barg = 1;
    }
    return(barg);
}

int getbarg(int *pac, char *av[], char *name)
{
    /*----------------------------------------------------*/
    /* Extracts a bolean parameter from the argument list */
    /*----------------------------------------------------*/
    int i,j,barg;
    barg = 0;
    for (i = 1; i < *pac; i++) {
        if (strcmp(av[i],name) == 0) {
            barg = 1;
            for (j = i; j < *pac-1; j++) av[j] = av[j+1];
            *pac -= 1;
        }
    }
    return(barg);
}

int getcarg(int *pac, char *av[], char *name, char cdef, char *pcarg)
{
    /*-------------------------------------------------------*/
    /* Extracts a character parameter from the argument list */
    /*-------------------------------------------------------*/
    int i,j;
    *pcarg = cdef;
    for (i = 1; i < *pac-1; i++) {
        if (strcmp(av[i],name) == 0) {
            if (sscanf(av[i+1],"%c",pcarg) != 1) {
                fprintf(stderr,"ERROR: %s option expects a char value\n",name);
                return(FAILURE);
            }
            for (j = i; j < *pac-2; j++) av[j] = av[j+2];
            *pac -= 2;
        }
    }
    return(SUCCESS);
}

int getiarg(int *pac, char *av[], char *name, int idef, int *piarg)
{
    /*------------------------------------------------------*/
    /* Extracts an integer parameter from the argument list */
    /*------------------------------------------------------*/
    int i,j;
    *piarg = idef;
    for (i = 1; i < *pac-1; i++) {
        if (strcmp(av[i],name) == 0) {
            if (sscanf(av[i+1],"%d",piarg) != 1) {
                fprintf(stderr,"ERROR: %s option expects a valid integer\n",name);
                return(FAILURE);
            }
            for (j = i; j < *pac-2; j++) av[j] = av[j+2];
            *pac -= 2;
        }
    }
    return(SUCCESS);
}

int getfarg(int *pac, char *av[], char *name, fdata fdef, fdata *pfarg)
{
    /*---------------------------------------------------*/
    /* Extracts a float parameter from the argument list */
    /*---------------------------------------------------*/
    int i,j;
    float f;
    f = (float) fdef;
    for (i = 1;i < *pac-1; i++) {
        if (strcmp(av[i],name) == 0) {
            if (sscanf(av[i+1],"%f",&f) != 1) {
                fprintf(stderr,"ERROR: %s option expects a valid float\n",name);
                return(FAILURE);
            }
            for (j = i; j < *pac-2; j++) av[j] = av[j+2];
            *pac -= 2;
        }
    }
    *pfarg = (fdata) f;
    return(SUCCESS);
}

int getsarg(int *pac, char *av[], char *name, char *sarg)
{
    /*----------------------------------------------------*/
    /* Extracts a string parameter from the argument list */
    /*----------------------------------------------------*/
    int i,j;
    *sarg = 0;
    for (i = 1;i < *pac-1; i++) {
        if (strcmp(av[i],name) == 0) {
            if (sscanf(av[i+1],"%s",sarg) != 1) {
                fprintf(stderr,"ERROR: %s option expects a valid filename\n",name);
                return(FAILURE);
            }
            for (j = i; j < *pac-2; j++) av[j] = av[j+2];
            *pac -= 2;
        }
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                             CHARACTER IMAGE FILE READ                              */
/*------------------------------------------------------------------------------------*/

int rdpgmsize(FILE *fp, int *pnx, int *pny, int *pnb)
{
    int c,nc;
    if (getc(fp) != 'P') return(FAILURE);
    if (((c = getc(fp)) != '5') && (c!= '2'))  return(FAILURE);
    *pnb = 1;
    if ((c = getc(fp)) == EOF) return(FAILURE);
    while ((((char) c) == ' ') || (((char) c) == '\n') || (((char) c) == '\t')) {
        if ((c = getc(fp)) == EOF) return(FAILURE);
    }
    while (((char) c) == '#') {
        while (((char) c) != '\n') {
            if ((c = getc(fp)) == EOF) return(FAILURE);
        }
        if ((c = getc(fp)) == EOF) return(FAILURE);
    }
    ungetc((char) c,fp);
    if (fscanf(fp,"%d",pnx) != 1) return(FAILURE);
    if (fscanf(fp,"%d",pny) != 1) return(FAILURE);
    if (fscanf(fp,"%d",&nc) != 1) return(FAILURE);
    getc(fp);
    return(SUCCESS);
}

int rdppmsize(FILE *fp, int *pnx, int *pny, int *pnb)
{
    int c,nc;
    if (getc(fp) != 'P') return(FAILURE);
    if (((c = getc(fp)) != '6') && (c!= '3'))  return(FAILURE);
    *pnb = 3;
    if ((c = getc(fp)) == EOF) return(FAILURE);
    while ((((char) c) == ' ') || (((char) c) == '\n') || (((char) c) == '\t')) {
        if ((c = getc(fp)) == EOF) return(FAILURE);
    }
    while (((char) c) == '#') {
        while (((char) c) != '\n') {
            if ((c = getc(fp)) == EOF) return(FAILURE);
        }
        if ((c = getc(fp)) == EOF) return(FAILURE);
    }
    ungetc((char) c,fp);
    if (fscanf(fp,"%d",pnx) != 1) return(FAILURE);
    if (fscanf(fp,"%d",pny) != 1) return(FAILURE);
    if (fscanf(fp,"%d",&nc) != 1) return(FAILURE);
    getc(fp);
    return(SUCCESS);
}

int rdlumsize(FILE *fp, int *pnx, int *pny, int *pnb)
{
    int i;
    for (i = 0; i < 240; i++) getc(fp);
    *pnx = getc(fp);
    *pnx = *pnx+256*getc(fp);
    getc(fp);
    getc(fp);
    *pny = getc(fp);
    *pny = *pny+256*getc(fp);
    getc(fp);
    getc(fp);
    *pnb = getc(fp);
    *pnb = *pnb+256*getc(fp);
    return(SUCCESS);
}

int rdcalsize(FILE *fp, int *pnx, int *pny, int *pnb)
{
    uchar nn[4];
    if (fread(nn,sizeof(uchar),4,fp) != 4) {
	fprintf(stderr,"Error rdcalsize: can't read image size\n");
	return(FAILURE);
    }
    if ((*pnx = (((int) nn[0])<<8) + ((int) nn[1])) == 0) {
	fprintf(stderr,"Error rdcalsize: null image width\n");
	return(FAILURE);
    }
    if ((*pny = (((int) nn[2])<<8) + ((int) nn[3])) == 0) {
	fprintf(stderr,"Error rdcalsize: null image height\n");
	return(FAILURE);
    }
    if (fseek(fp,0,2) != 0) {
        fprintf(stderr,"ERROR: rdcalsize: can't fseek to end of file\n");
        return(FAILURE);
    }
    *pnb = (ftell(fp)-4)/((*pnx)*(*pny));
    return(SUCCESS);
}

int rdrassize(FILE *fp, int *pnx, int *pny, int *pnb)
{
    /*--------------------------------------------*/
    /* get size parameters from a sun raster file */
    /*--------------------------------------------*/
    int header[8];
    if (fread(header,sizeof(int),8,fp) != 8) {
        fprintf(stderr,"ERROR: rdrassize: fail to read header\n");
        return(FAILURE);
    }
    bswap((void *) header, 8);
    *pnx = header[1];
    *pny = header[2];
    *pnb = 1;
    return(SUCCESS);
}

int getrawsize(int *pac, char *av[], int *pnx, int *pny, int *pnb, int *pfo)
{
    /*-----------------------------------------------------------*/
    /* get raw images size and number of bands from command line */
    /*-----------------------------------------------------------*/
    int i,j;
    for (i = 1; i< *pac-2; i++) {
        if (strcmp(av[i],"-B") == 0) {
            *pfo = RAW;
            sscanf(av[i+1],"%d",pnx);
            if (*pnx <= 0) {
                fprintf(stderr,"ERROR: invalid image width: %d\n",*pnx);
                return(FAILURE);
            }
            sscanf(av[i+2],"%d",pny);
            if (*pny <= 0) {
                fprintf(stderr,"ERROR: invalid image height: %d\n",*pny);
                return(FAILURE);
            }
            if ((i+3 < *pac) && (av[i+3][0] != '-')) {
                sscanf(av[i+3],"%d",pnb);
                if (*pnb <= 0) {
                    fprintf(stderr,"ERROR: invalid number of bands: %d\n",*pnb);
                    return(FAILURE);
                }
                for (j = i; j < *pac-4; j++) av[j] = av[j+4];
                *pac -= 4;
            }
            else {
                *pnb = 1;
                for (j = i; j < *pac-3; j++) av[j] = av[j+3];
                *pac -= 3;
            }
        }
    }
    return(SUCCESS);
}

int getrassize(char *fname, int *pnx, int *pny, int *pnb)
{
    /*--------------------------------------------*/
    /* get size parameters from a sun raster file */
    /*--------------------------------------------*/
    FILE *fp;
    int header[8];
    if ((fp=fopen(fname,"rb")) == NULL) {
        fprintf(stderr,"ERROR: getrastsize: can't open file: %s\n",fname);
        return(FAILURE);
    }
    if (fread(header,sizeof(int),8,fp) != 8) {
        fprintf(stderr,"ERROR: getrastsize: fail to read header: %s\n",fname);
        return(FAILURE);
    }
    fclose(fp);
    bswap((void *) header, 8);
    if (header[0] == RAS_MAGIC) {
        *pnx = header[1];
        *pny = header[2];
        *pnb = 1;
    } else {
        *pnx = *pny = *pnb = 0;
    }
    return(SUCCESS);
}

int getpgmsize(char *fname, int *pnx, int *pny, int *pnb)
{
    char fnamepgm[100];
    FILE *fp;
    sprintf(fnamepgm,"%s%s",fname,".pgm");
    if ((fp=fopen(fnamepgm,"rb")) == NULL) {
        fprintf(stderr,"ERROR: getpgmsize: can't open file: %s\n",fnamepgm);
        return(FAILURE);
    }
    if (rdpgmsize(fp,pnx,pny,pnb) == FAILURE) {
        fprintf(stderr,"ERROR: getpgmsize: can't read image size: %s\n",fnamepgm);
        return(FAILURE);
    }
    fclose(fp);
    return(SUCCESS);
}

int getppmsize(char *fname, int *pnx, int *pny, int *pnb)
{
    char fnameppm[100];
    FILE *fp;
    sprintf(fnameppm,"%s%s",fname,".ppm");
    if ((fp=fopen(fnameppm,"rb")) == NULL) {
        fprintf(stderr,"ERROR: getppmsize: can't open file: %s\n",fnameppm);
        return(FAILURE);
    }
    if (rdppmsize(fp,pnx,pny,pnb) == FAILURE) {
        fprintf(stderr,"ERROR: getppmsize: can't read image size: %s\n",fnameppm);
        return(FAILURE);
    }
    fclose(fp);
    return(SUCCESS);
}

int getlumsize(char *fname, int *pnx, int *pny, int *pnb)
{
    FILE *fp;
    if ((fp=fopen(fname,"rb")) == NULL) {
        fprintf(stderr,"ERROR: getlumsize: can't open file: %s\n",fname);
        return(FAILURE);
    }
    if (rdlumsize(fp,pnx,pny,pnb) == FAILURE) {
        fprintf(stderr,"ERROR: getlumsize: can't read image size: %s\n",fname);
        return(FAILURE);
    }
    fclose(fp);
    *pnb = 1;
    return(SUCCESS);
}

int getcalsize(char *fname, int *pnx, int *pny, int *pnb)
{
    FILE *fp;
    if ((fp=fopen(fname,"rb")) == NULL) {
        fprintf(stderr,"ERROR: getcalsize: can't open file: %s\n",fname);
        return(FAILURE);
    }
    if (rdcalsize(fp,pnx,pny,pnb) == FAILURE) {
        fprintf(stderr,"ERROR: getcalsize: can't read image size: %s\n",fname);
        return(FAILURE);
    }
    fclose(fp);
    return(SUCCESS);
}

int getimsize(int *pac, char *av[], char *fname, int *pnx, int *pny, int *pnb, int *pfo)
{
    /*----------------------------------*/
    /* get images size, number of bands */
    /*----------------------------------*/
    *pfo = RAS;
    if (getrawsize(pac,av,pnx,pny,pnb,pfo) == FAILURE) return(FAILURE);
    if (getbarg(pac,av,"-PGM")) {
        *pfo = PGM;
        if (getpgmsize(fname,pnx,pny,pnb) == FAILURE) return(FAILURE);
    }
    if (getbarg(pac,av,"-PPM")) {
        *pfo = PPM;
        if (getppmsize(fname,pnx,pny,pnb) == FAILURE) return(FAILURE);
    }
    if (getbarg(pac,av,"-LUM")) {
        *pfo = LUM;
        if (getlumsize(fname,pnx,pny,pnb) == FAILURE) return(FAILURE);
    }
    /*----------------------------------------------------------------------*/
    /* if -B -PGM and -PPM option are not used sun raster format is assumed */
    /*----------------------------------------------------------------------*/
    if ((*pfo == RAS) && (getrassize(fname,pnx,pny,pnb) == FAILURE)) return(FAILURE);
    /*----------------------------------------------------------------------*/
    /* if sun raster format is not found calfonc format is assumed          */
    /*----------------------------------------------------------------------*/
    if (*pnb == 0) {
        if (getcalsize(fname,pnx,pny,pnb) == FAILURE) return(FAILURE);
        else *pfo = CCF;
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int alcimage(CIMAGE *cim)
{
    int nx,ny,nb;
    uchar *im;
    nx = cim->nx;
    ny = cim->ny;
    nb = cim->nb;
    if (malloc1((void **) &im,nx*ny*nb,sizeof(uchar)) == FAILURE) {
        fprintf(stderr,"ERROR: alcimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    cim->im = im;
    return(SUCCESS);
}

int rdcimage(char *name, CIMAGE *cim, int fo, int al)
{
    FILE *fp;
    int i,nx,ny,nb,asc,c;
    uchar *im,*imr,*img,*imb;
    char fname[256],ext[6];
    strcpy(ext,"");
    if (fo == PGM) strcpy(ext,".pgm");
    if (fo == PPM) strcpy(ext,".ppm");
    sprintf(fname,"%s%s",name,ext);
    if ((fp = fopen(fname,"rb")) == NULL) {
        fprintf(stderr,"ERROR: rdcimage: can't open file: %s\n",fname);
        return(FAILURE);
    }
    if (fo == PGM) {
        getc(fp);
        asc = (getc(fp) == '2');
        fseek(fp,0,0);
    } else if (fo == PPM) {
        getc(fp);
        asc = (getc(fp) == '3');
        fseek(fp,0,0);
    } else {
        asc = 0;
    }
    if (((fo == PGM) && (rdpgmsize(fp,&nx,&ny,&nb) == FAILURE))
     || ((fo == PPM) && (rdppmsize(fp,&nx,&ny,&nb) == FAILURE))
     || ((fo == CCF) && (rdcalsize(fp,&nx,&ny,&nb) == FAILURE))
     || ((fo == RAS) && (rdrassize(fp,&nx,&ny,&nb) == FAILURE))) {
        fprintf(stderr,"ERROR: rdcimage: can't read image size: %s\n",fname);
        return(FAILURE);
    }
    if (fo == RAW) {
        nx = cim->nx;
        ny = cim->ny;
        nb = cim->nb;
    } else {
        cim->nx = nx;
        cim->ny = ny;
        cim->nb = nb;
    }
    if (al && (alcimage(cim) == FAILURE)) {
        fprintf(stderr,"ERROR: rdcimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    im = cim->im;
    if ((!asc) && (fseek(fp,-nx*ny*nb,2) != 0)) {
        fprintf(stderr,"ERROR: rdcimage: can't fseek to data: %s\n",fname);
        return(FAILURE);
    }
    if ((asc) && (fo == PGM)) {
        for (i = 0; i < nx*ny; i++) {
            fscanf(fp,"%d",&c); *im++ = c;
        }
    } else if ((asc) && (fo == PPM)) {
        imr = im;
        img = imr+nx*ny;
        imb = img+nx*ny;
        for (i = 0; i < nx*ny; i++) {
            fscanf(fp,"%d",&c); *imr++ = c;
            fscanf(fp,"%d",&c); *img++ = c;
            fscanf(fp,"%d",&c); *imb++ = c;
        }
    } else if (fo == PPM) {
        imr = im;
        img = imr+nx*ny;
        imb = img+nx*ny;
        for (i = 0; i < nx*ny; i++) {
            *imr++ = getc(fp);
            *img++ = getc(fp);
            *imb++ = getc(fp);
        }
    } else if (fread(im,sizeof(uchar),nx*ny*nb,fp) != nx*ny*nb) {
        fprintf(stderr,"ERROR: rdcimage: fail to read data: %s\n",fname);
        return(FAILURE);
    }
    fclose(fp);
    return(SUCCESS);
}

int rdlimage(char *name, CIMAGE *cim, int f, int il, int al)
{
    FILE *fp;
    int j,j0,nx,ny,ni;
    uchar *im;
    if ((fp = fopen(name,"rb")) == NULL) {
        fprintf(stderr,"ERROR: rdlimage: can't open file: %s\n",name);
        return(FAILURE);
    }
    if (rdlumsize(fp,&nx,&ny,&ni) == FAILURE) {
        fprintf(stderr,"ERROR: rdlimage: can't read image size: %s\n",name);
        return(FAILURE);
    }
    cim->nb = 1;
    if (al && (alcimage(cim) == FAILURE)) {
        fprintf(stderr,"ERROR: rdlimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    im = cim->im;
    ny *= (2-il);
    ni /= (2-il);
    if (fseek(fp,-nx*ny*(ni-f),2) != 0) {
        fprintf(stderr,"ERROR: rdlimage: can't fseek to data: %s\n",name);
        return(FAILURE);
    }
    for (j = 0; j < ny; j++) {
        j0 = il ? j : ((2*j)%ny)+((2*j)/ny);
        if (fread(im+j0*nx,sizeof(uchar),nx,fp) != nx) {
            fprintf(stderr,"ERROR: rdlimage: fail to read data: %s\n",name);
            return(FAILURE);
        }
    }
    fclose(fp);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int rdsimage(char *sname, SIMAGE *sim, int fo)
{
    /*----------------------------------------------*/
    /* read sequence image data from a set of files */
    /*----------------------------------------------*/
    int f,nx,ny,nb,nf,nl,*vf;
    uchar **im;
    char fname[256],ext[6];
    CIMAGE cim;
    strcpy(ext,"");
    if (fo == PGM) strcpy(ext,".pgm");
    if (fo == PPM) strcpy(ext,".ppm");
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    if (altens1((void **) (&vf),nf,nl,sizeof(int)) == FAILURE) {
        fprintf(stderr,"ERROR: rdsimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    sim->vf = vf;
    if (altens1((void **) (&im),nf,nl,sizeof(uchar *)) == FAILURE) {
        fprintf(stderr,"ERROR: rdsimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    sim->im = im;
    if (fo == RAW) {
        cim.nx = nx;
        cim.ny = ny;
        cim.nb = nb;
    }
    for (f = nf; f <= nl; f++) {
        sprintf(fname,"%s%d",sname,f);
        if (rdcimage(fname,&cim,fo,1) == FAILURE) {
            sprintf(fname,"%s%d%s",sname,f,ext);
            fprintf(stderr,"ERROR: rdsimage: fail to read data: %s\n",fname);
            return(FAILURE);
        }
        im[f] = cim.im;
        vf[f] = 1;
    }
    if (malloc1((void **) (&sim->val),nb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: rdsimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int frsimage(SIMAGE *sim)
{
    int f,nx,ny,nb,nf,nl,*vf;
    uchar **im;
    fdata *val;
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    vf = sim->vf;
    im = sim->im;
    val = sim->val;
    for (f = nf; f <= nl; f++) {
        if (free1((void **) (im+f),nx*ny*nb,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"ERROR: freesequence: fail to free buffer\n");
            return(FAILURE);
        }
    }
    if (free1((void **) (&val),nb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: freesequence: fail to free buffer\n");
        return(FAILURE);
    }
    if (frtens1((void **) (&im),nf,nl,sizeof(uchar *)) == FAILURE) {
        fprintf(stderr,"ERROR: freesequence: fail to free buffer\n");
        return(FAILURE);
    }
    if (frtens1((void **) (&vf),nf,nl,sizeof(int)) == FAILURE) {
        fprintf(stderr,"ERROR: freesequence: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                     Horizontal shrink/expand of image sequences                    */
/*------------------------------------------------------------------------------------*/

int simxshrink(SIMAGE *sim, uchar ***pwim)
{
    int i,f,nx,ny,nb,nf,nl,*vf;
    uchar **wim,*wimf,*himf;
    nx = sim->nx /= 2;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    vf = sim->vf;
    *pwim = wim = sim->im;
    if (altens1((void **) (&sim->im),nf,nl,sizeof(uchar *)) == FAILURE) {
        fprintf(stderr,"ERROR: simxshrink: fail to allocate buffer\n");
        return(FAILURE);
    }
    for (f = nf; f <= nl; f++) {
        if (vf[f]) {
            if (malloc1((void **) (&sim->im[f]),nx*ny*nb,sizeof(uchar)) == FAILURE) {
                fprintf(stderr,"ERROR: simxshrink: fail to allocate buffer\n");
                return(FAILURE);
            }
            wimf = wim[f];
            himf = sim->im[f];
            for (i = 0; i < nx*ny*nb; i++) {
                himf[i] = ((((int) (wimf[2*i]))+((int) (wimf[2*i+1]))))/2;
            }
        }
    }
    return(SUCCESS);
}

int simxexpand(SIMAGE *sim, uchar **wim)
{
    int f,nx,ny,nb,nf,nl,*vf;
    nx = sim->nx;
    ny = sim->ny;
    nb = sim->nb;
    nf = sim->nf;
    nl = sim->nl;
    vf = sim->vf;
    for (f = nf; f <= nl; f++) {
        if (vf[f]) {
            if (free1((void **) (&sim->im[f]),nx*ny*nb,sizeof(uchar)) == FAILURE) {
                fprintf(stderr,"ERROR: simxexpand: fail to allocate buffer\n");
                return(FAILURE);
            }
        }
    }
    if (frtens1((void **) (&sim->im),nf,nl,sizeof(uchar *)) == FAILURE) {
        fprintf(stderr,"ERROR: simxexpand: fail to allocate buffer\n");
        return(FAILURE);
    }
    sim->nx *= 2;
    sim->im = wim;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                   Add / remove smoothed bands to image sequences                   */
/*------------------------------------------------------------------------------------*/

int smocim(CIMAGE *cim0, CIMAGE *cim1)
{
    int i,j,k,nx,ny,nb,imm,icm,ipm,imc,icc,ipc,imp,icp,ipp,jm,jc,jp;
    uchar *im,*tm;
    nx = cim0->nx;
    ny = cim0->ny;
    nb = cim0->nb;
    tm = cim0->im;
    if (cim0->im == cim1->im) {
        if (malloc1((void **) &im,nx*ny*nb,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"Error smocim: Failed to allocate buffer\n");
            return(FAILURE);
        }
    } else {
        im = cim1->im;
    }
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny; j++) {
            if ((jm = j-1) < 0) jm = 0;
            jc = j;
            if ((jp = j+1) > (ny-1)) jp = ny-1;
            jm *= nx;
            jc *= nx;
            jp *= nx;
            for (i = 0; i < nx; i++) {
                if ((imc = i-1) < 0) imc = 0;
                icc = i;
                if ((ipc = i+1) > (nx-1)) ipc = nx-1;
                imm = imc+jm;
                icm = icc+jm;
                ipm = ipc+jm;
                imp = imc+jp;
                icp = icc+jp;
                ipp = ipc+jp;
                imc = imc+jc;
                icc = icc+jc;
                ipc = ipc+jc;
                im[icc] =  ((((int)(tm[imm]))+2*((int)(tm[imc]))+((int)(tm[imp])))
                        + 2*(((int)(tm[icm]))+2*((int)(tm[icc]))+((int)(tm[icp])))
                        +   (((int)(tm[ipm]))+2*((int)(tm[ipc]))+((int)(tm[ipp]))))/16;
            }
        }
        tm += nx*ny;
        im += nx*ny;
    }
    tm -= nx*ny*nb;
    im -= nx*ny*nb;
    if (cim0->im == cim1->im) {
        if (free1((void **) &tm,nx*ny*nb,sizeof(uchar)) == FAILURE) {
            fprintf(stderr,"Error smocim: Failed to free buffer\n");
            return(FAILURE);
        }
        cim1->im = im;
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                             CHARACTER IMAGE FILE WRITE                             */
/*------------------------------------------------------------------------------------*/

int wrpgmsize(FILE *fp, int nx, int ny, int nb)
{
    if (fprintf(fp,"P5\n%d %d\n255\n",nx,ny) == EOF) return(FAILURE);
    return(SUCCESS);
}

int wrppmsize(FILE *fp, int nx, int ny, int nb)
{
    if (fprintf(fp,"P6\n%d %d\n255\n",nx,ny) == EOF) return(FAILURE);
    return(SUCCESS);
}

int wrcalsize(FILE *fp, int nx, int ny, int nb)
{
    uchar nn[4];
    nn[0] = (uchar) (nx>>8); nn[1] = (uchar) (nx&0xFF);
    nn[2] = (uchar) (ny>>8); nn[3] = (uchar) (ny&0xFF);
    if (fwrite(nn,sizeof(uchar),4,fp) != 4) {
	fprintf(stderr,"Error wrcalsize: can't write image size\n");
	return(FAILURE);
    }
    return(SUCCESS);
}

int wrrassize(FILE *fp, int nx, int ny, int nb)
{
    int header[8];
    header[0] = 0x59a66a95;
    header[1] = nx;
    header[2] = ny;
    header[3] = 0x00000008;
    header[4] = 0x00013710;
    header[5] = 0x00000001;
    header[6] = 0x00000000;
    header[7] = 0x00000000;
    bswap((void *) header, 8);
    if (fwrite(header,sizeof(int),8,fp) != 8) {
        fprintf(stderr,"ERROR: wrrassize: fail to write header\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int frcimage(CIMAGE *cim)
{
    int nx,ny,nb;
    uchar *im;
    nx = cim->nx;
    ny = cim->ny;
    nb = cim->nb;
    im = cim->im;
    if (free1((void **) &im,nx*ny*nb,sizeof(uchar)) == FAILURE) {
        fprintf(stderr,"ERROR: frcimage: fail to free buffer\n");
        return(FAILURE);
    }
    cim->im = NULL;
    return(SUCCESS);
}

int wrcimage(char *name, CIMAGE *cim, int fo, int fr)
{
    FILE *fp;
    int i,nx,ny,nb;
    uchar *im,*imr,*img,*imb;
    char fname[256],ext[6];
    strcpy(ext,"");
    if (fo == PGM) strcpy(ext,".pgm");
    if (fo == PPM) strcpy(ext,".ppm");
    sprintf(fname,"%s%s",name,ext);
    nx = cim->nx;
    ny = cim->ny;
    nb = cim->nb;
    im = cim->im;
    if ((fp = fopen(fname,"wb")) == NULL) {
	fprintf(stderr,"Error wrcimage: can't open file: %s\n",fname);
	return(FAILURE);
    }
    if (((fo == PGM) && (wrpgmsize(fp,nx,ny,nb) == FAILURE))
     || ((fo == PPM) && (wrppmsize(fp,nx,ny,nb) == FAILURE))
     || ((fo == CCF) && (wrcalsize(fp,nx,ny,nb) == FAILURE))
     || ((fo == RAS) && (wrrassize(fp,nx,ny,nb) == FAILURE))) {
        fprintf(stderr,"ERROR: wrcimage: can't write image size: %s\n",fname);
        return(FAILURE);
    }
    if (fo == PPM) {
        imr = im;
        img = imr+nx*ny;
        imb = img+nx*ny;
        for (i = 0; i < nx*ny; i++) {
            putc((char) *imr++,fp);
            putc((char) *img++,fp);
            putc((char) *imb++,fp);
        }
    } else if (fwrite(im,sizeof(uchar),nx*ny*nb,fp) != nx*ny*nb) {
        fprintf(stderr,"ERROR: wrcimage: fail to write data: %s\n",fname);
        return(FAILURE);
    }
    fclose(fp);
    if (fr && (frcimage(cim) == FAILURE)) {
	fprintf(stderr,"Error wrcimage: fail to free buffer: %s\n",fname);
	return(FAILURE);
    }
    return(SUCCESS);
}

int wrcimagec(char *namec, char *name, CIMAGE *cim, int fo, int fr)
{
    FILE *fpc,*fp;
    int i,nx,ny,nb,len;
    uchar *im,*imr,*img,*imb;
    char fnamec[256],fname[256],ext[6];
    strcpy(ext,"");
    if (fo == PGM) strcpy(ext,".pgm");
    if (fo == PPM) strcpy(ext,".ppm");
    sprintf(fnamec,"%s%s",namec,ext);
    sprintf(fname,"%s%s",name,ext);
    nx = cim->nx;
    ny = cim->ny;
    nb = cim->nb;
    im = cim->im;
    if (filelen(fnamec,&len) == FAILURE) {
        fprintf(stderr,"ERROR: wrcimagec: can't get file size: %s\n",fnamec);
        return(FAILURE);
    }
    if ((fpc = fopen(fnamec,"rb")) == NULL) {
	fprintf(stderr,"Error wrcimagec: can't open file: %s\n",fnamec);
	return(FAILURE);
    }
    if ((fp = fopen(fname,"wb")) == NULL) {
	fprintf(stderr,"Error wrcimagec: can't open file: %s\n",fname);
	return(FAILURE);
    }
    for (len -= nx*ny*nb; len > 0; len--) putc(getc(fpc),fp);
    fclose(fpc);
    if (fo == PPM) {
        imr = im;
        img = imr+nx*ny;
        imb = img+nx*ny;
        for (i = 0; i < nx*ny; i++) {
            putc((char) *imr++,fp);
            putc((char) *img++,fp);
            putc((char) *imb++,fp);
        }
    } else if (fwrite(im,sizeof(uchar),nx*ny*nb,fp) != nx*ny*nb) {
        fprintf(stderr,"ERROR: wrcimage: fail to write data: %s\n",fname);
        return(FAILURE);
    }
    fclose(fp);
    if (fr && (frcimage(cim) == FAILURE)) {
	fprintf(stderr,"Error wrcimage: fail to free buffer: %s\n",fname);
	return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                              INTEGER IMAGE FILE READ                               */
/*------------------------------------------------------------------------------------*/

int aliimage(IIMAGE *iim)
{
    int nx,ny,nb,*im;
    nx = iim->nx;
    ny = iim->ny;
    nb = iim->nb;
    if (malloc1((void **) &im,nx*ny*nb,sizeof(int)) == FAILURE) {
        fprintf(stderr,"ERROR: aliimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    iim->im = im;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                              INTEGER IMAGE FILE WRITE                              */
/*------------------------------------------------------------------------------------*/

int friimage(IIMAGE *iim)
{
    int nx,ny,nb,*im;
    nx = iim->nx;
    ny = iim->ny;
    nb = iim->nb;
    im = iim->im;
    if (free1((void **) &im,nx*ny*nb,sizeof(int)) == FAILURE) {
        fprintf(stderr,"ERROR: friimage: fail to free buffer\n");
        return(FAILURE);
    }
    iim->im = NULL;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                          COLOR to GRAY and GRAY to COLOR                           */
/*------------------------------------------------------------------------------------*/

int coltogra(CIMAGE *cim3, int fr, CIMAGE *cim1, int al)
{
    int i,nx,ny;
    uchar *im,*imr,*img,*imb;
    cim1->nx = nx = cim3->nx;
    cim1->ny = ny = cim3->ny;
    if (cim3->nb != 3) {
        fprintf(stderr,"ERROR: coltogra: input image must have 3 band\n");
        return(FAILURE);
    }
    cim1->nb = 1;
    imr = cim3->im;
    img = imr+nx*ny;
    imb = img+nx*ny;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: coltogra: fail to allocate buffer\n");
        return(FAILURE);
    }
    im  = cim1->im;
    for (i = 0; i < nx*ny; i++)
        *im++ = (uchar)(0.5+0.30*(*imr++)+0.51*(*img++)+0.19*(*imb++));
    if (fr && (frcimage(cim3) == FAILURE)) {
        fprintf(stderr,"ERROR: coltogra: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int gratocol(CIMAGE *cim1, int fr, CIMAGE *cim3, int al)
{
    int i,nx,ny;
    uchar *im,*imr,*img,*imb;
    cim3->nx = nx = cim1->nx;
    cim3->ny = ny = cim1->ny;
    if (cim1->nb != 1) {
        fprintf(stderr,"ERROR: gratocol: input image must have 1 band\n");
        return(FAILURE);
    }
    cim3->nb = 3;
    im  = cim1->im;
    if (al && (alcimage(cim3) == FAILURE)) {
        fprintf(stderr,"ERROR: gratocol: fail to allocate buffer\n");
        return(FAILURE);
    }
    imr = cim3->im;
    img = imr+nx*ny;
    imb = img+nx*ny;
    for (i = 0; i < nx*ny; i++) *imr++ = *img++ = *imb++ = *im++;
    if (fr && (frcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: gratocol: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

void hpscim(CIMAGE *cim0, CIMAGE *cim1, fdata v)
{
    int i,j,k,nx,ny,nb,imm,icm,ipm,imc,icc,ipc,imp,icp,ipp,jm,jc,jp;
    uchar *im,*tm;
    nx = cim0->nx;
    ny = cim0->ny;
    nb = cim0->nb;
    tm = cim0->im;
    im = cim1->im;
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny; j++) {
            if ((jm = j-1) < 0) jm = 0;
            jc = j;
            if ((jp = j+1) > (ny-1)) jp = ny-1;
            jm *= nx;
            jc *= nx;
            jp *= nx;
            for (i = 0; i < nx; i++) {
                if ((imc = i-1) < 0) imc = 0;
                icc = i;
                if ((ipc = i+1) > (nx-1)) ipc = nx-1;
                imm = imc+jm;
                icm = icc+jm;
                ipm = ipc+jm;
                imp = imc+jp;
                icp = icc+jp;
                ipp = ipc+jp;
                imc = imc+jc;
                icc = icc+jc;
                ipc = ipc+jc;
                im[icc] =  ((1.0*((fdata)(tm[imm]))+2.0*((fdata)(tm[imc]))+1.0*((fdata)(tm[imp])))
                        +   (2.0*((fdata)(tm[icm]))+0.0*((fdata)(tm[icc]))+2.0*((fdata)(tm[icp])))
                        +   (1.0*((fdata)(tm[ipm]))+2.0*((fdata)(tm[ipc]))+1.0*((fdata)(tm[ipp]))))/12.0;
                im[icc] = (((fdata)(tm[icc]))+v*(255.0-((fdata)(im[icc]))))/(1+v);
            }
        }
        tm += nx*ny;
        im += nx*ny;
    }
}

int grasmcol(CIMAGE *cim1, int fr, CIMAGE *cim3, int al)
{
    int i,nx,ny;
    uchar *im,*imr,*img,*imb;
    cim3->nx = nx = cim1->nx;
    cim3->ny = ny = cim1->ny;
    if (cim1->nb != 1) {
        fprintf(stderr,"ERROR: grasmcol: input image must have 1 band\n");
        return(FAILURE);
    }
    cim3->nb = 3;
    im  = cim1->im;
    if (al && (alcimage(cim3) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to allocate buffer\n");
        return(FAILURE);
    }
    imr = cim3->im;
    img = imr+nx*ny;
    imb = img+nx*ny;
    for (i = 0; i < nx*ny; i++) img[i] = im[i];
    /* smocim(cim1); */
    for (i = 0; i < nx*ny; i++) imr[i] = im[i];
    /* smocim(cim1); */
    for (i = 0; i < nx*ny; i++) imb[i] = im[i];
    if (fr && (frcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int graoocol(CIMAGE *cim1, int fr, CIMAGE *cim3, int al, fdata th1, fdata th2)
{
    int i,j,nx,ny,imm,icm,ipm,imc,icc,ipc,imp,icp,ipp,jm,jc,jp;
    uchar *im,*imr,*img,*imb,cmp;
    fdata tmp;
    cim3->nx = nx = cim1->nx;
    cim3->ny = ny = cim1->ny;
    if (cim1->nb != 1) {
        fprintf(stderr,"ERROR: grasmcol: input image must have 1 band\n");
        return(FAILURE);
    }
    cim3->nb = 3;
    im  = cim1->im;
    if (al && (alcimage(cim3) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to allocate buffer\n");
        return(FAILURE);
    }
    imr = cim3->im;
    img = imr+nx*ny;
    imb = img+nx*ny;
    for (j = 0; j < ny; j++) {
        if ((jm = j-1) < 0) jm = 0;
        jc = j;
        if ((jp = j+1) > (ny-1)) jp = ny-1;
        jm *= nx;
        jc *= nx;
        jp *= nx;
        for (i = 0; i < nx; i++) {
            if ((imc = i-1) < 0) imc = 0;
            icc = i;
            if ((ipc = i+1) > (nx-1)) ipc = nx-1;
            imm = imc+jm;
            icm = icc+jm;
            ipm = ipc+jm;
            imp = imc+jp;
            icp = icc+jp;
            ipp = ipc+jp;
            imc = imc+jc;
            icc = icc+jc;
            ipc = ipc+jc;
            tmp = ((1.0*((fdata)(im[imm]))+ 2.0*((fdata)(im[imc]))+1.0*((fdata)(im[imp])))
                +  (2.0*((fdata)(im[icm]))-12.0*((fdata)(im[icc]))+2.0*((fdata)(im[icp])))
                +  (1.0*((fdata)(im[ipm]))+ 2.0*((fdata)(im[ipc]))+1.0*((fdata)(im[ipp]))))/3060.0;
            if (tmp <= -th2) {
                imr[icc] = 255;
                imb[icc] = 0;
            } else if (tmp < -th1) {
                imr[icc] = floor(0.5+127.5*(1-cos(M_PI*(tmp+th1)/(th1-th2))));
                imb[icc] = 0;
            } else if (tmp <= th1) {
                imr[icc] = 0;
                imb[icc] = 0;
            } else if (tmp < th2) {
                imr[icc] = 0;
                imb[icc] = floor(0.5+127.5*(1-cos(M_PI*(th1-tmp)/(th1-th2))));
            } else {
                imr[icc] = 0;
                imb[icc] = 255;
            }
            img[icc] = im[icc];
        }
    }
    for (j = 0; j < ny; j++) {
        if ((jm = j-1) < 0) jm = 0;
        jc = j;
        if ((jp = j+1) > (ny-1)) jp = ny-1;
        jm *= nx;
        jc *= nx;
        jp *= nx;
        for (i = 0; i < nx; i++) {
            if ((imc = i-1) < 0) imc = 0;
            icc = i;
            if ((ipc = i+1) > (nx-1)) ipc = nx-1;
            imm = imc+jm;
            icm = icc+jm;
            ipm = ipc+jm;
            imp = imc+jp;
            icp = icc+jp;
            ipp = ipc+jp;
            imc = imc+jc;
            icc = icc+jc;
            ipc = ipc+jc;
            cmp = 0;
            if (imr[imm] > cmp) cmp = imr[imm];
            if (imr[imc] > cmp) cmp = imr[imc];
            if (imr[imp] > cmp) cmp = imr[imp];
            if (imr[icm] > cmp) cmp = imr[icm];
            if (imr[icp] > cmp) cmp = imr[icp];
            if (imr[ipm] > cmp) cmp = imr[ipm];
            if (imr[ipc] > cmp) cmp = imr[ipc];
            if (imr[ipp] > cmp) cmp = imr[ipp];
            if (imr[icc] > cmp) imr[icc] = cmp;
            if (imb[imm] > cmp) cmp = imb[imm];
            if (imb[imc] > cmp) cmp = imb[imc];
            if (imb[imp] > cmp) cmp = imb[imp];
            if (imb[icm] > cmp) cmp = imb[icm];
            if (imb[icp] > cmp) cmp = imb[icp];
            if (imb[ipm] > cmp) cmp = imb[ipm];
            if (imb[ipc] > cmp) cmp = imb[ipc];
            if (imb[ipp] > cmp) cmp = imb[ipp];
            if (imb[icc] > cmp) imb[icc] = cmp;
        }
    }
    if (fr && (frcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimcolin(CIMAGE *cim0, CIMAGE *cim1, int fr, fdata l0, fdata l1, CIMAGE *cim2, int al)
{
    int i,nx,ny,nb,tmp;
    uchar *im0,*im1,*im2;
    im0 = cim0->im;
    im1 = cim1->im;
    cim2->nx = nx = cim0->nx;
    cim2->ny = ny = cim0->ny;
    cim2->nb = nb = cim0->nb;
    if (al && (alcimage(cim2) == FAILURE)) {
        fprintf(stderr,"ERROR: cimcolin: fail to allocate buffer\n");
        return(FAILURE);
    }
    im2 = cim2->im;
    for (i = 0; i < nx*ny*nb; i++) {
        tmp = 0.5+l0*im0[i]+l1*im1[i];
        if (tmp > 255) tmp = 255;
        if (tmp < 0) tmp = 0;
        im2[i] = tmp;
    }
    if (fr && ((frcimage(cim0) == FAILURE) || (frcimage(cim1) == FAILURE))) {
        fprintf(stderr,"ERROR: cimcolin: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimcrop(CIMAGE *cim0, int fr, CIMAGE *cim1, int al, int x0, int y0, int x1, int y1)
{
    int i,j,k,nx0,ny0,nx1,ny1,tm,nb;
    uchar *im0,*im1;
    im0 = cim0->im;
    nx0 = cim0->nx;
    ny0 = cim0->ny;
    nb = cim0->nb;
    if (x0 > x1) {
        tm = x0;
        x0 = x1;
        x1 = tm;
    }
    if (x0 < 0) x0 = 0;
    if (x1 > (nx0-1)) x1 = nx0-1;
    if (y0 > y1) {
        tm = y0;
        y0 = y1;
        y1 = tm;
    }
    if (y0 < 0) y0 = 0;
    if (y1 > (ny0-1)) y1 = ny0-1;
    cim1->nx = nx1 = (x1-x0)+1;
    cim1->ny = ny1 = (y1-y0)+1;
    cim1->nb = nb;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: cimcrop: fail to allocate buffer\n");
        return(FAILURE);
    }
    im1 = cim1->im;
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny0; j++) {
            for (i = 0; i < nx0; i++) {
                tm = *im0++;
                if ((i >= x0) && (i <= x1) && (j >= y0) && (j <= y1)) *im1++ = tm;
            }
        }
    }
    if (fr && (frcimage(cim0) == FAILURE)) {
        fprintf(stderr,"ERROR: cimcrop: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimrsize(CIMAGE *cim0, int fr, CIMAGE *cim1, int al, int nx1, int ny1)
{
    int i,j,k,nx0,ny0,nb,ip,jp,iq,jq;
    fdata x0,y0,x1,y1,sx,sy,px,qx,py,qy;
    uchar *im0,*im1;
    im0 = cim0->im;
    nx0 = cim0->nx;
    ny0 = cim0->ny;
    nb = cim0->nb;
    cim1->nx = nx1;
    cim1->ny = ny1;
    cim1->nb = nb;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: cimcrop: fail to allocate buffer\n");
        return(FAILURE);
    }
    im1 = cim1->im;
    sx = ((fdata) nx0)/((fdata) nx1);
    sy = ((fdata) ny0)/((fdata) ny1);
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny1; j++) {
            y0 = sy*(y1 = 0.5+j)-0.5;
            jq = (jp = floor(y0))+1;
            py = 1-(qy = y0-jp);
            if (jp < 0) jp = 0;
            if (jq > (ny0-1)) jq = ny0-1;
            for (i = 0; i < nx1; i++) {
                x0 = sx*(x1 = 0.5+i)-0.5;
                iq = (ip = floor(x0))+1;
                px = 1-(qx = x0-ip);
                if (ip < 0) ip = 0;
                if (iq > (nx0-1)) iq = nx0-1;
                im1[i+j*nx1+k*nx1*ny1] = (uchar)(0.5
                +py*px*((fdata)im0[ip+jp*nx0+k*nx0*ny0])+py*qx*((fdata)im0[iq+jp*nx0+k*nx0*ny0])
                +qy*px*((fdata)im0[ip+jq*nx0+k*nx0*ny0])+qy*qx*((fdata)im0[iq+jq*nx0+k*nx0*ny0]));
            }
        }
    }
    if (fr && (frcimage(cim0) == FAILURE)) {
        fprintf(stderr,"ERROR: cimcrop: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimswap(CIMAGE *cim0, int fr, CIMAGE *cim1, int al)
{
    int i,j,k,nx,ny,nb;
    uchar *im0,*im1;
    nx = cim0->nx;
    ny = cim0->ny;
    cim1->ny = nx;
    cim1->nx = ny;
    cim1->nb = nb = cim0->nb;
    im0 = cim0->im;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to allocate buffer\n");
        return(FAILURE);
    }
    im1 = cim1->im;
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                im1[i*ny+j] = im0[j*nx+i];
            }
        }
        im0 += nx*ny;
        im1 += nx*ny;
    }
    im0 -= nx*ny*nb;
    im1 -= nx*ny*nb;
    cim0->im = im0;
    if (fr && (frcimage(cim0) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to free buffer\n");
        return(FAILURE);
    }
    cim1->im = im1;
    return(SUCCESS);
}

int cimxshrink(CIMAGE *cim0, int fr, CIMAGE *cim1, int al)
{
    int i,nx,ny,nb;
    uchar *im0,*im1;
    cim1->nx = nx = (cim0->nx)/2;
    cim1->ny = ny = cim0->ny;
    cim1->nb = nb = cim0->nb;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to allocate buffer\n");
        return(FAILURE);
    }
    im0 = cim0->im;
    im1 = cim1->im;
    for (i = 0; i < nx*ny*nb; i++) *im1++ = ((((int) (*im0++))+((int) (*im0++))))/2;
    if (fr && (frcimage(cim0) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimxexpand(CIMAGE *cim0, int fr, CIMAGE *cim1, int al)
{
    int i,nx,ny,nb;
    uchar *im0,*im1;
    cim1->nx = 2*(nx = cim0->nx);
    cim1->ny = ny = cim0->ny;
    cim1->nb = nb = cim0->nb;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to allocate buffer\n");
        return(FAILURE);
    }
    im0 = cim0->im;
    im1 = cim1->im;
    for (i = 0; i < nx*ny*nb; i++) im1[2*i+1] = im1[2*i] = im0[i];
    if (fr && (frcimage(cim0) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimbshrink(CIMAGE *cim0, int fr, CIMAGE *cim1, int al)
{
    int i,j,nx,ny,nb;
    uchar *im00,*im01,*im1;
    cim1->nx = nx = (cim0->nx)/2;
    cim1->ny = ny = (cim0->ny)/2;
    cim1->nb = nb = cim0->nb;
    if (al && (alcimage(cim1) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to allocate buffer\n");
        return(FAILURE);
    }
    im00 = cim0->im;
    im01 = im00+2*nx;
    im1 = cim1->im;
    for (j = 0; j < ny*nb; j++) {
        for (i = 0; i < nx; i++) {
            im1[i] = (((int) (im00[2*i]))+((int) (im00[2*i+1]))
                     +((int) (im01[2*i]))+((int) (im01[2*i+1])))/4;
        }
        im1  += nx;
        im00 += 4*nx;
        im01 += 4*nx;
    }
    if (fr && (frcimage(cim0) == FAILURE)) {
        fprintf(stderr,"ERROR: grasmcol: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimsp2(CIMAGE *cim, int fr, CIMAGE *cim0, CIMAGE *cim1, int al)
{
    int i,j,k,nx,ny,nb,jm,jp;
    uchar *im,*im0,*im1,*imm,*imc,*imp;
    im = cim->im;
    cim1->nx = cim0->nx = nx = cim->nx;
    cim1->ny = cim0->ny = ny = cim->ny;
    cim1->nb = cim0->nb = nb = cim->nb;
    if (al && ((alcimage(cim0) == FAILURE) || (alcimage(cim1) == FAILURE))) {
        fprintf(stderr,"ERROR: cimspf: fail to allocate buffer\n");
        return(FAILURE);
    }
    im0 = cim0->im;
    im1 = cim1->im;
    for (i = 0; i < nx*ny*nb; i++) im1[i] = im0[i] = im[i];
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny; j++) {
            jm = j-1;
            if (jm < 0) jm += 2;
            jp = j+1;
            if (jp > (ny-1)) jp -= 2;
            for (i = 0; i < nx; i++) {
                if ((j > 22) && (j < 45) && 
                (((i > 138) && (i < 247)) || ((i > 309) && (i < 622)))) {
                    if (j&1) {
                        imm = im1+k*ny*nx+jp*nx;
                        imc = im1+k*ny*nx+j*nx;
                        imp = im1+k*ny*nx+jp*nx;
                    } else {
                        imm = im0+k*ny*nx+jm*nx;
                        imc = im0+k*ny*nx+j*nx;
                        imp = im0+k*ny*nx+jm*nx;
                    }
                } else {
                    if (j&1) {
                        imm = im1+k*ny*nx+jm*nx;
                        imc = im1+k*ny*nx+j*nx;
                        imp = im1+k*ny*nx+jp*nx;
                    } else {
                        imm = im0+k*ny*nx+jm*nx;
                        imc = im0+k*ny*nx+j*nx;
                        imp = im0+k*ny*nx+jp*nx;
                    }
                }
                imc[i] = (((int)(imm[i]))+((int)(imp[i])))/2;
            }
        }
    }
    if (fr && (frcimage(cim) == FAILURE)) {
        fprintf(stderr,"ERROR: cimspf: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimsp1(CIMAGE *cim, int fr, CIMAGE *cim2, int al)
{
    int i,j,k,nx,ny,nb,jm,jp;
    uchar *im,*im0,*im1,*imm,*imc,*imp;
    im = cim->im;
    cim2->nx = nx = cim->nx;
    cim2->ny = ny = cim->ny;
    cim2->nb = 2*(nb = cim->nb);
    if (al && (alcimage(cim2) == FAILURE)) {
        fprintf(stderr,"ERROR: cimspf: fail to allocate buffer\n");
        return(FAILURE);
    }
    im0 = cim2->im;
    im1 = im0+nx*ny*nb;
    for (i = 0; i < nx*ny*nb; i++) im1[i] = im0[i] = im[i];
    for (k = 0; k < nb; k++) {
        for (j = 0; j < ny; j++) {
            jm = j-1;
            if (jm < 0) jm += 2;
            jp = j+1;
            if (jp > (ny-1)) jp -= 2;
            if (j&1) {
                imm = im1+k*ny*nx+jm*nx;
                imc = im1+k*ny*nx+j*nx;
                imp = im1+k*ny*nx+jp*nx;
            } else {
                imm = im0+k*ny*nx+jm*nx;
                imc = im0+k*ny*nx+j*nx;
                imp = im0+k*ny*nx+jp*nx;
            }
            for (i = 0; i < nx; i++) {
                imc[i] = (((int)(imm[i]))+((int)(imp[i])))/2;
            }
        }
    }
    if (fr && (frcimage(cim) == FAILURE)) {
        fprintf(stderr,"ERROR: cimspf: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                         FDATA to UCHAR and UCHAR to FDATA                          */
/*------------------------------------------------------------------------------------*/

int fimtocim(FIMAGE *fim, int fr, CIMAGE *cim, int al)
{
    int i,nx,ny,nb,tmp;
    uchar *im;
    fdata *rim;
    cim->nx = nx = fim->nx;
    cim->ny = ny = fim->ny;
    cim->nb = nb = fim->nb;
    if (al && (alcimage(cim) == FAILURE)) {
        fprintf(stderr,"ERROR: fimtocim: fail to allocate buffer\n");
        return(FAILURE);
    }
    im = cim->im;
    rim = fim->im;
    for (i = 0; i < nx*ny*nb; i++) {
        tmp = (int) ((*rim++)+0.5);
        if (tmp > 255) tmp = 255;
        if (tmp < 0) tmp = 0;
        *im++ = tmp;
    }
    if (fr && (frfimage(fim) == FAILURE)) {
        fprintf(stderr,"ERROR: fimtocim: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

int cimtofim(CIMAGE *cim, int fr, FIMAGE *fim, int al)
{
    int i,nx,ny,nb;
    uchar *im;
    fdata *rim;
    fim->nx = nx = cim->nx;
    fim->ny = ny = cim->ny;
    fim->nb = nb = cim->nb;
    if (al && (alfimage(fim) == FAILURE)) {
        fprintf(stderr,"ERROR: cimtofim: fail to allocate buffer\n");
        return(FAILURE);
    }
    im = cim->im;
    rim = fim->im;
    for (i = 0; i < nx*ny*nb; i++) *rim++ = *im++;
    if (fr && (frcimage(cim) == FAILURE)) {
        fprintf(stderr,"ERROR: cimtofim: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                              VELOCITY IMAGE FILE READ                              */
/*------------------------------------------------------------------------------------*/

int alvimage(VIMAGE *vim)
{
    int nx,ny,nb;
    fdata ***d;
    nx = vim->nx;
    ny = vim->ny;
    nb = vim->nb;
    if (malloc3((void ****) &d,nb,nx,ny,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: alvimage: fail to allocate buffers\n");
        return(FAILURE);
    }
    vim->d = d;
    if (nb > 0) vim->dx = d[0];
    if (nb > 1) vim->dy = d[1];
    if (nb > 2) vim->dz = d[2];
    return(SUCCESS);
}

int rdvimage(char *fname, VIMAGE *vim, int fo, int al)
{
    /*--------------------------------*/
    /* read velocity data from a file */
    /*--------------------------------*/
    FILE *fp;
    int i,j,nx,ny,nb,nx0,ny0,nx1,ny1,n;
    char s[256];
    fdata **dx,**dy,**dz;
    float head[6],dd[2],cm,t,x0,x,y,z,vx,vy,vz,xs,ys;
    if (fo == UWO) {
        if ((fp = fopen(fname,"rb")) == NULL) {
            fprintf(stderr,"ERROR: rdvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        if (fread(head,sizeof(float),6,fp) != 6) {
            fprintf(stderr,"ERROR: rdvimage: can't read header: %s\n",fname);
            return(FAILURE);
        }
        bswap((void *) head, 6);
        vim->nb = nb = 2;
        vim->nx = nx = head[0];
        vim->ny = ny = head[1];
        vim->nx0 = nx0 = head[4];
        vim->ny0 = ny0 = head[5];
        vim->nx1 = nx1 = nx-(head[2]+nx0);
        vim->ny1 = ny1 = ny-(head[3]+ny0);
        if (al && (alvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: rdvimage: fail to allocate buffers\n");
            return(FAILURE);
        }
        dx = vim->dx;
        dy = vim->dy;
        for (j = 0; j < ny0; j++) {
            for (i = 0; i < nx; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
        for (j = ny0; j < ny-ny1; j++) {
            for (i = 0; i < nx0; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
            for (i = nx0; i < nx-nx1; i++) {
                if (fread(dd,sizeof(float),2,fp) != 2) {
                    fprintf(stderr,"ERROR: rdvimage: can't read data: %s\n",fname);
                    return(FAILURE);
                }
                bswap((void *) dd, 2);
                dx[i][j] = dd[0];
                dy[i][j] = -dd[1];
            }
            for (i = nx-nx1; i < nx; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
        for (j = ny-ny1; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
        fclose(fp);
        return(SUCCESS);
    } else if (fo == TEC) {
        nb = vim->nb = 2;
        nx = vim->nx;
        ny = vim->ny;
        vim->nx0 = vim->ny0 = vim->nx1 = vim->ny1 = 0;
        if (al && (alvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: rdvimage: fail to allocate buffers\n");
            return(FAILURE);
        }
        dx = vim->dx;
        dy = vim->dy;
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
        if ((fp = fopen(fname,"rb")) == NULL) {
            fprintf(stderr,"ERROR: rdvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        for (i = 0; i < 3; i++) fgets(s,256,fp);
        while (fscanf(fp,"%d %d %f %f %f",&i,&j,&vx,&vy,&cm) == 5) {
            if ((cm > 1.1) && (i >= 0) && (i < nx) && (j >= 0) && (j < ny)) {
                dx[i][ny-1-j] = vx;
                dy[i][ny-1-j] = -vy;
            }
        }
        fclose(fp);
        return(SUCCESS);
    } else if (fo == DAT) {
        nb = vim->nb = 2;
        nx = vim->nx;
        ny = vim->ny;
        vim->nx0 = vim->ny0 = vim->nx1 = vim->ny1 = 0;
        if (al && (alvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: rdvimage: fail to allocate buffers\n");
            return(FAILURE);
        }
        dx = vim->dx;
        dy = vim->dy;
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
        if ((fp = fopen(fname,"rb")) == NULL) {
            fprintf(stderr,"ERROR: rdvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        for (i = 0; i < 13; i++) fgets(s,256,fp);
        while (fscanf(fp,"%d %d %f %f %f",&i,&j,&vx,&vy,&vz) == 5) {
            if ((i >= 0) && (i < nx) && (j >= 0) && (j < ny)) {
                dx[i][ny-1-j] = vx;
                dy[i][ny-1-j] = -vy;
            }
        }
        fclose(fp);
        return(SUCCESS);
    } else if (fo == DA2) {
        nb = vim->nb = 3;
        nx = vim->nx;
        ny = vim->ny;
        vim->nx0 = vim->ny0 = vim->nx1 = vim->ny1 = 0;
        if (al && (alvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: rdvimage: fail to allocate buffers\n");
            return(FAILURE);
        }
        dx = vim->dx;
        dy = vim->dy;
        dz = vim->dz;
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                dx[i][j] = 100.0;
                dy[i][j] = -100.0;
            }
        }
        if ((fp = fopen(fname,"rb")) == NULL) {
            fprintf(stderr,"ERROR: rdvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        while (fscanf(fp,"%f %f %f %f %f %f",&x,&y,&t,&vx,&vy,&vz) == 6) {
            i = x;
            j = y;
            if ((i >= 0) && (i < nx) && (j >= 0) && (j < ny)) {
                dx[i][j] = vx;
                dy[i][j] = vy;
                dz[i][j] = vz;
            }
        }
        fclose(fp);
        return(SUCCESS);
    } else if (fo == DA3) {
        if ((fp = fopen(fname,"rb")) == NULL) {
            fprintf(stderr,"ERROR: rdvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        if (fscanf(fp,"%f %f %f %f %f %f",&x,&y,&z,&vx,&vy,&vz) != 6) return(FAILURE);
        vim->x0 = x;
        vim->y0 = y;
        x0 = x;
        ny = n = 1;
        xs = ys = 0.0;
        while (fscanf(fp,"%f %f %f %f %f %f",&x,&y,&z,&vx,&vy,&vz) == 6) {
            ny += (x == x0);
            xs = x;
            ys = y;
            n++;
        }
        nx = n/ny;
        vim->nx = nx = n/ny;
        vim->ny = ny;
        vim->nb = nb = 3;
        vim->nx0 = vim->ny0 = vim->nx1 = vim->ny1 = 0;
        vim->xs = (xs-vim->x0)/(nx-1);
        vim->ys = (ys-vim->y0)/(ny-1);
        if (al && (alvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: rdvimage: fail to allocate buffers\n");
            return(FAILURE);
        }
        dx = vim->dx;
        dy = vim->dy;
        dz = vim->dz;
        rewind(fp);
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                fscanf(fp,"%f %f %f %f %f %f",&x,&y,&z,&vx,&vy,&vz);
                dx[i][j] = vx;
                dy[i][j] = vy;
                dz[i][j] = vz;
            }
        }
        fclose(fp);
        return(SUCCESS);
    }
    return(FAILURE);
}

/*------------------------------------------------------------------------------------*/
/*                             VELOCITY IMAGE FILE WRITE                              */
/*------------------------------------------------------------------------------------*/

int frvimage(VIMAGE *vim)
{
    int nx,ny,nb;
    fdata ***d;
    nx = vim->nx;
    ny = vim->ny;
    nb = vim->nb;
    d  = vim->d;
    if ((d != NULL) && (free3((void ****) &d,nb,nx,ny,sizeof(fdata)) == FAILURE)) {
        fprintf(stderr,"ERROR: frvimage: fail to free buffers\n");
        return(FAILURE);
    }
    vim->d = NULL;
    return(SUCCESS);
}

int wrvimage(char *fname, VIMAGE *vim, int fo, int fr)
{
    /*-------------------------------*/
    /* write velocity data in a file */
    /*-------------------------------*/
    FILE *fp;
    int i,j,nx,ny,nx0,ny0,nx1,ny1;
    fdata **dx,**dy,**dz;
    float head[6],dd[2],x0,y0,xs,ys;
    nx = vim->nx;
    ny = vim->ny;
    dx = vim->dx;
    dy = vim->dy;
    dz = vim->dz;
    nx0 = vim->nx0;
    ny0 = vim->ny0;
    nx1 = vim->nx1;
    ny1 = vim->ny1;
    if (fo == UWO) {
        if ((fp = fopen(fname,"wb")) == NULL) {
            fprintf(stderr,"ERROR: wrvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        head[0] = nx;
        head[1] = ny;
        head[2] = nx-(nx0+nx1);
        head[3] = ny-(ny0+ny1);
        head[4] = nx0;
        head[5] = ny0;
        bswap((void *) head, 6);
        if (fwrite(head,sizeof(float),6,fp) != 6) {
            fprintf(stderr,"ERROR: wrvimage: can't write header: %s\n",fname);
            return(FAILURE);
        }
        for (j = ny0; j < ny-ny1; j++) {
            for (i = nx0; i < nx-nx1; i++) {
                dd[0] = dx[i][j];
                dd[1] = -dy[i][j];
                bswap((void *) dd,2);
                if (fwrite(dd,sizeof(float),2,fp) != 2) {
                    fprintf(stderr,"ERROR: wrvimage: can't write data: %s\n",fname);
                    return(FAILURE);
                }
            }
        }
        if (fr && (frvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: wrvimage: fail to free buffers\n");
            return(FAILURE);
        }
        fclose(fp);
        return(SUCCESS);
    } else if (fo == DA3) {
        if ((fp = fopen(fname,"wb")) == NULL) {
            fprintf(stderr,"ERROR: wrvimage: can't open file: %s\n",fname);
            return(FAILURE);
        }
        x0 = vim->x0;
        y0 = vim->y0;
        xs = vim->xs;
        ys = vim->ys;
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                fprintf(fp,"%.2f  %.2f   0.00  %.5f  %.5f  %.5f\n",
                x0+i*xs,y0+j*ys,dx[i][j],dy[i][j],dz[i][j]);
            }
        }
        if (fr && (frvimage(vim) == FAILURE)) {
            fprintf(stderr,"ERROR: wrvimage: fail to free buffers\n");
            return(FAILURE);
        }
        fclose(fp);
        return(SUCCESS);
    }
    return(FAILURE);
}

int vimxexpand(VIMAGE *vim0, int fr, VIMAGE *vim1, int al)
{
    int i,j,nx,ny,nx0,ny0,nx1,ny1;
    fdata **dx0,**dy0,**dx1,**dy1,***d1;
    VIMAGE vim[1];
    vim->nx  = 2*(nx = vim0->nx);
    vim->nx0 = 2*(nx0 = vim0->nx0);
    vim->nx1 = 2*(nx1 = vim0->nx1);
    vim->ny  = ny = vim0->ny;
    vim->ny0 = ny0 = vim0->ny0;
    vim->ny1 = ny1 = vim0->ny1;
    vim->nb  = 2;
    if (al && (alvimage(vim) == FAILURE)) {
        fprintf(stderr,"ERROR: vimxexpand: fail to allocate buffer\n");
        return(FAILURE);
    }
    dx0 = vim0->dx;
    dy0 = vim0->dy;
    d1  = al ? vim->d  : vim1->d;
    dx1 = al ? vim->dx : vim1->dx;
    dy1 = al ? vim->dy : vim1->dy;
    for (i = nx0; i < nx-nx1; i++) {
        for (j = ny0; j < ny-ny1; j++) {
            dx1[2*i+1][j] = dx1[2*i][j] = 2.0*dx0[i][j];
            dy1[2*i+1][j] = dy1[2*i][j] = dy0[i][j];
        }
    }
    if (fr && (frvimage(vim0) == FAILURE)) {
        fprintf(stderr,"ERROR: vimxexpand: fail to free buffer\n");
        return(FAILURE);
    }
    vim1->d   = d1;
    vim1->dx  = dx1;
    vim1->dy  = dy1;
    vim1->nx  = vim->nx;
    vim1->nx0 = vim->nx0;
    vim1->nx1 = vim->nx1;
    vim1->ny  = vim->ny;
    vim1->ny0 = vim->ny0;
    vim1->ny1 = vim->ny1;    
    vim1->nb  = 2;
    return(SUCCESS);
}

int vimcrop(VIMAGE *vim0, int fr, VIMAGE *vim1, int al, int x0, int y0, int x1, int y1)
{
    int i,j,nx0,ny0,nx1,ny1,tm;
    fdata **dx0,**dy0,**dx1,**dy1,*dx0i,*dy0i,*dx1i,*dy1i;
    dx0 = vim0->dx;
    dy0 = vim0->dy;
    nx0 = vim0->nx;
    ny0 = vim0->ny;
    if (x0 > x1) {
        tm = x0;
        x0 = x1;
        x1 = tm;
    }
    if (x0 < 0) x0 = 0;
    if (x1 > (nx0-1)) x1 = nx0-1;
    if (y0 > y1) {
        tm = y0;
        y0 = y1;
        y1 = tm;
    }
    if (y0 < 0) y0 = 0;
    if (y1 > (ny0-1)) y1 = ny0-1;
    vim1->nb  = 2;
    vim1->nx = nx1 = (x1-x0)+1;
    vim1->ny = ny1 = (y1-y0)+1;
    vim1->nx0 = 0;
    vim1->ny0 = 0;
    vim1->nx1 = 0;
    vim1->ny1 = 0;
    if (al && (alvimage(vim1) == FAILURE)) {
        fprintf(stderr,"ERROR: vimcrop: fail to allocate buffer\n");
        return(FAILURE);
    }
    dx1 = vim1->dx;
    dy1 = vim1->dy;
    for (i = 0; i < nx0; i++) {
        dx0i = *dx0++;
        dy0i = *dy0++;
        if ((i >= x0) && (i <= x1)) {
            dx1i = *dx1++;
            dy1i = *dy1++;
            for (j = 0; j < ny0; j++) {
                if ((j >= y0) && (j <= y1)) {
                    *dx1i++ = *dx0i++;
                    *dy1i++ = *dy0i++;
                } else {
                    dx0i++;
                    dy0i++;
                }
            }
        }
    }
    if (fr && (frvimage(vim0) == FAILURE)) {
        fprintf(stderr,"ERROR: vimcrop: fail to free buffer\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                               FLOAT IMAGE FILE READ                                */
/*------------------------------------------------------------------------------------*/

int alfimage(FIMAGE *fim)
{
    int nx,ny,nb;
    fdata *im;
    nx = fim->nx;
    ny = fim->ny;
    nb = fim->nb;
    if (malloc1((void **) &im,nx*ny*nb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: alfimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    fim->im = im;
    return(SUCCESS);
}

int rdfimage(char *name, FIMAGE *fim, int fo, int al)
{
    FILE *fp;
    int nx,ny,nb;
    fdata *im;
    if ((fp = fopen(name,"rb")) == NULL) {
        fprintf(stderr,"ERROR: rdfimage: can't open file: %s\n",name);
        return(FAILURE);
    }
    if ((fo == FCF) && (rdcalsize(fp,&nx,&ny,&nb) == FAILURE)) {
        fprintf(stderr,"ERROR: rdfimage: can't read image size: %s\n",name);
        return(FAILURE);
    }
    nb /= sizeof(fdata);
    if (fo == RAW) {
        nx = fim->nx;
        ny = fim->ny;
        nb = fim->nb;
    } else {
        fim->nx = nx;
        fim->ny = ny;
        fim->nb = nb;
    }
    if (al && (alfimage(fim) == FAILURE)) {
        fprintf(stderr,"ERROR: rdfimage: fail to allocate buffer\n");
        return(FAILURE);
    }
    im = fim->im;
    if (fseek(fp,-nx*ny*nb*sizeof(fdata),2) != 0) {
        fprintf(stderr,"ERROR: rdfimage: can't fseek to data: %s\n",name);
        return(FAILURE);
    }
    if (fread(im,sizeof(fdata),nx*ny*nb,fp) != nx*ny*nb) {
        fprintf(stderr,"ERROR: rdfimage: fail to read data: %s\n",name);
        return(FAILURE);
    }
    fclose(fp);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*                               FLOAT IMAGE FILE WRITE                               */
/*------------------------------------------------------------------------------------*/

int frfimage(FIMAGE *fim)
{
    int nx,ny,nb;
    fdata *im;
    nx = fim->nx;
    ny = fim->ny;
    nb = fim->nb;
    im = fim->im;
    if (free1((void **) &im,nx*ny*nb,sizeof(fdata)) == FAILURE) {
        fprintf(stderr,"ERROR: frcimage: fail to free buffer\n");
        return(FAILURE);
    }
    fim->im = NULL;
    return(SUCCESS);
}

int wrfimage(char *name, FIMAGE *fim, int fo, int fr)
{
    FILE *fp;
    int nx,ny,nb;
    fdata *im;
    nx = fim->nx;
    ny = fim->ny;
    nb = fim->nb;
    im = fim->im;
    if ((fp = fopen(name,"wb")) == NULL) {
	fprintf(stderr,"Error wrfimage: can't open file: %s\n",name);
	return(FAILURE);
    }
    if ((fo == FCF) && (wrcalsize(fp,nx,ny,nb) == FAILURE)) {
        fprintf(stderr,"ERROR: wrfimage: can't write image size: %s\n",name);
        return(FAILURE);
    }
    if (fwrite(im,sizeof(fdata),nx*ny*nb,fp) != nx*ny*nb) {
        fprintf(stderr,"ERROR: wrfimage: fail to write data: %s\n",name);
        return(FAILURE);
    }
    fclose(fp);
    if (fr && (frfimage(fim) == FAILURE)) {
	fprintf(stderr,"Error wrfimage: fail to free buffer: %s\n",name);
	return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/

int rdcalcam(char *fname, CALCAM *cal)
{
    FILE *fp;
    char str[256];
    if ((fp=fopen(fname,"rb")) == NULL) {
	fprintf(stderr,"Error rdcalcam: can't open file: %s\n",fname);
	return(FAILURE);
    }
    fscanf(fp,"%s %s %s %s %s %s %s",str,str,str,str,str,str,str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->u0);
    fscanf(fp,"%s %s %s",str,str,str);
    sscanf(str,"%lf",&cal->v0);
    fscanf(fp,"%s %s %s %s %s",str,str,str,str,str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->au);
    fscanf(fp,"%s %s %s",str,str,str);
    sscanf(str,"%lf",&cal->av);
    fscanf(fp,"%s %s %s %s %s",str,str,str,str,str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%d",&cal->nx);
    fscanf(fp,"%s %s %s",str,str,str);
    sscanf(str,"%d",&cal->ny);
    fscanf(fp,"%s %s %s %s %s %s",str,str,str,str,str,str);
    str[strlen(str)-1] = 0;
    sscanf(str+2,"%lf",&cal->rxx);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->rxy);
    fscanf(fp,"%s",str);
    str[strlen(str)-2] = 0;
    sscanf(str,"%lf",&cal->rxz);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str+1,"%lf",&cal->ryx);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->ryy);
    fscanf(fp,"%s",str);
    str[strlen(str)-2] = 0;
    sscanf(str,"%lf",&cal->ryz);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str+1,"%lf",&cal->rzx);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->rzy);
    fscanf(fp,"%s",str);
    str[strlen(str)-2] = 0;
    sscanf(str,"%lf",&cal->rzz);
    fscanf(fp,"%s %s",str,str);
    str[strlen(str)-1] = 0;
    sscanf(str+1,"%lf",&cal->tx);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->ty);
    fscanf(fp,"%s",str);
    str[strlen(str)-1] = 0;
    sscanf(str,"%lf",&cal->tz);
    fclose(fp);
    return(SUCCESS);
}

int wrcalcam(char *fname, CALCAM *cal)
{
    FILE *fp;
    if ((fp=fopen(fname,"wb")) == NULL) {
	fprintf(stderr,"Error rdcalcam: can't open file: %s\n",fname);
	return(FAILURE);
    }
    fprintf(fp,"Intrinsic parameters:\n");
    fprintf(fp,"Image center: u0 = %f, v0 = %f\n",(float) cal->u0,(float) cal->v0);
    fprintf(fp,"Scale factor: \n");
    fprintf(fp,"au = %f, av = %f\n",(float) cal->au,(float) cal->av);
    fprintf(fp,"Image size: dimx = %d, dimy = %d\n",cal->nx,cal->ny);
    fprintf(fp,"Focale: 1.0\n");
    fprintf(fp,"Extrinsic parameters:\n");
    fprintf(fp,"Rotation:\n");
    fprintf(fp,"{{%f, %f, %f},\n",(float) cal->rxx,(float) cal->rxy,(float)  cal->rxz);
    fprintf(fp," {%f, %f, %f},\n",(float) cal->ryx,(float) cal->ryy,(float)  cal->ryz);
    fprintf(fp," {%f, %f, %f}}\n",(float) cal->rzx,(float) cal->rzy,(float)  cal->rzz);
    fprintf(fp,"Translation:\n");
    fprintf(fp,"{%f, %f, %f}\n",(float) cal->tx,(float) cal->ty,(float)  cal->tz);
    fclose(fp);
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
/*---------------------------- LINEAR EQUATION ROUTINES ------------------------------*/
/*------------------------------------------------------------------------------------*/

fdata det3(fdata a00, fdata a10, fdata a20,
           fdata a01, fdata a11, fdata a21,
           fdata a02, fdata a12, fdata a22)
{
    return((a00*a11*a22+a01*a12*a20+a02*a10*a21)-(a00*a12*a21+a01*a10*a22+a02*a11*a20));
}

void eqlin3(fdata a00, fdata a10, fdata a20, fdata b0, fdata *px0,
            fdata a01, fdata a11, fdata a21, fdata b1, fdata *px1,
            fdata a02, fdata a12, fdata a22, fdata b2, fdata *px2)
{
    fdata d;
    d = det3(a00,a10,a20,a01,a11,a21,a02,a12,a22);
    *px0 = det3(b0,a10,a20,b1,a11,a21,b2,a12,a22)/d;
    *px1 = det3(a00,b0,a20,a01,b1,a21,a02,b2,a22)/d;
    *px2 = det3(a00,a10,b0,a01,a11,b1,a02,a12,b2)/d;
}

/*------------------------------------------------------------------------------------*/

double ddet3(double a00, double a10, double a20,
             double a01, double a11, double a21,
             double a02, double a12, double a22)
{
    return((a00*a11*a22+a01*a12*a20+a02*a10*a21)-(a00*a12*a21+a01*a10*a22+a02*a11*a20));
}

void deqlin3(double a00, double a10, double a20, double b0, double *px0,
             double a01, double a11, double a21, double b1, double *px1,
             double a02, double a12, double a22, double b2, double *px2)
{
    double d;
    d = ddet3(a00,a10,a20,a01,a11,a21,a02,a12,a22);
    *px0 = ddet3(b0,a10,a20,b1,a11,a21,b2,a12,a22)/d;
    *px1 = ddet3(a00,b0,a20,a01,b1,a21,a02,b2,a22)/d;
    *px2 = ddet3(a00,a10,b0,a01,a11,b1,a02,a12,b2)/d;
}

/*------------------------------------------------------------------------------------*/
