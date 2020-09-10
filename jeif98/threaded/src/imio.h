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

#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#ifdef slrs
#include <netinet/in.h>
#endif
#ifdef linux
#include <netinet/in.h>
#endif

/*------------------------------------------------------------------------------------*/

#define RAW     0                 /* raw file format, no filename extension           */
#define PGM     1                 /* pgm file format, .pgm filename extension         */
#define PPM     2                 /* ppm file format, .ppm filename extension         */
#define CCF     3                 /* calf image file format, no filename extension    */
#define RAS     4                 /* sun rasterfile format, no filename extension     */
#define UWO     5                 /* uwo velocity file format, no filename extension  */
#define FCF     6                 /* calf velocity file format, no filename extension */
#define LUM     7                 /* HHI AP video file format, no filename extension  */
#define TEC     8                 /* tec velocity file format, no filename extension  */
#define DAT     9                 /* dat velocity file format, vsj/std001-8 series    */
#define DA2     10                /* dat velocity file format, vsj/std301-2 series    */
#define DA3     11                /* dat velocity file format, vsj/std331-7 series    */
#define RAS_MAGIC 0x59a66a95

/*------------------------------------------------------------------------------------*/

typedef struct {
    /*------------------------------------------------*/
    /* For multiband raw image data (uchar)           */
    /*------------------------------------------------*/
    int nx;             /* image width                */
    int ny;             /* image height               */
    int nb;             /* image number of bands      */
    uchar *im;          /* multiband image data       */
    /*------------------------------------------------*/
} CIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* For multiband raw image data (int)             */
    /*------------------------------------------------*/
    int nx;             /* image width                */
    int ny;             /* image height               */
    int nb;             /* image number of bands      */
    int *im;            /* multiband image data       */
    /*------------------------------------------------*/
} IIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* For multiband raw image sequences (uchar)      */
    /*------------------------------------------------*/
    int nx;             /* image width                */
    int ny;             /* image height               */
    int nb;             /* image number of bands      */
    int nf;             /* first image number         */
    int nl;             /* last image number          */
    int *vf;            /* frame validation           */
    uchar **im;         /* multiband image sequence   */
    fdata *val;         /* sequence valuation         */
    /*------------------------------------------------*/
} SIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* For multiband float image data (fdata)         */
    /*------------------------------------------------*/
    int nx;             /* image width                */
    int ny;             /* image height               */
    int nb;             /* image number of bands      */
    fdata *im;          /* multiband fdata image data */
    /*------------------------------------------------*/
} FIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* For multiband complex image data (fdata)       */
    /*------------------------------------------------*/
    int nx;             /* image width                */
    int ny;             /* image height               */
    int nb;             /* image number of bands      */
    fdata *rim;         /* multiband fdata image data */
    fdata *iim;         /* multiband fdata image data */
    /*------------------------------------------------*/
} MIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* For transformed image data with a boolean flag */
    /* specifying if the tramsformation is valid      */
    /* (fdata) for the multiband image data           */
    /* (uchar) for the boolean flag                   */
    /*------------------------------------------------*/
    int nx;             /* image width                */
    int ny;             /* image height               */
    int nb;             /* image number of bands      */
    fdata ***xm;        /* multiband fdata image data */
    uchar **in;         /* image validity flag        */
    /*------------------------------------------------*/
} DIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* For velocity data (fdata)                      */
    /*------------------------------------------------*/
    int nx;          /* image width                   */
    int ny;          /* image height                  */
    int nb;          /* number of components          */
    int nx0;         /* image width left shrink       */
    int ny0;         /* image height top shrink       */
    int nx1;         /* image width right shrink      */
    int ny1;         /* image height bottom shrink    */
    fdata ***d;      /* all velocity components       */
    fdata **dx;      /* horizontal velocity component */
    fdata **dy;      /* vertical velocity component   */
    fdata **dz;      /* depth velocity component      */
    fdata x0;        /* x coordinate of point (0,0)   */
    fdata y0;        /* y coordinate of point (0,0)   */
    fdata xs;        /* x step                        */
    fdata ys;        /* y step                        */
    /*------------------------------------------------*/
} VIMAGE;

typedef struct {
    /*------------------------------------------------*/
    /* Camera calibration parameters                  */
    /*------------------------------------------------*/
    int nx;          /* image width                   */
    int ny;          /* image height                  */
    double u0;       /* image center horizontal       */
    double v0;       /* image center vertical         */
    double au;       /* scale factor horizontal       */
    double av;       /* scale factor vertical         */
    double rxx;      /* rotation parameter            */
    double rxy;      /* rotation parameter            */
    double rxz;      /* rotation parameter            */
    double tx;       /* translation parameter         */
    double ryx;      /* rotation parameter            */
    double ryy;      /* rotation parameter            */
    double ryz;      /* rotation parameter            */
    double ty;       /* translation parameter         */
    double rzx;      /* rotation parameter            */
    double rzy;      /* rotation parameter            */
    double rzz;      /* rotation parameter            */
    double tz;       /* translation parameter         */
    /*------------------------------------------------*/
} CALCAM;

/*------------------------------------------------------------------------------------*/

int filelen(char *, int *);
void bswap(void *, size_t);
void prflp(fdata, char *);
int testbarg(int *, char *[], char *);
int getbarg(int *, char *[], char *);
int getcarg(int *, char *[], char *, char, char *);
int getiarg(int *, char *[], char *, int, int *);
int getfarg(int *, char *[], char *, fdata, fdata *);
int getsarg(int *, char *[], char *, char *);
int getimsize(int *, char *[], char *, int *, int *, int *, int *);
int rdlumsize(FILE *, int *, int *, int *);
int alcimage(CIMAGE *);
int frcimage(CIMAGE *);
int rdcimage(char *, CIMAGE *, int, int);
int wrcimage(char *, CIMAGE *, int, int);
int wrcimagec(char *, char *, CIMAGE *, int, int);
int rdlimage(char *, CIMAGE *, int, int, int);
int rdsimage(char *, SIMAGE *, int);
int frsimage(SIMAGE *);
int aliimage(IIMAGE *);
int friimage(IIMAGE *);
int coltogra(CIMAGE *, int, CIMAGE *, int);
int gratocol(CIMAGE *, int, CIMAGE *, int);
int grasmcol(CIMAGE *, int, CIMAGE *, int);
int graoocol(CIMAGE *, int, CIMAGE *, int, fdata, fdata);
int alvimage(VIMAGE *);
int frvimage(VIMAGE *);
int rdvimage(char *, VIMAGE *, int, int);
int wrvimage(char *, VIMAGE *, int, int);
int alfimage(FIMAGE *);
int frfimage(FIMAGE *);
int rdfimage(char *, FIMAGE *, int, int);
int wrfimage(char *, FIMAGE *, int, int);
int rdcalcam(char *, CALCAM *);
int wrcalcam(char *, CALCAM *);

int cimswap(CIMAGE *, int, CIMAGE *, int);
int cimxshrink(CIMAGE *, int, CIMAGE *, int);
int cimbshrink(CIMAGE *, int, CIMAGE *, int);
int cimxexpand(CIMAGE *, int, CIMAGE *, int);
int simxshrink(SIMAGE *, uchar ***);
int simxexpand(SIMAGE *, uchar **);
int vimxexpand(VIMAGE *, int, VIMAGE *, int);
int smocim(CIMAGE *, CIMAGE *);
void hpscim(CIMAGE *, CIMAGE *, fdata);
int cimtofim(CIMAGE *, int, FIMAGE *, int);
int fimtocim(FIMAGE *, int, CIMAGE *, int);
int cimcolin(CIMAGE *, CIMAGE *, int, fdata, fdata, CIMAGE *, int);
int cimrsize(CIMAGE *, int, CIMAGE *, int, int, int);
int cimcrop(CIMAGE *, int, CIMAGE *, int, int, int, int, int);
int vimcrop(VIMAGE *, int, VIMAGE *, int, int, int, int, int);
int cimsp2(CIMAGE *, int, CIMAGE *, CIMAGE *, int);
int cimsp1(CIMAGE *, int, CIMAGE *, int);

fdata det3(fdata, fdata, fdata,
           fdata, fdata, fdata,
           fdata, fdata, fdata);
void eqlin3(fdata, fdata, fdata, fdata, fdata *,
            fdata, fdata, fdata, fdata, fdata *,
            fdata, fdata, fdata, fdata, fdata *);
double ddet3(double, double, double,
             double, double, double,
             double, double, double);
void deqlin3(double, double, double, double, double *,
             double, double, double, double, double *,
             double, double, double, double, double *);

/*------------------------------------------------------------------------------------*/
