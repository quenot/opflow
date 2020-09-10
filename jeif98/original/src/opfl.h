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

#define INFINI  1000000000
#define DXL     0.01
#define DXS     0.01
#define DYL     0.01
#define DYS     0.01
#define DLF     1.0
#define ID      'h'
#define NP      1
#define TAU     0.35
#define SPW     1.0
#define WIN     0.125
#define SMO     2.6
#define SMF     0
#define ENW     0
#define EIT     10
#define SC      3
#define MD      0.03
#define MF      2
#define MB      1

/*------------------------------------------------------------------------------------*/

typedef struct {
    /*----------------------------------------------------*/
    /* Orthogonal Dynamic Programming control parameters  */
    /*----------------------------------------------------*/
    fdata m;        /* computed frame number              */
    fdata dlin[2];  /* dispersion constraint linear part  */
    fdata dsig[2];  /* dispersion constraint sigmoid part */
    fdata eps;      /* streching path penality            */
    fdata dlf;      /* global distance specification      */
    fdata tau;      /* decreasing control parameter       */
    fdata spwfac;   /* strip spacing and width factor     */
    fdata winfac;   /* search window width factor         */
    fdata smofac;   /* smoothing factor                   */
    fdata md;       /* maximum searched displacement      */
    int vb;         /* verbose flag                       */
    int hx;         /* horizontal shrink of images        */
    int mb;         /* multiply bands with smoothing      */
    int goa;        /* gain and offset correction         */
    int dni;        /* do not interpolate distance matrix */
    int bsi;        /* use bicubic spline interpolation   */
    int id;         /* initial slicing direction          */
    int np;         /* number of passes                   */
    int enw;        /* strip extra width                  */
    int eit;        /* number of extra iterations         */
    int smf;        /* smootihng function option          */
    int mf;         /* multi-frame computation            */
    int sum;        /* use sum instead of max in distance */
    int sc;         /* slope constraint                   */
    int vs;         /* variable scale strategy            */
    /*----------------------------------------------------*/
} DPCONTROL;

/*------------------------------------------------------------------------------------*/

int imin(int, int);
int imax(int, int);
void dpusage(FILE *, int, fdata);
int getcontrol(int *, char **, DPCONTROL *, fdata, int, int);
void dpinfo(DPCONTROL, int);
int opflow(SIMAGE *, VIMAGE *, DPCONTROL *);
int opflow3(SIMAGE *, VIMAGE *, DPCONTROL *);

/*------------------------------------------------------------------------------------*/
