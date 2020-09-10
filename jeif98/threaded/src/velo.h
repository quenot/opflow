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

#define EPS_BAR 0.01

/*------------------------------------------------------------------------------------*/

void fillquad(int, int, fdata **, fdata **, fdata **, fdata **,
fdata **, fdata **, int, int, int **, fdata, int);
int vshift(VIMAGE *, int, VIMAGE *, int, fdata, fdata, fdata);
int invdep(VIMAGE *, int, VIMAGE *, int);
void vscale(VIMAGE *, fdata);
void vcolin(VIMAGE *, VIMAGE *, VIMAGE *, fdata, fdata);
int vsplit(VIMAGE *, VIMAGE *, VIMAGE *);
int vsplit2(VIMAGE *);
int vsplit3(VIMAGE *);
int ishift(CIMAGE *, CIMAGE *, VIMAGE *, int, int, fdata, fdata, fdata);
float rim(float, float, float **);
float rimc(float, float, float **, int, int);
int compdep(VIMAGE *, int, VIMAGE *, int, VIMAGE *, int);
void lininfo(fdata, fdata, fdata, int, fdata, fdata, fdata, int, fdata *, fdata *);
void autonorm(VIMAGE *vim);

/*------------------------------------------------------------------------------------*/

