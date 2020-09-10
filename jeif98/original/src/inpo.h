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

int inter1(fdata, fdata, int, int, int, uchar *, fdata *);
int inter2(fdata, fdata, int, int, int, uchar *, fdata *);
int apply(VIMAGE *, SIMAGE *, int, DIMAGE *, fdata, int, int, int);
void goadj(DIMAGE *, DIMAGE *, fdata);
int getimage(int, SIMAGE *, VIMAGE *, fdata, DIMAGE *, DIMAGE *, int, int, int);
int rerrorf(SIMAGE *, VIMAGE *, fdata, fdata *, fdata *);
int iminterp(SIMAGE *, CIMAGE *, VIMAGE *, int, fdata, int);
int triinterp(SIMAGE *, CIMAGE *, VIMAGE *, fdata *, int);
int steinterp(SIMAGE *, CIMAGE *, VIMAGE *, fdata, fdata, fdata, fdata, fdata,
fdata, int);

/*------------------------------------------------------------------------------------*/
