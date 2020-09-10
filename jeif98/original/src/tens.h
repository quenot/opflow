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

#include <malloc.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifdef sun
#ifndef slrs
#include <sys/stdtypes.h>
#endif
#endif

/*------------------------------------------------------------------------------------*/

int UsedMemory();
int MaxUsedMemory();
int malloc1(void **, int, size_t);
int free1(void **, int, size_t);
int malloc2(void ***, int, int, size_t);
int free2(void ***, int, int, size_t);
int malloc3(void ****, int, int, int, size_t);
int free3(void ****, int, int, int, size_t);
int altens1(void **, int, int, size_t);
int frtens1(void **, int, int, size_t);
int altens2(void ***, int, int, int, int, size_t);
int frtens2(void ***, int, int, int, int, size_t);
int altens3(void ****, int, int, int, int, int, int, size_t);
int frtens3(void ****, int, int, int, int, int, int, size_t);

/*------------------------------------------------------------------------------------*/
