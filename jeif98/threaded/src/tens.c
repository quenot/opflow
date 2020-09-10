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

#include <omp.h>

/*------------------------------------------------------------------------------------*/
/*------------------------ ARRAY ALLOCATE AND FREE PROCEDURES ------------------------*/
/*------------------------------------------------------------------------------------*/
/* mallock allocates a k-dimensional array of elements with parameterizable size.     */
/* freek frees a k-dimensional array of elements with parameterizable size.           */
/* altensk and frtensk are similar to mallock and freek except that for each          */
/* dimension, the size is defined by a lower and upper bound. k may be 1, 2 or 3.     */
/* The array pointer address must be passed casted to (void **), (void ***) and       */
/* (void ****) respectively. A SUCCESS/FAILURE status is returned.                    */
/*------------------------------------------------------------------------------------*/

int UsedMem = 0;
int MaxUsedMem = 0;

/*------------------------------------------------------------------------------------*/

int UsedMemory() {return(UsedMem);}
int MaxUsedMemory() {return(MaxUsedMem);}

#ifdef dec
void *memalign(size_t alignment, size_t size) {return malloc(size);}
#endif
#ifdef linux
void *memalign(size_t alignment, size_t size) {return calloc(size,1);}
#endif

/*------------------------------------------------------------------------------------*/

int malloc1(void **v,int n,size_t elsize)
{
    if (n <= 0) {
        fprintf(stderr,"Error in malloc1: n <= 0\n");
        return(FAILURE);
    }
    if (!(*v = (void *) memalign(4,n*elsize))) {
        fprintf(stderr,"Error in malloc1: malloc() failed\n");
        return(FAILURE);
    }
    #pragma omp atomic
    UsedMem += n*elsize;
    if (UsedMem > MaxUsedMem) MaxUsedMem = UsedMem;
    return(SUCCESS);
}

int free1(void **v,int n,size_t elsize)
{
    if (n <= 0) {
        fprintf(stderr,"Error in free1: n <= 0\n");
        return(FAILURE);
    }
    free(*v);
    if (errno == EINVAL) {
        fprintf(stderr,"Error in free1: free() failed\n");
        return(FAILURE);
    }
    #pragma omp atomic
    UsedMem -= n*elsize;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int malloc2(void ***m,int nr,int nc,size_t elsize)
{
    int i;
    if (malloc1((void **) m,nr,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in malloc2: malloc1 failed\n");
        return(FAILURE);
    }
    for (i = 0; i < nr; i++) {
        if (malloc1(*m+i,nc,elsize) == FAILURE) {
            fprintf(stderr,"Error in malloc2: malloc1 failed\n");
            return(FAILURE);
        }
    }
    return(SUCCESS);
}

int free2(void ***m,int nr,int nc,size_t elsize)
{
    int i;
    for (i = 0; i < nr; i++) {
        if (free1(*m+i,nc,elsize) == FAILURE) {
            fprintf(stderr,"Error in free2: free1 failed\n");
            return(FAILURE);
        }
    }
    if (free1((void **) m,nr,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in free2: free1 failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int malloc3(void ****t,int np,int nr,int nc,size_t elsize)
{
    int i;
    if (malloc1((void **) t,np,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in malloc3: malloc1 failed\n");
        return(FAILURE);
    }
    for (i = 0; i < np; i++) {
        if (malloc2(*t+i,nr,nc,elsize) == FAILURE) {
            fprintf(stderr,"Error in malloc3: malloc2 failed\n");
            return(FAILURE);
        }
    }
    return(SUCCESS);
}

int free3(void ****t,int np,int nr,int nc,size_t elsize)
{
    int i;
    for (i = 0; i < np; i++) {
        if (free2(*t+i,nr,nc,elsize) == FAILURE) {
            fprintf(stderr,"Error in free3: free2 failed\n");
            return(FAILURE);
        }
    }
    if (free1((void **) t,np,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in free3: free1 failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int altens1(void **v,int nl,int nh,size_t elsize)
{
    char *p;
    if (nh < nl) {
        fprintf(stderr,"Error in altens1: nh < nl\n");
        return(FAILURE);
    }
    p = (char *) memalign(4,(nh+1-nl)*elsize);
    if (!p) {
        fprintf(stderr,"Error in altens1: malloc() failed\n");
        return(FAILURE);
    }
    *v = (void *) (p-nl*elsize);
    #pragma omp atomic
    UsedMem += (nh+1-nl)*elsize;
    if (UsedMem > MaxUsedMem) MaxUsedMem = UsedMem;
    return(SUCCESS);
}

int frtens1(void **v,int nl,int nh,size_t elsize)
{
    char *p;
    if (nh < nl) {
        fprintf(stderr,"Error in frtens1: nh < nl\n");
        return(FAILURE);
    }
    p = ((char *) (*v))+nl*elsize;
    free((void *) p);
    if (errno == EINVAL) {
        fprintf(stderr,"Error in frtens1: free() failed\n");
        return(FAILURE);
    }
    #pragma omp atomic
    UsedMem -= (nh+1-nl)*elsize;
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int altens2(void ***m,int nrl,int nrh,int ncl,int nch,size_t elsize)
{
    int i;
    if (altens1((void **) m,nrl,nrh,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in altens2: altens1 failed\n");
        return(FAILURE);
    }
    for (i = nrl; i <= nrh; i++) {
        if (altens1(*m+i,ncl,nch,elsize) == FAILURE) {
            fprintf(stderr,"Error in altens2: altens1 failed\n");
            return(FAILURE);
        }
    }
    return(SUCCESS);
}

int frtens2(void ***m,int nrl,int nrh,int ncl,int nch,size_t elsize)
{
    int i;
    for (i = nrh; i >= nrl; i--) {
        if (frtens1(*m+i,ncl,nch,elsize) == FAILURE) {
            fprintf(stderr,"Error in frtens2: frtens1 failed\n");
            return(FAILURE);
        }
    }
    if (frtens1((void **) m,nrl,nrh,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in frtens2: frtens1 failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/

int altens3(void ****t,int npl,int nph,int nrl,int nrh,int ncl,int nch,size_t elsize)
{
    int i;
    if (altens1((void **) t,npl,nph,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in altens3: altens1 failed\n");
        return(FAILURE);
    }
    for (i = npl; i <= nph; i++) {
        if (altens2(*t+i,nrl,nrh,ncl,nch,elsize) == FAILURE) {
            fprintf(stderr,"Error in altens3: altens2 failed\n");
            return(FAILURE);
        }
    }
    return(SUCCESS);
}

int frtens3(void ****t,int npl,int nph,int nrl,int nrh,int ncl,int nch,size_t elsize)
{
    int i;
    for (i = nph; i >= npl; i--) {
        if (frtens2(*t+i,nrl,nrh,ncl,nch,elsize) == FAILURE) {
            fprintf(stderr,"Error in frtens3: frtens2 failed\n");
            return(FAILURE);
        }
    }
    if (frtens1((void **) t,npl,nph,sizeof(void *)) == FAILURE) {
        fprintf(stderr,"Error in frtens3: frtens1 failed\n");
        return(FAILURE);
    }
    return(SUCCESS);
}

/*------------------------------------------------------------------------------------*/
