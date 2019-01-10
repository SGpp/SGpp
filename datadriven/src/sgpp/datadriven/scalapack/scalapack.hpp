#pragma once

extern "C" {

/**
 *  ictxt   (global output) int
 *          ICTXT specifies the BLACS context handle identifying the
 *          created process grid.  The context itself is global.
 *
 *  nprow   (global input) int
 *          NPROW specifies the number of process rows in the grid
 *          to be created.
 *
 *  npcol   (global input) int
 *          NPCOL specifies the number of process columns in the grid
 *          to be created.
 */
void sl_init_(int &ictxt, int &nprow, int &npcol);

void pdgsev(int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, int *ipiv, double *b,
            int *ib, int *jb, int *descb, int *info);

void pdgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const double *alpha, const double *a, const int *ia, const int *ja, const int *desca,
            const double *b, const int *ib, const int *jb, const int *descb, const double *beta,
            double *c, const int *ic, const int *jc, const int *descc);

void pdlaset_(char &uplo, int &m, int &n, double &alpha, double &beta, double *a, int &ia, int &ja,
              int *desca);

void pdlaprnt_(int &m, int &n, double *a, int &ia, int &ja, int *desca, int &irprnt, int &icprnt,
               const char *cmatnm, int &nout, double *work);

void descinit_(int *desc, int m, int n, int mb, int nb, int irsrc, int icsrc, int ictxt, int lld,
               int info);
}