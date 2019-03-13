#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* GAMS API */
#include "gmomcc.h"
#include "gevmcc.h"
#include "optcc.h"

/* NuOpt API */
#include "nuoIf.h"

int solveMIQCP(
   gmoHandle_t gmo,
   gevHandle_t gev
)
{
   nuoptParam options;
   int n;
   int m;
   int minimize;
   double* x0;
   double* bL;
   double* bU;
   int* ibL;
   int* ibU;
   double* cL;
   double* cU;
   int* icL;
   int* icU;
   double* objL;
   int nAelem;
   int* rowstart;
   int* irowA;
   int* jcolA;
   double* a;
   int nQelem;
   int* irowQ = NULL;
   int* jcolQ = NULL;
   double* q = NULL;
   int nQCelem;
   int* ifunQC = NULL;
   int* irowQC = NULL;
   int* jcolQC = NULL;
   double* qc = NULL;
   int qcpos;
   int* ivtype;
   nuoptResult* res;

   int i;

   assert(gmo != NULL);
   assert(gev != NULL);

   if( gmoModelType(gmo) >= gmoProc_qcp )
      gmoUseQSet(gmo, 1);

   n = gmoN(gmo);
   m = gmoM(gmo);

   minimize = (gmoSense(gmo) == gmoObj_Min);

   x0 = (double*) malloc(n * sizeof(double));
   gmoGetVarL(gmo, x0);

   bL = (double*) malloc(n * sizeof(double));
   bU = (double*) malloc(n * sizeof(double));
   ibL = (int*) malloc(n * sizeof(int));
   ibU = (int*) malloc(n * sizeof(int));

   gmoGetVarLower(gmo, bL);
   gmoGetVarUpper(gmo, bU);
   for( i = 0; i < n; ++i )
   {
      ibL[i] = (bL[i] != gmoMinf(gmo));
      ibU[i] = (bU[i] != gmoPinf(gmo));
   }

   cL = (double*) malloc(m * sizeof(double));
   cU = cL;
   icL = (int*) malloc(m * sizeof(double));
   icU = (int*) malloc(m * sizeof(double));

   gmoGetRhs(gmo, cL);
   for( i = 0; i < m; ++i )
   {
      switch( gmoGetEquTypeOne(gmo, i+1) )
      {
         case gmoequ_E :  /* equality constraint */
            icL[i] = 1;
            icU[i] = 1;
            break;

         case gmoequ_G :  /* ax >= rhs constraint */
            icL[i] = 1;
            icU[i] = 0;
            break;

         case gmoequ_L :  /* ax <= rhs constraint */
            icL[i] = 0;
            icU[i] = 1;
            break;

         case gmoequ_N:  /* free constraint - unlikely */
            icL[i] = 0;
            icU[i] = 0;
            break;

         case gmoequ_X:
         case gmoequ_C:
         case gmoequ_B:
            /* these should not occur */
            goto TERMINATE;
      }
   }

   objL = (double*) malloc(n * sizeof(double));
   gmoGetObjVector(gmo, objL, NULL);

   nAelem = gmoNZ(gmo);
   rowstart = (int*) malloc((m+1) * sizeof(int));
   irowA = (int*) malloc(nAelem * sizeof(int));
   jcolA = (int*) malloc(nAelem * sizeof(int));
   a = (double*) malloc(nAelem * sizeof(double));
   gmoGetMatrixRow(gmo, rowstart, jcolA, a, NULL);
   for( i = 0; i < m; ++i )
   {
      int j;
      for( j = rowstart[i]; j < rowstart[i+1]; ++j )
         irowA[j-1] = i+1;
   }
   // somehow we might get less elements from gmoGetMatrixRow than gmoNZ(), probably because of Q-stuff
   nAelem = rowstart[m]-1;
   free(rowstart);

   nQelem = gmoObjQNZ(gmo);
   if( nQelem > 0 )
   {
      irowQ = (int*) malloc(nQelem * sizeof(int));
      jcolQ = (int*) malloc(nQelem * sizeof(int));
      q = (double*) malloc(nQelem * sizeof(double));
      gmoGetObjQ(gmo, jcolQ, irowQ, q);
      /* ??? NuOpt multiplies by 1/2; GAMS diagonal elements are already multiplied by 2 */
      for( i = 0; i < nQelem; ++i )
      {
         if( irowQ[i] != jcolQ[i] )
            q[i] *= 2.0;
         // gmoGetObjQ does not seem to have taken indexbase=1 into account!
         irowQ[i]++;
         jcolQ[i]++;
      }
   }

   nQCelem = 0;
   if( gmoModelType(gmo) >= gmoProc_qcp )
      for( i = 0; i < m; ++i )
         nQCelem += gmoGetRowQNZOne(gmo, i+1);
   if( nQCelem > 0 )
   {
      ifunQC = (int*) malloc(nQCelem * sizeof(int));
      irowQC = (int*) malloc(nQCelem * sizeof(int));
      jcolQC = (int*) malloc(nQCelem * sizeof(int));
      qc = (double*) malloc(nQCelem * sizeof(double));
      qcpos = 0;
      for( i = 0; i < m; ++i )
      {
         int qnz;
         int j;

         qnz = gmoGetRowQNZOne(gmo, i+1);
         if( qnz == 0 )
            continue;

         gmoGetRowQ(gmo, i+1, jcolQC + qcpos, irowQC + qcpos, qc + qcpos);
         for( j = 0; j < qnz; ++j )
         {
            ifunQC[qcpos + j] = i;
            if( jcolQC[qcpos + j] != irowQC[qcpos + j] )
               qc[qcpos + j] *= 2.0;
//            // gmoGetRowQ does not seem to have taken indexbase=1 into account!
//            jcolQC[qcpos+j]++;
//            irowQC[qcpos+j]++;
         }

         qcpos += qnz;
      }
   }

   ivtype = NULL;
   if( gmoNDisc(gmo) > 0 )
   {
      ivtype = (int*) malloc(n * sizeof(int));
      for( i = 0; i < n; ++i )
      {
         switch( gmoGetVarTypeOne(gmo, i+1) )
         {
            case gmovar_X :
               ivtype[i] = 0;
               break;
            case gmovar_B :
               ivtype[i] = 2;
               break;
            case gmovar_I :
               ivtype[i] = 1;
               break;

            /* not supported: */
            case gmovar_S1:  /* SOS1 */
            case gmovar_S2:  /* SOS2 */
            case gmovar_SC:  /* semi-continous */
               ivtype[i] = 0;
               goto TERMINATE;
            case gmovar_SI:  /* semi-integer */
               ivtype[i] = 1;
               goto TERMINATE;
         }
      }
   }

   //options.outputMode = "debug";

   res = solveQP(&options, n, m, minimize, x0,
      bL, bU, ibL, ibU,
      cL, cU, icL, icU,
      objL, nAelem, irowA, jcolA, a,
      nQelem, irowQ, jcolQ, q,
      nQCelem, ifunQC, irowQC, jcolQC, qc,
      ivtype, NULL, NULL, NULL, NULL, NULL);

   printf("result: %p\n", (void*)res);
   printf("error: %s\n", res->errorMessage());

   delete res;

TERMINATE:

   return 0;
}
