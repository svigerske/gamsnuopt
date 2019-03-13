#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* GAMS API */
#include "gmomcc.h"
#include "gevmcc.h"

/* NuOpt API */
#include "nuoIf.h"

#include "DPManager.h"
extern char currentDPM[1024];
extern char currentSolverDPM[1024];

#include "nuopt_integer.h"
using::nuopt::integer;
integer pscctf_(size_t len, const char* str);

// 以下はビルド用ダミー
//
//int SimpleMessageOutput(int){return 0;}
//void  wcspinfo(int,int,double,double,double,double,int){}
//int wcspDebug=0;
//void wcspinfo0(int,int,int){}

//int cpSetVar_para( const int* Vi, double* V ){return 0;}
//extern double cpGetCurrentValue_para( int i, double* V, double* &arg_ValidFlag, int* &arg_ValidFlagNum ){ return 0.0;}

int solveMINLP(
   gmoHandle_t gmo,
   gevHandle_t gev
)
{
   int nfnc, nvar, probType;

   strcpy(currentDPM, "basic_dpm");
   strcpy(currentSolverDPM, "solver_dpm");

   DPMMAP_CREATE(currentDPM);

   //nuoptDispVersion();
   //nuoptOst(stdout,"             (AMPL module)\n");

   nfnc = gmoM(gmo) + 1;
   nvar = gmoN(gmo);
   probType = gmoSense(gmo) == gmoObj_Max ? 1 : 0;

   nuoptParam p;
   p.iisDetect = "off";

   int retVal = setCommon(&p,&nvar,&nfnc,&probType);
   if ( !retVal )
     return 1;

   if( true /* nonlinear */ )
   {
      // TOOD if option "quadra" set, then do adjMethod(0) (line-search) instead
      adjMethod(1); // trust-region
   }
   else if( gmoNDisc(gmo) > 0 )
      adjMethod(3);
   else
      adjMethod(4);

   nuoptOpen22("gamsnuoptlog");
   pscctf_(strlen("gamsnuoptlog"),"gamsnuoptlog");

   nuoptResult* r = nuoptKernel();

   nuoptClose22();

//   if ( amplflag ) {
//     amplpt(r->errorMessage(),r->errorFlag()
//       ,r->optResidual()
//       ,r->optIteration()
//       ,r->nPProbs(),r->nSxPivots(),r->nSxPivotsForMip()
//       ,r->algorithm(),r->optValue(),r->getX(),r->getY());
//   } else {
//     nuoptOutput(r,outputRootName);
//   }

   delete r;

   return 0;
}
