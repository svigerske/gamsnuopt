#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>

/* GAMS API */
#include "gmomcc.h"
#include "gevmcc.h"

/* NuOpt API */
#include "nuoIf.h"
#include "argUtils.h"

#include "DPManager.h"
extern char currentDPM[1024];
extern char currentSolverDPM[1024];

#include "nuopt_integer.h"
using::nuopt::integer;
integer pscctf_(size_t len, const char* str);
//integer pscftc_(integer* len, char* fchar);

gmoHandle_t cur_gmo = NULL;
gevHandle_t cur_gev = NULL;

#define gamsDebug 0

// 以下はビルド用ダミー
//
//int SimpleMessageOutput(int){return 0;}
void  wcspinfo(int,int,double,double,double,double,int){}
//int wcspDebug=0;
void wcspinfo0(int,int,int){}

int cpSetVar_para( const int* Vi, double* V ){return 0;}
//extern double cpGetCurrentValue_para( int i, double* V, double* &arg_ValidFlag, int* &arg_ValidFlagNum ){ return 0.0;}

int nuopqpUsed=0;
int _nuopt_runningIn_VMS = 0; /* MPS 版ではダミー  */
int _SimpleSystem_fpCheck = 0;
int rangeAnalysisProbe = -1;
//extern int wcspResultBestSeed;
//int wcspResultBestSeed = 0;

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

   cur_gmo = gmo;
   cur_gev = gev;

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

extern "C"
{

static
int boundType(
   double* lower,
   double* upper
   )
{
   int lowerExist = 0;
   int upperExist = 0;

   if( *lower == *upper )
      return 0;

   if( *lower != gmoMinf(cur_gmo) )
      lowerExist = 1;
   else
      *lower = 0.0;

   if( *upper != gmoPinf(cur_gmo) )
      upperExist = 1;
   else
      *upper = 0.0;

   if( lowerExist && upperExist )
      return 3;

   if( lowerExist )
      return 1;

   if( upperExist )
      return 2;

   return 4;
}

void attrvrbc_(
   int*    i,
   int*    idir,
   double* dpc,
   double* upc,
   double* prio,
   int*    iunt
   )
{
   if( gamsDebug >= 1 )
      printf("CALL attrvrbc_(%d)\n", *i);

   intC_set(-1, idir);
   doubleC_set(-1.0, dpc);
   doubleC_set(-1.0, upc);
   doubleC_set(0.0, prio); /* TODO priority */
   intC_set(0, iunt);
}

void defvarbc_(
   int*    nvmax,
   int*    nvar,
   int*    kvar,
   double* bl,
   double* bu,
   char*   vrnam,
   int*    ir
   )
{
   char varname[GMS_SSSIZE];
   int i;

   if( gamsDebug >= 1 )
      printf("CALL defvarbc_\n");

   assert(cur_gmo != NULL);

   *nvar = gmoN(cur_gmo);

   if( *nvar > *nvmax )
   {
      if ( gamsDebug >= 1 )
         fprintf(stderr, "defvarbc: Too many variables. nvar = %d\n", *nvar);
      *ir = 10;
      return;
   }

   gmoGetVarLower(cur_gmo, bl);
   gmoGetVarUpper(cur_gmo, bu);
   memset(vrnam, ' ', 8 * *nvar);

   for( i = 0; i < *nvar; ++i)
   {
      gmoGetVarNameOne(cur_gmo, i+1, varname);
      memcpy(vrnam + 8*i, varname, std::min((size_t)8, strlen(varname)));

      kvar[i] = boundType(bl+i, bu+i);

      if( gamsDebug >= 2 )
         printf("DEBUG: bl = %15.8e bu = %15.8e kvar = %d\n", bl[i], bu[i], kvar[i]);
   }
}

void deffunbc_(
   char*   pbnam,
   int*    nfmax,
   int*    nfnc,
   int*    kfun,
   double* cl,
   double* cu,
   char*   funam,
   int*    khess,
   int*    ir
   )
{
   char equname[GMS_SSSIZE];
   int i;

   if( gamsDebug >= 1 )
      printf("CALL deffunbc_\n");

   //???
   //integer len = 40;
   //pscftc_(&len, pbnam);

   *nfnc = gmoM(cur_gmo) + 1;
   if( *nfnc > *nfmax )
   {
      if( gamsDebug >= 1 )
         fprintf(stderr, "deffunbc: Too many functions. nfnc = %d\n", *nfnc);
      *ir = 10;
      return;
   }

   for( i = 0 ; i < *nfnc-1; ++i )
   {
      switch( gmoGetEquTypeOne(cur_gmo, i+1) )
      {
         case gmoequ_E :
         {
            cl[i] = cu[i] = gmoGetRhsOne(cur_gmo, i+1);
            kfun[i] = 0;
            break;
         }

         case gmoequ_L :
         {
            cl[i] = 0.0;  /* nuopt-specific */
            cu[i] = gmoGetRhsOne(cur_gmo, i+1);
            kfun[i] = 2;
            break;
         }

         case gmoequ_G :
         {
            cl[i] = gmoGetRhsOne(cur_gmo, i+1);
            cu[i] = 0.0;  /* nuopt-specific */
            kfun[i] = 1;
            break;
         }

         case gmoequ_N :
         {
            cl[i] = cu[i] = 0.0;
            kfun[i] = 4;
            break;
         }

         default:
         {
            if( gamsDebug >= 1 )
               fprintf(stderr, "deffunbc: Unsupported equation type %d for equ %d\n", gmoGetEquTypeOne(cur_gmo, i), i);
            *ir = 10;
            break;
         }
      }

      if ( gamsDebug >= 2 )
         printf("DEBUG: cl = %15.8e cu = %15.8e kfun = %d\n", cl[i], cu[i], kfun[i]);

      if( gmoGetEquOrderOne(cur_gmo, i+1) != gmoorder_L )
         khess[i] = 3;  /* nonlinear */
      else
         khess[i] = 0;

      gmoGetEquNameOne(cur_gmo, i+1, equname);
      memcpy(funam + 8*i, equname, std::min((size_t)8, strlen(equname)));
   }

   /* objective function */
   kfun[*nfnc-1] = -1;
   if( gmoGetObjOrder(cur_gmo) != gmoorder_L )
      khess[*nfnc-1] = 3;
   else
      khess[*nfnc-1] = 0;

   gmoGetObjName(cur_gmo, equname);
   memcpy(funam + 8*(*nfnc-1), equname, std::min((size_t)8, strlen(equname)));
}

void funlbc_(
   int* nemax,
   int* ne,
   int* ifun,
   int* jvar,
   double* a,
   int* ir
   )
{
   int i, j;

   if( gamsDebug >= 1 )
      printf("CALL funlbc_\n");

   int nzjac = gmoNZ(cur_gmo);
   int nzobj = gmoObjNZ(cur_gmo);
   *ne = nzjac + nzobj;
   if( *nemax < *ne )
   {
      *ir = 10;
      return;
   }

   int* rowstart = (int*) malloc((gmoM(cur_gmo)+1) * sizeof(int));
   gmoGetMatrixRow(cur_gmo, rowstart, jvar, a, NULL);
   for( i = 0; i < gmoM(cur_gmo); ++i )
      for( j = rowstart[i]; j < rowstart[i+1]; ++j )
         ifun[j-1] = i+1;
   assert(nzjac == rowstart[gmoM(cur_gmo)]-1);
   free(rowstart);

   gmoGetObjSparse(cur_gmo, jvar+nzjac, a+nzjac, NULL, &nzobj, &i);
   for( j = 0; j < nzobj; ++j )
      ifun[nzjac+j] = gmoM(cur_gmo)+1;

   *ir = 0;
}

void funnlbc_(
   int* nvar,
   int* nfnc,
   double* x,
   double* f,
   int* ir
   )
{
   int numerr = 0;
   int i;

   if( gamsDebug >= 1 )
      printf("CALL funnlbc_\n");

   assert(*nvar == gmoN(cur_gmo));
   assert(*nfnc == gmoM(cur_gmo)+1);

   *ir = 0;

   for( i = 0; i < *nfnc-1; ++i )
   {
      if( gmoGetEquOrderOne(cur_gmo, i+1) != gmoorder_L )
      {
         gmoEvalFunc(cur_gmo, i+1, x, f+i, &numerr);
         if( numerr > 0 )
            *ir = 10;
      }
      else
         f[i] = 0.0;
   }

   if( gmoGetObjOrder(cur_gmo) != gmoorder_L )
   {
      gmoEvalFuncObj(cur_gmo, x, f+*nfnc-1, &numerr);
      if( numerr > 0 )
         *ir = 10;
   }
   else
      f[*nfnc-1] = 0.0;
}

void funqbc_(
   int* nemax,
   int* ne,
   int* ifun,
   int* j1var,
   int* j2var,
   double* h,
   int* ir
   )
{
   if( gamsDebug >= 1 )
      printf("CALL funqbc_\n");

   /* TODO */

   *ne = 0;
}

void gradnlbc_(
   int*    nvar,
   double* x,
   int*    nemax,
   int*    ne,
   int*    ifun,
   int*    jvar,
   double* a,
   int*    ir
   )
{
   int i;
   int pos;
   int numerr;
   double f;
   double gx;
   double* grad;

   if( gamsDebug >= 1 )
      printf("CALL gradnlbc_\n");

   assert(*nvar == gmoN(cur_gmo));

   int nzjac = gmoNZ(cur_gmo);
   int nzobj = gmoObjNZ(cur_gmo);
   *ne = nzjac + nzobj;
   if( *nemax < *ne )
   {
      *ir = 10;
      return;
   }

   // TODO it seems that we could skip over linear constraints/obj here

   grad = (double*)malloc(*nvar * sizeof(double));

   pos = 0;
   for( i = 0; i < gmoM(cur_gmo); ++i )
   {
      gmoEvalGrad(cur_gmo, i+1, x, &f, grad, &gx, &numerr);
      /* TODO numerr */

      void* jacptr = NULL;
      double jacval;
      int colidx;
      int nlflag;

      do
      {
         gmoGetRowJacInfoOne(cur_gmo, i+1, &jacptr, &jacval, &colidx, &nlflag);

         assert(pos < *ne);
         ifun[pos] = i+1;
         jvar[pos] = colidx;
         a[pos] = grad[colidx-1];

         if( gamsDebug >= 2 )
            printf("pos %d: ifun %d jvar %d a %g\n", pos, ifun[pos], jvar[pos], a[pos]);

         ++pos;
      }
      while( jacptr != NULL );
   }
   assert(pos == nzjac);

   gmoEvalGradObj(cur_gmo, x, &f, grad, &gx, &numerr);
   /* TODO numerr */

   int* objcolidx = new int[nzobj];
   gmoGetObjSparse(cur_gmo, objcolidx, a+pos, NULL, &nzobj, &i);
   for( i = 0; i < nzobj; ++i )
   {
      assert(pos < *ne);
      ifun[pos] = gmoM(cur_gmo)+1;
      jvar[pos] = objcolidx[i];
      a[pos] = grad[objcolidx[i]-1];
      ++pos;
   }
   assert(pos == nzjac+nzobj);
   delete[] objcolidx;

   free(grad);

   *ir = 0;
}

void hessnlbc_(
   int*    nvar,
   int*    nfnc,
   double* x,
   double* y,
   int*    nemax,
   int*    ne,
   int*    j1var,
   int*    j2var,
   double* h,
   int*    ir
   )
{
   if( gamsDebug >= 1 )
      printf("CALL hessnlbc_\n");

   /* TODO nonlinear */
   *ne = 0;

   *ir = 0;
}

void initvrbc_(
   int*    nvar,
   double* xini,
   int*    ir
   )
{
   if( gamsDebug >= 1 )
      printf("CALL initvrbc_\n");

   assert(*nvar == gmoN(cur_gmo));

   gmoGetVarL(cur_gmo, xini);

   *ir = 0;
}

// probably only called through typevr_C
void typevrbc_(
   int* nvar,
   int* ivtyp
   )
{
   int i;

   if( gamsDebug >= 1 )
      printf("CALL typevrbc_\n");

   assert(*nvar == gmoN(cur_gmo));

   for( i = 0; i < *nvar; ++i )
   {
      switch( gmoGetVarTypeOne(cur_gmo, i+1) )
      {
         case gmovar_X :
            ivtyp[i] = 0;
            break;
         case gmovar_B :
            ivtyp[i] = 2;
            break;
         case gmovar_I :
            ivtyp[i] = 1;
            break;
      }
   }
}

void typevr_C(
   int* nvar,
   int* ivtyp
   )
{
   if( gamsDebug >= 10 )
      printf("CALL typevr_C\n");

   typevrbc_(nvar, ivtyp);
}

}
