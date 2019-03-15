#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>

/* GAMS API */
#include "gmomcc.h"
#include "gevmcc.h"
#include "optcc.h"

#if defined(_WIN32)
#if !defined(STDCALL)
#define STDCALL __stdcall
#endif
#if !defined(DllExport)
#define DllExport __declspec(dllexport)
#endif
#else
#if !defined(STDCALL)
#define STDCALL
#endif
#if !defined(DllExport)
#define DllExport
#endif
#endif

/* NuOpt API */
#include "nuoIf.h"
#include "nuopt_exception.h"
#include "argUtils.h"

#include "DPManager.h"
extern char currentDPM[1024];
extern char currentSolverDPM[1024];

#include "nuopt_integer.h"
using::nuopt::integer;

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
integer pscctf_(size_t len, const char* str);
extern "C" integer pscftc_(integer* len, char* fchar);


#define gamsDebug 0
#define nuoptDebug 0

gmoHandle_t cur_gmo = NULL;
gevHandle_t cur_gev = NULL;

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
#if gamsDebug >= 1
   printf("CALL attrvrbc_(%d)\n", *i);
#endif

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

#if gamsDebug >= 1
   printf("CALL defvarbc_\n");
#endif

   assert(cur_gmo != NULL);

   *nvar = gmoN(cur_gmo);

   if( *nvar > *nvmax )
   {
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

#if gamsDebug >= 2
      printf("  bl = %15.8e bu = %15.8e kvar = %d\n", bl[i], bu[i], kvar[i]);
#endif
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

#if gamsDebug >= 1
   printf("CALL deffunbc_\n");
#endif

   integer len = 40;
   memset(pbnam, ' ', len);
   gmoNameModel(cur_gmo, equname);
   memcpy(pbnam, equname, std::min((size_t)len, strlen(equname)));

   pscftc_(&len, pbnam); // why do I call this?

   *nfnc = gmoM(cur_gmo) + 1;
   if( *nfnc > *nfmax )
   {
      *ir = 10;
      return;
   }

   memset(funam, ' ', 8 * *nfnc);

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
            fprintf(stderr, "deffunbc: Unsupported equation type %d for equ %d\n", gmoGetEquTypeOne(cur_gmo, i), i);
            *ir = 1;
            break;
         }
      }

      if( gmoGetEquOrderOne(cur_gmo, i+1) != gmoorder_L )
         khess[i] = 3;  /* nonlinear */
      else
         khess[i] = 0;

#if gamsDebug >= 2
      printf("  cl = %15.8e cu = %15.8e kfun = %d khess = %d\n", cl[i], cu[i], kfun[i], khess[i]);
#endif

      gmoGetEquNameOne(cur_gmo, i+1, equname);
      memcpy(funam + 8*i, equname, std::min((size_t)8, strlen(equname)));
   }

   /* objective function */
   cl[*nfnc-1] = cu[*nfnc-1] = 0.0;
   kfun[*nfnc-1] = -1;
   if( gmoGetObjOrder(cur_gmo) != gmoorder_L )
      khess[*nfnc-1] = 3;
   else
      khess[*nfnc-1] = 0;

#if gamsDebug >= 2
   printf("  cl = %15.8e cu = %15.8e kfun = %d khess = %d\n", cl[*nfnc-1], cu[*nfnc-1], kfun[*nfnc-1], khess[*nfnc-1]);
#endif

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
   int rownz;
   int i, j;

#if gamsDebug >= 1
   printf("CALL funlbc_\n");
#endif

   /* we switch to *ir == 10 if *nemax is too small
    * if that happens, we continue to update *ne only
    */
   *ir = 0;

   *ne = 0;
   for( i = 0; i < gmoM(cur_gmo); ++i )
   {
      if( gmoGetEquOrderOne(cur_gmo, i+1) != gmoorder_L )
         continue;

      gmoGetRowStat(cur_gmo, i+1, &rownz, &j, &j);
      if( *ne + rownz > *nemax )
         *ir = 10;

      if( *ir == 0 )
      {
         gmoGetRowSparse(cur_gmo, i+1, jvar + *ne, a + *ne, NULL, &rownz, &j);
         for( j = 0; j < rownz; ++j )
            ifun[*ne + j] = i+1;
      }
      *ne += rownz;
   }

   if( gmoGetObjOrder(cur_gmo) == gmoorder_L )
   {
      if( *ne + gmoObjNZ(cur_gmo) > *nemax )
         *ir = 10;

      if( *ir == 0 )
      {
         gmoGetObjSparse(cur_gmo, jvar + *ne, a + *ne, NULL, &rownz, &i);
         assert(rownz == gmoObjNZ(cur_gmo));
         for( j = 0; j < rownz; ++j )
            ifun[*ne + j] = gmoM(cur_gmo) + 1;
      }
      *ne += gmoObjNZ(cur_gmo);
   }

#if gamsDebug >= 2
   if( *ir == 0 )
      for( i = 0; i < *ne; ++i )
         printf("  ifun %d jvar %d a %g\n", ifun[i], jvar[i], a[i]);
#endif
#if gamsDebug >= 1
   if( *ir == 10 )
      printf("  *nemax = %d too small, requesting %d\n", *nemax, *ne);
#endif
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

#if gamsDebug >= 1
   printf("CALL funnlbc_\n");
#endif

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

#if gamsDebug >= 2
      printf("  equorder %d f = %g\n", gmoGetEquOrderOne(cur_gmo, i+1), f[i]);
#endif
   }

   if( gmoGetObjOrder(cur_gmo) != gmoorder_L )
   {
      gmoEvalFuncObj(cur_gmo, x, f+*nfnc-1, &numerr);
      if( numerr > 0 )
         *ir = 10;
   }
   else
      f[*nfnc-1] = 0.0;

#if gamsDebug >= 2
   printf("  objorder %d f = %g\n", gmoGetObjOrder(cur_gmo), f[*nfnc-1]);
#endif
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
#if gamsDebug >= 1
   printf("CALL funqbc_\n");
#endif

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
   int i, j;
   int rownz;
   int numerr;
   int evalerror = 0;
   double f;
   double gx;
   double* grad;

#if gamsDebug >= 1
   printf("CALL gradnlbc_\n");
#endif

   assert(*nvar == gmoN(cur_gmo));

   // OMG! ne と nemax に同じアドレスが亘る場合がある -> ne and nemax may have the same address
   int nemaxSaved = *nemax;

   /* we switch to *ir == 10 if *nemax is too small
    * if that happens, we continue to update *ne only
    */
   *ir = 0;
   *ne = 0;

   /* if everything linear, then just return */
   if( gmoObjNLNZ(cur_gmo) == 0 && gmoNLNZ(cur_gmo) == 0 )
      return;

   grad = new double[*nvar];

   for( i = 0; i < gmoM(cur_gmo); ++i )
   {
      if( gmoGetEquOrderOne(cur_gmo, i+1) == gmoorder_L )
         continue;

      gmoGetRowStat(cur_gmo, i+1, &rownz, &j, &j);

      if( *ir == 10 || *ne + rownz > nemaxSaved )
      {
         *ir = 10;
         *ne += rownz;
         continue;
      }

      gmoEvalGrad(cur_gmo, i+1, x, &f, grad, &gx, &numerr);
      if( numerr )
         evalerror = 1;

      void* jacptr = NULL;
      double jacval;
      int colidx;
      int nlflag;

      do
      {
         gmoGetRowJacInfoOne(cur_gmo, i+1, &jacptr, &jacval, &colidx, &nlflag);

         ifun[*ne] = i+1;
         jvar[*ne] = colidx;
         a[*ne] = grad[colidx-1];

         ++*ne;
      }
      while( jacptr != NULL );
   }

   if( gmoGetObjOrder(cur_gmo) != gmoorder_L )
   {
      if( *ir == 10 || *ne + gmoObjNZ(cur_gmo) > nemaxSaved )
      {
         *ir = 10;
         *ne += gmoObjNZ(cur_gmo);
      }
      else
      {
         gmoEvalGradObj(cur_gmo, x, &f, grad, &gx, &numerr);
         if( numerr )
            evalerror = 1;

         gmoGetObjSparse(cur_gmo, jvar + *ne, a + *ne, NULL, &rownz, &i);
         assert(rownz == gmoObjNZ(cur_gmo));
         for( j = 0; j < rownz; ++j )
         {
            ifun[*ne + j] = gmoM(cur_gmo) + 1;
            a[*ne + j] = grad[jvar[*ne + j] - 1];
         }
         *ne += rownz;
      }
   }

   delete[] grad;

#if gamsDebug >= 2
   if( *ir == 0 )
      for( i = 0; i < *ne; ++i )
         printf("  ifun %d jvar %d a %g\n", ifun[i], jvar[i], a[i]);
#endif
#if gamsDebug >= 1
   if( *ir == 10 )
      printf("  *nemax = %d too small, require %d\n", nemaxSaved, *ne);
#endif

   if( evalerror )
      *ir = 10;
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
   int numerr;

#if gamsDebug >= 1
   printf("CALL hessnlbc_\n");
#endif

   assert(*nvar == gmoN(cur_gmo));
   assert(*nfnc == gmoM(cur_gmo) + 1);

   /* if everything linear, then just return */
   if( gmoObjNLNZ(cur_gmo) == 0 && gmoNLNZ(cur_gmo) == 0 )
   {
      *ne = 0;
      return;
   }

   *ne = gmoHessLagNz(cur_gmo);
   if( *ne > *nemax )
   {
      *ir = 10;
      return;
   }

   gmoHessLagStruct(cur_gmo, j1var, j2var);
   gmoHessLagValue(cur_gmo, x, y, h, 1.0, 1.0, &numerr);

#if gamsDebug >= 2
   for( int i = 0; i < *ne; ++i )
      printf("  j1var %d j2var %d h %g\n", j1var[i], j2var[i], h[i]);
#endif

   *ir = numerr ? 10 : 0;
}

void initvrbc_(
   int*    nvar,
   double* xini,
   int*    ir
   )
{
#if gamsDebug >= 1
   printf("CALL initvrbc_\n");
#endif

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

#if gamsDebug >= 1
   printf("CALL typevrbc_\n");
#endif

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
         default:
            gevLogStat(cur_gev, "Unsupported variable type.");
            // how to signal an error here?
            ivtyp[i] = 0;
            break;
      }
#if gamsDebug >= 2
      printf("  ivtype = %d\n", ivtyp[i]);
#endif
   }
}

void typevr_C(
   int* nvar,
   int* ivtyp
   )
{
#if gamsDebug >= 1
   printf("CALL typevr_C\n");
#endif

   typevrbc_(nvar, ivtyp);
}


DllExport void STDCALL nuoXCreate(void** Cptr)
{
   assert(Cptr != NULL);

   // GAMS needs a non-Null pointer as return
   *Cptr = (void*)0x1;
}

DllExport int STDCALL nuocreate(void** Cptr, char* msgBuf, int msgBufLen)
{
   assert(Cptr != NULL);
   assert(msgBufLen > 0);
   assert(msgBuf != NULL);

   // GAMS needs a non-Null pointer as return
   *Cptr = (void*)0x1;

   msgBuf[0] = 0;

   return 1;
}

DllExport void STDCALL nuoXFree(void** Cptr)
{
   assert(Cptr != NULL);
   assert(*Cptr == (void*)0x1);

   cur_gmo = NULL;
   cur_gev = NULL;

   gmoLibraryUnload();
   gevLibraryUnload();
}

DllExport int STDCALL nuofree(void** Cptr)
{
   nuoXFree(Cptr);

   return 1;
}

/* comp returns the compatibility mode:
           0: client is too old for the DLL, no compatibility
           1: client version and DLL version are the same, full compatibility
           2: client is older than DLL, but defined as compatible, backward compatibility
           3: client is newer than DLL, forward compatibility
           FIXME: for now, we just claim full compatibility
 */
DllExport int STDCALL C__nuoXAPIVersion(int api, char* Msg, int* comp)
{
   *comp = 1;
   return 1;
}

DllExport int STDCALL D__nuoXAPIVersion(int api, char* Msg, int* comp)
{
   *comp = 1;
   return 1;
}

DllExport int STDCALL C__nuoXCheck(const char* funcn, int ClNrArg, int Clsign[], char* Msg) { return 1; }

DllExport int STDCALL D__nuoXCheck(const char* funcn, int ClNrArg, int Clsign[], char* Msg) { return 1; }

DllExport int STDCALL C__nuoReadyAPI(void* Cptr, gmoHandle_t Gptr, optHandle_t Optr)
{
   assert(Cptr == (void*)0x1);
   assert(Gptr != NULL);
   assert(Optr == NULL);

   char msg[256];
   if( !gmoGetReady(msg, sizeof(msg)) )
      return 1;
   if( !gevGetReady(msg, sizeof(msg)) )
      return 1;

   cur_gmo = Gptr;
   cur_gev = (gevHandle_t)gmoEnvironment(cur_gmo);

   strcpy(currentDPM, "basic_dpm");
   strcpy(currentSolverDPM, "solver_dpm");

   DPMMAP_CREATE(currentDPM);

   nuoptDispVersion();
   //nuoptOst(stdout,"             (GAMS module)\n");

   return 0;
}

DllExport int STDCALL C__nuoCallSolver(void* Cptr)
{
   assert(Cptr == (void*)0x1);
   assert(cur_gmo != NULL);
   assert(cur_gev != NULL);

   int nfnc, nvar, probType;

   /* reformulate objective variable out of model, if possible */
   gmoObjStyleSet(cur_gmo, gmoObjType_Fun);
   gmoObjReformSet(cur_gmo, 1);

   /* NuOpt uses 1-based indexing */
   gmoIndexBaseSet(cur_gmo, 1);

   nfnc = gmoM(cur_gmo) + 1;
   nvar = gmoN(cur_gmo);
   probType = gmoSense(cur_gmo) == gmoObj_Max ? 1 : 0;

   nuoptParam p;
   p.iisDetect = "off";

   if( gmoOptFile(cur_gmo) )
   {
      char optfilename[GMS_SSSIZE];
      gmoNameOptFile(cur_gmo, optfilename);
      printf("Read optfile %s\n", optfilename);
      p.read(optfilename);
      p.outputMode = "normal";
   }
#if nuoptDebug > 1
   p.outputMode = "debug";
#endif

   int retVal = setCommon(&p, &nvar, &nfnc, &probType);
   if( !retVal )
     return 1;

   if( gmoNLNZ(cur_gmo) > 0 || gmoObjNLNZ(cur_gmo) > 0 )  /* nonlinear */
   {
      // TODO if option "quadra" set, then do adjMethod(0) (line-search) instead
      adjMethod(1); // trust-region

      int do2dir = 0;
      int dohess = 1;
      gmoHessLoad(cur_gmo, 0, &do2dir, &dohess);
      if( !dohess )
      {
         gevLogStat(cur_gev, "Failed to initialize Hessian structure!");
         return 1;
      }
   }
   else if( gmoNDisc(cur_gmo) > 0 )
      adjMethod(3);
   else
      adjMethod(4);

#if nuoptDebug > 0
   nuoptOpen22("gamsnuoptlog");
   pscctf_(strlen("gamsnuoptlog"),"gamsnuoptlog");
#endif

   nuoptResult* r = NULL;
   try
   {
      r = nuoptKernel();
   }
   catch( const std::exception& e )
   {
      gevLogStat(cur_gev, e.what());
      return 1;
   }

#if nuoptDebug > 0
   nuoptClose22();

   nuoptOutput(r, "nuopt.out");
#endif

#if gamsDebug > 0
   if( r->errorMessage()[0] == '\0' )
      printf("Optimal\n");
   else
      printf("Non-Optimal: %s\n", r->errorMessage());

   if( r->errorFlag() != 5 )
   {
      printf("Objective: %g\n", r->optValue());
      int i;
#if gamsDebug > 1
      for( i = 0; i < gmoN(cur_gmo); ++i )
         printf("X[%d] = %g\n", i, r->getX()[i]);
#endif
   }
#endif

   // https://translate.googleusercontent.com/translate_c?depth=1&hl=en&prev=search&rurl=translate.google.com&sl=ja&sp=nmt4&u=http://www.msi.co.jp/nuopt/docs/v20/manual/html/0A-02-01.html&xid=17259,15700023,15700186,15700191,15700248,15700253&usg=ALkJrhjy7-fngwOOYutfKNgGqtO6PURY9w
   int havesol = 0;
   switch( r->errorFlag() )
   {
      case 0:
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         if( gmoNDisc(cur_gmo) > 0 )
            gmoModelStatSet(cur_gmo, gmoModelStat_OptimalGlobal);
         else
            gmoModelStatSet(cur_gmo, gmoModelStat_OptimalLocal);
         havesol = 1;
         break;
      }

      case 2 : // infeasible (linear constraints and variable bounds)
      case 5 : // infeasible (inefficient (integer) variable bounds)
      case 16: // infeasible MIP (relaxation solution is provided)
      case 49: // Variable fixed to infeasible value
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleNoSolution);
         break;
      }

      case 6: // unbounded (linear constraints and variable bounds)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         gmoModelStatSet(cur_gmo, gmoModelStat_UnboundedNoSolution);
         break;
      }

      case 10: // IPM iteration limit exceeded
      case 40: // SQP iteration limit exceeded
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Iteration);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleIntermed);
         havesol = 1;
         break;
      }

      case 11: // infeasible with solution
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleLocal);
         havesol = 1;
         break;
      }

      case 13: // unbounded with solution
      case 29: // iteration diverged with solution
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         gmoModelStatSet(cur_gmo, gmoModelStat_Unbounded);
         havesol = 1;
         break;
      }

      case 17: // B & B node limit reached (with feas.sol.)
      case 37: // B & B terminated with given # of feas.sol. (solution limit)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_User);
         gmoModelStatSet(cur_gmo, gmoModelStat_Feasible);
         havesol = 1;
         break;
      }

      case 18: // MIP iteration failed (with feas.sol.)
      case 83: // Feas.sol found (numerical error in B & B)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Solver);
         gmoModelStatSet(cur_gmo, gmoModelStat_Feasible);
         havesol = 1;
         break;
      }

      case 19: // B & B node limit reached (no feas.sol.)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Iteration);
         gmoModelStatSet(cur_gmo, gmoModelStat_NoSolutionReturned);
         break;
      }

      case 20: // MIP iter. Failed (no feas.sol.)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Solver);
         gmoModelStatSet(cur_gmo, gmoModelStat_NoSolutionReturned);
         break;
      }

      case 21: // B & B itr. Timeout (with feas.sol.)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Resource);
         gmoModelStatSet(cur_gmo, gmoModelStat_Feasible);
         havesol = 1;
         break;
      }

      case 22: // B & B itr. Timeout (no feas.sol.)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Resource);
         gmoModelStatSet(cur_gmo, gmoModelStat_NoSolutionReturned);
         break;
      }

      case 27: // SIMPLEX iteration limit exceeded
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Iteration);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleIntermed);
         havesol = 1;
         break;
      }

      case 30: // terminated by user with solution
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_User);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleIntermed);
         havesol = 1;
         break;
      }

      case 31: // B & B terminated by user (with feas.sol.)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_User);
         gmoModelStatSet(cur_gmo, gmoModelStat_Feasible);
         havesol = 1;
         break;
      }

      case 32: // B & B terminated by user (no feas.sol.)
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_User);
         gmoModelStatSet(cur_gmo, gmoModelStat_NoSolutionReturned);
         break;
      }

      case 33: // Bound violated
      case 34: // Bound and Constraint violated
      case 35: // Constraint violated
      case 36: // Equality constraint violated
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Solver);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleIntermed);
         havesol = 1;
         break;
      }

      case 38: // dual infeasible with solution
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         gmoModelStatSet(cur_gmo, gmoModelStat_Unbounded);
         havesol = 1;
         break;
      }

      case 39: // IPM iteration timeout
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Resource);
         gmoModelStatSet(cur_gmo, gmoModelStat_InfeasibleIntermed);
         havesol = 1;
         break;
      }

      case 45: // B & B gap reaches under the limit
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_Normal);
         gmoModelStatSet(cur_gmo, gmoModelStat_Feasible);
         havesol = 1;
         break;
      }

      case 100: // Can not open NUOPT License file
      case 101: // Invalid License file
      case 102: // Machine key (YY) not consistent to License file
      case 103: // License expired on XX days ago
      case 104: // Invalid License limit
      case 105: // Can not get license information check machine configuration
      {
         gmoSolveStatSet(cur_gmo, gmoSolveStat_License);
         gmoModelStatSet(cur_gmo, gmoModelStat_LicenseError);
         break;
      }

      default:
      {
         gevLogStat(cur_gev, r->errorMessage());
         gmoSolveStatSet(cur_gmo, gmoSolveStat_InternalErr);
         gmoModelStatSet(cur_gmo, gmoModelStat_ErrorNoSolution);
         break;
      }
   }

   if( havesol )
   {
      // TODO duals, basis, etc
      gmoSetSolutionPrimal(cur_gmo, r->getX());
      gmoCompleteSolution(cur_gmo);
   }

   gmoSetHeadnTail(cur_gmo, gmoHresused, r->optTime());
   // TODO: return more like itercount, nodecount, dualbound

   delete r;

   return 0;
}

DllExport int STDCALL C__nuoHaveModifyProblem(void* Cptr)
{
   return 0;
}

DllExport int STDCALL C__nuoModifyProblem(void* Cptr)
{
   assert(Cptr == (void*)0x1);
   return 1;
}

}
