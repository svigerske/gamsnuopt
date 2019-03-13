#include "gmomcc.h"
#include "gevmcc.h"

#include "nuopt_exception.h"

extern "C"
{
#include "loadgms.h"
}

int solveMIQCP(
   gmoHandle_t gmo,
   gevHandle_t gev
);

int solveMINLP(
   gmoHandle_t gmo,
   gevHandle_t gev
);

int main(int argc, char** argv)
{
   gmoHandle_t gmo;
   gevHandle_t gev;

//   nuoptDispVersion();

   loadGMS(&gmo, &gev, argv[1]);

   gmoIndexBaseSet(gmo, 1);

   try
   {
      if( gmoModelType(gmo) <= gmoProc_rmip || gmoModelType(gmo) >= gmoProc_qcp )
         solveMIQCP(gmo, gev);
      else
         solveMINLP(gmo, gev);
   }
   catch( const NuoptException& e )
   {
      gevLogStat(gev, e.what());
   }

   freeGMS(&gmo, &gev);
}