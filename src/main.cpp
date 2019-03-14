#include "gmomcc.h"
#include "gevmcc.h"

#include "nuopt_exception.h"

extern "C"
{
#include "loadgms.h"
}

#if 0
int solveMIQCP(
   gmoHandle_t gmo,
   gevHandle_t gev
);
#endif

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
      solveMINLP(gmo, gev);
   }
   catch( const NuoptException& e )
   {
      gevLogStat(gev, e.what());
   }

   freeGMS(&gmo, &gev);
}