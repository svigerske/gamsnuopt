#ifndef LOADGMS_H
#define LOADGMS_H

typedef enum
{
    RETURN_OK = EXIT_SUCCESS,
    RETURN_ERROR = EXIT_FAILURE
} RETURN;

extern
RETURN loadGMS(
   struct gmoRec** gmo,
   struct gevRec** gev,
   const char* gmsfile
);


extern
void freeGMS(
   struct gmoRec** gmo,
   struct gevRec** gev
);

#endif
