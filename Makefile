FC=gfortran

IFLAGS = -Igams/apifiles/C/api -Inuopt/libnuopt -Inuopt/arch -Inuopt/ampl -Inuopt/dp -Inuopt/libsimple
WFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-ignored-qualifiers
CFLAGS = $(IFLAGS) -g -std=c99 -fPIC
CXXFLAGS = $(IFLAGS) -g -std=c++11 -fPIC
FFLAGS = -fPIC

LDFLAGS = -shared
LDFLAGS += -ldl -Wl,-rpath,\$$ORIGIN
#LDFLAGS += -Lnuopt/userapp/lib -lnuopt_unix
#LDFLAGS += -Lnuopt/libnuopt -lnuopt_nolapack -llapack -lblas
LDFLAGS += -Lnuopt/libnuopt -lnuopt
LDFLAGS += -Lnuopt/f2c -lI77 -lF77
LDFLAGS += -Lnuopt/dp -ldp

DUMMY_OBJECTS = \
  myAD.o \
  myArgmmEF.o \
  myCountEF.o \
  myEF.o \
  myMaxNEF.o \
  myMinNEF.o \
  myPGVar.o \
  myProdEF.o \
  mySumEF.o \
  mycpnuopt.o \
  mywcspSolve.o \
  myrcpspSolve.o \
  myrcpspnuopt.o \
  myoptseqSolve.o \
  mysimpleenv.o \
  mymtxfree.o

OBJS = $(addprefix obj/, gamsnuopt.o callbackwrap.o gmomcc.o gevmcc.o $(DUMMY_OBJECTS))

libgamsnuopt.so : $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf obj
	rm -f libgamsnuopt.so

obj :
	mkdir -p obj

obj/%.o : src/%.cpp | obj
	$(CXX) -o $@ -c $^ $(CXXFLAGS) $(WFLAGS)

obj/%.o : src/%.f | obj
	$(FC) -o $@ -c $^ $(FFLAGS)

obj/%.o : gams/apifiles/C/api/%.c | obj
	$(CC) -o $@ -c $^ $(CFLAGS)

obj/%.o : nuopt/ampl/%.cc | obj
	$(CXX) -o $@ -c $^ $(CXXFLAGS)
