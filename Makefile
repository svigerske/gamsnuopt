all : gamsnuopt

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

gamsnuopt : src/main.o src/minlp.o src/minlp_f.o src/loadgms.o gmomcc.o gevmcc.o optcc.o $(DUMMY_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o src/*.o gamsnuopt

%.c : gams/apifiles/C/api/%.c
	cp $< $@

%.cc : nuopt/ampl/%.cc
	cp $< $@


FC=gfortran

IFLAGS = -Igams/apifiles/C/api -Inuopt/libnuopt -Inuopt/arch -Inuopt/ampl -Inuopt/dp -Inuopt/libsimple
#WFLAGS = -Wall -Wextra -Wno-unused-parameter
CFLAGS = $(IFLAGS) $(WFLAGS) -g -std=c99 -DGAMSDIR=\"gams\" -fPIC
CXXFLAGS = $(IFLAGS) $(WFLAGS) -g -std=c++11 -fPIC
FFLAGS = -fPIC

LDFLAGS = -ldl -Wl,-rpath,\$$ORIGIN -Wl,-rpath,$(realpath gams)
#LDFLAGS += -Lnuopt/userapp/lib -lnuopt_unix
#LDFLAGS += -Lnuopt/libnuopt -lnuopt_nolapack -llapack -lblas
LDFLAGS += -Lnuopt/libnuopt -lnuopt
LDFLAGS += -Lnuopt/f2c -lI77 -lF77
LDFLAGS += -Lnuopt/dp -ldp
