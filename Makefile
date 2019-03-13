all : gamsnuopt

gamsnuopt : src/main.o src/miqcp.o src/loadgms.o gmomcc.o gevmcc.o optcc.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o src/*.o gamsnuopt

%.c : gams/apifiles/C/api/%.c
	cp $< $@

IFLAGS = -Igams/apifiles/C/api -Inuopt/userapp/include
#WFLAGS = -Wall -Wextra -Wno-unused-parameter
CFLAGS = $(IFLAGS) $(WFLAGS) -g -std=c99 -DGAMSDIR=\"gams\"
CXXFLAGS = $(IFLAGS) $(WFLAGS) -g -std=c++11

LDFLAGS = -ldl -Wl,-rpath,\$$ORIGIN -Wl,-rpath,$(realpath gams)
LDFLAGS += -Lnuopt/userapp/lib -lnuopt_unix
#LDFLAGS += -Lnuopt/libnuopt -lnuopt
