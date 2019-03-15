# GAMS/NuOpt interface

Interface to allow using [NTT Data NuOpt](https://www.msi.co.jp/nuopt/)
as solver within [GAMS](https://www.gams.com).

## Getting Started

You need a GAMS system and NuOpt libraries including their AMPL interface source.

1. Ensure that a directory `gams` is available that contains the GAMS
   distribution.
2. Ensure that a directory `nuopt` is available that contains the NuOpt
   libraries, header files, and AMPL link sources.
3. Have a look into Makefile on what actually happens during the build.
4. Run `make`.
5. Edit `gmscmpun.txt` in the GAMS system directory to make the solver link
   known to GAMS. Add the following section:

       NUOPT 11 5 0001020304 1 0 2 LP RMIP MIP QCP MIQCP RMIQCP NLP DNLP RMINLP MINLP
       gmsgenus.run
       gmsgenux.out
       /path/to/libgamsnuopt.so nuo 1 1

References:
- [NuOpt Connection Manual](https://translate.google.com/translate?hl=en&sl=ja&u=http://www.msi.co.jp/nuopt/docs/v20/connection/&prev=search)
- [GAMS GMO API](https://www.gams.com/latest/docs/apis/expert-level/gmoqdrep.html)
