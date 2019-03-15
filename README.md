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

       NUOPT 11 5 0001020304 1 0 2 LP RMIP MIP QCP RMIQCP NLP DNLP RMINLP
       gmsgenus.run
       gmsgenux.out
       /path/to/libgamsnuopt.so nuo 1 1

## Usage

Use GAMS option `SOLVER=NUOPT` to call NuOpt when installed.

Use GAMS option `optfile=1` to read `nuopt.opt`, which is expected to be
in [NuOpt parameter file format](https://translate.googleusercontent.com/translate_c?depth=1&hl=en&prev=search&rurl=translate.google.com&sl=ja&sp=nmt4&u=http://www.msi.co.jp/nuopt/docs/v20/manual/html/15-01-00.html&xid=17259,15700023,15700186,15700191,15700248,15700253&usg=ALkJrhhguuoumGmmZv9YKM83WaWj1YKcxg), e.g.,
```
begin
method: trust
scaling: on
end
```

## TODO

- implement special treatment for quadratics
- redirect NuOpt output
- pass on GAMS options

## References

- [NuOpt Connection Manual](https://translate.google.com/translate?hl=en&sl=ja&u=http://www.msi.co.jp/nuopt/docs/v20/connection/&prev=search)
- [GAMS GMO API](https://www.gams.com/latest/docs/apis/expert-level/gmoqdrep.html)
