# FYS3150 Prosjekt 3
Dette repoet inneholder alt materiale brukt av Lyder Rumohr Blingsmo og Bendik Samseth i prosjekt 3 i kurset FYS3150 ved UiO. 

## Rapport
Den fulle rapporten genereres ved å skrive:
``` bash
$ cd report/ && make
```
[Makefile-en](report/Makefile) bruker `pdflatex` og `biber` for å lage rapporten.


## Reproduksjon
Alle data og figurer som brukes i rapporten kan reproduseres ved å skrive

``` bash
$ ./benchmark.sh
```
Merk at dette forutsetter at `bash` finnes i `/bin/bash`. Dersom dette ikke er tilfellet må noe annet spesifiseres.

Dette burde reprodusere innholdet i [benchmarks/](benchmarks/).
