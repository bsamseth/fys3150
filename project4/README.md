# FYS3150 Prosjekt 4
Dette repoet inneholder alt materiale brukt av Lyder Rumohr Blingsmo og Bendik Samseth i prosjekt 4 i kurset FYS3150 ved UiO. 

## Rapport
Den fulle rapporten genereres ved å skrive:
``` bash
$ cd report/ && make
```
[Makefile-en](report/Makefile) bruker `pdflatex` og `biber` for å lage rapporten.


## Reproduksjon
Alle data og figurer som brukes i rapporten kan reproduseres ved å skrive

``` 
$ bash benchmark.sh [--no-compute]
```
Legg ved valgfritt argument `--no-compute` for å *ikke* generere datafiler på nytt. Dette vil bruke de filene som ligger i mappen fra før (som svarer til dataen bruk i rapporten). Alle datafiler kan reproduseres ved å kjøre skriptet uten ekstra argument. Merk at dette vil ta *lang* tid.

Dette burde reprodusere innholdet i [benchmarks/](benchmarks/).
