# FYS3150 Prosjekt 2
Dette repoet inneholder alt materiale brukt av Lyder Rumohr Blingsmo og Bendik Samseth i prosjekt 2 i kurset FYS3150 ved UiO. 

Den fulle rapporten genereres ved å skrive:
``` bash
$ pdflatex report.tex
$ biber report
$ pdflatex report.tex
```

Alle data og figurer som brukes i rapporten kan reproduseres ved å skrive

``` bash
$ ./benchmark.sh
```
Merk at dette forutsetter at `bash` finnes i `/bin/bash`. Dersom dette ikk er tilfellet må noe annet spesifiseres.

Dette vil reprodusere innholdet i [benchmarks/](benchmarks/).
