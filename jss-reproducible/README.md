Reproducing the analyses in the paper

There are several different levels of reproducibility.
- For the most basic level of reproducibility you can re-make the
  figures using the jss-figure-*.R files provided. Type "make clean"
  to remove the figure files, then type "make" to re-make the
  figures. On my computer it took about 2 minutes to re-make the
  figure files.
- For the next level of reproducibility you can type "rm *.rds" to
  remove the RDS files which store the results of the model fitting
  computations used to make those figures. Re-running
  jss.disk.memory.R, jss.variable.peaks.R, and
  jss.evaluations.R to create the correponding rds files should
  take several hours/days. They will download the 35GB UCI chipseq
  data set (https://github.com/tdhock/feature-learning-benchmark),
  then run the algo with various parameters on some of them.
- The final level of reproducibility consists of re-running the algo
  with a grid of penalty values in (log N, N) for each of the chipseq
  data sets, in order to re-do the (fp,fn,errors) computations in
  jss.bench.models.csv. I initially did this on a compute cluster, and
  it took 1.7 years of computation time. I provided the R script
  jss.bench.models.one.R which can re-do the computation of one line
  of the jss.bench.models.csv file. There are 115925 lines total in
  this file. The script should be invoked on the command line via
  "Rscript jss.bench.models.one.R N" where N is an integer between 1
  and 115925, which specifies which line of the sorted file to
  re-compute. (smaller N for smaller data/penalty values, and quicker
  computation time) The script will print the line in
  jss.bench.models.csv and the newly computed data, and then it will
  save the newly computed data to a 1-line csv file,
  jss.bench.models.one/N.csv. So the entire CSV file could be
  re-created by running the script for all N, and then concatenating
  the jss.bench.models.one/*.csv files. For example one run yields the
  following output (note there are insignificant differences in
  megabytes/seconds/mean.intervals, which is expected due to
  differences in hardware).

```
tdhock@recycled:~/projects/PeakSegFPOP-paper/jss-reproducible(master*)$ Rscript jss.bench.models.one.R 1
...
       data
1:   stored
2: computed
...
   bedGraph.lines  penalty megabytes seconds peaks   bases mean.pen.cost
1:           2173 277.5592  2.187500   0.867   391 1413146     0.1902719
2:           2173 277.5592  1.992188   0.396   391 1413146     0.1902719
   total.loss mean.intervals max.intervals fn fp errors
1:   160356.4       4.439254            11  0  4      4
2:   160356.4       4.439715            11  0  4      4
tdhock@recycled:~/projects/PeakSegFPOP-paper/jss-reproducible(master*)$ 
```

Note the following timings as a function of the row number:

```
> bench.models[, row := 1:.N];bench.models[seq(1, .N, l=10), list(row, minutes=seconds/60)]
       row   minutes
 1:      1  0.014450
 2:  12881  1.208983
 3:  25761  1.389383
 4:  38642  2.593617
 5:  51522  3.409533
 6:  64403  6.896133
 7:  77283  7.883450
 8:  90164 10.041050
 9: 103044 14.364517
10: 115925 52.412700
> 
```

If you want to re-do all the computations, the best way would be to
use a compute cluster. The jss.bench.models.several.R script can be
used to create jss.bench.models.several.csv, which assigns a job ID to
each row. For a job time limit of 24 hours there are 645 jobs. Each
job can be run via `Rscript jss.bench.models.several.one.R JOB_ID`
where `JOB_ID` is an integer from 1 to 645. Each runs several
rows/models of different sizes:

```
> rand.models[, list(
+   rows=.N,
+   min.minutes=min(minutes),
+   max.minutes=max(minutes),
+   total.minutes=sum(minutes)
+ ), by=list(job)]
     job rows min.minutes max.minutes total.minutes
  1:   1  182  0.01741667    56.53437      1429.750
  2:   2  195  0.02666667    50.84578      1434.791
  3:   3  176  0.02061667   111.52672      1454.279
  4:   4  163  0.01795000   109.11963      1438.732
  5:   5  181  0.11725000   104.61587      1432.241
 ---                                               
641: 641  192  0.05641667    81.40693      1433.808
642: 642  185  0.01698333    60.08692      1444.709
643: 643  197  0.14301667    45.55493      1440.653
644: 644  172  0.05743333   137.72985      1435.945
645: 645  170  0.01878333    81.14487      1326.871
> 
```

The jss-packages.R script should take care of installing all the
necessary R packages for you. 
