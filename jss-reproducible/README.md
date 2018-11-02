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
necessary R packages for you. You need to install R and Berkeley DB
STL Development libraries (which is required to install PeakSegDisk --
it is used to write a large temporary file to disk).

Installing BerkeleyDB STL

If you see the following, you need to install BerkeleyDB STL headers, or tell R where to find them.

Ubuntu/g++
```
* installing *source* package ‘PeakSegPipeline’ ...
** libs
g++ -std=c++0x -I/home/thocking/lib64/R/include -DNDEBUG  -I/usr/local/include    -fpic  -g -O2 -c PeakSegFPOPLog.cpp -o PeakSegFPOPLog.o
PeakSegFPOPLog.cpp:10:26: error: dbstl_vector.h: No such file or directory
PeakSegFPOPLog.cpp: In function ‘int PeakSegFPOP_disk(char*, char*)’:
```

MacOS/clang++
```
* installing to library '/Library/Frameworks/R.framework/Versions/3.2/Resources/library'
* installing *source* package 'PeakSegPipeline' ...
** libs
clang++ -std=c++11 -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include -I/usr/local/include/freetype2 -I/opt/X11/include    -fPIC  -Wall -mtune=core2 -g -O2 -c PeakSegFPOPLog.cpp -o PeakSegFPOPLog.o
PeakSegFPOPLog.cpp:10:10: fatal error: 'dbstl_vector.h' file not found
#include <dbstl_vector.h>
         ^
1 error generated.
make: *** [PeakSegFPOPLog.o] Error 1
ERROR: compilation failed for package 'PeakSegPipeline'
```

If you install BerkeleyDB STL to a standard directory, e.g. /usr/local, then you should be able to redo the package installation/compilation and everything should work fine. The easiest way to do this is with your system package manager, for example on Ubuntu the following commands install a C++ compiler and the BerkeleyDB STL development library.

```
sudo aptitude install build-essential libdb6.0++-dev libdb6.0-stl-dev

```

However if you don't have root, you may need to install BerkeleyDB STL to a non-standard directory, e.g. your home directory:
```
curl -OL http://download.oracle.com/berkeley-db/db-6.2.23.NC.tar.gz
tar xf db-6.2.23.NC.tar.gz
cd db-6.2.23.NC/build_unix
../dist/configure --prefix=$HOME --enable-stl
make
make install
```

In that case, you need to tell R where to find BerkeleyDB STL. The [R Installation and Administration Manual, section Customizing Package Compilation](https://cran.r-project.org/doc/manuals/R-admin.html#Customizing-package-compilation) explains that you can do this by putting the following code in your ~/.R/Makevars file.

MacOS/clang++
```
CPPFLAGS=-I${HOME}/include
LDFLAGS=-L${HOME}/lib
```

Ubuntu/g++
```
CPPFLAGS=-I${HOME}/include
LDFLAGS=-L${HOME}/lib -Wl,-rpath=${HOME}/lib
```

After you do that, you should see the following when you install PeakSegPipeline.

MacOS/clang++
```
* installing to library '/Library/Frameworks/R.framework/Versions/3.2/Resources/library'
* installing *source* package 'PeakSegPipeline' ...
** libs
clang++ -std=c++11 -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/Users/wireless/include    -fPIC  -Wall -mtune=core2 -g -O2 -c PeakSegFPOPLog.cpp -o PeakSegFPOPLog.o
clang++ -std=c++11 -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/Users/wireless/include    -fPIC  -Wall -mtune=core2 -g -O2 -c funPieceListLog.cpp -o funPieceListLog.o
clang++ -std=c++11 -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/Users/wireless/include    -fPIC  -Wall -mtune=core2 -g -O2 -c interface.cpp -o interface.o
clang++ -std=c++11 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/Users/wireless/lib -o PeakSegPipeline.so PeakSegFPOPLog.o funPieceListLog.o interface.o -ldb_stl -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
installing to /Library/Frameworks/R.framework/Versions/3.2/Resources/library/PeakSegPipeline/libs
** R
** preparing package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded
* DONE (PeakSegPipeline)
```

Ubuntu/g++
```
* installing to library ‘/gs/home/thocking/lib64/R/library’
* installing *source* package ‘PeakSegPipeline’ ...
** libs
g++ -std=c++0x -I/home/thocking/lib64/R/include -DNDEBUG  "-I/home/thocking/include"    -fpic  -g -O2 -c PeakSegFPOPLog.cpp -o PeakSegFPOPLog.o
g++ -std=c++0x -I/home/thocking/lib64/R/include -DNDEBUG  "-I/home/thocking/include"    -fpic  -g -O2 -c funPieceListLog.cpp -o funPieceListLog.o
g++ -std=c++0x -I/home/thocking/lib64/R/include -DNDEBUG  "-I/home/thocking/include"    -fpic  -g -O2 -c interface.cpp -o interface.o
g++ -std=c++0x -shared -L/home/thocking/lib -Wl,-rpath=/home/thocking/lib -o PeakSegPipeline.so PeakSegFPOPLog.o funPieceListLog.o interface.o -ldb_stl
installing to /gs/home/thocking/lib64/R/library/PeakSegPipeline/libs
** R
** preparing package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded
* DONE (PeakSegPipeline)
```

