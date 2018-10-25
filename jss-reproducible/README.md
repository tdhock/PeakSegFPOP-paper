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
  data sets. I initially did this on a compute cluster and saved the
  results to target.intervals.models.csv. However I don't provide a
  reproducibility script, because re-doing that using one processor
  would probably take several months/years.

The script should take care of installing all the necessary R packages
for you. You need to install R and Berkeley DB STL Development
libraries (which is required to install PeakSegDisk -- it is used to
write a large temporary file to disk).

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

