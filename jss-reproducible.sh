#!/bin/bash
set -o errexit
cp -r jss*.R jss*.rds jss-figure-more-likely-models jss-figure-label-error jss-reproducible
rm -f jss-reproducible/*~
##rm -rf jss-reproducible/jss-data jss-reproducible/jss-figure-label-error
tar czf jss-reproducible.tgz jss-reproducible
cp jss-reproducible.tgz /tmp
cd /tmp
rm -rf jss-reproducible
tar xf jss-reproducible.tgz
cd jss-reproducible
make
