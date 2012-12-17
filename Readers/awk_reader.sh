#! /bin/bash

cat h2o.evp | gawk '{printf"%d  %f \t \n", $1, $8}' > outfile
