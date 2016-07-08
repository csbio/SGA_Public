#!/usr/bin/env bash

# checks the entire score file for differences
diff -q Sci_Batch40_standard.txt Sci_Batch40_output.txt

if [ $? -eq 0 ]
then
   echo "Global check passed. Results are identical."
   exit 0
else
   echo "Global check failed. Run further tests."
   exit 1
fi


