#!/usr/bin/env bash

# removes all of the files in this directory except 
# those in the repo and input files.

files_to_remove=$(find . -type f \
   -not -name Sci_Batch40.txt \
   -not -name Sci_Batch40.txt.gz \
   -not -name 'Sci_Batch40_standard.*' \
   -not -name cleanup.sh \
   -not -name 'test_*')

echo
echo "You are about to remove:"
echo "$files_to_remove"
echo
while true; do
    read -p "Proceed?  " yn
    case $yn in
        [Yy]* ) rm $files_to_remove; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

