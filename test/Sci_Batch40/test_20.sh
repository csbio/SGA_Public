#!/usr/bin/env bash
# checks the strain manifests for differences

# I don't remember exactly what determines the order of strains in the final
# output. So we will sort them just in case before calling diff.
sort Sci_Batch40_standard.orf > tmp_left
sort Sci_Batch40_output.orf   > tmp_right

echo "--------------------------------------------------------------"
diff --side-by-side --width=80 --suppress-common-lines tmp_left tmp_right

if [ $? -eq 0 ]
then
   echo "Strain check passed."
   rm tmp_left tmp_right
   exit 0
else
   echo "--------------------------------------------------------------"
   echo "< Sci_Batch40_standard.orf            > Sci_Batch40_output.orf"
   echo
   echo "Strain check failed."
   rm tmp_left tmp_right
   exit 1
fi


