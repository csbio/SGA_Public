Data Description
================

40 randomly selected batches from the 2010 dataset
(supplementary data file 4)
With \_strainIDs added later

Sci_40.txt.gz summary:
type   sn   tsq   damp   dma   tsa   y   unann   total
query  122  8     7      0      0    0   0       137
array   0   0     0      4293   0    0   0       4293


To test the pipeline:
=====================
(I think this requires about <3 GB of memory, and runs end to end in 20 minutes
on my laptop.)

prep the input data:
```
# gzip file under revision control, don't gunzip in place...
zcat Sci_Batch40.txt.gz > Sci_Batch40.txt
```

run the job:
```
% MATLAB
addpath('test');
Sci_Batch40
```

Compare the output to established data:
=======================================
Tests for Sci_Batch40 should compare a newly generated `_output` version with
the established `_standard` version in the repo.

Test files should be called `test_##.ext` where ## is a number, and ext can be
whatever (.sh, .m, .py etc). They go from global to specific. So if you pass
the first test, you can stop, but if you fail, you can run more tests to find
the problem. 

The repo will ignore files that begin with `test` so if you want to spit out a
pdf or a log, you should name it `result` so that git will ignore it. Assume they
are run in order, so you can write state to result files or the matlab environment
for use in subsequent tests.

```
test_10.sh diffs the txt file to check all scores
test_20.sh diffs the orf file to check for strain census changes
```

*You may also want to examine a diff of the log files, but as file paths are subject to change
without affecting scores, this is tricky to do automatically.*
 
