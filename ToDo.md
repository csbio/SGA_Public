SGA tasks
==========

Rectangularize
--------------
The code keeps query strains and array strains in a common name space. This means array strains
simply occupy a higher range of the "address space," which is harmless. However, near the end
of the script, the data are reorganized into a preallocated "square" space, which will be only 
about 1/4 full. Logic for this part of the script should be rewritten to avoid allocating 
(and subsequently saving) huge tracts of unused memory.

Condense moudules
-----------------
Move the following sections to external functions:
   * Array-means (batch block 1: prep)
   * batch correction (batch blocks 2,3)
   * triple-specific array removel (to filt_colsize)
   * model-fitting

Hard-coded Options
-------
Several scripts (e.g. clustergram generators) still have several modes which are 
selected by commenting code. These should be fixed with input parameters and default
options and /or prompts. The goal is to keep spurious changes out of the GIT log.

Green-Block removal
-------------------
Green-Block mitigation is now a stardard part of the pipeline. It should be moved inside the script
with an option to disable. 
