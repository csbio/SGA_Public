SGA
===

SGA Scoring script and related utilities



Additional directories (outside of GIT). These are customary, and I think not hardcoded anywhere.
Thus, you can change them as long as you set your parameters up accordingly. 
Sub-directories and symbolic links also work fine. 

  - `rawdata/`: for input files
  - `scored/`: for output files
  - `refdata/: for supplemental input files

  
Parameters:
  - Mandatory
    -  `inputfile` database dump to score (9 columns)
    -  `outputfile` output file pattern (no extension)
    -  `smfitnessfile` path to query / array fitness file
    -  `linkagefile` path to linkage definition file (ignored if linkage is skipped, see below)
    -  `coord_file` path to orf coordinate file (ignored if linkage is skipped)
    -  `removearraylist` path to list of known bad arrays (to be removed)
  
  -  Overrideable (have a hard-coded default, but may drop you into the debugger for confirmation if not supplied):
    -  `border_strain_orf` ORF_strainid found on the border
    -  `wild_type` ORF_strainid to use as wild-type
    -  `skip_perl_step` skips a preprocessing step. If you're scoring an `inputfile` you've already scored, you can skip this to save some time
    
  - Optional:
    - `skip_linkage_detection` set to `true` to skip linkage detection altogether
    - `skip_linkage_mask` set the to `true` to replace linkage colonies after corrections
    - `skip_wt_remove` set this to `true` to include WT query strain/s in the output
    - `eps_qnorm_ref` path to mat-file containing a quantile normalization table

  
  
Example scoring job (using the script):
You can paste this into matlab, (% are comments)

```
% -------------------------------------------------------------------------------------
% October 9, 2014 
% nolink score of fg30 for adrian verster
% in which linkage data are held out from corrections, then replaced for analysis
% -------------------------------------------------------------------------------------

    inputfile = 'rawdata/Collab/AdrianVerster/raw_sga_fg_t30_131130_adrian.txt';
    outputfile= 'scored/Collab/AdrianVerster/nolink_sga_fg_t30_131130_adrian_scored_141009';
    skip_perl_step = false;
    skip_linkage_mask = true;
    skip_wt_remove = false;
    wild_type = 'URA3control_sn4757';
    border_strain_orf = 'YOR202W_dma1';
    smfitnessfile = 'refdata/smf_t30_130417.txt';
    linkagefile = 'refdata/linkage-est_sga_merged_131122.txt';
    coord_file = 'refdata/chrom_coordinates_111220.txt';
    removearraylist = 'refdata/bad_array_strains_140526.csv';
    compute_sgascore
```




  
  
