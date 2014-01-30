function[sga] = mask_collab_and_linkage(sga, bad_linkage_orfs, fix_linkage_orfs)
%function[sga] = mask_collab_and_linkage(sga, bad_linkage_orfs, fix_linkage_orfs)
% sga is expected to be either ts_merge or fg_merge
% Does several things:
% 1) unset isQuery for any strain in bad_linkage_orfs
% 2) extend linkage region for anything in fix_linkage_orfs 
%       (column 1 is affected orf, column 2 is neighbor to emulate)
% 3) unset some is Array for new bad arrays
% 4) put a comment in the struct


	% 1) -----------------------------------------------------------------------
	unset_ix = ismember(sga.Cannon.Orf, bad_linkage_orfs);
	fprintf('found %d / %d orfs to remove\n', sum(unset_ix), length(bad_linkage_orfs));
	sga.Cannon.isQuery(unset_ix) = false;


	% 2) -----------------------------------------------------------------------
	linkagefile = '~/SGA/Main/refdata/linkage-est_sga_merged_131122.txt';
	coord_file = '~/SGA/Main/refdata/chrom_coordinates_111220.txt';

	fid = fopen(linkagefile, 'r');
	linkage = textscan(fid, '%s%d%d', 'Delimiter', '\t', 'ReturnOnError', false);
	fclose(fid);
	
	% coords go opposite ways for C/W so we sort them
	fid = fopen(coord_file, 'r');
	coords = textscan(fid, '%s%d%d%d', 'Delimiter', '\t', 'ReturnOnError', false);
	fclose(fid);
	CD = sort([coords{3} coords{4}], 2);

	for i=1:size(fix_linkage_orfs, 1)
		% This round these were all SN strains...
		afflicted_ix = strmatch(fix_linkage_orfs{i,1}, sga.Cannon.Orf, 'exact');
		emulate_ix = strmatch(fix_linkage_orfs{i,2}, linkage{1});
		mask_region = [linkage{2}(emulate_ix) linkage{3}(emulate_ix)];

		% find any gene on this chromosome which begins or ends in this region
		this_chr = strmatch(['Y' fix_linkage_orfs{i,2}(2)], coords{1});
		begins_within = CD(:,1) > mask_region(1) & CD(:,1) < mask_region(2);
		ends_within   = CD(:,2) > mask_region(1) & CD(:,2) < mask_region(2);

		arrays_to_mask = coords{1}(intersect(this_chr, find(begins_within | ends_within)));
		array_ix = sga.Cannon.isArray & ismember(StripOrfs(sga.Cannon.Orf), arrays_to_mask)';
	
		sga.eps(afflicted_ix, array_ix) = NaN;
		sga.pvl(afflicted_ix, array_ix) = NaN;
		sga.dbl(afflicted_ix, array_ix) = NaN;
		sga.dbl_std(afflicted_ix, array_ix) = NaN;
		
	end
		
	% 3) -----------------------------------------------------------------------
	removearraylist = '~/SGA/Main/refdata/bad_array_strains_130108.csv';
	fid = fopen(removearraylist, 'r');
	bad_arrays = textscan(fid, '%s');
	fclose(fid);
	array_ix = ismember(sga.Cannon.Orf', bad_arrays{1}) & sga.Cannon.isArray;
	fprintf('Removing %d additional bad arrays\n', sum(array_ix));
	sga.Cannon.isArray(array_ix) = false;

	% 4) -----------------------------------------------------------------------
	sga.comment = ['post processing linkage fixes ' date()];

