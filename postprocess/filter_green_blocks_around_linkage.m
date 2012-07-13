function [sga] = filter_green_blocks_around_linkage(sga, linkage_file, coord_file, layout_file, wild_type, border_strain)
%function [sga] = filter_green_blocks_around_linkage(sga, linkage_file, coord_file, layout_file, wild_type)

	% convert my style struct to anastasia's style, then use her original code
	dataset = struct();
	dataset.queries = StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery));
	dataset.arrays  = StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray));
	dataset.scores = sga.escore(sga.Cannon.isQuery,sga.Cannon.isArray);
	dataset = helper_dataset_wrapper(dataset, linkage_file, coord_file, layout_file, wild_type, border_strain);


	% convert back - unpack
	sga.escore(sga.Cannon.isQuery,sga.Cannon.isArray) = dataset.scores;

	% do our own sync step
	removed = find(sga.escore == 0);
	sga.escore(removed) =  NaN;
	sga.eps(removed) =     NaN;
	sga.dbl(removed) =     NaN;
	sga.pvl(removed) =     NaN;
	sga.dbl_std(removed) = NaN;

end

function[dataset] = helper_dataset_wrapper(dataset, linkage_file, coord_file, layout_file, wild_type, border_strain)

   inds_q = 1:length(dataset.queries);

	l_fid = fopen(linkage_file, 'r');
	A = textscan(l_fid, '%s%d%d', 'Delimiter', '\t');
	fclose(l_fid);
	qsl = struct();
	qsl.orf=A{1};
	qsl.start=A{2};
	qsl.end=A{3};

	for i = 1 : length(qsl.orf)
		 qsl.chr(i,1) = findstr(qsl.orf{i}(2),'ABCDEFGHIJKLMNOP');
	end

	[chrcoord.orf, chrcoord.chrom, chrcoord.start, chrcoord.end] = textread(coord_file,'%s %f %f %f');
	ind = find(chrcoord.start > chrcoord.end);
	tmp = chrcoord.start;
	chrcoord.start(ind) = chrcoord.end(ind);
	chrcoord.end(ind) = tmp(ind);
	clear tmp;

	[fgmap.plate, fgmap.row, fgmap.col, fgmap.orf] = textread(layout_file,'%f %f %f %s');

	block_filter_count = 0;
	ARRAY = split_by_delim('_', layout_file);
	ARRAY = ARRAY{3};
	if(~ismember(ARRAY, {'fg', 'ts'}))
		error('unrecognized array format!');
	end
	if(strcmp(ARRAY, 'fg'))
		FG_ARRAY = true;
	else
		FG_ARRAY = false;
	end

	for q = 1 : length(inds_q)
		 
		orfs_this_query = split_by_delimiter('+', dataset.queries{inds_q(q)});
		for i=1:length(orfs_this_query)
			if(strcmp(orfs_this_query{i}, wild_type))
				continue
			else
				[chr, strt, endd] = GetLnkCoord(orfs_this_query{i}, qsl, chrcoord);
				linked_orfs = chrcoord.orf(find(chrcoord.chrom == chr & ...
			  		((chrcoord.start >= strt & chrcoord.start <= endd) | ...
			  		(chrcoord.end >= strt & chrcoord.end <= endd))));

				% only Do this part if working with the FG array
				
				if(FG_ARRAY)				
					[t,ind1,ind2] = intersect(linked_orfs, fgmap.orf);
					for a = 1 : length(ind2)
						affected_orfs = fgmap.orf(find(fgmap.plate == fgmap.plate(ind2(a)) & fgmap.row >= ...
								fgmap.row(ind2(a))-1 & fgmap.row <= fgmap.row(ind2(a))+1 & ...
								fgmap.col >= fgmap.col(ind2(a))-1 & fgmap.col <= fgmap.col(ind2(a))+1));
						[t,ind3,ind4] = intersect(dataset.arrays, affected_orfs);
						ii = find(dataset.scores(inds_q(q),ind3) > 0);
						dataset.scores(inds_q(q),ind3(ii)) = 0;
						block_filter_count = block_filter_count + length(ii);

					end
				end

				%% Removing CONTROL green blocks; do this for ALL arrays
				% when the border is linked, scan the next row and col
				if ismember(border_strain, linked_orfs)
					affected_orfs = fgmap.orf(find(fgmap.row == 2 | fgmap.row == 15 | fgmap.col == 2 | fgmap.col == 23));
					[t,ind3,ind4] = intersect(dataset.arrays, affected_orfs);
					ii = find(dataset.scores(inds_q(q),ind3) > 0);
					dataset.scores(inds_q(q),ind3(ii)) = 0;
					block_filter_count = block_filter_count + length(ii);
				end
			end
		end
			  
	end
	fprintf('green-block-filter removed scores: %d\n', block_filter_count)
end


function[chr, strt, endd] = GetLnkCoord(orf, qsl, chrcoord)

	% check for exact linkage match
	ind = strmatch(orf, qsl.orf, 'exact');
	if(~isempty(ind))
		chr  = qsl.chr(ind);
		strt = qsl.start(ind);
		endd = qsl.end(ind);
		return
	end

	% Strip any annotation, look for approx match
	ix = strfind(orf, '_');
	if(~isempty(ix))
		orf = orf(1:ix-1);
	end
	ind = strmatch(orf, qsl.orf, 'exact');
	if(~isempty(ind))
		chr  = qsl.chr(ind);
		strt = qsl.start(ind);
		endd = qsl.end(ind);
		return
	end
	
	% look for coordinate match
	ind = strmatch(orf, chrcoord.orf, 'exact');
	if(~isempty(ind))
		chr  = chrcoord.chrom(ind);
		strt = chrcoord.start(ind) - 200000;
		endd = chrcoord.end(ind) + 200000;
		return
	end

	fprintf('error getting linkage coords for %s\n', orf);
end
