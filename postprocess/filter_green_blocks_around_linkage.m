function [sga] = filter_green_blocks_around_linkage(sga, linkage_file, coord_file, layout_file)
%function [sga] = filter_green_blocks_around_linkage(sga, linkage_file, coord_file, layout_file)

	% convert my style struct to anastasia's style, then use her original code
	dataset = struct();
	dataset.queries = StripOrfs(sga.Cannon.Orf(sga.Cannon.isQuery));
	dataset.arrays  = StripOrfs(sga.Cannon.Orf(sga.Cannon.isArray));
	dataset.scores = sga.escore(sga.Cannon.isQuery,sga.Cannon.isArray);
	dataset = helper_dataset_wrapper(dataset, linkage_file, coord_file, layout_file);


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

function[dataset] = helper_dataset_wrapper(dataset, linkage_file, coord_file, layout_file)

	 inds_q = 1:length(dataset.queries);



	[qsl.orf, qsl.start, qsl.end] = ...
		 textread(linkage_file,'%s %d %d %d %d','delimiter','\t');
	for i = 1 : length(qsl.orf)
		 qsl.chr(i,1) = findstr(qsl.orf{i}(2),'ABCDEFGHIJKLMNOP');
	end

	% Ignore this when dealing with double mutants, not triples
	qsl.start2 = zeros(length(qsl.orf),1);
	qsl.end2 = zeros(length(qsl.orf),1);
	qsl.chr2 = zeros(length(qsl.orf),1);

	% Triples
	% [qsl.orf, qsl.chr, qsl.start, qsl.end, qsl.chr2, qsl.start2, qsl.end2] = ...
	%     textread('/Users/Anastasia/Laboratory/Datasets/Utils_SGA/Linkage_query_specific/query_specific_linkage_101101.txt', ...
	%     '%s %d %d %d %d %d %d','delimiter','\t');

	% [qsl.orf, qsl.chr, qsl.start, qsl.end, qsl.chr2, qsl.start2, qsl.end2] = ...
	%     textread('/Users/Anastasia/Laboratory/People/_Collab/Brandi_Mahaney/10-10-21_Scoring/query_specific_linkage.txt', ...
	%     '%s %d %d %d %d %d %d','delimiter','\t');




	[chrcoord.orf, chrcoord.chrom, chrcoord.start, chrcoord.end] = textread(coord_file,'%s %f %f %f');
	ind = find(chrcoord.start > chrcoord.end);
	tmp = chrcoord.start;
	chrcoord.start(ind) = chrcoord.end(ind);
	chrcoord.end(ind) = tmp(ind);
	clear tmp;

	chr1 = zeros(length(dataset.queries),1);
	chr2 = zeros(length(dataset.queries),1);
	st1 = zeros(length(dataset.queries),1);
	st2 = zeros(length(dataset.queries),1);
	en1 = zeros(length(dataset.queries),1);
	en2 = zeros(length(dataset.queries),1);

	for i = 1 : length(inds_q)
		 i
		 ind = strmatch(dataset.queries{inds_q(i)}, qsl.orf,'exact');
		 if(length(ind)>1), ind = ind(1); end
		 if ~isempty(ind)
			  chr1(inds_q(i),1) = qsl.chr(ind);
			  st1(inds_q(i),1) = qsl.start(ind);
			  en1(inds_q(i),1) = qsl.end(ind);
			  chr2(inds_q(i),1) = qsl.chr2(ind);
			  st2(inds_q(i),1) = qsl.start2(ind);
			  en2(inds_q(i),1) = qsl.end2(ind);
		 else
			  fields = split_by_delimiter('_', dataset.queries{inds_q(i)});
			  x = strmatch(fields{1}, chrcoord.orf, 'exact');
			  chr1(inds_q(i),1) = findstr(fields{1}(2),'ABCDEFGHIJKLMNOP');
			  st1(inds_q(i),1) = chrcoord.start(x)-200000;
			  en1(inds_q(i),1) = chrcoord.end(x)+200000;
			  
			  if length(fields)>1
					x = strmatch(fields{2},chrcoord.orf,'exact');
					if ~isempty(x)
						 chr2(inds_q(i),1) = findstr(fields{2}(2),'ABCDEFGHIJKLMNOP');
						 st2(inds_q(i),1) = chrcoord.start(x)-200000;
						 en2(inds_q(i),1) = chrcoord.end(x)+200000;
					end
			  end
		 end
	end


	[fgmap.plate, fgmap.row, fgmap.col, fgmap.orf] = textread(layout_file,'%f %f %f %s');

	fprintf(['|', blanks(100), '|\n']);
	fprintf('|');
	y = 0;

	for q = 1 : length(inds_q)
		 
		 x = fix(q * 100 / length(inds_q));
		 if x > y
			  fprintf('*');
			  y = x;
		 end
		 
		 % First ORF
		 linked_orfs = chrcoord.orf(find(chrcoord.chrom == chr1(inds_q(q)) & ...
			  ((chrcoord.start >= st1(inds_q(q)) & chrcoord.start <= en1(inds_q(q))) | ...
			  (chrcoord.end >= st1(inds_q(q)) & chrcoord.end <= en1(inds_q(q))))));

		 [t,ind1,ind2] = intersect(linked_orfs, fgmap.orf);

		 for a = 1 : length(ind1)

			  affected_orfs = fgmap.orf(find(fgmap.plate == fgmap.plate(ind2(a)) & fgmap.row >= fgmap.row(ind2(a))-1 & fgmap.row <= fgmap.row(ind2(a))+1 & ...
					fgmap.col >= fgmap.col(ind2(a))-1 & fgmap.col <= fgmap.col(ind2(a))+1));
			  [t,ind3,ind4] = intersect(dataset.arrays, affected_orfs);
			  ii = find(dataset.scores(inds_q(q),ind3) > 0);
			  dataset.scores(inds_q(q),ind3(ii)) = 0;

		 end
		 
		 % Second ORF
		 linked_orfs = chrcoord.orf(find(chrcoord.chrom == chr2(inds_q(q)) & ...
			  ((chrcoord.start >= st2(inds_q(q)) & chrcoord.start <= en2(inds_q(q))) | ...
			  (chrcoord.end >= st2(inds_q(q)) & chrcoord.end <= en2(inds_q(q))))));

		 [t,ind1,ind2] = intersect(linked_orfs, fgmap.orf);

		 for a = 1 : length(ind1)

			  affected_orfs = fgmap.orf(find(fgmap.plate == fgmap.plate(ind2(a)) & fgmap.row >= fgmap.row(ind2(a))-1 & fgmap.row <= fgmap.row(ind2(a))+1 & ...
					fgmap.col >= fgmap.col(ind2(a))-1 & fgmap.col <= fgmap.col(ind2(a))+1));
			  [t,ind3,ind4] = intersect(dataset.arrays, affected_orfs);
			  ii = find(dataset.scores(inds_q(q),ind3) > 0);
			  dataset.scores(inds_q(q),ind3(ii)) = 0;

		 end
			  
			  
	end

	fprintf('|\n');

	%% Removing HIS3 green blocks
	for q = 1 : length(inds_q)
			  
		 linked_orfs1 = chrcoord.orf(find(chrcoord.chrom == chr1(inds_q(q)) & ...
			  ((chrcoord.start >= st1(inds_q(q)) & chrcoord.start <= en1(inds_q(q))) | ...
			  (chrcoord.end >= st1(inds_q(q)) & chrcoord.end <= en1(inds_q(q))))));
		 
		 linked_orfs2 = chrcoord.orf(find(chrcoord.chrom == chr2(inds_q(q)) & ...
			  ((chrcoord.start >= st2(inds_q(q)) & chrcoord.start <= en2(inds_q(q))) | ...
			  (chrcoord.end >= st2(inds_q(q)) & chrcoord.end <= en2(inds_q(q))))));
			  
		 if ismember('YOR202W',linked_orfs1) | ismember('YOR202W', linked_orfs2)      
			  affected_orfs = fgmap.orf(find(fgmap.row == 2 | fgmap.row == 15 | fgmap.col == 2 | fgmap.col == 23));
			  [t,ind3,ind4] = intersect(dataset.arrays, affected_orfs);
			  ii = find(dataset.scores(inds_q(q),ind3) > 0);
			  dataset.scores(inds_q(q),ind3(ii)) = 0;
		 end
	end

	%% Sync all the matrices
	% MOVED OUT OF WRAPPER
	%{

	inds = find(dataset.scores == 0);
	dataset.scores(inds) = NaN;
	%dataset.scores_sd(inds) = NaN;
	dataset.scores_eps(inds) = NaN;
	dataset.scores_eps_sd(inds) = NaN;
	dataset.scores_pvalue(inds) = NaN;
	%}

end
