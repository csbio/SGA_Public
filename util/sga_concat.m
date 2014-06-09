function[sga] = sga_concat(sga1, sga2, cat_type, inter)
%function[sga] = sga_concat(sga1, sga2, cat_type, inter)
% cat_type = 'horizontal' | horz | vertical | vert]
% members of 1 overwrite 2
% output struct will be abbr
% common dimension will be intersected if inter==true (default)
% if set to false will be unioned, potentially leaving a big whole
% respects isQuery and isArray


	if ~exist('cat_type', 'var')
		cat_type = 'horz';
		fprintf('default cat: side by side\n');
	end
	assert(ismember(cat_type, {'horizontal', 'horz', 'vertical', 'vert'}));

	if ~exist('inter', 'var')
		inter = true;
		fprintf('default intersect: common dimension restricted\n');
	end
	assert(islogical(inter))

	% HORIZONTAL CAT - SIDE BY SIDE
	if strcmp(cat_type(1:3), 'hor')
		% side by side, intersect the queries, union the arrays
		% remove common arrays from set2
		sga2.Cannon.isArray(ismember(sga2.Cannon.Orf, sga1.Cannon.Orf)) = false;

		Qorf1 = sga1.Cannon.Orf(sga1.Cannon.isQuery);
		Qorf2 = sga2.Cannon.Orf(sga2.Cannon.isQuery);
		Aorf1 = sga1.Cannon.Orf(sga1.Cannon.isArray);
		Aorf2 = sga2.Cannon.Orf(sga2.Cannon.isArray);
		if(inter)
			[ComQ, ix1, ix2] = intersect(Qorf1, Qorf2);
		else
			ComQ = union(Qorf1, Qorf2);
			map1 = Hash([], Qorf1);
			map2 = Hash([], Qorf2);
			ix1 = apply_map(map1, ComQ); % leaves 0's
			ix2 = apply_map(map2, ComQ);
		end
		ComA = [Aorf1; Aorf2];

	% VERTICAL CAT - TOP AND BOTTOM
	else
		% top and bottom, intersect the arrays and union the queries
		% remove common queries from set2
		sga2.Cannon.isQuery(ismember(sga2.Cannon.Orf, sga1.Cannon.Orf)) = false;

		Qorf1 = sga1.Cannon.Orf(sga1.Cannon.isQuery);
		Qorf2 = sga2.Cannon.Orf(sga2.Cannon.isQuery);
		Aorf1 = sga1.Cannon.Orf(sga1.Cannon.isArray);
		Aorf2 = sga2.Cannon.Orf(sga2.Cannon.isArray);
		ComQ = [Qorf1; Qorf2];
		if(inter)
			[ComA, iy1, iy2] = intersect(Aorf1, Aorf2);
		else
			ComA = union(Aorf1, Aorf2);
			map1 = Hash([], Aorf1);
			map2 = Hash([], Aorf2);
			iy1 = apply_map(map1, ComA); % leaves 0's
			iy2 = apply_map(map2, ComA);
		end

	end


	sga = struct();
	sga.Cannon = struct();
	sga.Cannon.Orf = [ComQ; ComA];
	sga.Cannon.Common = OrfToCommon(sga.Cannon.Orf);
	sga.Cannon.Map = Hash([], sga.Cannon.Orf);
	sga.Cannon.GENES = length(sga.Cannon.Orf);
	sga.Cannon.isQuery = logical(zeros(sga.Cannon.GENES,1));
	sga.Cannon.isQuery(1:length(ComQ)) = true;
	sga.Cannon.isArray = logical(zeros(1,sga.Cannon.GENES));
	sga.Cannon.isArray(length(ComQ)+1:end) = true;

	% do the data
	% get all available fields (confined to this list)
	fields = {'eps', 'pvl', 'dbl', 'escore', 'dbl_std'};
	fields1 = fieldnames(sga1);
	fields2 = fieldnames(sga2);
	fields = intersect(fields, intersect(fields1, fields2));

	for f = 1:length(fields)
		sga.(fields{f}) = NaN(sga.Cannon.GENES);
		E1 = sga1.(fields{f})(sga1.Cannon.isQuery, sga1.Cannon.isArray);
		E2 = sga2.(fields{f})(sga2.Cannon.isQuery, sga2.Cannon.isArray);

		% add a "not found" row/col of NaNs and use its index in place of union-misses
		% allows us to recycle the intersect code
		if strcmp(cat_type(1:3), 'hor')
			E1 = [E1; NaN(1, size(E1,2))]; % add a nan row
			E2 = [E2; NaN(1, size(E2,2))];
			ix1(ix1 == 0) = size(E1,1);
			ix2(ix2 == 0) = size(E2,1);
			sga.(fields{f})(sga.Cannon.isQuery, sga.Cannon.isArray) = [E1(ix1,:) E2(ix2,:)];
		else
			E1 = [E1 NaN(size(E1,1),1)]; % add a nan col
			E2 = [E2 NaN(size(E2,1),1)];
			iy1(iy1 == 0) = size(E1,2);
			iy2(iy2 == 0) = size(E2,2);
			sga.(fields{f})(sga.Cannon.isQuery, sga.Cannon.isArray) = [E1(:,iy1); E2(:,iy2)];
		end
	end


end
