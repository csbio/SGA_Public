function[sga] = sga_concat(sga1, sga2, cat_type)
%function[sga] = sga_concat(sga1, sga2, cat_type)
% cat_type = 'horizontal' | horz | vertical | vert]
% members of 1 overwrite 2
% output struct will be abbr
% common dimension will be intersected

	if ~exist('cat_type', 'var')
		cat_type = 'horz';
		fprintf('default cat: side by side\n');
	end
	assert(ismember(cat_type, {'horizontal', 'horz', 'vertical', 'vert'}));


	if strcmp(cat_type(1:3), 'hor')
		% side by side, intersect the queries, union the arrays
		% remove common arrays from set2
		sga2.Cannon.isArray(ismember(sga2.Cannon.Orf, sga1.Cannon.Orf)) = false;

		Qorf1 = sga1.Cannon.Orf(sga1.Cannon.isQuery);
		Qorf2 = sga2.Cannon.Orf(sga2.Cannon.isQuery);
		Aorf1 = sga1.Cannon.Orf(sga1.Cannon.isArray);
		Aorf2 = sga2.Cannon.Orf(sga2.Cannon.isArray);
		[ComQ, ix1, ix2] = intersect(Qorf1, Qorf2);
		ComA = [Aorf1; Aorf2];

	else
		% top and bottom, intersect the arrays and union the queries
		% remove common queries from set2
		sga2.Cannon.isQuery(ismember(sga2.Cannon.Orf, sga1.Cannon.Orf)) = false;

		Qorf1 = sga1.Cannon.Orf(sga1.Cannon.isQuery);
		Qorf2 = sga2.Cannon.Orf(sga2.Cannon.isQuery);
		Aorf1 = sga1.Cannon.Orf(sga1.Cannon.isArray);
		Aorf2 = sga2.Cannon.Orf(sga2.Cannon.isArray);
		ComQ = [Qorf1; Qorf2];
		[ComA, iy1, iy2] = intersect(Aorf1, Aorf2);

	end

	sga = struct();

	%do the cannon
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
	fields = {'eps', 'pvl'};
	%fields = {'eps', 'pvl', 'dbl', 'escore', 'dbl_std'};
	for f = 1:length(fields)
		sga.(fields{f}) = NaN(sga.Cannon.GENES);
		E1 = sga1.(fields{f})(sga1.Cannon.isQuery, sga1.Cannon.isArray);
		E2 = sga2.(fields{f})(sga2.Cannon.isQuery, sga2.Cannon.isArray);
		if strcmp(cat_type(1:3), 'hor')
			sga.(fields{f})(sga.Cannon.isQuery, sga.Cannon.isArray) = [E1(ix1,:) E2(ix2,:)];
		else
			sga.(fields{f})(sga.Cannon.isQuery, sga.Cannon.isArray) = [E1(:,iy1); E2(:,iy2)];
		end
	end
end


% function [common, ixa, ixb] = internal_union(SetA, SetB)
% 	% asserts disjoint sets
% 	% returns ixa, ixb S.T. common = union(A,B) && common = [seta(ix); setb(ix);]

% 	assert(isempty(intersect(SetA, SetB)));

% 	common = set_A