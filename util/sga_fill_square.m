function [sga] = sga_fill_square(query_list, sga_fg, sga_ts)
%function [sga] = sga_fill_square(query_list, sga_fg, sga_ts)
% given a list of querys, fill out a square with available epsilon and pvalue data
% queries are matched exactly, arrays are mapped by known equivalence

	% load the equiv and keep only tsq
    equiv = Csv2Cell('~/SGA/refdata/array_query_map_equiv_good_only_121002.csv');

	sga = struct();
	sga.Cannon = struct();
	sga.Cannon.Orf = query_list;
	sga.Cannon.Common = OrfToCommon(query_list);
	sga.Cannon.GENES = length(query_list);
	sga.Cannon.isArray = logical(ones(1,sga.Cannon.GENES));
	sga.Cannon.isQuery = logical(ones(sga.Cannon.GENES,1));

	sga.eps = NaN(length(query_list));
	sga.pvl = NaN(length(query_list));
	sga.dbl = NaN(length(query_list));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fill in anything from TS first, will over write later... %%%%%%%%%%%%%%%%%
	% find the exact queries
	[Qint, qixa, qixb] = intersect(query_list, sga_ts.Cannon.Orf(sga_ts.Cannon.isQuery));

	% rename the arrays to queries and then exact match...
	sga_ts.Cannon.Orf(sga_ts.Cannon.isArray) = ...
		convert_array_to_query(sga_ts.Cannon.Orf(sga_ts.Cannon.isArray), equiv);

	[Aint, aixa, aixb] = intersect(query_list, sga_ts.Cannon.Orf(sga_ts.Cannon.isArray));

	EPS = sga_ts.eps(sga_ts.Cannon.isQuery, sga_ts.Cannon.isArray);
	PVL = sga_ts.pvl(sga_ts.Cannon.isQuery, sga_ts.Cannon.isArray);
	DBL = sga_ts.dbl(sga_ts.Cannon.isQuery, sga_ts.Cannon.isArray);

	%fill them in 
	sga.eps(qixa, aixa) = EPS(qixb, aixb);
	sga.pvl(qixa, aixa) = PVL(qixb, aixb);
	sga.dbl(qixa, aixa) = DBL(qixb, aixb);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Now fill in anything we can find from FG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% there will be some overlap so this can overwrite previously found values
	% find the exact queries
	[Qint, qixa, qixb] = intersect(query_list, sga_fg.Cannon.Orf(sga_fg.Cannon.isQuery));

	% rename the arrays to queries and then exact match...
	sga_fg.Cannon.Orf(sga_fg.Cannon.isArray) = ...
		convert_array_to_query(sga_fg.Cannon.Orf(sga_fg.Cannon.isArray), equiv);

	[Aint, aixa, aixb] = intersect(query_list, sga_fg.Cannon.Orf(sga_fg.Cannon.isArray));

	EPS = sga_fg.eps(sga_fg.Cannon.isQuery, sga_fg.Cannon.isArray);
	PVL = sga_fg.pvl(sga_fg.Cannon.isQuery, sga_fg.Cannon.isArray);
	DBL = sga_fg.dbl(sga_fg.Cannon.isQuery, sga_fg.Cannon.isArray);

	%fill them in 
	sga.eps(qixa, aixa) = EPS(qixb, aixb);
	sga.pvl(qixa, aixa) = PVL(qixb, aixb);
	sga.dbl(qixa, aixa) = DBL(qixb, aixb);

end


function [orfs] = convert_array_to_query(orfs, equiv)
	% if missing, they won't get used?

	[orfs, suffix] = StripOrfs(orfs);
	for i=1:length(orfs)
		ix = strmatch(suffix{i}, equiv(:,1), 'exact');
		if ~isempty(ix)
			orfs{i} = [orfs{i} '_' equiv{ix,2}];
		else
			orfs{i} = [orfs{i} '_' suffix{i}];
		end
	end
end
