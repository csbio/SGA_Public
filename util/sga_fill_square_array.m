function [sga] = sga_fill_square_array(array_list, sga_fg, sga_ts)
%function [sga] = sga_fill_square_array(array_list, sga_fg, sga_ts)
% given a list of arrays, fill out a square with available epsilon and pvalue data
% arrays are matched exactly, queries are mapped by known equivalence
% Order in is order out?

	% load the equiv and keep only tsq
    % equiv = Csv2Cell('~/SGA/refdata/array_query_map_equiv_good_only_140105.csv');
    % equiv = Csv2Cell('/heap/data/stage/array_query_map_equiv_good_only_141106.csv');
    equiv = Csv2Cell('~/stage/array_query_map_equiv_good_only_141106.csv');
	query_equiv = convert_array_to_query(array_list, equiv);

	sga = struct();
	sga.Cannon = struct();
	sga.Cannon.Orf = array_list;
	sga.Cannon.Common = OrfToCommon(array_list);
	sga.Cannon.GENES = length(array_list);
	sga.Cannon.isArray = logical(ones(1,sga.Cannon.GENES));
	sga.Cannon.isArray = logical(ones(sga.Cannon.GENES,1));

	% unlike the usual case, these will collide, which should work well...
	sga.Cannon.Map = Hash([], array_list);
	sga.Cannon.Map = Hash(sga.Cannon.Map, query_equiv);

	sga.eps = NaN(length(array_list));
	sga.pvl = NaN(length(array_list));
	sga.dbl = NaN(length(array_list));

	% allocate temp versions, all the same size
	EPS_ts = sga.eps;
	PVL_ts = sga.eps;
	DBL_ts = sga.eps;
	EPS_fg = sga.eps;
	PVL_fg = sga.eps;
	DBL_fg = sga.eps;	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fill in anything from TS first, will over write later... %%%%%%%%%%%%%%%%%
	% find the exact arrays. They should be in order and we need to mark which ones
	% are found vs not, then remove misses so we can use as an index.
	a_id = apply_map(sga_ts.Cannon.Map, array_list);
	has_a_id = a_id > 0;
	a_id = nonzeros(a_id);

	% convert these arrays to query equiv, and do a similar exact lookup
	q_id = apply_map(sga_ts.Cannon.Map, query_equiv);
	has_q_id = q_id > 0;
	q_id = nonzeros(q_id);

	%fill them in 
	EPS_ts(has_q_id, has_a_id) = sga_ts.eps(q_id, a_id);
	PVL_ts(has_q_id, has_a_id) = sga_ts.pvl(q_id, a_id);
	DBL_ts(has_q_id, has_a_id) = sga_ts.dbl(q_id, a_id);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Now fill in anything we can find from FG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% there will be some overlap so this can overwrite previously found values

	a_id = apply_map(sga_fg.Cannon.Map, array_list);
	has_a_id = a_id > 0;
	a_id = nonzeros(a_id);

	% convert these arrays to query equiv, and do a similar exact lookup
	q_id = apply_map(sga_fg.Cannon.Map, query_equiv);
	has_q_id = q_id > 0;
	q_id = nonzeros(q_id);

	%fill them in 
	EPS_fg(has_q_id, has_a_id) = sga_fg.eps(q_id, a_id);
	PVL_fg(has_q_id, has_a_id) = sga_fg.pvl(q_id, a_id);
	DBL_fg(has_q_id, has_a_id) = sga_fg.dbl(q_id, a_id);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The final square will combine valid values, with a preference for FG
	sga.eps(~isnan(EPS_ts)) = EPS_ts(~isnan(EPS_ts));
	sga.pvl(~isnan(PVL_ts)) = PVL_ts(~isnan(PVL_ts));
	sga.dbl(~isnan(DBL_ts)) = DBL_ts(~isnan(DBL_ts));

	sga.eps(~isnan(EPS_fg)) = EPS_fg(~isnan(EPS_fg));
	sga.pvl(~isnan(PVL_fg)) = PVL_fg(~isnan(PVL_fg));
	sga.dbl(~isnan(DBL_fg)) = DBL_fg(~isnan(DBL_fg));
end

function [orfs] = convert_array_to_query(orfs, equiv)
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
