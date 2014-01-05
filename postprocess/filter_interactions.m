function[filt, fitness_struct] = filter_interactions(sga, fitness_file, sga_inputfile, array_query_eqiv_file, skip_filter)
%function[filt, fitness_struct] = filter_interactions(sga, fitness_file, sga_inputfile, array_query_equiv_file)

	% load the fitness data
	smf_fid = fopen(fitness_file, 'r');
	A = textscan(smf_fid, '%s%f%f', 'Delimiter', '\t', 'ReturnOnError', false);
	fitness_struct = [A{1} num2cell(A{2}) num2cell(A{3})];
	fclose(smf_fid);

	% load the equivalence file
	equiv_fid = fopen(array_query_eqiv_file, 'r');
	equiv = textscan(equiv_fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
	fclose(equiv_fid);
	equiv = lower([equiv{1} equiv{2}]);

	% load / calculate the cobatch_scores
	cobatch_scores = calculate_cobatch_by_query(sga, sga_inputfile);


	% remove _collab screens
	% and trigenic screens
	sga.Cannon.isQuery(substrmatch('_collab', sga.Cannon.Orf))=false;
	sga.Cannon.isQuery(substrmatch('+',       sga.Cannon.Orf))=false;
	

	% filter 
	filt = help_filter_by_cobatch(sga, cobatch_scores);
	filt = help_filter_by_Agreement(filt, equiv);

	% this step is optional

	if ~exist('skip_filter', 'var')
		skip_filter = false;
	end

	if(~skip_filter)
		filt = help_filter_by_ABBA(filt, fitness_struct, equiv);
	end

end

function[filt] = help_filter_by_cobatch(sga, cobatch_scores)
% sets scores to NaN and removes from query index

	filt = sga;
	cobatch_target = filt.Cannon.isQuery;
	for i=find(cobatch_target')
		ix = strmatch(filt.Cannon.Orf{i}, cobatch_scores(:,1), 'exact');
		if(cobatch_scores{ix,2} > 0.2)
			cobatch_target(i) = true;
			filt.eps(i,:) = NaN;
			filt.pvl(i,:) = NaN;
		else
			cobatch_target(i) = false;
		end
	end
	
	fprintf('%d queries targeted for CoBatch filtering\n', sum(cobatch_target));
	filt.cobatch_target = cobatch_target;
	filt.Cannon.isQuery(cobatch_target) = false;

end

function[sga] = help_filter_by_Agreement(sga, equiv)
% remove any ABBA that is significant in 2 directions:

	ABBA_disagree = 0;

	Cannon_strains = sga.Cannon.Orf;
	for i=1:length(Cannon_strains)
		ix = strfind(Cannon_strains{i}, '_');
		Cannon_strains{i} = Cannon_strains{i}(ix+1:end);
	end

	for i=find(sga.Cannon.isQuery')
		equiv_array = help_equiv_set(sga.Cannon.Orf{i}, Cannon_strains, equiv);
		if(~isempty(equiv_array))
			% NEG -- they're opposite so if there's a positive here we can skip it , we'll catch the negative on the flip side
			for j=find(sga.eps(i,:) < -0.08 & sga.pvl(i,:)<0.05)
				equiv_query = help_equiv_set(sga.Cannon.Orf{j}, Cannon_strains, equiv);
				if(~isempty(equiv_query))
					subnet = sga.eps(equiv_query, equiv_array) > 0.08 & sga.pvl(equiv_query, equiv_array)<0.05; 
					if(sum(sum(subnet)) > 0)
						sga.eps(i,j) = 0;
						sga.pvl(i,j) = 1;
						sga.pvl(equiv_query, equiv_array) = 1;
						ABBA_disagree = ABBA_disagree +1;
					end
				end
			end
		end
	end

	fprintf('ABBA_disagree (x2): %d\n', ABBA_disagree);
end

function[filt] = help_filter_by_ABBA(sga, fitness_struct, equiv)
	% tosses intermediate negatives for high fitness queries
	% unless rescued by ABBA

	% abs(eps) > 0.12 & pvl || AB & BA

	% isolate all intermediate interactions in query
	% fitness range. Target for removal, check in a loop...
	
	% always read from sga and write to filt
	filt = sga;
	filt.tossed = sparse(zeros(size(sga.eps)));
	fitness_target = sga.Cannon.isQuery;
	ABBA_rescue = 0;
	ABBA_refuse = 0;
	ABBA_unavail= 0;
	for i=find(fitness_target')
		ix = strmatch(sga.Cannon.Orf{i}, fitness_struct(:,1), 'exact');
		if(length(ix)>1)
			fprintf('warning: multiple fitness values for %s. Taking first\n', sga.Cannon.Orf{i});
			ix = ix(1);
		end
		if(~isempty(ix) && fitness_struct{ix,2}>0.97) % target aquired
			fitness_target(i) = true;
		else
			fitness_target(i) = false;
		end
	end

	Cannon_strains = sga.Cannon.Orf;
	for i=1:length(Cannon_strains)
		ix = strfind(Cannon_strains{i}, '_');
		Cannon_strains{i} = Cannon_strains{i}(ix+1:end);
	end

	% for the sick queries, we need to isolate the interactions
	
	fprintf('%d queries targeted for ABBA filtering\n', sum(fitness_target));
	filt.fitness_target = fitness_target;
	for i=find(fitness_target')

		% intermediate negatives
		interactions = find(sga.eps(i,:) > -0.12 & sga.eps(i,:) < -0.08 & ...
                          sga.pvl(i,:) < 0.05 );
		for j=interactions
			% is this array also a (BA) query?
			ixBAq = help_equiv_set(filt.Cannon.Orf{j}, Cannon_strains, equiv);
			% is this (AB) query on the (BA) array?
			ixBAa = help_equiv_set(filt.Cannon.Orf{i}, Cannon_strains, equiv);

			% the mapping isn't quite 1 - 1; This affects maybe 20 queries and/or arrays
			if(length(ixBAq) > 1)
				ixBAq = RandomSubset(ixBAq,1);
			end
			if(length(ixBAa) > 1)
				ixBAa = RandomSubset(ixBAa,1);
			end

			if(isempty(ixBAq) || isempty(ixBAa)) % unavailable BA -> TOSS

				filt.eps(i,j) = 0;
				filt.pvl(i,j) = 1;
				filt.tossed(i,j) = 1;
				ABBA_unavail = ABBA_unavail+1;

			elseif(sga.eps(i,j) < 0 && ...
                ~(sga.eps(ixBAq,ixBAa)<-0.08 && ...
					   sga.pvl(ixBAq,ixBAa)<0.05)) % BA negative failed test -> TOSS

				filt.eps(i,j) = 0;
				filt.pvl(i,j) = 1;
				filt.tossed(i,j) = 1;
				ABBA_refuse = ABBA_refuse+1;

			else
				ABBA_rescue = ABBA_rescue+1;
			end
		end

		% ALL positives

		% record the "significant" positives
		sig_pos = sga.eps(i,:) > 0.08 & sga.pvl(i,:)<0.05;
		filt.tossed(i,sig_pos) = 1;

		% remove anything > 0
		filt.pvl(i,filt.eps(i,:)>0) = 1;
		filt.eps(i,filt.eps(i,:)>0) = 0;
	end
	fprintf('ABBA_unavail: %d\nABBA_rescue:  %d\nABBA_refuse: %d\n\n', ABBA_unavail,ABBA_rescue,ABBA_refuse);
end

function[array_ix] = help_equiv_set(query, Cannon_strains, equiv)
% match tsa with tsq and sn with dma
% and vice versa!
% damp with ?? (tsa)

	ix = strfind(query, '_');
	strain_id = query(ix+1:end);
	array_ix = [];
	
	ix = strmatch(strain_id, equiv(:,1), 'exact');
	if(~isempty(ix))
		array = equiv{ix,2};
		array_ix = strmatch(array, Cannon_strains, 'exact');
		
		return
	else

		ix = strmatch(strain_id, equiv(:,2), 'exact');
		if(~isempty(ix))
			array = equiv{ix,1};
			array_ix = strmatch(array, Cannon_strains, 'exact');
			return
		end
	end
end

function[cobatch_scores] = calculate_cobatch_by_query(sga, inputfile)
	% rip through the input file and get bach assignments
	% in python
	cobatch_file = [inputfile(1:end-4) '.bch'];
	exec_string = ['GenerateCoBatchStandard.py ' inputfile ' > ' cobatch_file];
	[status] = system(exec_string);
	if(status > 0)
		fprintf('Warning Python CoBatch process returned nonzero exit status');
	end

	c_fid = fopen(cobatch_file, 'r');
	cobatch = textscan(c_fid, '%s%s'); % orf_strain orf_strain
	fclose(c_fid);
	cobatch = [cobatch{1} cobatch{2}];

	Queries = sga.Cannon.Orf(sga.Cannon.isQuery);
	results = zeros(length(Queries),1); % co_batch IP

	the_map = Hash(java.util.HashMap(length(Queries)), Queries);
	CoBatch = boolean(sparse(length(cobatch), length(cobatch)));
	for i=1:length(cobatch)
		ida = the_map.get(cobatch{i,1});
		idb = the_map.get(cobatch{i,2});
		if(~isempty(ida) && ~isempty(idb))
			CoBatch(ida,idb) = true;
		end 
	end
	CoBatch = CoBatch | CoBatch';
	CoBatch(boolean(eye(size(CoBatch,1)))) = false;

	EPS = sga.eps(sga.Cannon.isQuery, sga.Cannon.isArray);
	EPS(isnan(EPS)) = 0;
	%IP = EPS * EPS';
	IP = corrcoef(EPS', 'rows', 'pairwise'); % PEARSON is BETTER
	for i=1:length(Queries)
		results(i,1) = median(IP(CoBatch(:,i),i));
	end
		
	cobatch_scores = [Queries num2cell(results)];
end

function[vec] = substrmatch(str, cellary)
	%function[vec] = substrmatch(str, cellary)
	% see also: Code/cellgrep

	vec = find(~cellfun(@isempty, strfind(cellary, str)));

end


