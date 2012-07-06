function[filt, fitness_struct] = filter_interactions(sga, fitness_file, sga_inputfile)
%function[filt, fitness_struct] = filter_interactions(sga, fitness_file, cobatch_scores)

	% load the fitness data
	smf_fid = fopen(fitness_file, 'r');
	A = textscan(smf_fid, '%s%s%s', 'Delimiter', '\t', 'ReturnOnError', 'False');
	fitness_struct = [A{1} num2cell(A{2}) num2cell(A{3})];
	fclose(smf_fid);

	% load / calculate the cobatch_scores
	cobatch_scores = calculate_cobatch_by_query(sga, sga_inputfile)


	% filter 
	filt = help_filter_by_cobatch(sga, cobatch_scores);
	filt = help_filter_by_ABBA(filt, fitness_struct);

end

function[filt] = help_filter_by_cobatch(sga, cobatch_scores)
% sets scores to NaN and removes from query index

	filt = sga;
	cobatch_target = filt.Cannon.isQuery;
	for i=find(cobatch_target')
		ix = strmatch(filt.Cannon.Orf{i}, cobatch_scores(:,1));
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

function[filt] = help_filter_by_ABBA(sga, fitness_struct)
	% abs(eps) > 0.12 & pvl || AB & BA
	% sets ALL non interactions to ZERO

	% isolate all intermediate interactions in query
	% fitness range. Target for removal, check in a loop...
	
	filt = sga;
	filt.tossed = sparse(zeros(size(sga.eps)));
	fitness_target = filt.Cannon.isQuery;
	for i=find(fitness_target')
		ix = strmatch(filt.Cannon.Orf{i}, fitness_struct(:,1));
		if(~isempty(ix) && fitness_struct{ix,2}>0.97) % target aquired
			fitness_target(i) = true;
		else
			fitness_target(i) = false;
		end
	end

	% for the sick queries, we need to isolate the interactions
	
	fprintf('%d queries targeted for ABBA filtering\n', sum(fitness_target));
	filt.fitness_target = fitness_target;
	for i=find(fitness_target')



		% intermediate negatives
		interactions = find(filt.eps(i,:) > -0.12 & filt.eps(i,:) < -0.08 & ...
                          filt.pvl(i,:) < 0.05 );
		for j=interactions
			% is this array also a (BA) query?
			ixBAq = help_equiv_set(filt.Cannon.Orf{j}, filt.Cannon);
			% is this (AB) query on the (BA) array?
			ixBAa = help_equiv_set(filt.Cannon.Orf{i}, filt.Cannon);

			% for now we pick one random allele; affects maybe 25% of the data
			if(length(ixBAq) > 1)
				ixBAq = RandomSubset(ixBAq, 1);
			end
			if(length(ixBAa) > 1)
				ixBAa = RandomSubset(ixBAa, 1);
			end

			if(isempty(ixBAq) || isempty(ixBAa)) % unavailable BA -> TOSS

				filt.eps(i,j) = 0;
				filt.pvl(i,j) = 0;
				filt.tossed(i,j) = 1;

			elseif(filt.eps(i,j) < 0 && ...
                ~(filt.eps(ixBAq,ixBAa)<-0.08 && ...
					   filt.pvl(ixBAq,ixBAa)<0.05)) % BA negative failed test -> TOSS

				filt.eps(i,j) = 0;
				filt.pvl(i,j) = 0;
				filt.tossed(i,j) = 1;

			end
		end

		% ALL positives

		% record the "significant" positives
		sig_pos = filt.eps(i,:) > 0.08 & filt.pvl(i,:)<0.05;
		filt.tossed(i,sig_pos) = 1;

		% remove anything > 0
		filt.pvl(i,filt.eps(i,:)>0) = 0;
		filt.eps(i,filt.eps(i,:)>0) = 0;



	end
end

function[array_set] = help_equiv_set(query, cannon)
% match tsa with tsq and sn with dma
% and vice versa!
% damp with ?? (tsa)

	if(~isempty(strfind(query, '_tsq')))
		array_set = strmatch([StripOrfs(query) '_tsa'], cannon.Orf); 
	elseif(~isempty(strfind(query, '_tsa')))
		array_set = strmatch([StripOrfs(query) '_tsq'], cannon.Orf); 
	elseif(~isempty(strfind(query, '_sn')))
		array_set = strmatch([StripOrfs(query) '_dma'], cannon.Orf); 
	elseif(~isempty(strfind(query, '_dma')))
		array_set = strmatch([StripOrfs(query) '_sn' ], cannon.Orf); 
	elseif(~isempty(strfind(query, '_damp')))
		array_set = strmatch([StripOrfs(query) '_tsa' ], cannon.Orf); 
	else
		error('Unexpected StrainID suffix');
	end
end

function[cobatch_scores] = calculate_cobatch_by_query(sga, inputfile)
	% rip through the input file and get bach assignments
	% in python
	batch_file = [inputfile '.bch'];
	exec_string = ['GenerateCoBatchStandard.py ' intputfile '.txt > ' batch_file];
	[status] = system(exec_string);
	if(status > 0)
		fprintf('Warning Python CoBatch process returned nonzero exit status');
	end



	% load and return the cobatch scores
	fid = fopen(batch_file, 'r');
	A = textscan(fid, '%s%f');
	fclose(fid);
	cobatch_scores = [A{1} num2cell(A{2})];
end
