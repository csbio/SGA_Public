function[duplicates, counts, ix] = FindDuplicates(cellarr)
% function[duplicates, counts, cell_ix] = FindDuplicates(cellarr)
% returns duplicates (unique list), 
% the number of times each appears
% and a cell array of their instances in the original cellarr

	srt = sort(cellarr);
	duplicates = {};
	counts = [];

	ptr = 1;
	while ptr < length(srt)
		count = 1;
		ptrr = ptr+1;
		while( ptrr <= length(srt) && strcmp(srt{ptr}, srt{ptrr}))
			count = count +1;
			ptrr = ptrr +1;
		end

		if(count > 1)
			duplicates = [duplicates; srt{ptr}];
			counts = [counts; count];
		end

		ptr = ptrr;
	end


	ix = cell(size(counts));
	for i=1:length(ix)
		ix{i} = strmatch(duplicates{i}, cellarr, 'exact');
	end
end
