function[sga] = reduce_alleles(sga)
%function[sga] = reduce_alleles(sga)
% Choose one representitive ORF for each (Q, A) at random...
% and mask out the rest.

	All_Qs = strip_annotation(sga.Cannon.Orf(sga.Cannon.isQuery), 'first');
	Qs = unique(All_Qs);
	All_As = strip_annotation(sga.Cannon.Orf(sga.Cannon.isArray), 'first');
	As = unique(All_As);

	QQ=find(sga.Cannon.isQuery);
	AA=find(sga.Cannon.isArray);


	for i=1:length(Qs)
		ix = strmatch(Qs{i}, All_Qs, 'exact');
		if(length(ix)>1)
			ix = random_subset(ix, length(ix)-1);
			sga.Cannon.isQuery(QQ(ix)) = false;
		end
	end
		
		
	for i=1:length(As)
		ix = strmatch(As{i}, All_As, 'exact');
		if(length(ix)>1)
			ix = random_subset(ix, length(ix)-1);
			sga.Cannon.isArray(AA(ix)) = false;
		end
	end

end
