function[v1, v2] = verify_scores(sga1, sga2)
%function[] = verify_scores(sga1, sga2)
% grab a bunch of random epsilons and cross-reference

	NUM = 10000;

	EPS1 = sga1.eps(sga1.Cannon.isQuery, sga1.Cannon.isArray);
	Q1 = sga1.Cannon.Orf(sga1.Cannon.isQuery);
	A1 = sga1.Cannon.Orf(sga1.Cannon.isArray);

	EPS2 = sga2.eps(sga2.Cannon.isQuery, sga2.Cannon.isArray);
	Q2 = sga2.Cannon.Orf(sga2.Cannon.isQuery);
	A2 = sga2.Cannon.Orf(sga2.Cannon.isArray);

	EPS1(isnan(EPS1)) = 0;

	[r,c,v1] = find(EPS1);
	length(r)

	rand_ix = RandomSubset(1:length(r), NUM);

	r = r(rand_ix);
	c = c(rand_ix);
	v1 = v1(rand_ix);

	v2 = nan(size(v1));
	for i=1:length(v2)
		ixq = strmatch(Q1{r(i)}, Q2, 'exact');
		ixa = strmatch(A1{c(i)}, A2, 'exact');
		v2(i) = EPS2(ixq, ixa);
		
	end


	[sum(isnan(v1)) sum(v1 ~= v2)]


	v1 = EPS1 * EPS1';
	v1 = sort(v1(:));

	v2 = EPS2 * EPS2';
	v2 = sort(v2(:));


	keyboard

	sum(v1 ~= v2);
end
