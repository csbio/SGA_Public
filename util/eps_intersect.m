function[E1, E2, P1, P2, ComQ, ComA, D1, D2] = eps_intersect(sga1, sga2)
%function[E1, E2, P1, P2, ComQ, ComA] = eps_intersect(sga1, sga2)
% takes in two structs and returns "matched" matricies
% see also: sga_intersect

	% old epsilon code
	% E1 and E2 will be "matched"
	E1 = sga1.eps(sga1.Cannon.isQuery, sga1.Cannon.isArray);
	E2 = sga2.eps(sga2.Cannon.isQuery, sga2.Cannon.isArray);

	P1 = sga1.pvl(sga1.Cannon.isQuery, sga1.Cannon.isArray);
	P2 = sga2.pvl(sga2.Cannon.isQuery, sga2.Cannon.isArray);

	D1 = sga1.dbl(sga1.Cannon.isQuery, sga1.Cannon.isArray);
	D2 = sga2.dbl(sga2.Cannon.isQuery, sga2.Cannon.isArray);

	[ComQ, ix1, ix2] = intersect(sga1.Cannon.Orf(sga1.Cannon.isQuery), ...
		                         sga2.Cannon.Orf(sga2.Cannon.isQuery));
	[ComA, iy1, iy2] = intersect(sga1.Cannon.Orf(sga1.Cannon.isArray), ...
		                         sga2.Cannon.Orf(sga2.Cannon.isArray));

	E1 = E1(ix1,iy1);
	E2 = E2(ix2,iy2);

	P1 = P1(ix1,iy1);
	P2 = P2(ix2,iy2);

	D1 = D1(ix1,iy1);
	D2 = D2(ix2,iy2);

end
