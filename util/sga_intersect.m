function[sga1, sga2] = sga_intersect(sga1, sga2)
%function[sga1, sga2] = sga_intersect(sga1, sga2)
% takes in two structs and returns two structs with is* 
% variables manipulated to make the Q&A sets identical
% see also: eps_intersect

	% new struct code
	% epsilon mats will not be matched but "sets" will be the same

	sga1.Cannon.isQuery(~ismember(sga1.Cannon.Orf, sga2.Cannon.Orf(sga2.Cannon.isQuery))) = false;
	sga2.Cannon.isQuery(~ismember(sga2.Cannon.Orf, sga1.Cannon.Orf(sga1.Cannon.isQuery))) = false;

	sga1.Cannon.isArray(~ismember(sga1.Cannon.Orf, sga2.Cannon.Orf(sga2.Cannon.isArray))) = false;
	sga2.Cannon.isArray(~ismember(sga2.Cannon.Orf, sga1.Cannon.Orf(sga1.Cannon.isArray))) = false;


%{
	% old epsilon code
	% E1 and E2 will be "matched"
	E1 = sga1.eps(sga1.Cannon.isQuery, sga1.Cannon.isArray);
	E2 = sga2.eps(sga2.Cannon.isQuery, sga2.Cannon.isArray);

	P1 = sga1.pvl(sga1.Cannon.isQuery, sga1.Cannon.isArray);
	P2 = sga2.pvl(sga2.Cannon.isQuery, sga2.Cannon.isArray);


	[ComQ, ix1, ix2] = intersect(sga1.Cannon.Orf(sga1.Cannon.isQuery), ...
		                         sga2.Cannon.Orf(sga2.Cannon.isQuery));
	[ComA, iy1, iy2] = intersect(sga1.Cannon.Orf(sga1.Cannon.isArray), ...
		                         sga2.Cannon.Orf(sga2.Cannon.isArray));


	E1 = E1(ix1,iy1);
	E2 = E2(ix2,iy2);

	P1 = P1(ix1,iy1);
	P2 = P2(ix2,iy2);
%}

end
