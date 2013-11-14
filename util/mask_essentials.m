function[sga] = mask_essentials(sga, select)
%function[sga] = mask_essentials(sga, select)
% select = e.g. 'EE' or 'AN'
%
% Query Options:
% E Essential
% N Non-Essential
% D Damp
% d !Damp
% A All
% T Trigenic (DM queries)
% C Trigenic (DM controls)
% M Misc (_y... && _u...)
%
% Array Options:
% E Essential
% N Non-Essential
% A All

	assert(length(select) == 2, 'rtfm');

	Qess = logical(zeros(sga.Cannon.GENES,1));
	Qess(substrmatch('_tsq', sga.Cannon.Orf)) = true;

	Qdmp = logical(zeros(sga.Cannon.GENES,1));
	Qdmp(substrmatch('_damp', sga.Cannon.Orf)) = true;

	Qnss = logical(zeros(sga.Cannon.GENES,1));
	Qnss(substrmatch('_sn', sga.Cannon.Orf)) = true;

	Qdmq = logical(zeros(sga.Cannon.GENES,1));
	Qdmq(substrmatch('+', sga.Cannon.Orf)) = true;
	Qdmq(substrmatch('YDL227C', sga.Cannon.Orf)) = false;

	Qdmc = logical(zeros(sga.Cannon.GENES,1));
	Qdmc(substrmatch('YDL227C+', sga.Cannon.Orf)) = true;
	Qdmc(substrmatch('+YDL227C+', sga.Cannon.Orf)) = true;

	Qmisc = logical(zeros(sga.Cannon.GENES,1));
	Qmisc(substrmatch('_y', sga.Cannon.Orf)) = true;
	Qmisc(substrmatch('_u', sga.Cannon.Orf)) = true;
	Qmisc(substrmatch('_col', sga.Cannon.Orf)) = true;

	Qall = sga.Cannon.isQuery;

	Aess = logical(zeros(sga.Cannon.GENES,1));
	Aess(substrmatch('_tsa', sga.Cannon.Orf)) = true;

	Anss = logical(zeros(sga.Cannon.GENES,1));
	Anss(substrmatch('_dma', sga.Cannon.Orf)) = true;

	Aall = sga.Cannon.isArray;


	if select(1) == 'E'
		sga.Cannon.isQuery(~Qess) = false;
	elseif select(1) == 'N'
		sga.Cannon.isQuery(~Qnss) = false;
	elseif select(1) == 'D'
		sga.Cannon.isQuery(~Qdmp) = false;
	elseif select(1) == 'd'
		sga.Cannon.isQuery(Qdmp) = false;
	elseif select(1) == 'A'
		sga.Cannon.isQuery(~Qall) = false; % tautology
	elseif select(1) == 'T'
		sga.Cannon.isQuery(~Qdmq) = false;
	elseif select(1) == 'C'
		sga.Cannon.isQuery(~Qdmc) = false;
	elseif select(1) == 'M'
		sga.Cannon.isQuery(~Qmisc) = false;
	else
		error('unrecognized option [ENDdTCMA]')
	end
		
	if select(2) == 'E'
		sga.Cannon.isArray(~Aess) = false;
	elseif select(2) == 'N'
		sga.Cannon.isArray(~Anss) = false;
	elseif select(2) == 'A'
		sga.Cannon.isArray(~Aall) = false; % tautology
	else
		error('unrecognized option [ENA]')
	end
end
