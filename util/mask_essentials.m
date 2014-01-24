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

	if select(1) == 'E'
		Qess = logical(zeros(sga.Cannon.GENES,1));
		Qess(substrmatch('_tsq', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(~Qess) = false;
	elseif select(1) == 'N'
		Qnss = logical(zeros(sga.Cannon.GENES,1));
		Qnss(substrmatch('_sn', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(~Qnss) = false;
	elseif select(1) == 'D'
		Qdmp = logical(zeros(sga.Cannon.GENES,1));
		Qdmp(substrmatch('_damp', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(~Qdmp) = false;
	elseif select(1) == 'd'
		Qdmp = logical(zeros(sga.Cannon.GENES,1));
		Qdmp(substrmatch('_damp', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(Qdmp) = false;
	elseif select(1) == 'A'
		Qall = sga.Cannon.isQuery;
		sga.Cannon.isQuery(~Qall) = false; % tautology
	elseif select(1) == 'T'
		Qdmq = logical(zeros(sga.Cannon.GENES,1));
		Qdmq(substrmatch('+', sga.Cannon.Orf)) = true;
		Qdmq(substrmatch('YDL227C', sga.Cannon.Orf)) = false;
		sga.Cannon.isQuery(~Qdmq) = false;
	elseif select(1) == 'C'
		Qdmc = logical(zeros(sga.Cannon.GENES,1));
		Qdmc(substrmatch('YDL227C+', sga.Cannon.Orf)) = true;
		Qdmc(substrmatch('+YDL227C+', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(~Qdmc) = false;
	elseif select(1) == 'M'
		Qmisc = logical(zeros(sga.Cannon.GENES,1));
		Qmisc(substrmatch('_y', sga.Cannon.Orf)) = true;
		Qmisc(substrmatch('_u', sga.Cannon.Orf)) = true;
		Qmisc(substrmatch('_col', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(~Qmisc) = false;
	else
		error('unrecognized option [ENDdTCMA]')
	end
		
	if select(2) == 'E'
		sga.Cannon.isArray(~Aess) = false;
		Aess = logical(zeros(sga.Cannon.GENES,1));
		Aess(substrmatch('_tsa', sga.Cannon.Orf)) = true;
	elseif select(2) == 'N'
		Anss = logical(zeros(sga.Cannon.GENES,1));
		Anss(substrmatch('_dma', sga.Cannon.Orf)) = true;
		sga.Cannon.isArray(~Anss) = false;
	elseif select(2) == 'A'
		Aall = sga.Cannon.isArray;
		sga.Cannon.isArray(~Aall) = false; % tautology
	else
		error('unrecognized option [ENA]')
	end
end
