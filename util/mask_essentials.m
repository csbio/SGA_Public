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
% t (toss any with a +)
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
		[EsuppQ, ~] = E_supp(sga);
		sga.Cannon.isQuery(~Qess) = false;
		sga.Cannon.isQuery(EsuppQ) = true; % reset essential suppressors
	elseif select(1) == 'N'
		Qnss = logical(zeros(sga.Cannon.GENES,1));
		Qnss(substrmatch('_sn', sga.Cannon.Orf)) = true;
		[NsuppQ, ~] = N_supp(sga);
		sga.Cannon.isQuery(~Qnss) = false;
		sga.Cannon.isQuery(NsuppQ) = true; % reset N suppressors
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
		Qdmc(substrmatch('+YDL227C', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(~Qdmc) = false;
	elseif select(1) == 't'
		Qdm =logical(zeros(sga.Cannon.GENES,1));
		Qdm(substrmatch('+', sga.Cannon.Orf)) = true;
		sga.Cannon.isQuery(Qdm) = false;
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
		Aess = logical(zeros(sga.Cannon.GENES,1));
		Aess(substrmatch('_tsa', sga.Cannon.Orf)) = true;
		[~, EsuppA] = E_supp(sga);
		sga.Cannon.isArray(~Aess) = false;
		sga.Cannon.isArray(EsuppA) = true; % reset essential supp arrays
	elseif select(2) == 'N'
		Anss = logical(zeros(sga.Cannon.GENES,1));
		Anss(substrmatch('_dma', sga.Cannon.Orf)) = true;
		sga.Cannon.isArray(~Anss) = false;
	elseif select(2) == 'A'
		Aall = sga.Cannon.isArray;
		[~, NsuppA] = E_supp(sga);
		sga.Cannon.isArray(~Aall) = false; % tautology
		sga.Cannon.isArray(NsuppA) = true; % reset N supp arrays
	else
		error('unrecognized option [ENA]')
	end
end


% for suppressor strain support in E and N where labels are _S
function [isQuery, isArray] = E_supp(sga)
   % return vectors: isQuery & is_essential & is_suppressor
   supp_def_file = '~/SGA/refdata/suppressor_strain_essentiality_160425.csv';
   supp_def = Csv2Cell(supp_def_file);
   E_ix = strcmp('E', supp_def(:,2));
   supp_ess = supp_def(E_ix,1);

   [~, tails] = StripOrfs(sga.Cannon.Orf);
   is_ess = ismember(tails, supp_ess); %is_supp implied
   isQuery = sga.Cannon.isQuery & is_ess;
   isArray = sga.Cannon.isArray & is_ess';
end
function [isQuery, isArray] = N_supp(sga)
   % return vectors: isQuery & is_essential & is_suppressor
   supp_def_file = '~/SGA/refdata/suppressor_strain_essentiality_160425.csv';
   supp_def = Csv2Cell(supp_def_file);
   N_ix = strcmp('N', supp_def(:,2));
   supp_non = supp_def(N_ix,1);

   [~, tails] = StripOrfs(sga.Cannon.Orf);
   is_non = ismember(tails, supp_non); %is_supp implied
   isQuery = sga.Cannon.isQuery & is_non;
   isArray = sga.Cannon.isArray & is_non';
end
