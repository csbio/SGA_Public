function [sga1] = sga_set_diff(sga1, sga2)
%function [sga1] = sga_set_diff(sga1, sga2)
% remove any sig interactions in 1 from 2
% ASSUMES MATCHED STRUCTS

	assert(sga1.Cannon.GENES == sga2.Cannon.GENES);

	sga1.eps(sga2.eps <-0.08 & sga2.pvl < 0.05) = NaN;
	sga1.eps(sga2.eps > 0.08 & sga2.pvl < 0.05) = NaN;

end

