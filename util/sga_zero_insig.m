function [sga] = sga_zero_insig(sga, val, pos)
%function [sga] = sga_zero_insig(sga, val, pos)
% set all insignif to val (i.e. zero or NaN)
% Doesn't pay attn to isQuery or isArray
% if pos = true (default false) then only >0



	if ~exist('val', 'var')
		val = NaN;
	end

	if(pos)
		mask = sga.eps > 0 & sga.eps < 0.08 ;
		mask = mask | (sga.eps > 0 & sga.pvl > 0.05);
	else
		mask = sga.eps > -0.08 & sga.eps < 0.08 ;
		mask = mask | sga.pvl > 0.05;
	end

	sga.eps(mask) = val;
end