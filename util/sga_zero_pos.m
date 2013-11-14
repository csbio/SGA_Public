function [sga] = sga_zero_pos(sga, val)
%function [sga] = sga_zero_pos(sga, val)

% set all pos to val (i.e. zero or NaN)
% Doesn't pay attn to isQuery or isArray
% if pos = true (default false) then only >0



	if ~exist('val', 'var')
		val = NaN;
	end

	mask = sga.eps > 0;
	mask = mask | (sga.eps > 0 & sga.pvl > 0.05);
	sga.eps(mask) = val;
end
