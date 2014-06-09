function [sga] = sga_zero_pos(sga, val)
%function [sga] = sga_zero_pos(sga, [val=NaN])
% set all pos to val (e.g. zero or NaN)
% Doesn't pay attn to isQuery or isArray
% if pos = true (default false) then only >0

	if ~exist('val', 'var')
		val = NaN;
	end

	mask = sga.eps > 0;
	sga.eps(mask) = val;
end
