function [ cellarr ] = StripHO(cellarr)
%function [ cellarr ] = StripHO(cellarr)
% removes HO or YDL227C from strainIDs
% also the leading or trailing + (by assumption)
% replaces strain IDs (use StripOrfs() to modify)

	for i=1:numel(cellarr)
		% check for common name
		ix_ho = findstr('HO', cellarr{i});
		if ~isempty(ix_ho)
			if ix_ho == 1
				cellarr{i} = cellarr{i}(4:end); % HO+xxx -> xxx
			else
				cellarr{i} = cellarr{i}([1:ix_ho-2, ix_ho+2:end]); % will this blow up with no _xxx?
			end
		end
	
		% check for ORFS
		ix_ydl = findstr('YDL227C', cellarr{i});
		if ~isempty(ix_ydl)
			if ix_ydl == 1
				cellarr{i} = cellarr{i}(9:end); % HO+xxx -> xxx
			else
				cellarr{i} = cellarr{i}([1:ix_ydl-2, ix_ydl+7:end]); % will this blow up with no _xxx?
			end
		end
	end
end
