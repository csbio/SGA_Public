function[joined] = JoinOrfs(a,b, delim)
%function[joined] = JoinOrfs(a, b, [delim=+])
	% joins cellarrays into a single col
	% a may by Nx2

	% maybe a is N,2 and there is no B?
	if size(a,2)==2 && ~exist('delim', 'var')
		if exist('b', 'var')
			delim = b; 
		end
		b = a(:,2);
		a = a(:,1);
	end

	if ~exist('delim', 'var')
		delim = '+';
	end

	% support single strings also
	if isstr(a)
		[joined] = single_join(a, b, delim);
	else
		joined = cell(size(a));
		for i=1:length(a)
			joined{i} = single_join(a{i}, b{i}, delim);
		end
	end
end

function[string] = single_join(a,b,d)
	string = [a d b];
end

