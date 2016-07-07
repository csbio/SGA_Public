function [ NVec ] = random_subset( Vector, n )
%function [ NVec ] = random_subset( Vector, n )
%RANDOMSUBSET returns random subset of column vector (length n)
%   if Vector is boolean then "n" ones will be preserved at random
%   without n, random_subset returns a random permutation

if(~exist('n', 'var'))
	n = length(Vector);
end

if(islogical(Vector))
	if(n > sum(Vector))
		error('You are asking too much! sum(V)=%d n=%d\n',sum(Vector), n);
	else
		NVec = boolean(zeros(size(Vector)));
		NVec(MyChoose(find(Vector), n)) = true;
	end
else
	if(n>length(Vector))
		error('You are asking too much! length(V)=%d n=%d\n', length(Vector), n);
	else
		NVec = MyChoose(Vector, n);
	end
end
end

function[NVec] = MyChoose(Vector, n)

	if(iscell(Vector))
		NVec = cell(n,1);
	else
		NVec = zeros(n,1);
	end

	Vsize = length(Vector);

	A = 1;
	for i=1:n
		ind = ceil(rand()*Vsize);
		if(iscell(Vector))
			NVec{i} = Vector{ind};
			Vector{ind} = Vector{Vsize};
		else
			NVec(i) = Vector(ind);
			Vector(ind) = Vector(Vsize);
		end
		Vsize = Vsize-1;
	end

end

