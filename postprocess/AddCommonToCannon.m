function[Cannon] = AddCommonToCannon(Cannon, common_name_map_file)
%function[Cannon] = AddCommonToCannon(Cannon, common_name_map_file)
%

if(~exist('common_name_map_file', 'var'))
	common_name_map_file = '/project/csbio/benjamin/Data/Master_Common_Ref_SGD.txt';
end

% Hash the namemap file
map = java.util.HashMap(Cannon.GENES);
fid = fopen(common_name_map_file, 'r');
line = fgetl(fid);
while(ischar(line))
	A = textscan(line, '%s', 'Delimiter', '\t');
	if(length(A{1}) > 1)
		map.put(A{1}{1}, line);
	end
	line = fgetl(fid);
end
fclose(fid);

% try to look up every orf
Cannon.Common = cell(size(Cannon.Orf));
for i=1:Cannon.GENES
	orf = Cannon.Orf{i};
	ix = findstr('_', orf);
	suffix = [];
	if(~isempty(ix))
		suffix = orf(ix:end);
		orf = orf(1:ix-1);
	end
	line = map.get(orf);	
	if(~isempty(line))
		A = textscan(line, '%s', 'Delimiter', '\t');
		for j=2:length(A{1})
			Cannon.Map.put(java.lang.String(A{1}{j}), java.lang.Integer(i));
		end
		Cannon.Common{i} = [A{1}{2} suffix];
	else
		Cannon.Common{i} = Cannon.Orf{i};
	end
end
		
