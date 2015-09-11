function cell2csv(filename,cellArray,seperator)
%function cell2csv(filename,cellArray,seperator)
% Writes cell array content into a *.csv file.
% 
% CELL2TAB(filename,cellArray,separator)
%
% filename(fid)= Name of the file to save. [ i.e. 'text.csv' ]
%                can also be the descriptor (int) of an open file
% cellarray    = Name of the Cell Array where the data is in
% seperator    = seperating sign, normally:'\t' (it's default)

if(~exist('seperator', 'var'))
	seperator = '\t';
end

if(isnumeric(filename))
	fid = filename;
else
	fid = fopen(filename,'w');
end

for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)
        
        var = eval(['cellArray{z,s}']);
        
        if size(var,1) == 0
            var = '';
        end
        
        if isnumeric(var) == 1
            var = num2str(var);
        end
        
        if islogical(var) == 1
            if var == 1
                var = 'TRUE';
            else
                var = 'FALSE';
            end
        end

        
        % we need to escape '%'s  % -> %%
        ix = findstr('%', var);
        ix = ix(end:-1:1); % start at the end to preserve index
        for i=1:length(ix)
           var = [var(1:ix(i)) '%' var(ix(i)+1:end)];
        end

        fprintf(fid,var);
        
        if s ~= size(cellArray,2)
            fprintf(fid,seperator);
        end
    end
    fprintf(fid,'\n');
end

if(~isnumeric(filename))
	fclose(fid);
end
