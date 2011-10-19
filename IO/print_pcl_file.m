function print_pcl_file(data,orf_names,orf_annotations,exp_labels,filename),

fid = fopen(filename,'w');
fprintf(fid,'ORF\tNAME\tGWEIGHT\t');
for i=1:length(exp_labels),
    if(i < length(exp_labels))
        fprintf(fid,'%s\t',exp_labels{i});
    else
        fprintf(fid,'%s\n',exp_labels{i});        
    end
end

fprintf(fid,'EWEIGHT\t\t\t');
fprintf(fid,'%d\t',ones(1,size(data,2)-1));
fprintf(fid,'1\n');

for i=1:size(data,1),
   fprintf(fid,'%s\t%s\t1\t',orf_names{i},orf_annotations{i});
   for j=1:size(data,2)-1,
       if(isnan(data(i,j)))
           fprintf(fid,'\t');
       else
           fprintf(fid,'%f\t',data(i,j));
       end
   end
   if(isnan(data(i,end)))
       fprintf(fid,'\n');
   else
       fprintf(fid,'%f\n',data(i,end));
   end
end
fclose(fid);