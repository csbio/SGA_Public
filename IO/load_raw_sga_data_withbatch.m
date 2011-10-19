%%
% LOAD_RAW_SGA_DATA_WITHBATCH
%
% Description:
%   This function loads a text file with raw SGA colony size data.
%   The file is assumed to be in the following tab-delimited format:
%   query ORF, array ORF, arrayplateid, setID, plateID, batchID, row, column, colony
%   size (pixels)
%
% Syntax:
%   data = load_raw_sga_data_withbatch(inputfile)
%
% Inputs:
%   inputfile - input filename
%
% Outputs:
%   data - a structure with fields corresponding to each of the columns 
%
% Authors: Chad Myers (cmyers@cs.umn.edu), 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca),
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2011-02-09
%
%%
function sgadata = load_raw_sga_data_withbatch(inputfile, skip_perl_step)

% Print the name and path of this script
p = mfilename('fullpath');
fprintf('\nLoad raw SGA data with batch script: %s\n\n',p);

% First, run perl script to convert data to completely numeric (this speeds
% the input process up)
if(~skip_perl_step)
    fprintf('beginning perl preprocessing\n');
    tic
    perl('process_rawsga_dmdata.pl',inputfile);
    toc
else
    fprintf('skipping perl preprocessing...\n');
end

% Fields: 
% 1 query_id 
% 2 array_id 
% 3 arrayplateids 
% 4 setids 
% 5 plateids 
% 6 batchid	
% 7 rowid 
% 8 colid 
% 9 colsize 

field_string = '%n%n%n%n%n%n%n%n%n';

% Buffered textscan is much faster than importdata
input_fid = fopen([inputfile '_numeric'], 'r');
rawdata = textscan(input_fid, field_string, 'BufSize', 2^20, 'ReturnOnError', false);
fclose(input_fid);

sgadata.querys = rawdata{1};
sgadata.arrays = rawdata{2};
sgadata.arrayplateids = rawdata{3};
sgadata.setids = rawdata{4};
sgadata.plateids = rawdata{5};
sgadata.batch = rawdata{6};
sgadata.rows = rawdata{7};
sgadata.cols = rawdata{8};
sgadata.colsize = rawdata{9};

data=importdata([inputfile,'_orfidmap']);
[vals,ind] = sort(data.data);
sgadata.orfnames = data.textdata(ind);
