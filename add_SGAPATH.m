function add_SGAPATH()
% function add_SGAPATH()
% add project directories relative to path of this_file 

   base_dir = fileparts(mfilename('fullpath'));
   addpath(base_dir)
   addpath([base_dir '/IO']);
   addpath([base_dir '/corrections']);
   addpath([base_dir '/util']);
   addpath([base_dir '/postprocess']);
   addpath([base_dir '/linkage_est']);
   addpath([base_dir '/test']);
end
