function[ret_val] = log_printf(logfile_fid, varargin)
%function[ret_val] = log_printf(logfile_fid, varargin)
%
% log_fprintf behaves the same as printf, except it prints to BOTH
% standard out and to an open file (first argument like printf)
% The standard "failed fid" is -1, we'll check for a -11 as a way
% to suppress both writing to a file and pesky error messages
%
% Example 1: keep a log
%    log_fid = fopen('logfile' 'w');
%    log_printf(log_fid, 'standard messages');
%    fclose(log_fid)
%
% Example 2: log open fails, will throw warnings
%    log_fid = fopen('/invalid/path/or/no/permission', 'w'); % returns -1
%    log_printf(log_fid, 'standard messages')
%       standard messages
%       WARNING Could not print to log file
% 
% Revision: 111206
% Author:   Benjamin VanderSluis (bvander@cs.umn.edu)


return2 = fprintf(varargin{:});

if(logfile_fid ~= -11)
    try
       return1 = fprintf(logfile_fid, varargin{:});
    catch
       fprintf('WARNING Could not print to log file\n');
       return1 = -1; 
    end
else
    return1 = inf;
end

% return an error if either operation returned an error
% or the number of bytes written which will be the same
% for either operation
ret_val = min(return1, return2);
