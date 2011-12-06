%%
% PRINT_PROGRESS - prints a progress bar, given a total and a position in the progress
%
% Inputs:
%   total - the total number of steps
%   step - the position in the progress
%
% Authors: Chad Myers (cmyers@cs.umn.edu), Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%
% Last revision: 2010-07-19
%
%%

function print_progress(total, step, lfid)
    
    % Determine how much progress made since the last step of the loop
    x = fix(step / (total/50)) - fix((step-1)/(total/50));
    
    for i = 1 : x
        log_printf(lfid, '*');
    end
    
    % Fixes a little rounding issue at the end
    if step == total 
        df = 50 - fix(step / (total/50));
        for i = 1 : df
            log_printf(lfid, '*');
        end
    end
