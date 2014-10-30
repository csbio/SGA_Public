function [] = compute_sgascore_function(parameter_struct)
%function [] = compute_sgascore_function(parameter_struct)
% for when you really need a function instead of a script

    % unpack the parameters
    parameters = fieldnames(parameter_struct);
    for i = 1:length(parameters)
        eval(sprintf('%s = parameter_struct.%s;', parameters{i}, parameters{i}));
    end

    % clear up the unpack (just in case)
    clear parameters parameter_struct i

    % turn on the debugger to gain all the advantages of a script!
    dbstop if error

    % run the script
    compute_sgascore
end
