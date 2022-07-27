close all; clear; clc

% add all folders to your path
file_path = which(mfilename('fullpath'));
folder_name = fullfile(fileparts(file_path));
addpath(genpath(folder_name))

% list of examples to run
ex = {};
ex{end+1} = @ex_building_wire;
ex{end+1} = @ex_newton;
ex{end+1} = @ex_production_planning;
ex{end+1} = @ex_traffic_network;
ex{end+1} = @ex_assignment;
ex{end+1} = @ex_network;
ex{end+1} = @ex_null_space;
ex{end+1} = @ex_simplex;
ex{end+1} = @ex_newton_minimization;
ex{end+1} = @ex_nn;
ex{end+1} = @ex_quasi_newton;
ex{end+1} = @ex_steepest_descent;
ex{end+1} = @ex_trust_region;
ex{end+1} = @ex_active_set;
ex{end+1} = @ex_barrier;
ex{end+1} = @ex_sqp;
ex{end+1} = @ex_truss;
ex{end+1} = @ex_golden_section;
ex{end+1} = @ex_matlab_ps_ga_pso;
ex{end+1} = @ex_simple_ga;
ex{end+1} = @ex_simple_pso;
% ex{end+1} = @ex_surrogateopt;

% number of examples
n = length(ex);

% initialize progress bar
f = waitbar(0, 'Starting');

% run through each example
for idx = 1:n

    % run the example
    local_function(ex{idx})

    % increment progress bar
    waitbar(idx/n, f, sprintf('Progress: %d %%', floor(idx/n*100)));

end

% delete the progress bar
delete(f)

% display inuse licenses
license('inuse')

% local function (to control the script workspace)
function local_function(fun)

fun()

end