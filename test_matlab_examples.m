close all; clear; clc

% add all folders to your path
file_path = which(mfilename('fullpath'));
folder_name = fullfile(fileparts(file_path));
addpath(genpath(folder_name))

% list of examples to run
ex = {};
ex{end+1} = @ex_building_wire; % 1
ex{end+1} = @ex_newton; % 1
ex{end+1} = @ex_production_planning; % 1
ex{end+1} = @ex_traffic_network; % 1
ex{end+1} = @ex_assignment; % 2
ex{end+1} = @ex_maximum_flow; % 2
ex{end+1} = @ex_network; % 2
ex{end+1} = @ex_null_space; % 2
ex{end+1} = @ex_simplex; % 2
ex{end+1} = @ex_line_search_cubic; % 3
ex{end+1} = @ex_newton_minimization; % 3
ex{end+1} = @ex_nn; % 3
ex{end+1} = @ex_quasi_newton; % 3
ex{end+1} = @ex_steepest_descent; % 3
ex{end+1} = @ex_trust_region; % 3
ex{end+1} = @ex_active_set; % 4
ex{end+1} = @ex_barrier; % 4
ex{end+1} = @ex_sqp; % 4
ex{end+1} = @ex_termination_conditions; % 4
ex{end+1} = @ex_truss; % 4
ex{end+1} = @ex_direct_1d; % 5
ex{end+1} = @ex_golden_section; % 5
ex{end+1} = @ex_matlab_ps_ga_pso; % 5
ex{end+1} = @ex_simple_ga; % 5
ex{end+1} = @ex_simple_pso; % 5

% number of examples
n = length(ex);

% initialize progress bar
f = waitbar(0, 'Starting');

% run through each example
for idx = 1:n

    % run the example
    local_function(ex{idx})

    % increment progress bar
    waitbar(idx/n, f, strcat("Progress: ",string(floor(idx/n*100))," \%"));
       
end

% delete the progress bar
delete(f)

% get a list of licenses in use
license('inuse')

% local function (to control the script workspace)
function local_function(fun)

fun()

end