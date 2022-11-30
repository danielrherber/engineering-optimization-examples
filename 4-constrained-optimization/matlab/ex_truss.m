close all; clear; clc

% problem data
stress_limit = 25e3;
min_area = 0.1;

% compute minimum mass (for scaling)
[mass0, ~] = truss(min_area*ones(10,1));

% problem functions and constraints
FUN = @(x) objective(x,mass0); % objective function
X0 = ones(10,1); % initial guess
A = []; % A matrix in linear inequality constraints
B = []; % b vector in linear inequality constraints
Aeq = []; % A matrix in linear equality constraints
Beq = []; % b vector in linear equality constraints
LB = min_area*ones(10,1); % simple lower bounds
UB = []; % simple upper bounds
NONLCON = @(x) constraints(x,stress_limit); % nonlinear constraint function

% options
OPTIONS = optimoptions('fmincon');
OPTIONS.Display = 'iter';
% OPTIONS.FiniteDifferenceType = 'central';
OPTIONS.MaxFunctionEvaluations = 1e4;
OPTIONS.OptimalityTolerance = 1e-16;
OPTIONS.ConstraintTolerance = 1e-16;
OPTIONS.StepTolerance = 1e-16;

% solve the optimization problem
X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS);

% calculate the objective and constraints at minimizer
[mass, stress] = truss(X);

% display solution
disp_helper("X",X,4)
disp_helper("stress",stress,4)

% plot of the truss members
start = [5, 3, 6, 4, 4, 2, 5, 6, 3, 4];
finish = [3, 1, 4, 2, 3, 1, 4, 3, 2, 1];
weights = X;
names = mat2str(round(stress,1));
hf = figure; hf.Color = 'w';
G = graph(start,finish,weights);
H = plot(G,'EdgeLabel',round(G.Edges.Weight,1),'LineWidth',G.Edges.Weight);
H.XData = [2 2 1 1 0 0];
H.YData = [1 0 1 0 1 0];

%% helper functions
% objective function
function f = objective(x,mass0)

% compute mass
[mass, ~] = truss(x);

% assign objective value
% f = mass;
f = mass/mass0; % scaled version

end

% nonlinear constraint function
function [c,ceq] = constraints(x,stress_limit)

% compute stress
[~, stress] = truss(x);

% inequality constraints
% c = abs(stress) - stress_limit;
% c = [stress - stress_limit; -stress - stress_limit]; % differentiable version
c = [stress/stress_limit - 1; -stress/stress_limit - 1]; % scaled version

% equality constraints
ceq = []; % empty

end

% function to make it easier to display things in the command window
function disp_helper(name,number,n)

% default value of the number of digits
if isempty(n)
    n = 5;
end

% form string
str = strcat(string(name)," = ",mat2str(round(number,n)));

% display string
disp(str)

end