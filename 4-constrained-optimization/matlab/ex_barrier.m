close all; clear; clc

% problem functions
f = @(x) x(1) - 2*x(2);
g = @(x) [1 + x(1) - x(2)^2; x(2)];

% initial strictly feasible point
x = [0.5;0.5];

% initial barrier parameter
mu = 1;

% barrier term
phi = @(x) -sum(log(g(x)));

% barrier function
beta = @(x,mu) f(x) + mu*phi(x);

% fminunc options
OPTIONS = optimoptions('fminunc');
OPTIONS.Display = 'none';
OPTIONS.StepTolerance = 0;
OPTIONS.OptimalityTolerance = 0;
OPTIONS.FunctionTolerance = 0;

% go through each iteration
for iter = 1:9

    % solve the unconstrained optimization problem
    x = fminunc(@(x) beta(x,mu),x,OPTIONS);

    % display stuff
    disp_helper("--- iteration",iter,[])
    disp_helper("mu",mu,10)
    disp_helper("x",x,10)

    % update mu
    mu = mu/10;

end

%% helper functions
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