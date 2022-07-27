close all; clear; clc

method = 1; % see cases below

% create objective function
syms x y real
fs = x^4 + y^4 - 4*x^3 - 3*y^3 + 2*x^2 + 2*x*y;
obj = matlabFunction(fs);

% bounds
ub = [4 3]; % upper bounds
lb = [-2 -2]; % lower bounds

% number of variables
n = 2;

% solve the problem based on the selected method
switch method
    %----------------------------------------------------------------------
    case 1 % pattern search

    % inputs
    FUN = @(x) obj(x(1),x(2));
    X0 = [1;1];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    LB = lb;
    UB = ub;
    NONLCON = [];

    % options
    OPTIONS = optimoptions('patternsearch');
    OPTIONS.Display = 'iter';
    OPTIONS.PlotFcns = {'psplotbestf','psplotmeshsize','psplotbestx'};

    % solve the problem
    x = patternsearch(FUN,X0,A,b,Aeq,beq,LB,UB,NONLCON,OPTIONS);

    %----------------------------------------------------------------------
    case 2 % genetic algorithm

    % inputs
    FITNESSFCN = @(x) obj(x(1),x(2));
    NVARS = n;
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    LB = lb;
    UB = ub;
    NONLCON = [];

    % options
    OPTIONS = optimoptions('ga');
    OPTIONS.Display = 'iter';
    OPTIONS.PlotFcns = {'gaplotdistance','gaplotgenealogy',...
        'gaplotselection','gaplotbestf','gaplotexpectation','gaplotrange'};
    OPTIONS.PopulationSize = 30;

    % solve the problem
    x = ga(FITNESSFCN,NVARS,A,b,Aeq,beq,LB,UB,NONLCON,OPTIONS);

    %----------------------------------------------------------------------
    case 3 % particle swarm optimization

    % inputs
    FUN = @(x) obj(x(1),x(2));
    NVARS = n;
    LB = lb;
    UB = ub;

    % options
    OPTIONS = optimoptions('particleswarm');
    OPTIONS.Display = 'iter';
    OPTIONS.SwarmSize = 10;

    % solve the problem
    x = particleswarm(FUN,NVARS,LB,UB,OPTIONS);

end