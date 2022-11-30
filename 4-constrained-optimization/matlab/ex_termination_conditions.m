close all; clear; clc

% fmincon options and tolerances
OPTIONS = optimoptions('fmincon');
OPTIONS.Display = 'iter';
OPTIONS.OptimalityTolerance = 1e-8;
OPTIONS.ConstraintTolerance = 1e-8;
OPTIONS.StepTolerance = 1e-8;

% problem setup
fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
x0 = [1;0];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @circlecon;

% solve the problem
[x,FVAL,EXITFLAG,OUTPUT] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,OPTIONS)

% nonlinear constraint definition
function [c,ceq] = circlecon(x)
ceq = 3*x(1)^2 + x(2)^2 - 1;
c = [];
end