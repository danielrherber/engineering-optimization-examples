% ex_production_planning.m
% example of a linear program (LP) representing a production planning
% optimization problem
% [reference] Section 1.4 in LNO
% [course] Session 1 - Introduction to Linear Optimization
close all; clear; clc

% problem definition
f = [100 150 200 400]; % profit coefficients
a1 = [10 12 25 20]; % wood resource coefficients
b1 = 5000; % amount of wood resource
a2 = [2 4 8 12]; % labor resource coefficients
b2 = 1500; % amount of labor resource
LB = [0 0 0 0]; % nonnegativity constraint

OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';

% solve the optimization problem
X = linprog(-f,[a1;a2],[b1;b2],[],[],LB,[],OPTIONS);

% display the optimal solution
names = ["bookshelves";"cabinets with doors";"tall cabinets with doors";"fancy cabinets"];
disp(strcat(names,"=",string(X)))