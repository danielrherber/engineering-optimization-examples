% ex_production_planning_integer.m
% example of a linear program (LP) representing a production planning
% optimization problem
% solved both as a continuous problem and integer-valued one
% [reference] Section 1.4 in LNO
% [course] Session 14 - Additional Topics
close all; clear; clc

% problem definition
f = [100 150 200 400]; % profit coefficients
a1 = [10 12 25 20]; % wood resource coefficients
b1 = 5000; % amount of wood resource
a2 = [2 4 8 12]; % labor resource coefficients
b2 = 1500; % amount of labor resource
LB = [0 0 0 0]; % nonnegativity constraint

% solve the optimization problem (continuous)
OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';
[X,F] = linprog(-f,[a1;a2],[b1;b2],[],[],LB,[],OPTIONS);

% solve the optimization problem (integer)
intcon = 1:4; % all variables are required to be integers
OPTIONS = optimoptions('intlinprog');
OPTIONS.Display = 'iter';
[Xint,Fint] = intlinprog(-f,intcon,[a1;a2],[b1;b2],[],[],LB,[],[],OPTIONS);

% display the optimal solutions
names = ["bookshelves";"cabinets with doors";"tall cabinets with doors";"fancy cabinets"];
disp("--- linprog ---")
disp(strcat("profit = ",string(-F)))
disp(strcat(names,"=",string(X)))

disp("--- intlinprog ---")
disp(strcat("profit = ",string(-Fint)))
disp(strcat(names,"=",string(Xint)))