close all; clear; clc

% data
A = [1 -1 0 0; 0 0 1 1];
b = [2; 2];

% orthonormal basis for the null space of A
% Z = null(A)
% "rational" basis for the null space of A
Z = null(A,'r')

% feasible point
xbar = [3; 1; 0; 2]

% check that it is feasible
A*xbar - b

% add null space
v = rand(2,1)*100; % random point
xnew = xbar + Z*v

% check that it is still feasible
A*xnew - b