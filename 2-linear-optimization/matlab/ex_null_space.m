% ex_null_space.m
% example of the null-space matrix and moving in feasible directions for a
% linear constraints Ax = b
% [reference] Section 3.2 in LNO
% [course] Session 4 - Linear Optimization (1)
close all; clear; clc

%--------------------------------------------------------------------------
% Example of the Basis Matrix for the Null Space
% data
A = [1 -1 0 0; 0 0 1 1];
b = [2; 2];

% feasible point
xbar = [3; 1; 0; 2]

% orthonormal basis for the null space of A, Z'*Z = I
% Z = null(A)
% "rational" basis for the null space of A
Z = null(A,'r')

% check that the point is feasible
A*xbar - b

% add a null-space component
v = rand(2,1)*100; % random point in R2
xnew = xbar + Z*v

% check that it is still feasible
A*xnew - b

%--------------------------------------------------------------------------
% Relationship between Null and Range Space
% solve a linear system for null-space and range-space component of xbar
pq = [Z,A']\xbar;

% compute null-space and range-space components
p = Z*pq(1:2)
q = A'*pq(3:4)

% confirm that p + q = xbar
isequal(p+q,xbar)

% confirm the vectors are orthogonal
isequal(p'*q,0)