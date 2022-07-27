close all; clear; clc

% test number (see below)
test = 2;

switch test
    case 1 % Session 5 Class Example
        A = [-2 1 1 0 0; -1 2 0 1 0; 1 0 0 0 1];
        b = [2; 7; 3];
        c = [-1; -2; 0; 0; 0];
        x = [0; 0; 2; 7; 3];
    case 2 % Example for the Two-phase Method Continued (Phase 1)
        A = [3 2 0 0 1 0; 2 -4 -1 0 0 1; 4 3 0 1 0 0];
        b = [14 2 19]';
        c = [0 0 0 0 1 1]';
        x = [0 0 0 19 14 2]';
    case 3 % Example for the Two-phase Method Continued (Phase 2)
        A = [3 2 0 0; 2 -4 -1 0; 4 3 0 1];
        b = [14 2 19]';
        c = [2 3 0 0]';
        x = [4 1 2 0]';
    case 4
        A = [2 0 0 1 0; 1 1 2 0 1];
        b = [4 2]';
        c = [0 0 0 1 1]';
        x = [0 0 0 4 2]';
end

% initialize as not optimal and is not unbounded
isOptimal = false; isUnbounded = false;

% continue iterating until optimal point is found
while ~isOptimal && ~isUnbounded

[x,isOptimal,isUnbounded] = simplex_iteration(A,b,c,x);

end

%% check with linprog
X = linprog(c,[],[],A,b,zeros(length(x),1),[],[],[]);
disp("x with linprog = ")
disp(X)

% simplex iteration
% NOTE: this assumes that you start at a basic feasible solution
function [x,isOptimal,isUnbounded] = simplex_iteration(A,b,c,x)

%% initial stuff
% display current value of x (optional)
disp("x = ")
disp(x)

% initial values
isUnbounded = false;
isOptimal = false;

% number of variables, basic variables, and nonbasic variables
n = size(A,2);
nb = size(A,1);
nn = n - nb;

% determine indices nonbasic variables (smallest values, equal to 0)
[~,In] = mink(x,nn); % NOTE: this is not generally a good approach because of degeneracy

% determine basic variable indices in x
Ib = (1:n)';
Ib(In) = [];

% extract N and B matrices
N = A(:,In);
B = A; B(:,In) = [];

%% optimality test
% extract cN and cB
cN = c(In);
cB = c; cB(In) = [];

% compute simplex multipliers
y = (B'\cB);

% compute reduced cost vector
chat = cN - N'*y;

% check if current basis is optimal
if all(chat >= 0)
    disp("optimal!")
    isOptimal = true;
    return
end

% select entering variable
[~,It_] = min(chat);
It = In(It_); % get actual variable index

%% step
% compute constraint coefficients corresponding to the entering variable
At = B\A(:,It);
bt = B\b;

% compute ratios
Alpha = bt./At;

% replace negative numbers with infinite step size
Alpha(Alpha<0) = inf;

% check if the problem is unbounded
if all(isinf(Alpha))
    disp("unbounded!")
    isUnbounded = true;
    return
end

% find smallest positive number for value of the entering variable
xt = min(Alpha);

%% update
% compute null space matrix (using canonical ordering)
Z = [eye(nn);-B\N];

% restore Z back to original ordering
Z([In;Ib],:) = Z;

% step using correct column of Z
x = x + xt*Z(:,It_);

end