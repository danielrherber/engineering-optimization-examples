% ex_simplex.m
% naive simplex method implementation for nondegenerate linear problems
% that start with a basic feasible solution
% [reference] Section 5.2 in LNO
% [course] Session 5 - Linear Optimization (2)
% [reference] Sections 5.3.1, 5.4, and 5.5 in LNO
% [course] Session 6 - Linear Optimization (3)
close all; clear; clc

% test number (see below)
test = 1;

switch test
    case 1 % session 5 class example, Section 5.2 in LNO
        A = [-2 1 1 0 0; -1 2 0 1 0; 1 0 0 0 1];
        b = [2; 7; 3];
        c = [-1; -2; 0; 0; 0];
        basis = [3,4,5];
    case 2 % session 6 multiple solutions example, Section 5.3.1 in LNO
        A = [-2 1 1 0 0; -1 1 0 1 0; 1 0 0 0 1];
        b = [2 3 3]';
        c = [-1 0 0 0 0]';
        basis = [3,4,5];
    case 3 % session 6 two-phase method (phase 1), Section 5.4 in LNO
        A = [3 2 0 0 1 0; 2 -4 -1 0 0 1; 4 3 0 1 0 0];
        b = [14 2 19]';
        c = [0 0 0 0 1 1]';
        basis = [4,5,6];
    case 4 % session 6 two-phase method (phase 2), Section 5.4 in LNO
        A = [3 2 0 0; 2 -4 -1 0; 4 3 0 1];
        b = [14 2 19]';
        c = [2 3 0 0]';
        basis = [1,2,3];
    case 5 % session 6 cycling degeneracy example, Section 5.5 in LNO
        A = [1/4 -60 -1/25 9 1 0 0; 1/2 -90 -1/50 3 0 1 0; 0 0 1 0 0 0 1];
        b = [0 0 1]';
        c = [-3/4 150 -1/50 6 0 0 0]';
        basis = [5,6,7];
end

% initialize as not optimal and is not unbounded
isOptimal = false; isUnbounded = false;

% continue iterating until optimal point is found
while ~isOptimal && ~isUnbounded

    [x,basis,isOptimal,isUnbounded] = simplex_iteration(A,b,c,basis);

end

% check with linprog
options = optimoptions('linprog','Display','none');
X = linprog(c,[],[],A,b,zeros(length(x),1),[],options);
disp("x with linprog = ")
disp(X)

%--------------------------------------------------------------------------
% simplex iteration
% NOTE: this assumes that you start at a basic feasible solution
function [x,Ib,isOptimal,isUnbounded] = simplex_iteration(A,b,c,Ib)

%--- initial stuff
% initial values
isUnbounded = false;
isOptimal = false;

% number of variables, basic variables, and nonbasic variables
n = size(A,2);
nb = numel(Ib);
nn = n - nb;

% determine nonbasic variable indices in x
In = 1:n;
In(Ib) = [];

% extract N and B matrices
N = A(:,In);
B = A(:,Ib);

% determine current basic feasible solution (using canonical ordering)
xN = zeros(nn,1);
xB = B\b;
x = [xN;xB];

% restore x back to original ordering
x([In,Ib],:) = x;

% display current value of x and basis (optional)
disp(strcat("x = ",mat2str(x)))
disp(strcat("basis = ",mat2str(Ib)))

%--- (i) check optimality of xk
% extract cN and cB vectors
cN = c(In);
cB = c(Ib);

% compute simplex multipliers
y = (B'\cB);

% compute reduced cost vector
chat = cN - N'*y;

% check if current basis is optimal
if all(chat >= 0)
    disp("optimal!")
    isOptimal = true;
    return % stop
end

%--- (ii) select nonbasic variable satisfying chat < 0 as entering variable
% select entering variable
[~,Ie_] = min(chat); % most negative
Ie = In(Ie_); % get actual variable index

%--- (iii) determine xk+1, a new estimate of the solution
% compute constraint coefficients corresponding to the entering variable
ae = B\A(:,Ie);
be = B\b;

% compute ratios
Alphas = be./ae;

% replace negative numbers with infinite step size (never selected)
Alphas(Alphas<0) = inf;

% check if the problem is unbounded
if all(isinf(Alphas))
    disp("unbounded!")
    isUnbounded = true;
    return
end

% find smallest positive number for value of the entering variable and
% leaving variable index
[xe,Il_] = min(Alphas);
Il = Ib(Il_);

% compute null space matrix (using canonical ordering)
Z = [eye(nn);-B\N];

% restore Z back to original ordering
Z([In,Ib],:) = Z;

% step using correct column of Z
x = x + xe*Z(:,Ie_);

% remove leaving variable from basis
Ib(Ib==Il) = [];

% add entering variable to basis
Ib = [Ib,Ie];

end