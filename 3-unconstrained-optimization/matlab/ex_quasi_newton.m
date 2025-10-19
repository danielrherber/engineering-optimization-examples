% ex_quasi_newton.m
% examples of the update formulas used in quasi-Newton methods including
% symmetric rank-one update and BFGS
% [reference] Section 12.3 in LNO
% [course] Session 8 - Unconstrained Optimization (2)
close all; clear; clc

example = 4;

updateformulaflag = 2; % 1: symmetric rank-one update, 2:BFGS

switch example
    %----------------------------------------------------------------------
    case 1 % quadratic function Example 12.9 in LNO
        quadraticFlag = true; % is this a quadratic function?
        Q = diag([2,3,4]); % x'*Q*x/2
        c = [-8;-9;-8]; % c'*x
        x0 = zeros(size(Q,1),1); % initial point
        B = eye(length(Q)); % initial approximate hessian matrix
    %----------------------------------------------------------------------
    case 2 % randomized 3D quadratic function example
        quadraticFlag = true; % is this a quadratic function?
        rng(3253) % set random seed
        q = rand(3);
        Q = q+q'; % x'*Q*x/2
        c = [-8;-9;-8]; % c'*x
        x0 = zeros(size(Q,1),1); % initial point
        B = eye(length(Q)); % initial approximate hessian matrix
    %----------------------------------------------------------------------
    case 3 % randomized nD quadratic function example
        quadraticFlag = true; % is this a quadratic function?
        rng(5754) % set random seed
        n = 20;
        q = rand(n);
        Q = q + q'+ 8*eye(n); % x'*Q*x/2
        c = rand(n,1); % c'*x
        x0 = zeros(n,1); % initial point
        B = eye(n); % initial approximate hessian matrix
    %----------------------------------------------------------------------
    case 4 % nonlinear 2D example
        quadraticFlag = false; % is this a quadratic function?
        n = 2; x = sym('x',[n 1]); % create symbolic variables
        f = x(1)^4 + 2*x(1)^3 + 24*x(1)^2 + x(2)^4 + 12*x(2)^2; % objective function
        B = eye(n); % initial approximate hessian matrix
        x0 = [2; 1]; % initial point
end

% create problem functions
if quadraticFlag
    G = @(x) Q*x - c; % exact gradient
    H = @(x) Q; % exact hessian
else
    g = gradient(f,x); % exact gradient
    h = hessian(f,x); % exact hessian
    F = matlabFunction(f,'Vars',{x}); % objective function
    G = matlabFunction(g,'Vars',{x}); % gradient
    H = matlabFunction(h,'Vars',{x}); % hessian
end

% maximum number of iterations
max_iterations = 100;

% compute initial gradient
g = G(x0);

% assign initial value
x = x0;

% go through each iteration
for iteration = 0:max_iterations-1

    % (i) check if xk is optimal
    if norm(g,inf) <= 100*(eps)
        disp('Optimal!')
        break % stop iterating
    else
        disp_helper("iteration",iteration)
        disp_helper("norm(g)",norm(g,inf))
    end

    % (ii) solve for the search direction
    p = -B\g;
    disp_helper("p",p)

    % form string
    str = strcat("B =",mat2str(round(B,3)));

    % display string
    disp(str)

    % (iii) determine step using line search
    if quadraticFlag % quadratic objective function exact line search
        alpha = -(p'*g)/(p'*Q*p);
    else % general case
        alpha = fminbnd(@(alpha) F(x+alpha*p),0,100,...
            optimset('TolX',1e-10,'TolFun',1e-10));
    end
    disp_helper("alpha",alpha)

    x_new = x + alpha*p; % take the step
    disp_helper("x",x_new)

    % (iv) compute intermediate quantities
    g_new = G(x_new); % gradient
    s = x_new - x; % step
    y = g_new - g; % gradient difference

    disp_helper("g",g_new)
    disp_helper("s",s)
    disp_helper("y",y)

    % (v) use a quasi-Newton update formula
    switch updateformulaflag
        case 1 % use symmetric rank-one update formula
            z = y - B*s;
            B_new = B + z*z'/(z'*s);
        case 2 % BFGS update formula
            z = B*s;
            B_new = B - z*z'/(z'*s) + y*y'/(y'*s);
    end
    disp_helper("B",B_new)

    % update stored variables
    B = B_new;
    g = g_new;
    x = x_new;
    disp(" ")

end

% compute error compared to the exact hessian
e_Q = norm(B-H(x),inf);
disp_helper("error Q",e_Q)

% compute error compared to the "exact" solution
if quadraticFlag
    e_norm = norm(x - (Q\c),inf);
    disp_helper("max abs error",e_norm)
end

%--------------------------------------------------------------------------
% function to make it easier to display things in the command window
function disp_helper(name,number)

% form string
str = strcat(string(name)," = ",mat2str(round(number,16)));

% display string
disp(str)

end