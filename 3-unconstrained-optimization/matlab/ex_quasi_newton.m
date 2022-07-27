close all; clear; clc

test = 1;

updateformulaflag = 1;

switch test

    case 1 % Example 12.9 from LNO
        Q = diag([2,3,4]); % x'*Q*x/2
        c = [-8;-9;-8]; % c'*x
        B = eye(length(Q)); % initial approximate hessian matrix
        x = zeros(size(Q,1),1); % initial point
    case 2
        rng(3253) % set random seed
        q = rand(3);
        Q = q+q'; % x'*Q*x/2
        c = [-8;-9;-8]; % c'*x
        B = eye(length(Q)); % initial approximate hessian matrix
        x = zeros(size(Q,1),1); % initial point
    case 3
        rng(5754) % set random seed
        n = 20;
        q = rand(n);
        Q = q + q'+ 8*eye(n); % x'*Q*x/2
        c = rand(n,1); % c'*x
        B = eye(length(Q)); % initial approximate hessian matrix
        x = zeros(size(Q,1),1); % initial point
end

% exact gradient
G = @(x) Q*x - c;

% extract hessian
H = @(x) Q;

% maximum number of iterations
max_iterations = 100;

% compute initial gradient
g = G(x);

% go through each iteration
for iteration = 1:max_iterations

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

    % (iii) using exact line search
    alpha = -(p'*g)/(p'*Q*p);
    disp_helper("alpha",alpha)

    x_new = x + alpha*p;
    disp_helper("x",x_new)

    % (iv) compute intermediate quantities
    g_new = G(x_new); % gradient
    s = x_new - x; % step
    y = g_new - g; % gradient difference

    disp_helper("g",g_new)
    disp_helper("s",s)
    disp_helper("y",y)

    % (v) use symmetric rank-one update formula
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

% compute error compared to the "exact" solution
e_norm = norm(x - (Q\c),inf);
disp_helper("max abs error",e_norm)

% compute error compared to the exact hessian
e_Q = norm(B-Q,inf);
disp_helper("error Q",e_Q)

% function to make it easier to display things in the command window
function disp_helper(name,number)

% form string
str = strcat(string(name)," = ",mat2str(round(number,16)));

% display string
disp(str)

end