% ex_newton.m
% illustration of Newton's method for solving f(x) = 0 for several 1-d and
% 2-d examples
% [reference] Section 2.7 in LNO
% [course] Session 3 - Fundamental Concepts in Optimization (2)
close all; clear; clc

% select the example number (see cases below)
example_number = 1;

switch example_number
    case 1 % one-dimensional example
        f = @(x) [7 3 2 9 4]*[x^4 x^3 x^2 x^1 x^0]';
        fd = @(x) [28 9 4 9]*[x^3 x^2 x^1 x^0]';
        X(1) = 0; % initial point
        n = 8; % number of iterations
        xmin = -1.25; % minimum x value for plotting
        xmax = 0.25; % maximum x value for plotting
    case 2 % one-dimensional example when f'(x*) = 0
        f = @(x) [1 -7 17 -17 6]*[x^4 x^3 x^2 x^1 x^0]';
        fd = @(x) [4 -21 34 -17]*[x^3 x^2 x^1 x^0]';
        X(1) = 1.1; % initial point (try 1.1 and 2.1)
        n = 25; % number of iterations
        xmin = 0.5; % minimum x value for plotting
        xmax = 3; % maximum x value for plotting
    case 3 % failure of Newton's method
        f = @(x) (exp(x) - exp(-x))/(exp(x) + exp(-x));
        fd = @(x) 4*exp(2*x)/(exp(2*x) + 1)^2;
        X(1) = 1.1; % initial point (try 1 and 1.1)
        n = 10; % number of iterations
        xmin = -6; % minimum x value for plotting
        xmax = 6; % maximum x value for plotting
    case 4 % multidimensional case
        f = @(x) [3*x(1)*x(2) + 7*x(1) + 2*x(2) - 3;...
            5*x(1)*x(2) - 9*x(1) - 4*x(2) + 6];
        fd = @(x) [3*x(2) + 7, 3*x(1) + 2; 5*x(2) - 9, 5*x(1) - 4];
        X(:,1) = [1;2]; % initial point
        n = 10; % number of iterations
        xmin = -6; % minimum x value for plotting
        xmax = 6; % maximum x value for plotting
end

% create plot (see below for the functions)
if isscalar(X(:,1))
    plot_start_1d(X,f,xmin,xmax) % 1-d case
else
    plot_start_2d(X,f,xmin,xmax) % 2-d case
end

% go through each iteration
for k = 1:n

    % 1-d case (but not needed as newton_iteration also handles 1-d case)
    % X(:,k+1) = newton_iteration_1d(X(:,k),f,fd);

    % compute Newton iteration
    X(:,k+1) = newton_iteration(X(:,k),f,fd);

    % extract current value
    x_ = X(:,k+1);

    % plot iteration (see below for the functions)
    if isscalar(x_)
        plot_1d(x_,f,k,X(:,k),fd) % 1-d case
    else
        plot_2d(x_,f,k) % 2-d case
    end

end

%--------------------------------------------------------------------------
% 1-d Newton's method
function xk1 = newton_iteration_1d(xk,f,fd)

% Newton iteration
xk1 = xk - f(xk)/fd(xk);

end

%--------------------------------------------------------------------------
% multidimensional Newton's method
function xk1 = newton_iteration(xk,f,fd)

% compute step size from linear system: A*p = -b
A = fd(xk);
b = f(xk);
p = -A\b;

% Newton iteration
xk1 = xk + p;

end

%--------------------------------------------------------------------------
% plotting functions for 1d case
function plot_1d(x,f,k,xold,fd)

% plot tangent line
X = linspace(xold,x,2);
plot(X,f(xold) + fd(xold)*(X-xold),'m--','linewidth',1.5)
plot([x x],[0 f(x)],'k--','linewidth',1.5)

% plot next point
plot(x,f(x),'r.','markersize',24)

% display stuff
disp(strcat(string(k)," x:",string(vpa(x,16))," f(x):",string(vpa(f(x),16))))

end

%--------------------------------------------------------------------------
function plot_start_1d(X,f,xmin,xmax)

set(0,'defaultTextInterpreter','latex');
hf = figure; hold on
hf.Color = 'w';
ha = gca; ha.FontSize = 18;
xlim([xmin xmax])
ha.LineWidth = 1;
xlabel('$x$'); ylabel('$f(x)$'); % label axes

Xg = linspace(xmin,xmax,1e6);
for k = 1:length(Xg)
    Fg(k) = f(Xg(k));
end
plot([xmin xmax],[0 0],'b-','linewidth',1.5) % desired value
plot(Xg,Fg,'k','linewidth',1.5) % function value
plot(X(:,1),f(X(:,1)),'g.','markersize',24) % initial point

end

%--------------------------------------------------------------------------
% plotting functions for 2d case
function plot_2d(x,f,k)

% plot current point
plot3(x(1),x(2),norm(f(x)),'r.','markersize',24)

% display stuff
disp(strcat(string(k)," x1:",string(vpa(x(1),24)),...
    " x2:",string(vpa(x(2),16))," e:",string(vpa(norm(f(x)),16))))

end

%--------------------------------------------------------------------------
function plot_start_2d(X,f,xmin,xmax)

hf = figure; hold on
hf.Color = 'w';
ha = gca;
ha.LineWidth = 1;
ha.FontSize = 14;
xlabel('x_1')
ylabel('x_2')
zlabel('log(norm(f(x)))')

x = linspace(xmin,xmax,300);
[X1,X2] = meshgrid(x,x);
E = nan(size(X1));
for k = 1:numel(X1)
    E(k) = norm(f([X1(k),X2(k)]));
end

contourf(X1,X2,log(E),40) % plot errors
plot(X(1),X(2),'g.','markersize',24) % initial point

end