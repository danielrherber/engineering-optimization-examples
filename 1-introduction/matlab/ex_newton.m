close all; clear; clc

% select the example number (see cases below)
example_number = 3;

switch example_number
    case 1 % One-dimensional Example
        f = @(x) [7 3 2 9 4]*[x^4 x^3 x^2 x^1 x^0]';
        fd = @(x) [28 9 4 9]*[x^3 x^2 x^1 x^0]';
        X(1) = 0; % initial point
        n = 8; % number of iterations
        xmin = -1.5; % minimum x value for plotting
        xmax = 0.5; % maximum x value for plotting
    case 2 % One-dimensional Example when f'(x*) = 0
        f = @(x) [1 -7 17 -17 6]*[x^4 x^3 x^2 x^1 x^0]';
        fd = @(x) [4 -21 34 -17]*[x^3 x^2 x^1 x^0]';
        X(1) = 1.1; % initial point (try 1.1 and 2.1)
        n = 25; % number of iterations
        xmin = 0.5; % minimum x value for plotting
        xmax = 3; % maximum x value for plotting
    case 3 % Multidimensional Case
        f = @(x) [3*x(1)*x(2) + 7*x(1) + 2*x(2) - 3;...
            5*x(1)*x(2) - 9*x(1) - 4*x(2) + 6];
        fd = @(x) [3*x(2) + 7, 3*x(1) + 2; 5*x(2) - 9, 5*x(1) - 4];
        X(:,1) = [1;2]; % initial point
        n = 10; % number of iterations
        xmin = -6; % minimum x value for plotting
        xmax = 6; % maximum x value for plotting
    case 4 % Failure of Newton’s Method
        f = @(x) (exp(x) - exp(-x))/(exp(x) + exp(-x));
        fd = @(x) 4*exp(2*x)/(exp(2*x) + 1)^2;
        X(1) = 1.1; % initial point (try 1 and 1.1)
        n = 10; % number of iterations
        xmin = -6; % minimum x value for plotting
        xmax = 6; % maximum x value for plotting
end

% create plot (see below)
if isscalar(X(:,1)) % 1-d case
    plot_start_1d(X,f,xmin,xmax)
else
    plot_start_2d(X,f,xmin,xmax)
end

% go through each iteration
for k = 1:n

    % 1d case (but not needed as newton_iteration also handled 1d case)
    % X(:,k+1) = newton_iteration_1d(X(:,k),f,fd);

    % compute Newton iteration
    X(:,k+1) = newton_iteration(X(:,k),f,fd);

    % extract current value
    x_ = X(:,k+1);

    % check if 1-d or 2-d case
    if isscalar(x_) % 1-d case

        % plot current point
        plot(x_,f(x_),'r.','markersize',16)

        % display stuff
        disp(strcat(string(k)," x:",...
            string(vpa(x_,16))," f(x):",string(vpa(f(x_),16))))
    else % 2-d case

        hold on
        % plot current point
        plot3(x_(1),x_(2),norm(f(x_)),'r.','markersize',16)

        % display stuff
        disp(strcat(string(k)," x1:",string(vpa(x_(1),16)),...
            " x2:",string(vpa(x_(2),16))," e:",string(vpa(norm(f(x_)),16))))

    end

end

% 1-d Newton's method
function xk1 = newton_iteration_1d(xk,f,fd)

% Newton iteration
xk1 = xk - f(xk)/fd(xk);

end

% n-d Newton's method
function xk1 = newton_iteration(xk,f,fd)

% compute step size from linear system: A*p = -b
A = fd(xk);
b = f(xk);
p = -A\b;

% Newton iteration
xk1 = xk + p;

end

% plotting function for 1d case
function plot_start_1d(X,f,xmin,xmax)

hf = figure; hold on
hf.Color = 'w';
ha = gca;
ha.LineWidth = 1;
ha.FontSize = 14;
xlabel('x')
ylabel('f(x)')

Xg = linspace(xmin,xmax,1e6);
for k = 1:length(Xg)
    Fg(k) = f(Xg(k));
end
plot([xmin xmax],[0 0],'b-','linewidth',1.5) % desired value
plot(Xg,Fg,'k','linewidth',1.5) % function value
plot(X(:,1),f(X(:,1)),'g.','markersize',16) % initial point

end

% plotting function for 2d case
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
plot(X(1),X(2),'g.','markersize',16) % initial point

end