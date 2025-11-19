% ex_newton.m
% illustration of Newton's method for solving f(x) = 0 for several 1D and
% 2D examples
% [reference] Section 2.7 in LNO
% [course] Session 3 - Fundamental Concepts in Optimization (2)
close all; clear; clc

% select the example number (see cases below)
example = 1;

switch example
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
        n = 5; % number of iterations
        xmin = -6; % minimum x value for plotting
        xmax = 6; % maximum x value for plotting
end

% create plot (see below for the functions)
if isscalar(X(:,1))
    plot_1d(X(:,1),f,0,[],[],X,xmin,xmax,'init') % 1D case
else
    plot_2d(X(:,1),f,0,X,xmin,xmax,'init') % 2D case
end

% go through each iteration
for k = 1:n

    % 1D case (but not needed as newton_iteration also handles 1D case)
    % X(:,k+1) = newton_iteration_1d(X(:,k),f,fd);

    % compute Newton iteration
    X(:,k+1) = newton_iteration(X(:,k),f,fd);

    % extract current value
    x_ = X(:,k+1);

    % plot iteration (see below for the functions)
    if isscalar(x_)
        plot_1d(x_,f,k,X(:,k),fd,[],[],[],'iter') % 1D case
    else
        plot_2d(x_,f,k,[],[],[],'iter') % 2D case
    end

end

%--------------------------------------------------------------------------
% 1D Newton's method
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
% plotting functions for 1D case
% (not the main content)
function plot_1d(x,f,k,xold,fd,X,xmin,xmax,tasknum)

% colors and other parameters
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;
LineWidth = 1.5;
MarkerSize = 24;
FontSize = 12;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

switch tasknum
    %----------------------------------------------------------------------
    case 'init' % initial point

    % initialize figure
    hf = figure; hf.Color = 'w'; hold on

    % function values
    Xg = linspace(xmin,xmax,1e6);
    for idx = 1:length(Xg)
        Fg(idx) = f(Xg(idx));
    end
    plot([xmin xmax],[0 0],'-',plotOpts{:},'Color',niceblue) % desired value
    plot(Xg,Fg,plotOpts{:},'Color','k') % function value
    plot(X(:,1),f(X(:,1)),'.',plotOpts{:},'Color',nicegreen) % initial point

    xlim([xmin xmax])

    xlabel('$x$','Interpreter','latex');
    ylabel('$f(x)$','Interpreter','latex');

    ha = gca; ha.XColor = 'k'; ha.YColor = 'k'; ha.Color = 'none';
    ha.LineWidth = 1; ha.FontSize = FontSize;

    %----------------------------------------------------------------------
    case 'iter' % plot current iteration

    % plot tangent line
    X = linspace(xold,x,2);
    plot(X,f(xold) + fd(xold)*(X-xold),'--',plotOpts{:},'Color',nicered)
    plot([x x],[0 f(x)],'--',plotOpts{:},'Color','k')

    % plot next point
    plot(x,f(x),'.',plotOpts{:},'Color',nicered)

end

% display stuff
disp(strcat(string(k)," x:",string(vpa(x,16))," f(x):",string(vpa(f(x),16))))

end

%--------------------------------------------------------------------------
% plotting functions for 2D case
% (not the main content)
function plot_2d(x,f,k,X,xmin,xmax,tasknum)

% colors and other parameters
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;
LineWidth = 1.5;
MarkerSize = 24;
FontSize = 12;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

switch tasknum
    %----------------------------------------------------------------------
    case 'init' % initial point

    % initialize figure
    hf = figure; hf.Color = 'w'; hold on

    % calculate norm of errors
    x_ = linspace(xmin,xmax,800);
    [X1,X2] = meshgrid(x_,x_);
    E = nan(size(X1));
    for idk = 1:numel(X1)
        E(idk) = norm(f([X1(idk),X2(idk)]));
    end

    % contourf(X1,X2,log(E),40) % plot errors with contourf
    hf = surf(X1,X2,log(E)); % plot errors with surf
    hf.LineStyle = 'none'; hf.FaceAlpha = 0.75; % style surf
    plot3(X(1),X(2),log(norm(f(X))),'.',plotOpts{:},'Color',nicegreen) % initial point

    xlabel('$x_1$','Interpreter','latex');
    ylabel('$x_2$','Interpreter','latex');
    zlabel('$\log(\|f(x)\|)$','Interpreter','latex');

    colorbar('Color','k','LineWidth',1);

    ha = gca; ha.XColor = 'k'; ha.YColor = 'k'; ha.ZColor = 'k'; ha.Color = 'none';
    ha.LineWidth = 1; ha.FontSize = FontSize;

    %----------------------------------------------------------------------
    case 'iter' % plot current iteration

    % plot current point
    plot3(x(1),x(2),log(norm(f(x))),'.',plotOpts{:},'Color',nicered)

end

% display stuff
disp(strcat(string(k)," x1:",string(vpa(x(1),24)),...
    " x2:",string(vpa(x(2),16))," e:",string(vpa(norm(f(x)),16))))

end