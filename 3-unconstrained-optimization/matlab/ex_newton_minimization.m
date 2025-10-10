% ex_newton_minimization.m
% illustration of basic Newton's methods for a 1D function
% [reference] Section 11.3 in LNO
% [course] Session 7 - Unconstrained Optimization (1)
close all; clear; clc

example_number = 3;

switch example_number
    case 1
        xk = 1; % starting point
        xlimits = [-2.5 2.5]; ylimits = [-0.1 3]; % plotting limits
        niter = 6; % number of iterations
    case 2
        xk = -1.9; % starting point
        xlimits = [-3 3]; ylimits = [-0.3 6]; % plotting limits
        niter = 7; % number of iterations
    case 3
        xk = -2.8; % starting point
        xlimits = [-6 -0.5]; ylimits = [-0.05 0.55]; % plotting limits
        niter = 3; % number of iterations
end

% symbolic stuff (function, gradient, hessian)
syms x
f = exp(0.5*x-1)*(x+1)^2;
g = gradient(f);
h = hessian(f);

% create matlab functions
F = matlabFunction(f);
G = matlabFunction(g);
H = matlabFunction(h);

% quadratic approximation of f at point x
q = @(x,p) F(x) + G(x)*p + 1/2*p'*H(x)*p;

% create plot
plot_1d(xk,F,q,[],0,xlimits,ylimits,niter,'init')

% go through each iteration
for iter = 1:niter

    % compute step length
    p = newton_step(xk,G,H);

    % plotting
    plot_1d(xk,F,q,p,iter,xlimits,ylimits,niter,'iter-1')
    pause % so you can see the steps

    % take step
    xk = xk + p;
    disp(strcat("x",string(iter)," = ",string(xk)))

    % plotting
    plot_1d(xk,F,q,p,iter,xlimits,ylimits,niter,'iter-2')
    pause % so you can see the steps

end

% display final point
disp(strcat("x = ",string(xk)))

return

%--------------------------------------------------------------------------
% Newton's method step
function p = newton_step(x,G,H)

% solve Newton equations
p = -H(x)\G(x);

end

%--------------------------------------------------------------------------
% plotting functions
% (not the main content)
function plot_1d(xk,F,q,p,iter,xlimits,ylimits,niter,tasknum)

% colors and other parameters
nicered = [225, 86, 86]/255;
LineWidth = 1.5;
MarkerSize = 24;
FontSize = 12;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};
N = 1e5; % number of points to plot
xmin = -10; xmax = 10; % x limits
X = linspace(xmin,xmax,N)'; % vector of x points
pmin = -1; pmax = 1; % p limits
colors = parula(niter); % colors for plotting

switch tasknum
    %----------------------------------------------------------------------
    case 'init' % initial point

    % initialize figure
    hf = figure; hf.Color = 'w'; hold on

    % function values
    plot(X,F(X),'-',plotOpts{:},'Color','k') % plot f(x)
    plot(xk,F(xk),'.',plotOpts{:},'Color',nicered) % plot initial point

    % change limits
    xlim(xlimits);
    ylim(ylimits);

    % labels
    xlabel('$x$','Interpreter','latex');
    ylabel('$f(x)$','Interpreter','latex');

    % axis properties
    ha = gca; ha.XColor = 'k'; ha.YColor = 'k'; ha.Color = 'none';
    ha.LineWidth = 1; ha.FontSize = FontSize;

    %----------------------------------------------------------------------
    case 'iter-1' % plot quadratic approximation

    P = linspace(min(0,p)+pmin,max(0,p)+pmax,N);
    hp = plot(xk+P,eval_q(xk,P,q),'-',plotOpts{:},'color',colors(iter,:));
    hp.Color = [hp.Color 0.75]; % change transparency

    %----------------------------------------------------------------------
    case 'iter-2' % plot step and function value

    hp = plot(xk,q(xk-p,p),'.',plotOpts{:},'color',colors(iter,:));
    hp.Color = [hp.Color 0.75]; % change transparency
    hp = plot([xk,xk],[q(xk-p,p) F(xk)],'-',plotOpts{:},'color',colors(iter,:));
    hp.Color = [hp.Color 0.35]; % change transparency
    hp = plot(xk,F(xk),'.',plotOpts{:},'color',colors(iter,:));
    hp.Color = [hp.Color 0.75]; % change transparency

end

end

%--------------------------------------------------------------------------
% evaluate q with vector-valued P
function Q = eval_q(x,P,q)

% initialize
Q = zeros(length(P),1);

% go through each point
for k = 1:length(P)
    Q(k) = q(x,P(k)); % evaluate the quadratic function
end

end