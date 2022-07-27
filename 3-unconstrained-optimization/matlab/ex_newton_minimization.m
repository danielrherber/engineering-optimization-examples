close all; clear; clc

test = 1;

switch test
    case 1
        xk = 1; % starting point
        xlimits = [-2.5 2.5]; ylimits = [-0.1 3]; % plotting limits
        niter = 4; % number of iterations
    case 2
        xk = -1.9; % starting point
        xlimits = [-3 3]; ylimits = [-0.3 6]; % plotting limits
        niter = 5; % number of iterations
    case 3
        xk = -3; % starting point
        xlimits = [-6 -0.5]; ylimits = [-0.05 0.5]; % plotting limits
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

% plot setup
N = 1e5; % number of points to plot
xmin = -10; xmax = 10; % x limits
X = linspace(xmin,xmax,N)'; % vector of x points
pmin = -1; pmax = 1; % p limits
colors = parula(niter); % colors for plotting

% create plot
set(0,'defaultTextInterpreter','latex');
hf = figure; hf.Color = 'w'; hold on
xlabel('$x$'); ylabel('$f(x)$'); % label axes
xlim(xlimits); ylim(ylimits); % change limits
ha = gca; ha.FontSize = 24;
set(gca,'TickLabelInterpreter','latex')
plot(X,F(X),'k-','linewidth',3) % plot f(x)
plot(xk,F(xk),'.r','markersize',30) % plot initial point

% go through each iteration
for iter = 1:niter

% compute step length
p = newton_step(xk,G,H);

% plotting
P = linspace(min(0,p)+pmin,max(0,p)+pmax,N);
hp = plot(xk+P,eval_q(xk,P,q),'r-','linewidth',2,'color',colors(iter,:));
hp.Color = [hp.Color 0.75]; % change transparency

% take step
xk = xk + p;

% plotting
pause
hp = plot(xk,q(xk-p,p),'.','markersize',20,'color',colors(iter,:));
hp.Color = [hp.Color 0.75]; % change transparency
hp = plot([xk,xk],[q(xk-p,p) F(xk)],'-','color',colors(iter,:),'linewidth',1);
hp.Color = [hp.Color 0.35]; % change transparency
hp = plot(xk,F(xk),'.','markersize',20,'color',colors(iter,:));
hp.Color = [hp.Color 0.75]; % change transparency

% pause so that you can see the steps
pause

end

% display final point
disp(strcat("x = ",string(xk)))

return

% Newton's method step
function p = newton_step(x,G,H)
p = -H(x)\G(x); % solve Newton equations
end

% evaluate q with vector-valued P
function Q = eval_q(x,P,q)

% initialize
Q = zeros(length(P),1);

% go through each point
for k = 1:length(P)
    Q(k) = q(x,P(k)); % evaluate the quadratic function
end

end