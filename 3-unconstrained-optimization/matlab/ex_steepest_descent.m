% ex_steepest_descent.m
% examples of the steepest-descent method on quadratic functions
% [reference] Section 12.2 in LNO
% [course] Session 8 - Unconstrained Optimization (2)
close all; clear; clc

example = 4;

switch example
    %----------------------------------------------------------------------
    case 1 % Example 12.1 in LNO
        Q = diag([1 5 25]);
        c = [-1;-1;-1];
        x0 = [0;0;0];
        n = 217;
    %----------------------------------------------------------------------
    case 2 % cond(Q) is large
        d = [1 25];
        s = 3535;
        Q = create_rand_pd_matrix(d,s);
        c = [-1;-1];
        x0 = [0;0];
        n = 25;
    %----------------------------------------------------------------------
    case 3 % cond(Q) is smallish
        d = [1 2];
        s = 3243;
        Q = create_rand_pd_matrix(d,s);
        c = [-1;-1];
        x0 = [0;0];
        n = 20;
    %----------------------------------------------------------------------
    case 4 % cond(Q) = 1
        Q = diag([5 5]);
        c = [-1;-1];
        x0 = [0;0];
        n = 5;
end

% optimal solution
x_opt = Q\c;

% check if 2D problem that can be plotted
if length(Q) == 2
    n2flag = true;
else
    n2flag = false;
end

% initialize figure (2D problem only)
plot_2d([],[],[],[],[],[],n2flag,'init')

% assign initial point
x = x0;

% initialize function value storage
f = nan(1,n);

% go through each iteration
for k = 1:n

    % store old value
    xold = x;

    % compute current function value
    if k < 25
        f(k) = x'*Q*x/2 - c'*x;
    end

    % steepest-descent direction
    p = -(Q*x - c);

    % norm of the gradient
    norm_g = norm(Q*x - c);
    disp(strcat("norm(g) = ",string(norm_g)))

    % exact line search
    alpha = (p'*p)/(p'*Q*p);

    % next point
    x = x + alpha*p;

    % plot iteration (2D problem only)
    plot_2d(x,xold,x_opt,f,Q,c,n2flag,'iter')

end

% plot contours and optimal point (2D problem only)
plot_2d(x,xold,x_opt,f,Q,c,n2flag,'final')

%--------------------------------------------------------------------------
% create a random positive-definite matrix with given eigenvalues
function A = create_rand_pd_matrix(d,s)

% set random seed
if ~isempty(s)
    rng(s)
end

% size of the matrix
n = length(d);

% get a random rotation matrix
[Q, ~] = qr(rand(n));

% apply the transformation
A = Q * diag(d) * transpose(Q);

end

%--------------------------------------------------------------------------
% plotting functions
% (not the main content)
function plot_2d(x,xold,x_opt,f,Q,c,n2flag,tasknum)

% only do stuff if a 2D problem
if ~n2flag
    return
end

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

    % labels
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');

    % axis properties
    ha = gca; ha.XColor = 'k'; ha.YColor = 'k'; ha.Color = 'none';
    ha.LineWidth = 1; ha.FontSize = FontSize;
    axis equal

    %----------------------------------------------------------------------
    case 'iter'

    % plot line between old and new point
    plot([xold(1) x(1)],[xold(2) x(2)],'-',plotOpts{:},'Color',niceblue);


    % plot current point
    plot(xold(1),xold(2),'.',plotOpts{:},'Color',nicered);

    %----------------------------------------------------------------------
    case 'final'

    % create a good set of contour lines
    f(isnan(f)) = [];
    fmax = 15;
    f_ = linspace(max(f),fmax,40);
    f = unique(sort([f,f_]));

    % create grid
    N = 1000;
    x = linspace(-1.5,0.5,N);
    y = linspace(-1,0.5,N);
    [X,Y] = meshgrid(x,y);

    % evaluate quadratic function
    C = zeros(size(X));
    for k = 1:numel(X)
        x = [X(k);Y(k)];
        C(k) = x'*Q*x/2 - c'*x;
    end

    % create contour plot
    [~,hc] = contour(X,Y,C,f);
    hc.LineWidth = LineWidth;
    colormap(bone)

    % plot optimal point
    plot(x_opt(1),x_opt(2),'.',plotOpts{:},'Color',nicegreen);

    % axis equal to visualize condition number
    axis equal

end

end