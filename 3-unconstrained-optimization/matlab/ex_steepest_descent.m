close all; clear; clc

test = 2;

switch test
    case 1 % example from the book
        Q = diag([1 5 25]);
        c = [-1;-1;-1];
        x0 = [0;0;0];
        n = 217;
    case 2 % cond(Q) is large
        d = [1 25];
        s = 3535;
        Q = create_rand_pd_matrix(d,s);
        c = [-1;-1];
        x0 = [0;0];
        n = 20;
    case 3 % cond(Q) is smallish
        d = [1 2];
        s = 3243;
        Q = create_rand_pd_matrix(d,s);
        c = [-1;-1];
        x0 = [0;0];
        n = 20;
    case 4 % cond(Q) = 1
        Q = diag([5 5]);
        c = [-1;-1];
        x0 = [0;0];
        n = 5;
end

% optimal solution
x_opt = Q\c;

% check if this is a 2d problem for plotting
if length(Q) == 2
    n2flag = true;
else
    n2flag = false;
end

% setup figure
if n2flag
    hf = figure; hf.Color = 'w'; hold on
    xlabel('x')
    ylabel('y')
    axis equal
    ha = gca;
    ha.FontSize = 20;
end

% assign initial point
x = x0;

% initialize
f = nan(1,n);

% go through each iteration
for k = 1:n

    % store old value
    xold = x;

    % compute current function value
    if k < 20
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

    % plot if it is a two dimensional problem
    if n2flag

        % plot line between old and new point
        hp = plot([xold(1) x(1)],[xold(2) x(2)],'-b');
        hp.LineWidth = 1;

        % plot current point
        hp = plot(xold(1),xold(2),'.r');
        hp.MarkerSize = 12;

    end

end

% plot if it is a two dimensional problem
if n2flag

    % create a good set of contour lines
    f(isnan(f)) = [];
    fmax = 15;
    f_ = linspace(max(f),fmax,30);
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

    % labels
    xlabel('x'); ylabel('y');

    % plot optimal point
    hp = plot(x_opt(1),x_opt(2),'.g');
    hp.MarkerSize = 12;

    % axis equal to visualize condition number
    axis equal

end

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