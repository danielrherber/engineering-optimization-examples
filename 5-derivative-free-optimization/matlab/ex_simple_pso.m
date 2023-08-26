% ex_simple_pso.m
% illustration of particle swarm optimization (PSO) on 2d problems
% [reference] Section 7.7 in EDO
% [course] Session 13 - Derivative-free Optimization (2)
close all; clear; clc

test = 2; % see cases below

switch test
    case 1
        fun = 1; % objective function (see below)
        max_iterations = 40; % maximum number of iterations
        np = 20; % population size (even number for simplicity)
        p_norm_max = 2; % maximum step size norm
        rng(2153) % set random seed for repeatable results
    case 2
        fun = 2; % objective function (see below)
        max_iterations = 40; % maximum number of iterations
        np = 10; % population size (even number for simplicity)
        p_norm_max = 0.5; % maximum step size norm
        rng(1561) % set random seed for repeatable results
    case 3
        fun = 3; % objective function (see below)
        max_iterations = 30; % maximum number of iterations
        np = 30; % population size (even number for simplicity) try 30 and 40
        p_norm_max = 2; % maximum step size norm
        rng(46346) % set random seed for repeatable results
end

% algorithm parameters
alpha = 0.7; % inertia parameter
beta_max = 1; % self influence parameter
gamma_max = 1; % social influence parameter

% create objective function
syms x y real
con = []; pen = []; n = 2;
switch fun
    case 1 % spring example
        l1 = 12; l2 = 8; k1 = 1; k2 = 10; mg = 7; % problem data
        fs =  k1/2*(sqrt((l1 + x).^2 + y.^2) - l1).^2 + ...
            k2/2*(sqrt((l2 - x).^2 + y.^2) - l2).^2 + ...
            -mg*y;
        ub = [16 12]; % upper bounds
        lb = [-6 -9]; % lower bounds
    case 2 % Jones function (multiple minima)
        fs = x^4 + y^4 - 4*x^3 - 3*y^3 + 2*x^2 + 2*x*y;
        ub = [4 3]; % upper bounds
        lb = [-2 -2]; % lower bounds
    case 3 % spring example with linear constraint
        l1 = 12; l2 = 8; k1 = 1; k2 = 10; mg = 7; % problem data
        fs =  k1/2*(sqrt((l1 + x).^2 + y.^2) - l1).^2 + ...
            k2/2*(sqrt((l2 - x).^2 + y.^2) - l2).^2 + ...
            -mg*y;
        gs = 2*x + y - 8;
        pen = 100*gs^2;
        fs = fs + pen;
        con = matlabFunction(gs);
        pen = matlabFunction(pen);
        ub = [16 12]; % upper bounds
        lb = [-6 -9]; % lower bounds
end

% create matlab functions
obj = matlabFunction(fs);

% determine initial swarm
x1_ = linspace(lb(1),ub(1),round(sqrt(np)));
x2_ = linspace(lb(2),ub(2),round(sqrt(np)));
[x1_,x2_] = meshgrid(x1_,x2_); % stratified sample
% x1_ = lb(1) + rand(np,1)*(ub(1)-lb(1)); % randomize sample
% x2_ = lb(2) + rand(np,1)*(ub(2)-lb(2));

% combine
x = [x1_(:),x2_(:)];

% update population size
np = size(x,1);

% randomize initial directions (velocities)
p = 0.1*(1-2*rand(np,2)).*(ub-lb);

% limit step based on maximum step
p = step_norm_limiter(p,p_norm_max);

% create initial plot
[hp,parent_color,colors] = plot_helper1(obj,pen,con,lb,ub,x,p,max_iterations);

% get optimal solution using fminunc so that it can be plotted
switch fun
    case 1
        options = optimoptions('fminunc','Display','none');
        X_optimal = fminunc(@(x) obj(x(1),x(2)),[1,1],options);
    case 2 % multiple local optimum
        options = optimoptions('fminunc','Display','none');
        X_optimal(1,:) = fminunc(@(x) obj(x(1),x(2)),[2.5,-1],options);
        X_optimal(2,:) = [-0.4495,2.2928];
        X_optimal(3,:) = [2.4239,1.9219];
    case 3
        options = optimoptions('fmincon','Display','iter');
        X_optimal = fmincon(@(x) obj(x(1),x(2)),[0,0],[],[],[],[],[],[],@(x)mycon(x,con),options);
end

% initialize
X_best = nan(max_iterations,2); F_best = nan(max_iterations,1); F_mean = nan(max_iterations,1);
x_best_particles = inf(np,n); f_best_particles = inf(np,1); f_best_overall = inf; x_best_overall = [];

% go through each generation
for k = 1:max_iterations

    % evaluate current swarm
    f = obj(x(:,1),x(:,2));

    % display and calculate things related to this iteration
    [X_best,F_best,F_mean] = disp_iteration(f,x,X_best,F_best,F_mean,k);

    % update best individual points
    for i = 1:np
        if f_best_particles(i) > f(i) % update if current point is better
            f_best_particles(i) = f(i);
            x_best_particles(i,:) = x(i,:);
        end
    end

    % update best swarm point
    [F_best_,I] = min(f_best_particles); % find best point in this swarm
    if f_best_overall > F_best_
        f_best_overall  = F_best_;
        x_best_overall = x_best_particles(I,:);
    end

    % compute each particle's step
    beta = beta_max*rand(np,1);
    gamma = gamma_max*rand(np,1);
    p = alpha*p + beta.*(x_best_particles - x) + gamma.*(x_best_overall - x);

    % limit step based on maximum step
    p = step_norm_limiter(p,p_norm_max);

    % update each particle's position
    x = x + p;

    % plot new swarm
    hp = plot_helper2(hp,parent_color,x,p,colors,X_optimal,k,con,lb,ub);
    drawnow % draw the plot
    pause(0.5) % pause for the animation

end

% plot connections between best generation points
plot(X_best(:,1),X_best(:,2),'c-','LineWidth',2)

% go through each iteration a plot
for k = 1:max_iterations

    % plot best points in each generation
    plot(X_best(k,1),X_best(k,2),'.','markersize',30,'color',colors(k,:))

    % text for best points in each generation
    text(X_best(k,1),X_best(k,2),string(k),'HorizontalAlignment','center',...
        'FontSize',8)

end

% plot optimal point
plot(X_optimal(:,1),X_optimal(:,2),'m.','markersize',30);

% plot convergence behavior
hf = figure; hf.Color = 'w'; hold on
ha = gca; ha.FontSize = 18; ha.LineWidth = 1;
plot(0:max_iterations-1,F_best,'b.-','LineWidth',2);
plot(0:max_iterations-1,F_mean,'r.-','LineWidth',2);
legend('Best','Mean')
xlabel('Iteration Number'); ylabel('$f$'); % axis labels
figure(1); % bring the first figure to the front

%--------------------------------------------------------------------------
% function to make it easier to display things in the command window
function disp_helper(name,number,n)

% default value of the number of digits
if isempty(n)
    n = 5;
end

% form string
str = strcat(string(name)," = ",mat2str(round(number,n)));

% display string
disp(str)

end

%--------------------------------------------------------------------------
% display and calculate things related to this iteration (generation)
function [X_best,F_best,F_mean] = disp_iteration(f,x,X_best,F_best,F_mean,k)

% get best point
[f_best,I_best] = min(f);
x_best = x(I_best,:);
f_mean = mean(f);
X_best(k,:) = x_best;
F_best(k) = f_best;
F_mean(k) = f_mean;

% display stuff
disp_helper("*iteration",k,[])
disp_helper("mean(f)",f_mean,[])
disp_helper("min(f)",f_best,[])

end

%--------------------------------------------------------------------------
% create initial plot
function [hp,parent_color,colors] = plot_helper1(obj,pen,con,lb,ub,x,p,max_iterations)

% create grid of evaluation points
N = 200;
x1 = linspace(lb(1),ub(1),N);
x2 = linspace(lb(2),ub(2),N);
[X1,X2] = meshgrid(x1,x2);

% evaluate the function
F_ = obj(X1,X2);

% remove penalty
if ~isempty(pen)
    F_ = F_ - pen(X1,X2);
end

% create figure
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
hf = figure; hf.Color = 'w'; hold on
contourf(X1,X2,F_,30);
xlim([lb(1) ub(1)]); ylim([lb(2) ub(2)]); % axis limits
xlabel('$x_1$'); ylabel('$x_2$'); % axis labels
ha = gca; ha.FontSize = 18; ha.LineWidth = 1;

% plot initial population
hp(1) = plot(x(:,1),x(:,2),'g.','markersize',20);
hp(2) = quiver(x(:,1),x(:,2),p(:,1),p(:,2),'off','g');

% plot colors
parent_color = [0.5 0.5 0.5];
colors = autumn(max_iterations);

end

%--------------------------------------------------------------------------
% plot new population
function hp = plot_helper2(hp,parent_color,x,p,colors,X_optimal,k,con,lb,ub)

% change old population points to gray and smaller
hp(end-1).Color = parent_color;
hp(end-1).MarkerSize = 10;
hp(end).Color = parent_color;
hp(end).MarkerSize = 10;
hp(end).LineStyle = 'none';

% plot current population
hp(end+1) = plot(x(:,1),x(:,2),'.','markersize',20,'color',colors(k,:));

% plot current steps
hp(end+1) = quiver(x(:,1),x(:,2),p(:,1),p(:,2),'off','color',colors(k,:));

% plot constraint line
if ~isempty(con)
    y1 = fsolve(@(y) con(lb(1),y),lb(2),optimoptions('fsolve','Display','none'));
    y2 = fsolve(@(y) con(ub(1),y),lb(2),optimoptions('fsolve','Display','none'));
    plot([lb(1) ub(1)],[y1,y2],'g-','linewidth',2)
end

% plot optimal point
plot(X_optimal(:,1),X_optimal(:,2),'m.','markersize',30);

% customize
ha = gca;
ha.Layer = 'top'; % place the axes on top of the data
box on

end

%--------------------------------------------------------------------------
% constraint function for fmincon
function [c,ceq] = mycon(x,con)

c = [];
ceq = con(x(1),x(2));

end

%--------------------------------------------------------------------------
% potentially limit norm of the step
function p = step_norm_limiter(p,p_norm_max)

% compute norm of each step
p_norm = sqrt(sum(p.^2,2));

% get indices that are larger than the maximum
I = p_norm > p_norm_max;

% normalize
p_limited = p_norm_max*p./p_norm;

% modify only the large steps
p(I,:) = p_limited(I,:);

end