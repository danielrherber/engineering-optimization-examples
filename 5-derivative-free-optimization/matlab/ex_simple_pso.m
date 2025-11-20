% ex_simple_pso.m
% illustration of particle swarm optimization (PSO) on 2D problems
% [reference] Section 7.7 in EDO
% [course] Session 13 - Derivative-Free Optimization (2)
close all; clear; clc

example = 2; % see cases below

switch example
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

% create initial plot
[hp,parent_color,colors] = plot_helper1(obj,pen,lb,ub,x,p,max_iterations,X_optimal);

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
    hp = plot_helper2(hp,parent_color,x,p,colors,k,con,lb,ub,obj,pen);
    drawnow % draw the plot
    pause(0.5) % pause for the animation

end

% plot final results
plot_helper3(X_best,obj,pen,max_iterations,X_optimal,F_best,F_mean)

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
function [hp,parent_color,colors] = plot_helper1(obj,pen,lb,ub,x,p,max_iterations,X_optimal)

% colors and other parameters
nicegreen = [109, 195, 80]/255;
FontSize = 18;
LineWidth = 1.5;
MarkerSize = 24;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

% create grid of evaluation points
N = 200;
x1 = linspace(lb(1),ub(1),N);
x2 = linspace(lb(2),ub(2),N);
[X1,X2] = meshgrid(x1,x2);

% evaluate the function
F_contour = obj(X1,X2);
F = obj(x(:,1),x(:,2));
F_optimal = obj(X_optimal(:,1),X_optimal(:,2));

% remove penalty
if ~isempty(pen)
    F_contour = F_contour - pen(X1,X2);
    F = F - pen(x(:,1),x(:,2));
    F_optimal = F_optimal - pen(X_optimal(:,1),X_optimal(:,2));
end

% initialize figure
hf = figure; hf.Color = 'w'; hold on

% labels
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');

% axis properties
ha = gca; ha.XColor = 'k'; ha.YColor = 'k'; ha.ZColor = 'k';
ha.Color = 'none'; ha.LineWidth = 1; ha.FontSize = FontSize;

% plot objective function contours
contour3(X1,X2,F_contour,60,'LineWidth',LineWidth);

% axis limits
xlim([lb(1) ub(1)]); ylim([lb(2) ub(2)]); zlim([min(F_optimal) max(F_contour(:))]);

% plot initial population
hp(1) = plot3(x(:,1),x(:,2),F,'.',plotOpts{:},'Color',nicegreen);
hp(2) = quiver3(x(:,1),x(:,2),F,p(:,1),p(:,2),zeros(size(F)),'off','Color',nicegreen);

% plot colors
parent_color = [0.5 0.5 0.5];
colors = autumn(max_iterations);

end

%--------------------------------------------------------------------------
% plot new population
function hp = plot_helper2(hp,parent_color,x,p,colors,k,con,lb,ub,obj,pen)

% colors and other parameters
nicegreen = [109, 195, 80]/255;
LineWidth = 2*1.5;
MarkerSize = 24;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

% change old population points to gray and smaller
hp(end-1).Color = parent_color;
hp(end-1).MarkerSize = 10;
hp(end).Color = parent_color;
hp(end).MarkerSize = 10;
hp(end).LineStyle = 'none';

% current population function values
F = obj(x(:,1),x(:,2));

% remove penalty
if ~isempty(pen)
    F = F - pen(x(:,1),x(:,2));
end

% plot current population
hp(end+1) = plot3(x(:,1),x(:,2),F,...
    '.','markersize',MarkerSize,'color',colors(k,:));

% plot current steps
hp(end+1) = quiver3(x(:,1),x(:,2),F,p(:,1),p(:,2),zeros(size(F)),'off',...
    'Color',colors(k,:),'LineWidth',LineWidth);

% plot constraint line
if ~isempty(con)
    y1 = fsolve(@(y) con(lb(1),y),lb(2),optimoptions('fsolve','Display','none'));
    y2 = fsolve(@(y) con(ub(1),y),lb(2),optimoptions('fsolve','Display','none'));
    x1 = linspace(lb(1),ub(1),1e3)';
    x2 = linspace(y1,y2,1e3)';
    plot3(x1,x2,obj(x1,x2)-pen(x1,x2),'-',plotOpts{:},'Color',nicegreen)
    % plot([lb(1) ub(1)],[y1,y2],'g-','linewidth',2)
end

% customize
ha = gca;
ha.Layer = 'top'; % place the axes on top of the data
box on

end

%--------------------------------------------------------------------------
% plot best points and new plot of convergence behavior
function plot_helper3(X_best,obj,pen,max_iterations,X_optimal,F_best,F_mean)

% colors and other parameters
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicepink = [239, 142, 219]/255;
nicepurple = [150, 103, 126]/255;
FontSize = 18;
LineWidth = 1.5;
MarkerSize = 24;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

% current population function values
F = obj(X_best(:,1),X_best(:,2));
F_optimal = obj(X_optimal(:,1),X_optimal(:,2));

% remove penalty
if ~isempty(pen)
    F = F - pen(X_best(:,1),X_best(:,2));
    F_optimal = F_optimal - pen(X_optimal(:,1),X_optimal(:,2));
end

% plot connections between best generation points
plot3(X_best(:,1),X_best(:,2),F,...
    '.-','LineWidth',LineWidth,'MarkerSize',MarkerSize*2,'Color',nicepink)

% text for best points in each generation
for k = 1:max_iterations

    text(X_best(k,1),X_best(k,2),F(k),string(k),...
        'HorizontalAlignment','center','FontSize',FontSize/1.5)

end

% plot optimal point
plot3(X_optimal(:,1),X_optimal(:,2),F_optimal,'.',plotOpts{:},'Color',nicepurple);

% initialize figure
hf = figure; hf.Color = 'w'; hold on

% labels
xlabel('Iteration Number');
ylabel('$f$','Interpreter','latex');

% axis properties
ha = gca; ha.XColor = 'k'; ha.YColor = 'k';
ha.Color = 'none'; ha.LineWidth = 1; ha.FontSize = FontSize;

% plot convergence behavior
plot(0:max_iterations-1,F_best,'.-',plotOpts{:},'Color',niceblue);
plot(0:max_iterations-1,F_mean,'.-',plotOpts{:},'Color',nicered);

% legend
legend('Best','Mean')

% bring the first figure to the front
figure(1);

end