% ex_simple_ga.m
% illustration of a genetic algorithm (GA) on 2D problems using roulette
% wheel selection, linear crossover, and uniform random mutation
% [reference] Section 7.6 in EDO
% [course] Session 13 - Derivative-Free Optimization (2)
close all; clear; clc

example = 2; % see cases below

switch example
    case 1
        fun = 1;
        max_iterations = 25; % maximum number of iterations
        np = 30; % population size (even number for simplicity)
        rng(876044) % set random seed for repeatable results
    case 2
        fun = 2;
        max_iterations = 30; % maximum number of iterations
        np = 30; % population size (even number for simplicity)
        rng(876044) % set random seed for repeatable results
    case 3
        fun = 3;
        max_iterations = 30; % maximum number of iterations
        np = 30; % population size (even number for simplicity)
        rng(876044) % set random seed for repeatable results
end

% algorithm parameters
p = 0.01; % mutation probability
Delta = 0.5; % maximum mutation amount
alpha = 0.5; % linear crossover parameter

% symbolically create objective function and derivatives
syms x y real
con = []; pen = [];
switch fun
    case 1 % spring example
        l1 = 12; l2 = 8; k1 = 1; k2 = 10; mg = 7; % problem data
        fs =  k1/2*(sqrt((l1 + x).^2 + y.^2) - l1).^2 + ...
            k2/2*(sqrt((l2 - x).^2 + y.^2) - l2).^2 + ...
            -mg*y;
        ub = [16 12];
        lb = [-6 -9];
    case 2 % Jones function
        fs = x^4 + y^4 - 4*x^3 - 3*y^3 + 2*x^2 + 2*x*y;
        ub = [4 3];
        lb = [-2 -2];
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
        ub = [16 12];
        lb = [-6 -9];
end

% create matlab functions
obj = matlabFunction(fs);

%--- random initial population
x1_= lb(1) + rand(np,1)*(ub(1)-lb(1));
x2_= lb(2) + rand(np,1)*(ub(2)-lb(2));
x = [x1_,x2_];

% create initial plot
[hp,parent_color,colors] = plot_helper1(obj,pen,lb,ub,x,max_iterations);

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

% go through each generation
for k = 1:max_iterations

    % evaluate population
    f = obj(x(:,1),x(:,2));

    % display and calculate things related to this iteration (generation)
    [X_best,F_best,F_mean] = disp_iteration(f,x,X_best,F_best,F_mean,k);

    % selection step
    parent_indices = roulette_wheel_selection(f);

    % crossover step
    children = linear_crossover(f,x,parent_indices,alpha);

    % mutation step
    children = uniform_random_mutation(children,p,Delta);

    % fix points outside the bounds
    children(children(:,1)<lb(1),1) = lb(1);
    children(children(:,2)<lb(2),2) = lb(2);
    children(children(:,1)>ub(1),1) = ub(1);
    children(children(:,2)>ub(2),2) = ub(2);

    % plot new population
    hp = plot_helper2(hp,parent_color,children,colors,k,con,lb,ub,obj,pen);
    drawnow % draw the plot
    pause(0.5) % pause for the animation

    % set current children are the next population
    x = children;

end

% plot final results
plot_helper3(X_best,obj,pen,max_iterations,X_optimal,F_best,F_mean)

%--------------------------------------------------------------------------
% implementation of roulette wheel selection
function parent_indices = roulette_wheel_selection(f)

% determine number of points in a generation
np = length(f);

% minimum, maximum, and spread of the function values
fmin = min(f);
fmax = max(f);
Fdelta = 1.1*fmax - 0.1*fmin;

% scaled fitness values
F = (-f + Fdelta)/(max(1,Fdelta-fmin));

% normalized cumulative sum of the scaled fitness values
S = cumsum(F)/sum(F);

% random [0,1] matrix for creating the pairs
R = rand(np/2,2);

% determine the index of the parent
parent_indices = nan(np/2,2);
for i = 1:np/2

    % determine first parent
    parent_indices(i,1) = find(S>=R(i,1),1);

    % determine second parent
    parent_indices(i,2) = find(S>=R(i,2),1);

end

end

%--------------------------------------------------------------------------
% implementation of linear crossover on 2 variables
function children = linear_crossover(f,x,parent_indices,alpha)

% determine number of points in a generation
np = length(f);

% determine children
children = nan(np,2);

% go through each pair of parents
for i = 1:np/2

    % extract parent indices
    parent1_index = parent_indices(i,1);
    parent2_index = parent_indices(i,2);

    % flip points is ordering is wrong
    if f(parent2_index) > f(parent1_index)
        temp = parent1_index;
        parent1_index = parent2_index;
        parent2_index = temp;
    end

    % extract parents
    xp1 = x(parent1_index,:);
    xp2 = x(parent2_index,:);

    % child 1
    children(2*i-1,:) = 0.5*xp1 + 0.5*xp2;

    % child 2
    children(2*i,:) = (1+alpha)*xp2 - alpha*xp1;

end

end

%--------------------------------------------------------------------------
% implementation of uniform random mutation on 2 variables
function children = uniform_random_mutation(children,p,Delta)

% determine number of points in a generation
np = size(children,1);

% random [0,1] matrix for determining if mutation should occur
pr = rand(np,2);

% random [0,1] matrix for determining how much mutation should occur
rr = rand(np,2);

% mutate some children
children = children + (pr <= p).*(rr-0.5)*Delta;

end

%--------------------------------------------------------------------------
% constraint function for fmincon
function [c,ceq] = mycon(x,con)

c = [];
ceq = con(x(1),x(2));

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
function [hp,parent_color,colors] = plot_helper1(obj,pen,lb,ub,x,max_iterations)

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
F_ = obj(X1,X2);
F = obj(x(:,1),x(:,2));

% remove penalty
if ~isempty(pen)
    F_ = F_ - pen(X1,X2);
    F = F - pen(x(:,1),x(:,2));
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
contour3(X1,X2,F_,60,'LineWidth',LineWidth);

% axis limits
xlim([lb(1) ub(1)]); ylim([lb(2) ub(2)]); zlim([min(F_(:)) max(F_(:))]);

% plot initial population
hp(1) = plot3(x(:,1),x(:,2),F,'.',plotOpts{:},'Color',nicegreen);

% plot colors
parent_color = [0.5 0.5 0.5];
colors = autumn(max_iterations);

end

%--------------------------------------------------------------------------
% plot new population
function hp = plot_helper2(hp,parent_color,children,colors,k,con,lb,ub,obj,pen)

% colors and other parameters
nicegreen = [109, 195, 80]/255;
LineWidth = 2*1.5;
MarkerSize = 24;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

% change old population points to gray and smaller
hp(end).Color = parent_color;
hp(end).MarkerSize = MarkerSize/2;

% current population function values
F = obj(children(:,1),children(:,2));

% remove penalty
if ~isempty(pen)
    F = F - pen(children(:,1),children(:,2));
end

% plot current population
hp(end+1) = plot3(children(:,1),children(:,2),F,...
    '.','markersize',MarkerSize,'color',colors(k,:));

% plot constraint line
if ~isempty(con)
    y1 = fsolve(@(y) con(lb(1),y),lb(2),optimoptions('fsolve','Display','none'));
    y2 = fsolve(@(y) con(ub(1),y),lb(2),optimoptions('fsolve','Display','none'));
    x1 = linspace(lb(1),ub(1),1e3)';
    x2 = linspace(y1,y2,1e3)';
    plot3(x1,x2,obj(x1,x2)-pen(x1,x2),'-',plotOpts{:},'Color',nicegreen)
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