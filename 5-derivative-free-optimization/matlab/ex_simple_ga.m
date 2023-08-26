% ex_simple_ga.m
% illustration of a genetic algorithm (GA) on 2d problems using roulette
% wheel selection, linear crossover, and uniform random mutation
% [reference] Section 7.6 in EDO
% [course] Session 13 - Derivative-free Optimization (2)
close all; clear; clc

test = 2; % see cases below

switch test
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
[hp,parent_color,colors] = plot_helper1(obj,pen,con,lb,ub,x,max_iterations);

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
    hp = plot_helper2(hp,parent_color,children,colors,X_optimal,k,con,lb,ub);
    drawnow % draw the plot
    pause(0.5) % pause for the animation

    % set current children are the next population
    x = children;

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
function [hp,parent_color,colors] = plot_helper1(obj,pen,con,lb,ub,x,max_iterations)

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

% plot colors
parent_color = [0.5 0.5 0.5];
colors = autumn(max_iterations);

end

%--------------------------------------------------------------------------
% plot new population
function hp = plot_helper2(hp,parent_color,children,colors,X_optimal,k,con,lb,ub)

% change old population points to gray and smaller
hp(end).Color = parent_color;
hp(end).MarkerSize = 10;

% plot current population
hp(end+1) = plot(children(:,1),children(:,2),'.','markersize',20,'color',colors(k,:));

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