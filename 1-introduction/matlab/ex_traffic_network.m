% ex_traffic_network.m
% example of a nonlinear program (NLP) seeking to minimize the total travel
% time through the network for a volume of X cars per hour
% [reference] Section 1.6 in LNO
% [course] Session 1 - Introduction to Nonlinear Optimization
close all; clear; clc

% problem data
p.t = [0.1 0.14 0.13 0.25 0.12]; % constant travel time [hour]
p.a = [0.0001 0.0001 0.0002 0.0001 0.0001]; % travel rate density constant [hour^2/car]
p.c = [2200 2200 1500 1000 2000]; % road capacities [car/hour]
p.X = 2000; % volume of incoming/outgoing cars [car/hour]
p.lower = [0 0 0 0 0]; % no negative rates [car/hour]

% linear traffic network constraints
% x12 = x(1); x13 = x(2); x23 = x(3); x24 = x(4); x34 = x(5);
Aeq(1,:) = [1 1 0 0 0];
Aeq(2,:) = [-1 0 1 1 0];
Aeq(3,:) = [0 -1 -1 0 1];
Aeq(4,:) = [0 0 0 1 1];
Beq = [p.X;0;0;p.X];

% initial guess (happens to be feasible)
x0 = [p.X 0 0 p.X 0];

% optimization options
OPTIONS = optimoptions('fmincon');
OPTIONS.Display = "Iter";
OPTIONS.Algorithm = "sqp";
OPTIONS.SpecifyObjectiveGradient = true;

% small capacity offset to avoid inf values near maximum capacities
ep = 0.1;

% solve the optimization problem
[x,f] = fmincon(@(x) travel_time(x,p),x0,[],[],Aeq,Beq,p.lower,p.c-ep,[],OPTIONS);

% visualize the solution
plot_solution(x,f,p)

%--------------------------------------------------------------------------
% total travel time (objective) function
function [f,g] = travel_time(x,p)

% extract problem data
t = p.t; a = p.a; c = p.c;

% calculate total travel time
f = sum(x.*(t + a.*(x./(1-x./c))));

% compute gradient if needed
if nargout  > 1
    t0 = x.*a;
    t1 = ones(size(x)) - x./c;
    g = t + 2*t0./t1 + t0.*x./(t1.*t1)./c;
end

end

%--------------------------------------------------------------------------
% plot the optimal solution
function plot_solution(x,f,p)

% display optimal solution
disp(x)
disp(f)

% visualize the solution
s = []; t = []; w = [];
s(end+1) = 1; t(end+1) = 2; w(end+1) = p.X;
s(end+1) = 2; t(end+1) = 3; w(end+1) = x(1);
s(end+1) = 2; t(end+1) = 4; w(end+1) = x(2);
s(end+1) = 3; t(end+1) = 4; w(end+1) = x(3);
s(end+1) = 3; t(end+1) = 5; w(end+1) = x(4);
s(end+1) = 4; t(end+1) = 5; w(end+1) = x(5);
s(end+1) = 5; t(end+1) = 6; w(end+1) = p.X;

% create graph/network
G = digraph(s,t,w);

% edge colors
EdgeColor = repmat([0, 0, 0]/255,length(p.t)+2,1);
EdgeColor([false logical(x) false],:) = repmat([75, 120, 166]/255,nnz(x),1);

% edge labels
EdgeLabels = strcat(string(round(G.Edges.Weight,0)));
EdgeLabels(2:end-1) = strcat(string(round(G.Edges.Weight(2:end-1),0)),"/",string(p.c'));

% edge widths
EdgeWidths = (5*[p.X x p.X]/max([x p.X]));

% node colors
NodeColors = repmat([225, 87, 87]/255,6,1);
NodeColors(1,:) = [255 255 255]/255;
NodeColors(end,:) = [255 255 255]/255;

% node labels
NodeLabels = ["","1","2","3","4",""];

% initialize figure
hf = figure; hf.Color = 'w';

% create network
hp = plot(G,'EdgeLabel',EdgeLabels,'LineWidth',EdgeWidths+eps,...
    'EdgeColor',EdgeColor,'layout','force',...
    'NodeColor',NodeColors,'NodeLabel',NodeLabels);

% customize plot
axis off
axis equal
hp.XData = [13.918-20, 13.918, 41.573, 63.697, 100.564, 100.564+20];
hp.YData = [22.826, 22.826, 33.886, 4.760, 22.826, 22.826];

end