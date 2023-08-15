close all; clear; clc

% node names
NodeInfo = {'Seattle','Denver','Albuquerque','Chicago','Nashville',...
    'Columbus','Charlotte','Boston'};

% number of nodes
N = numel(NodeInfo);

% initialize adjacency matrix for the network
u = zeros(N);

% add arcs with capacities
u(1,2) = 7;
u(1,3) = 5;
u(1,4) = 12;
u(1,8) = 15;
u(2,4) = 11;
u(3,4) = 6;
u(3,5) = 4;
u(3,7) = 5;
u(4,2) = 10;
u(4,5) = 7;
u(4,6) = 21;
u(5,4) = 8;
u(5,7) = 6;
u(6,4) = 21;
u(6,7) = 15;
u(6,8) = 25;
u(7,5) = 7;
u(7,6) = 15;
u(7,8) = 18;
u(8,6) = 25;
u(8,7) = 15;

% find vertex pairs that define all edges in the network
[I,J,U] = find(u);

% number of arcs
n = numel(U);

% create primary coefficient matrix
Aeq_ = zeros(N,n);
Aeq_(sub2ind([N n],I',1:n)) = 1;
Aeq_(sub2ind([N n],J',1:n)) = -1;

% add additional variable f to the coefficient matrix
Aeq = [Aeq_,zeros(N,1)];
Aeq(1,end) = -1; % source node variable -f
Aeq(8,end) = 1; % sink node variable f

% other matrices
c = [zeros(n,1);1]; % max f
A = []; % not present
b = []; % not present
beq = zeros(N,1);
LB = zeros(length(c),1); % zero lower bounds
UB = [U;inf]; % nonzero upper bounds

% options
OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';

% solve using linprog
[X,FVAL,EXITFLAG,OUTPUT] = linprog(-c,A,b,Aeq,beq,LB,UB,OPTIONS)

%% visualize the solution
hf = figure; hf.Color = 'w'; hold on

% colors
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;

% create the graph
EdgeInfo = table(U,(1:n)','VariableNames',["Weight","Index"]);
G = digraph(I,J,EdgeInfo,NodeInfo);

% determine line widths
LWidths = 3*G.Edges.Weight/max(G.Edges.Weight);

% get original edge indices (because matlab sorts the edges)
Ir = G.Edges.Index;

% determine if edges/variables were nonzero (using the original ordering)
Xorig = logical(X(Ir));

% assign edge colors
EdgeColor = repmat([0 0 0],length(c)-1,1);
EdgeColor(Xorig(1:end-1),:) = repmat(nicegreen,sum(Xorig(1:end-1)),1);

% create edge labels in the format x(i)/max value of x(i)
EdgeLabel = strcat(string(X(Ir)),"/",string(G.Edges.Weight));

% create the figure
hp = plot(G,'layout','circle','LineWidth',LWidths,'EdgeColor',EdgeColor,...
    'EdgeLabel',EdgeLabel);

% set node color
hp.NodeColor = nicered;

% manually position nodes consistent with book picture
XData = hp.XData;
XData(1) = 14.924;
XData(2) = 37.337;
XData(3) = 34.570;
XData(4) = 53.662;
XData(5) = 64.869;
XData(6) = 76.076;
XData(7) = 84.377;
XData(8) = 95.445;
hp.XData = XData;

YData = hp.YData;
YData(1) = 27.978;
YData(2) = 14.281;
YData(3) = 3.213;
YData(4) = 28.116;
YData(5) = 11.376;
YData(6) = 25.210;
YData(7) = 8.746;
YData(8) = 28.116+5;
hp.YData = YData;

axis off
axis equal