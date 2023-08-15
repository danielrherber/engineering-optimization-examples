close all; clear; clc

% problem data
c = [12 15 9 8 7 6 8 4 5 3 11];
% c = [12 15 1 8 7 6 8 4 5 3 1]; % alternative c vector
A = []; % not present
b = []; % not present
Aeq = zeros(7,11);
Aeq(1,1:2) = [1 1];
Aeq(2,1:4) = [-1 0 1 1];
Aeq(3,2:6) = [-1 0 0 1 1];
Aeq(4,3:8) = [-1 0 -1 0 1 1];
Aeq(5,6:10) = [-1 -1 0 1 1];
Aeq(6,4:11) = [-1 0 0 0 -1 -1 0 1];
Aeq(7,10:11) = [-1 -1];
beq = [50 0 0 0 0 -20 -30];
LB = zeros(length(c),1);
UB = 30*ones(length(c),1);

% options
OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';

% solve using linprog
[X,FVAL,EXITFLAG,OUTPUT] = linprog(c,A,b,Aeq,beq,LB,UB,OPTIONS)

%% visualize the solution
hf = figure; hf.Color = 'w'; hold on

% colors
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;

% create the graph
s = [1 1 2 2 3 3 4 4 5 5 6];
t = [2 3 4 6 4 5 5 6 6 7 7];
weights = c;
G = digraph(s,t,weights);
LWidths = (5*X/max(X))+eps;

% assign edge colors
EdgeColor = repmat([0.7 0.7 0.7],length(c),1);
EdgeColor(logical(X),:) = repmat(nicegreen,nnz(X),1);

% create edge labels in the format x(i)/max value of x(i)
EdgeLabel = strcat("x=",string(X)," @ c=",string(G.Edges.Weight));

% create the figure
hp = plot(G,'LineWidth',LWidths,'EdgeColor',EdgeColor,...
    'EdgeLabel',EdgeLabel);

% set node color
hp.NodeColor = nicered;

% manually position nodes consistent with book picture
XData = hp.XData;
XData(1) = 8.506;
XData(2) = 24.216;
XData(3) = 24.217;
XData(4) = 45.989;
XData(5) = 45.989;
XData(6) = 76.076;
XData(7) = 89.689;
hp.XData = XData;

YData = hp.YData;
YData(1) = 10.586;
YData(2) = 18.360;
YData(3) = 2.806;
YData(4) = 18.359;
YData(5) = 2.806;
YData(6) = 18.320+10;
YData(7) = 2.808;
hp.YData = YData;

% labels sources and sinks
text(hp.XData(1)-3,hp.YData(1)+4,'50')
text(hp.XData(6)+3,hp.YData(6)+4,'-20')
text(hp.XData(7)+3,hp.YData(7)+4,'-30')

axis off
axis equal