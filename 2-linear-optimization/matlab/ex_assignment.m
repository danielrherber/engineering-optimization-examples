% ex_assignment.m
% example of an assignment network problem where a company is planning to
% assign three people to three jobs
% [reference] Example 8.4 in LNO
% [course] Session 6 - Linear Optimization (3)
close all; clear; clc

% problem data
c = [11 5 2 15 12 8 3 1 10];
% c = [11 5 2 20 12 8 3 1 10]; % alternative c vector
A = []; % not present
b = []; % not present
Aeq = [1 1 1 0 0 0 0 0 0;
    0 0 0 1 1 1 0 0 0;
    0 0 0 0 0 0 1 1 1;
    -1 0 0 -1 0 0 -1 0 0;
    0 -1  0 0 -1 0 0 -1 0;
    0 0 -1 0 0 -1 0 0 -1];
beq = [1 1 1 -1 -1 -1];
LB = zeros(length(c),1);
UB = ones(length(c),1);

% options
OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';

% solve using linprog
[X,FVAL,EXITFLAG,OUTPUT] = linprog(-c,A,b,Aeq,beq,LB,UB,OPTIONS)

%--------------------------------------------------------------------------
% visualize the solution
hf = figure; hf.Color = 'w'; hold on

% colors
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;

% create the graph
s = [1 1 1 2 2 2 3 3 3];
t = [4 5 6 4 5 6 4 5 6];
weights = c;
names = {'Person 1','Person 2','Person 3',...
    'Accountant','Budget Director','Personnel Manager'};
G = digraph(s,t,weights,names);

% determine line widths
LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);

% assign edge colors
EdgeColor = repmat([0 0 0],length(c),1);
EdgeColor(logical(X),:) = repmat(nicegreen,sum(X),1);

% create the figure
hp = plot(G,'LineWidth',LWidths,'Layout','layered',...
    'Direction','right','EdgeColor',EdgeColor);

% set node color
hp.NodeColor = nicered;
hp.NodeFontSize = 12;

% manually position nodes consistent with book picture
% XData = hp.XData;
% XData(1) = 18.709;
% XData(2) = 18.709;
% XData(3) = 18.709;
% XData(4) = 87.734;
% XData(5) = 87.734;
% XData(6) = 87.734;
% hp.XData = XData;
% 
% YData = hp.YData;
% YData(1) = 46.260;
% YData(2) = 25.103;
% YData(3) = 3.948;
% YData(4) = 46.260;
% YData(5) = 25.103;
% YData(6) = 3.948;
% hp.YData = YData;

axis off
axis equal