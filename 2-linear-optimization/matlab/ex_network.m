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

% initial point
X0 = []; % not needed

% options
OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';

% solve using linprog
[X,FVAL,EXITFLAG,OUTPUT] = linprog(c,A,b,Aeq,beq,LB,UB,X0,OPTIONS)

%% visualize the solution
hf = figure; hf.Color = 'w'; hold on

s = [1 1 2 2 3 3 4 4 5 5 6];
t = [2 3 4 6 4 5 5 6 6 7 7];
weights = c;
G = digraph(s,t,weights);
LWidths = (5*X/max(X))+eps;

EdgeColor = repmat([0.7 0.7 0.7],length(c),1);
EdgeColor(logical(X),:) = repmat([1 0 0],nnz(X),1);

hp = plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,...
    'EdgeColor',EdgeColor);
text(hp.XData(1)-0.1,hp.YData(1)+0.3,'50')
text(hp.XData(6)-0.1,hp.YData(6)-0.3,'-20')
text(hp.XData(7)-0.1,hp.YData(7)-0.3,'-30')

axis off