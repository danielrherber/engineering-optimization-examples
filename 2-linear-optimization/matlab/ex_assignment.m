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

% initial point
X0 = []; % not needed

% options
OPTIONS = optimoptions('linprog');
OPTIONS.Display = 'iter';

% solve using linprog
[X,FVAL,EXITFLAG,OUTPUT] = linprog(-c,A,b,Aeq,beq,LB,UB,X0,OPTIONS)

%% visualize the solution
hf = figure; hf.Color = 'w'; hold on

s = [1 1 1 2 2 2 3 3 3];
t = [4 5 6 4 5 6 4 5 6];
weights = c;
names = {'Person 1','Person 2','Person 3',...
    'Accountant','Budget Director','Personnel Manager'};
G = digraph(s,t,weights,names);
LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);

EdgeColor = repmat([0 0 0],length(c),1);
EdgeColor(logical(X),:) = repmat([1 0 0],sum(X),1);

hp = plot(G,'LineWidth',LWidths,'Layout','layered',...
    'Direction','right','EdgeColor',EdgeColor);

axis off