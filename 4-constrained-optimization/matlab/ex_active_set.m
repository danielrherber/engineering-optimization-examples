% ex_active_set.m
% illustration of an active-set method to solve a 2d quadratic program with
% inequality constraints
% [reference] Section 15.4 in LNO
% [course] Session 11 - Constrained Optimization (3)
close all; clear; clc

% problem data from Example 15.7 in LNO
% f(x) = 1/2*x'*Q*x + c'*x + d and Ax >= b
Q = [1 0; 0 2];
c = [-3; -4];
d = 17/2;
A = [2 -1; -1 -1; 0 1];
b = [0; -4; 0];

% create the plot
plot_helper(Q,c,d,A,b)

% initial feasible point
x = [0; 0];
plot(x(1),x(2),'.b','markersize',24)
text(x(1)-0.2,x(2),string(0),'color','w','FontSize',14)

% tolerances
ConstraintTolerance = 1e-14;
OptimalityTolerance = 1e-14;
MaxIterations = 100;

% determine initial working set
Iw = find(abs(A*x - b) <= ConstraintTolerance);

% update working set
Ab = A(Iw,:); % constraint matrix
Zb = null(Ab,'r'); % null space of the constraint matrix
Abr = pinv(Ab); % right inverse of the constraint matrix

% number of constraints
m = size(A,1);

% initial optimality flag
IsOptimal = false;

% go through each iteration
for iter = 1:MaxIterations

    %--- The Optimality Test
    % compute gradient
    g = Q*x + c;

    % check reduced gradient
    while norm(Zb'*g) <= OptimalityTolerance

        % compute (approximate) Lagrange multipliers
        Lb = Abr'*g;

        % find all negative multipliers
        In = Lb < 0;

        % stop if all multipliers are nonnegative
        if ~any(In)
            disp("Optimal!")
            IsOptimal = true;
            break
        end

        % select most negative multiplier
        [~,Ir] = min(Lb(In));

        % remove the constraint from the working set
        Iw(Ir) = [];

        % update working set
        Ab = A(Iw,:); % constraint matrix
        Zb = null(Ab,'r'); % null space of the constraint matrix
        Abr = pinv(Ab); % right inverse of the constraint matrix

    end

% check if optimal
if IsOptimal
    break
end

%--- The Search Direction
% compute reduced Newton search direction
p = -Zb*((Zb'*Q*Zb)\(Zb'*g));

%--- The Step
% initialize maximum step as the newton step
alpha = 1;

% go through each constraint
for k = 1:m

    % only if the constraint is not in the working set
    if ~any(k == Iw)

        % extract current constraint row
        a = A(k,:);

        % check if this is a descent direction
        if a*p < 0

            % compute maximum step length with ratio test
            alpha_ = -(a*x - b(k))/(a*p);

            % update maximum step length
            alpha = min(alpha_,alpha);

        end

    end

end

%--- The Update
% take step
x = x + alpha*p;

% update working set (this add all approximately active constraints)
Iw = find(abs(A*x - b) <= ConstraintTolerance);

% update working set
Ab = A(Iw,:); % constraint matrix
Zb = null(Ab,'r'); % null space of the constraint matrix
Abr = pinv(Ab); % right inverse of the constraint matrix

% display stuff
disp_helper("---Iteration",iter,[])
disp_helper("x",x,[])
disp_helper("Iw",Iw,[])
disp_helper("Zb",Zb,[])
plot(x(1),x(2),'.r','markersize',24)
text(x(1)-0.2,x(2),string(iter),'color','w','FontSize',14)

% compute Lagrange multipliers
g = Q*x + c; Lb = Abr'*g; E = Ab'*Lb - g;
disp_helper("Lb",Lb,[])
disp_helper("Lb Error",E,[])
if norm(E) > 1e-8
    text(x(1)+0.15,x(2),'inconsistent multipliers','color','w','FontSize',14)
end

end

%--------------------------------------------------------------------------
% function to create the initial plot
function plot_helper(Q,c,d,A,b)

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
hf = figure; hf.Color = 'w'; hold on
N = 1e3;
x1 = linspace(-0.5,4.5,N);
x2 = linspace(-0.5,3,N);
xlim([min(x1),max(x1)]); ylim([min(x2),max(x2)]); % axis limits
[X1,X2] = meshgrid(x1,x2);
F_ = nan(size(X1));
for k = 1:numel(X1)
    x = [X1(k);X2(k)];
    F_(k) = 1/2*x'*Q*x + c'*x + d;
end
contourf(X1,X2,F_,50)
xlabel('$x_1$'); ylabel('$x_2$'); % axis labels
ha = gca;
ha.LineWidth = 1; ha.FontSize = 18;

% hard-coded constraint boundary
hp = patch('XData',[4/3 0 4],'YData',[8/3 0 0]);
hp.EdgeColor = 'g';
hp.LineWidth = 1.5;
hp.FaceAlpha = 0.15;

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