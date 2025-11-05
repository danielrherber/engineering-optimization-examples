% ex_active_set.m
% illustration of an active-set method to solve a 2D quadratic program with
% inequality constraints
% [reference] Section 15.4 in LNO
% [course] Session 11 - Constrained Optimization (3)
close all; clear; clc

% problem data from Example 15.7 in LNO
% f(x) = 1/2*x'*Q*x + c'*x + d and Ax >= b
Q = [1 0; 0 2];
c = [-3; -4];
% c = [-2; -1]; % <- try this too
% c = [-2; 0.5]; % <- try this too
d = 17/2;
A = [2 -1; -1 -1; 0 1];
b = [0; -4; 0];

% initial feasible point
x = [0; 0];

% create the plot
plot_helper(Q,c,d,x,0,[],'init')

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

    % update based on current working set
    Ab = A(Iw,:); % constraint matrix
    Zb = null(Ab,'r'); % null space of the constraint matrix
    Abr = pinv(Ab); % right inverse of the constraint matrix

    % compute Lagrange multipliers
    g = Q*x + c; % objective gradient
    Lb = Abr'*g; % use right inverse
    E = Ab'*Lb - g; % error in FONC

    % display stuff
    disp_helper("---Iteration",iter,[])
    disp_helper("x",x,[])
    disp_helper("Iw",Iw,[])
    disp_helper("Zb",Zb,[])
    disp_helper("Lb",Lb,[])
    disp_helper("Lb Error",E,[])

    % plot iteration
    plot_helper([],[],[],x,iter,E,'iter')

end

%--------------------------------------------------------------------------
% function to create the initial plot
function plot_helper(Q,c,d,x,iter,E,tasknum)

% colors and other parameters
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;
LineWidth = 1.5;
MarkerSize = 24;
FontSize = 18;
plotOpts = {'LineWidth',LineWidth,'MarkerSize',MarkerSize};

switch tasknum
    %----------------------------------------------------------------------
    case 'init' % initial point

    % initialize figure
    hf = figure; hf.Color = 'w'; hold on

    % labels
    xlabel('$x_1$','Interpreter','latex');
    ylabel('$x_2$','Interpreter','latex');

    % axis properties
    ha = gca; ha.XColor = 'k'; ha.YColor = 'k'; ha.Color = 'none';
    ha.LineWidth = 1; ha.FontSize = FontSize;

    % plot objective function contours
    N = 1e3;
    x1 = linspace(-0.5,4.5,N);
    x2 = linspace(-0.5,3,N);
    xlim([min(x1),max(x1)]); ylim([min(x2),max(x2)]); % axis limits
    [X1,X2] = meshgrid(x1,x2);
    F_ = nan(size(X1));
    for k = 1:numel(X1)
        x_ = [X1(k);X2(k)];
        F_(k) = 1/2*x_'*Q*x_ + c'*x_ + d;
    end
    contour(X1,X2,F_,50)

    % hard-coded constraint boundary (does NOT depend on A and b)
    hp = patch('XData',[4/3 0 4],'YData',[8/3 0 0]);
    hp.EdgeColor = nicegreen;
    hp.LineWidth = LineWidth;
    hp.FaceAlpha = 0.15;

    % plot initial feasible point
    plot(x(1),x(2),'.',plotOpts{:},'Color',niceblue)
    text(x(1)-0.2,x(2),string(iter),'FontSize',FontSize,'Color','k')

    %----------------------------------------------------------------------
    case 'iter' % plot current iteration

    % plot current point
    plot(x(1),x(2),'.',plotOpts{:},'Color',nicered)
    text(x(1)-0.2,x(2),string(iter),'FontSize',FontSize,'Color','k')

    % check if inconsistent multipliers
    if norm(E) > 1e-8
        text(x(1)+0.15,x(2),'inconsistent multipliers',...
            'FontSize',FontSize,'Color','k')
    end

end

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