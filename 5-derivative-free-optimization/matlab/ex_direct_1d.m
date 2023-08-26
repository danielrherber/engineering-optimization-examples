% ex_direct_1d.m
% illustration of the dividing rectangles (DIRECT) algorithm in 1d
% [reference] Section 7.5 in EDO
% [course] Session 13 - Derivative-free Optimization (2)
close all; clear; clc

test = 1;
plotflag = true;

switch test
    %----------------------------------------------------------------------
    case 1
        fun = @(x) sin(x) + sin(10/3*x);
        a = -3; % lower limit
        b = 8; % upper limit
        YLimits = [-2.5 2]; % function value limits for plotting
    %----------------------------------------------------------------------
    case 2
        fun = @(x) round(-x.^2.*abs(1+sin(x)).*exp(-x/8)/3);
        a = 0; % lower limit
        b = 150; % upper limit
        YLimits = [-25 1]; % function value limits for plotting
    %----------------------------------------------------------------------
    case 3
        fun = @(x) -x.^2.*abs(1+sin(x)).*exp(-x/10);
        a = 0; % lower limit
        b = 150; % upper limit
        YLimits = [-115 1]; % function value limits for plotting
    %----------------------------------------------------------------------
end

% potential decrease condition parameter
epsilon = 1e-4;

% maximum number of iterations
maxiter = 30;

% maximum number of function calls
maxfunCount = 75;

% (potentially) create figure for visualization
if plotflag
    in1.fun = fun; in1.a = a; in1.b = b;
    plot_helper(plotflag,1,YLimits,in1,[],[]);
end

% initialize iteration counter
k = 0;

% initialize first center point
pt(1).a = a;
pt(1).b = b;
pt(1).c = (a+b)/2; % center value
pt(1).d = (b-a)/2; % half width
pt(1).f = fun(pt(1).c); % function value at the center
pt(1).DivisionLevel = 0;

% initialize
iter = 0; funCount = 1;
in2.output = []; A = a; B = b;

% while not converged
while (funCount < maxfunCount) && (iter < maxiter)

    % increment iteration counter
    iter = iter + 1;

    % create vectors with current point data
    D = [pt(:).d];
    F = [pt(:).f];
    DL = [pt(:).DivisionLevel];

    % determine set of potentially optimal segments
    [IsPotOpt,Ltest] = potentiallyoptimal(D,F,DL,epsilon);

    % (potentially) plot D vs. F
    if plotflag
        in2.D = D; in2.F = F; in2.epsilon = epsilon;
        in2.IsPotOpt = IsPotOpt; in2.Ltest = Ltest;
        in2.output = plot_helper(plotflag,2,YLimits,[],in2,[]);
    end

    % determine the current number of points
    n = length(pt);

    % preallocate
    pt(end+3*length(IsPotOpt)).f = nan;

    % for each segment, split into thirds
    for k = IsPotOpt(:)'

        % extract current point data
        currentPoint = pt(k);
        a = currentPoint.a;
        b = currentPoint.b;
        f = currentPoint.f;
        DL = currentPoint.DivisionLevel;

        % determine trisection boundaries
        a2 = a*2/3 + b/3;
        b2 = a/3 + b*2/3;

        % divide the segment into thirds as a---a2---b2---b
        pt(n+1) = updatePoint(a,a2,fun,DL); % left
        pt(n+2) = updatePoint(a2,b2,f,DL); % center (use previous f(c))
        pt(n+3) = updatePoint(b2,b,fun,DL); % right

        % increment current point and function call counters
        n = n + 3;
        funCount = funCount + 2; % only two new function calls

    end

% remove old points (they have been trisected)
pt(IsPotOpt) = [];

% display current iteration information
disp(strcat(string(iter),...
    " fmin = ",string(min([pt.f])),...
    " funCount = ",string(funCount)))

% plot newly added points
if plotflag; plot_helper(plotflag,3,YLimits,[],[],pt); end

end

% determine best point found
[F,I] = min([pt.f]);
X = pt(I).c;

% display final value
disp(strcat("x = ",string(X)))
disp(strcat("f = ",string(F)))

%--- optional comparisons
% display value found using equidistant sampling
[Feq,Xeq] = bestEquidistant(fun,A,B,funCount);
disp(strcat("x_eq = ",string(Xeq)))
disp(strcat("f_eq = ",string(Feq)))

% try to find a better minimizer in the final segment
[X_,F_] = refineMinimizer(fun,pt(I));
disp(strcat("x_best = ",string(X_)))
disp(strcat("f_best = ",string(F_)))

abs(X-X_)
abs(F-F_)

abs(Xeq-X_)
abs(Feq-F_)

%--------------------------------------------------------------------------
% update point for the new segment boundaries
function pt = updatePoint(a,b,f,DivisionLevel)

% set endpoints of the segment
pt.a = a;
pt.b = b;

% calculate midpoint of the segment
pt.c = (a+b)/2;

% calculate width of the segment
pt.d = (b-a)/2;

% determine midpoint value
if isnumeric(f)
    pt.f = f;
else
    pt.f = f(pt.c);
end

% increment division level
pt.DivisionLevel = DivisionLevel + 1;

end

%--------------------------------------------------------------------------
% determine the potentially optimal segments
function [K,Ltest] = potentiallyoptimal(D,F,DLevel,epsilon)

% update minimum function value found
fmin = min(F);

% current number of points
N = length(D);

% initialize all points as potentially optimal
IsPotOpt = true(N,1);

% determine potential decrease condition Lipschitz constant value
[Ltest,I] = min((F - fmin + epsilon*abs(fmin))./D);

% DL min
DLmin = DLevel(I);

% go through each point
for j = 1:length(D)

    % check if the point passes the potential decrease condition
    if DLevel(j) > DLmin
        IsPotOpt(j) = false; % not potentially optimal
        continue
    end

    % check if the point passes the domination condition
    if any( (F(j) - F > 0) & (DLevel(j) - DLevel >= 0) )
        IsPotOpt(j) = false; % not potentially optimal
        continue
    end

end

% determine the indices of the potentially optimal points
K = find(IsPotOpt);

end

%--------------------------------------------------------------------------
% plot helper (to separate the plotting code from the algorithm)
function output = plot_helper(plotflag,flag,YLimits,in1,in2,in3)

% check if we should plot things
if ~plotflag
    output = [];
    return
end

% example-specific parameters
dmin = 1e-8;
N = 1e4;

% formatting parameters
pauseTime = 1; % seconds
MarkerSize = 14;
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;

switch flag
    %----------------------------------------------------------------------
    case 1
        % initialize figure
        hf = figure; hf.Color = 'w';
        hf.Position(3:4) = [840 420];

        % initialize subplot 2
        subplot(1,2,2); hold on
        ha = gca; ha.LineWidth = 1; ha.XColor = 'k'; ha.YColor = 'k';
        ha.FontSize = 16; ha.XScale = 'log';
        xlabel('$d$','interpreter','latex')
        ylabel('$f$','interpreter','latex')
        ylim(YLimits)
        xlim([dmin in1.b-in1.a])

        % initialize subplot 1
        subplot(1,2,1); hold on
        ha = gca; ha.LineWidth = 1; ha.XColor = 'k'; ha.YColor = 'k';
        ha.FontSize = 16;
        xlabel('$x$','interpreter','latex')
        ylabel('$f$','interpreter','latex')
        ylim(YLimits)
        xlim([in1.a in1.b])

        % plot the function on a fine grid
        X = linspace(in1.a,in1.b,N);
        F = in1.fun(X);
        plot(X,F,'Color','k','LineWidth',1)
    %----------------------------------------------------------------------
    case 2
        % select subplot 2
        subplot(1,2,2);

        % clear old lines and change old points to light gray
        if ~isempty(in2.output)
            delete(in2.output(2:3))
            in2.output(1).Color = 0.9*[1 1 1];
        end

        % display all points
        hp(1) = plot(in2.D,in2.F,'.','Color','k','MarkerSize',MarkerSize);

        % display potentially optimal points
        hp(2) = plot(in2.D(in2.IsPotOpt),in2.F(in2.IsPotOpt),'.','Color',nicered,'MarkerSize',MarkerSize);

        %
        Dv = logspace(log10(dmin),log10(max(in2.D)),1000);
        hp(3) = plot(Dv,min(in2.F)-in2.epsilon*abs(min(in2.F)) + in2.Ltest*Dv,...
            '-','Color',niceblue,'LineWidth',1);

        % save lines (for deletion later)
        output = hp;

        % pause for some time
        pause(pauseTime)
    %----------------------------------------------------------------------
    case 3
        % extract data
        C = [in3.c];
        F = [in3.f];

        % select subplot 1
        subplot(1,2,1);

        % display currently evaluated points
        plot(C,F,'.','Color',nicered,'MarkerSize',14)

        % display line segments for currently evaluated points
        h = abs(diff(YLimits));
        plot([C;C],[YLimits(1),YLimits(1)+h*0.05]'*ones(1,length(in3)),...
            '-','Color',nicegreen,'LineWidth',0.75)
    %----------------------------------------------------------------------
end

end

%--------------------------------------------------------------------------
% find best function value using equidistant sampling
function [F,X] = bestEquidistant(fun,a,b,N)

% equidistant sampling
X = linspace(a,b,N);

% find minimum value
[F,I] = min(fun(X));

% get x value
X = X(I);

end

%--------------------------------------------------------------------------
% try to find a minimizer in the optimal segment
% [X_,F_] = refineMinimizer(fun,pt(I)) % where I is the best segment index
function [X,F] = refineMinimizer(fun,pt)

% set options
options = optimset();
options.Display = 'none';
options.TolX = 0;

% find minimum using fminbnd
[X,F] = fminbnd(@(x) fun(x),pt.a,pt.b,options);

end