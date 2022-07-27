close all; clear; clc

% minimize the function from Session 7
f = @(x) exp(0.5*x-1).*(x+1).^2;
% f = @(x) (x+1).^2; % same result

% initial bracket
a = -8;
b = 1;

% stopping tolerance (bracket size
tolerance = 1e-8;

% maximum number of iterations
max_iterations = 100;

% golden ratio
phi = (1 + sqrt(5))/2;

% create plot
plot_helper(f,a,b)

% go through each iteration
for k = 1:max_iterations

    % optimality check
    if (b-a) < tolerance
        disp("Optimal!")
        break
    end

    % compute interior points
    c = b - (b-a)/phi;
    d = a + (b-a)/phi;

    % plot interior points
    hc = plot(c,f(c),'r.','markersize',30);
    hd = plot(d,f(d),'g.','markersize',30);

    % check which one contains the minimum
    if f(c) < f(d)

        % point to remove
        q = b;

        % update
        b = d;

    else

        % point to remove
        q = a;

        % update
        a = c;

    end

    % update plot with removed point
    plot(q,f(q),'.','markersize',30,'Color',[0.5 0.5 0.5]);

    % pause
    pause

end

% midpoint is the optimal point
x = (a+b)/2;
disp(x)
disp(k)

%% helper functions
function plot_helper(f,a,b)

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
hf = figure; hf.Color = 'w'; hold on
ha = gca; ha.LineWidth = 1; ha.FontSize = 18;
x = linspace(a,b,1e5);
plot(x,f(x),'k-','LineWidth',2)
xlabel('$x$'); ylabel('$f(x)$');
plot(a,f(a),'b.','markersize',30)
plot(b,f(b),'b.','markersize',30)

end