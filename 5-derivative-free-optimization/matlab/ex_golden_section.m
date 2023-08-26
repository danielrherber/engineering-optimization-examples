% ex_golden_section.m
% illustration of golden-section search bracketing for minimizing a 1d
% function
% [course] Session 12 - Constrained Optimization (4) and Derivative-free
% Optimization (1)
close all; clear; clc

% minimize the function from Session 7
f = @(x) exp(0.5*x-1).*(x+1).^2;
% f = @(x) (x+1).^2;

% initial bracket
a = -8;
b = 1;

% stopping tolerance (bracket size)
tolerance = 1e-4;

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

    % compute interior function values
    % (note that this is not done efficiently)
    fc = f(c);
    fd = f(d);

    % check which one contains the minimum
    if fc < fd
        r = b; % point to remove
        b = d; % update
        plot_helper_update(1,f,c,d,r) % update figure
    else
        r = a; % point to remove
        a = c; % update
        plot_helper_update(2,f,c,d,r) % update figure
    end

    % display iteration information
    disp_helper("--- iteration",k,[])
    disp_helper("[a,b]",[a,b],[])
    disp_helper("best f(x)",min(fc,fd),[])

    pause(1)

end

%--------------------------------------------------------------------------
function plot_helper(f,a,b)

% colors
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;

% create plot
hf = figure; hf.Color = 'w'; hold on
ha = gca; ha.LineWidth = 1; ha.FontSize = 18;
xlabel('$x$','Interpreter','latex');
ylabel('$f(x)$','Interpreter','latex');

% plot function and points
x = linspace(a,b,1e5);
plot(x,f(x),'k-','LineWidth',2)
plot(a,f(a),'b.','markersize',30,'color',nicered)
plot(b,f(b),'b.','markersize',30,'color',niceblue)

end

%--------------------------------------------------------------------------
function plot_helper_update(flag,f,c,d,r)

% colors
niceblue = [77, 121, 167]/255;
nicered = [225, 86, 86]/255;
nicegreen = [109, 195, 80]/255;
nicegray = [110, 110, 110]/255;

% plot new interval
switch flag
    case 1
        plot(c,f(c),'.','markersize',30,'color',nicegreen);
        plot(d,f(d),'.','markersize',30,'color',niceblue);
    case 2
        plot(c,f(c),'.','markersize',30,'color',nicered);
        plot(d,f(d),'.','markersize',30,'color',nicegreen);
end

% gray out old point
plot(r,f(r),'.','markersize',30,'Color',nicegray);

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