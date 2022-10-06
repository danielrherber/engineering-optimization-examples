close all; clear; clc

test = 2;

switch test
    case 1 % Example 11.9 in Linear and Nonlinear Optimization
        syms x alpha % initialize symbolic variables
        f(x) = 5 - x - log(4.5-x); % function of interest
        xk = 0; % current point
        pk = 1; % search direction (descent only)
        l = 2; r = 4; % current bracket [l,r] should contain a minimum
        F(alpha) = f(xk + alpha*pk); % line search function
    case 2
        syms x1 x2 alpha % initialize symbolic variables
        x = [x1;x2];
        f(x) = 100*(x2-x1^2)^2 + (x1-1)^2; % function of interest
        xk = [-5;60]; % current point
        pk = [1;-10]; % search direction (descent only)
        l = 0; r = 8; % current bracket [l,r] should contain a minimum
        xk1 = xk + alpha*pk; % formula for the next point
        F(alpha) = f(xk1(1),xk1(2)); % line search function
end

% Wolfe condition parameter
eta = 0.1;

% calculate directional derivative
dF = diff(F,alpha,1);

% create Matlab functions
F_ = matlabFunction(F);
dF_ = matlabFunction(dF);

% store some things for later
info.xk = xk; info.pk = pk; info.n = length(xk);
info.f = f; info.F_ = F_; info.dF_ = dF_;
info.plotflag = true;

% evaluate initial endpoints
F_l = F_(l); dF_l = dF_(l); % left endpoint
F_r = F_(r); dF_r = dF_(r); % right endpoint

% save directional derivative value for the Wolfe condition
dF_0 = dF_l;

% initialize the line search algorithm
findingAlpha = true;
iter = 0;
disp(strcat("F(xk) = ",mat2str(F_(0))))
if info.plotflag, plot_initial, end

% search for a suitable alpha value
while findingAlpha

    % increment iteration number
    iter = iter + 1;
    disp(strcat("iter = ",string(iter)))

    % determine candidate alpha using the cubic polynomial
    alpha_hat = cubic_alpha(l,r,F_l,F_r,dF_l,dF_r,info,iter);
    disp(strcat("alpha = ",string(alpha_hat)))

    % evaluate the directional derivative at the candidate point
    dF_alpha_hat = dF_(alpha_hat);

    % check if the Wolfe condition is satisfied; otherwise, update bracket
    if abs(dF_alpha_hat) <= eta*abs(dF_0)

        findingAlpha = false;
        disp("Wolfe condition satisfied")
        break

    end

    % determine which part of the bracket should be updated
    if dF_alpha_hat < 0 % smaller value

        l = alpha_hat;
        F_l = F_(l); dF_l = dF_(l); % left endpoint

    else % larger value

        r = alpha_hat;
        F_r = F_(r); dF_r = dF_(r); % right endpoint

    end
    disp(strcat("bracket now [",string(l)," ",string(r),"]"))

end

% take the step
xk1 = xk + alpha_hat*pk;
disp(strcat("xk1 = ",mat2str(xk1)))
disp(strcat("F(xk1) = ",mat2str(F_(alpha_hat))))

% plot this point in the 2-d case
if info.plotflag, plot_final(info,xk1), end

% determine candidate step length using a specific cubic polynomial
function alpha = cubic_alpha(l,r,Fl,Fr,Fpl,Fpr,info,iter)

% cubic coefficients
a = (Fpl + Fpr - 2*(Fr - Fl)/(r - l))/(l - r)^2;
b = (Fpr - Fpl)/(r - l)/2 - 3/2*(l + r)*a;
c = Fpl - 3*l^2*a - 2*l*b;
d = Fl - l^3*a - l^2*b - l*c;

% minimum value of the cubic in [l,r]
alpha = (-b + sqrt(b^2 - 3*a*c))/(3*a);

% create a visualization of the cubic (1d or 2d case)
if info.plotflag, myplot(l,r,a,b,c,d,alpha,info,iter), end

end

% plotting function for this example;
function myplot(l,r,a,b,c,d,alpha_hat,info,iter)

% colors
niceblue = [77, 121, 167]/255;
nicegreen = [109, 195, 80]/255;
nicegray = [110, 110, 110]/255;

% vector of alpha values to plot
A = linspace(l,r,1000);

% cubic polynomial function
P3 = @(alpha) a*alpha.^3 + b*alpha.^2 + c*alpha.^1 + d*alpha.^0;

if info.n == 2

    % plot original function
    subplot(1,2,1); hold on
    ha = gca;
    ha.LineWidth = 1; ha.FontSize = 16;
    ha.XColor = 'k'; ha.YColor = 'k';
    xlabel('$x_1$')
    ylabel('$x_2$')

    % plot original function and current point
    if iter == 1
        f_ = matlabFunction(info.f);
        [X1,X2] = meshgrid(linspace(-7,7,60),linspace(-30,70,60));
        contourf(X1,X2,f_(X1,X2),50)
        plot(info.xk(1),info.xk(2),'.','MarkerSize',18,'Color','k')
    end

    % plot current bracket and line search line
    Y = info.xk + A.*info.pk;
    gcolors = gray(8);
    plot(Y(1,:),Y(2,:),'-','LineWidth',2,'Color',gcolors(iter+1,:))

    % select next axis
    subplot(1,2,2); hold on

end

% plot line search function F
ha = gca;
ha.LineWidth = 1; ha.FontSize = 16;
ha.XColor = 'k'; ha.YColor = 'k';
xlabel('$\alpha$')
ylabel('$F(\alpha)$')

% plot F(alpha)
if iter == 1
    plot(A,info.F_(A),'linewidth',2,'Color',nicegray)
    xlim([l r])
end

% plot current cubic polynomial
plot(A,P3(A),'linewidth',2,'Color',niceblue)

% plot minimum of the cubic, alpha_hat
plot(alpha_hat,P3(alpha_hat),'.','Color',nicegreen,'MarkerSize',24)

end

% create the initial figure and change some settings
function plot_initial

hf = figure(1); hold on
hf.Color = 'w';
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

end

% plot the final step taken
function plot_final(info,x)

% only supports the 2d case for now
if info.n == 2

    % colors
    nicered = [225, 86, 86]/255;

    % select the correct subplot
    subplot(1,2,1); hold on

    % plot the point
    plot(x(1),x(2),'.','MarkerSize',18,'Color',nicered)

end

end