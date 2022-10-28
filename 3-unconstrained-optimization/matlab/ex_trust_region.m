close all; clear; clc

test = 1;

switch test
    %----------------------------------------------------------------------
    case 1
        syms x1 x2
        x = [x1; x2];
        f = x1^4 + 2*x1^3 + 24*x1^2 + x2^4 + 12*x2^2;
        g = gradient(f);
        h = hessian(f);
        x0 = [2; 1]; % initial guess
        limits = [-1 3 -1 1.5]; % plot limits
    %----------------------------------------------------------------------
    case 2
        syms x1 x2
        x = [x1; x2];
        f = x1^4 + 2*x1^3 + 24*x1^2 + x2^4 + 12*x2^2 + x1^6*x2^4;
        g = gradient(f);
        h = hessian(f);
        x0 = [2.5; 1]; % initial guess
        limits = [-1 3 -1 1.5]; % plot limits
end

% number of iterations
n = 4;

% do you want to visualize the trust region model?
modelflag = true;

% initial trust region size
D0 = 1;

% algorithm parameters
eta = 0.75;
mu = 0.25;

% create functions
F_ = matlabFunction(f);
G_ = matlabFunction(g);
H_ = matlabFunction(h);

% plot stuff (see below)
output = plot_helper(1,{F_,modelflag,n,limits});
colors = output{1};

% display initial function value
f0 = F_(x0(1),x0(2));
disp_helper("f(x)",f0)

% initial handles
hs = []; hp_new = [];

% go through each iteration
for k = 1:n

    disp(" ")

    % current function value, gradient, and hessian
    f0 = F_(x0(1),x0(2));
    g0 = G_(x0(1),x0(2));
    h0 = H_(x0(1),x0(2));

    % plot stuff (see below)
    output = plot_helper(2,{x0,colors,k,D0,modelflag,f0,g0,h0,hs,hp_new});
    hs = output{1};

    % compute Newton direction (solving a set of linear equations)
    pN = -h0\g0;

    % check norm of Newton direction
    if norm(pN,2) < D0 % small enough?

        % assign step as the newton step
        p0 = pN;
        disp("Using Newton step")

    else

        % solve for l
        l = fminbnd(@(l) (norm((h0+l*eye(size(h0)))\-g0,2) - D0).^2,...
            0,1e10,optimset('TolX',1e-12));

        % compute trust-region step
        p0 = -(h0+l*eye(size(h0)))\g0;

    end

    % compute, for this step, the trust-region model value
    q0 = f0 + g0'*p0 + 0.5*p0'*h0*p0;

    % compute function value at x0 + p0
    f1 = F_(x0(1)+p0(1),x0(2)+p0(2));

    % compute ratio of actual to predicted reduction
    rho0 = (f0 - f1)/(f0 - q0);

    % display function value and prediction ratio
    disp_helper("f(x)",f1)
    disp_helper("rho",rho0)

    % determine if the step was successful or not
    if norm(f0 - q0) < eps
        disp("f and q the same")
        x0 = x0 + p0;
        disp("Stopping")
        break
    elseif rho0 <= mu
        D0 = 0.5*D0; % decrease
        disp("Trust region size decreased")
        disp("Step unsucessful")
    elseif rho0 >= eta
        D0 = 2*D0; % increase
        disp("Trust region size increased")
        x0 = x0 + p0;
        disp("Step sucessful")
    else
        disp("Trust region size the same")
        x0 = x0 + p0;
        disp("Step sucessful")
    end

    % plot stuff (see below)
    output = plot_helper(3,{modelflag,x0,q0,colors,k,p0,pN,g0});
    hp_new = output{1};

end

% function to make it easier to display things in the command window
function disp_helper(name,number)

% form string
str = strcat(string(name)," = ",mat2str(round(number,5)));

% display string
disp(str)

end

% plot circle function
function plotcircle(r,x,y,c)

    % parameter between 0 and 2pi
    th = 0:pi/100:2*pi;

    % complex-valued circle equation
    f = r * exp(1i*th) + x+1i*y;

    % plot
    plot(real(f), imag(f),'color',c,'linewidth',2);

end

% plot helper (to separate the plotting code from the algorithm)
function output = plot_helper(flag,data)

    switch flag

    %----------------------------------------------------------------------
    case 1

    % extract
    F_ = data{1};
    modelflag = data{2};
    n = data{3};
    limits = data{4};

    % change to latex
    set(groot,'defaulttextinterpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');

    % create a grid of points
    N = 600;
    x1 = linspace(limits(1),limits(2),N)';
    x2 = linspace(limits(3),limits(4),N)';
    [X1,X2] = meshgrid(x1,x2);

    % create figure
    hf = figure; hf.Color = 'w'; hold on

    % create contour plot
    F__ = F_(X1,X2);
    F__(F__ > 600) = nan; % limit the data
    [~,hc] = contourf(X1,X2,F__,50);
    % hc.LineStyle = 'none';

    % only if we are NOT visualizing the quadratic model
    if ~modelflag
        axis equal % so the circles look like circles
    end

    % limits
    xlim([min(X1,[],"all") max(X1,[],"all")]); % x limits
    ylim([min(X2,[],"all") max(X2,[],"all")]); % y limits

    % change font size
    ha = gca;
    ha.FontSize = 20;

    % colors for each iteration
    colors = flipud(summer(n));
    output{1} = colors;

    %----------------------------------------------------------------------
    case 2

    % extract
    x0 = data{1};
    colors = data{2};
    k = data{3};
    D0 = data{4};
    modelflag = data{5};
    f0 = data{6};
    g0 = data{7};
    h0 = data{8};
    hs = data{9};
    hp_new = data{10};

    % plot current point
    plot(x0(1),x0(2),'.','color',colors(k,:),'markersize',20)

    % plot current trust-region circle
    plotcircle(D0,x0(1),x0(2),colors(k,:))

    % only if we are visualizing the quadratic model
    if modelflag

        % create grid of points around current point
        D0_ = linspace(-D0,D0,600)';
        [D0_1,D0_2] = meshgrid(D0_,D0_);

        % values of the quadratic model
        Q0 = zeros(size(D0_1));
        for idx = 1:numel(D0_1)
            P0 = [D0_1(idx); D0_2(idx)];
            Q0(idx) = f0 + g0'*P0 + 0.5*P0'*h0*P0;
        end

        % ignore values outside trust region
        I = (D0_1).^2 + (D0_2).^2 > D0^2;
        Q0(I) = nan;

        % ignore values outside plot limits
        ha = gca;
        I = ha.XLim(1) > D0_1 + x0(1);
        Q0(I) = nan;
        I = ha.XLim(2) < D0_1 + x0(1);
        Q0(I) = nan;
        I = ha.YLim(1) > D0_2 + x0(2);
        Q0(I) = nan;
        I = ha.YLim(2) < D0_2 + x0(2);
        Q0(I) = nan;

        % limits based on  maximum value
        zlim([min(0,min(Q0(:))) max(Q0(:))])

        % delete old model
        delete(hs)
        delete(hp_new)

        % create surface plot
        hs = surf(x0(1) + D0_1, x0(2) + D0_2,Q0,Q0,...
            'LineStyle','none','FaceAlpha',0.75);
        output{1} = hs;

    else
        output{1} = [];
    end

    %----------------------------------------------------------------------
    case 3

    % extract
    modelflag = data{1};
    x0 = data{2};
    q0 = data{3};
    colors = data{4};
    k = data{5};
    p0 = data{6};
    pN = data{7};
    g0 = data{8};

    % only if we are visualizing the quadratic model
    if modelflag

        % plot next point
        hp_new = plot3(x0(1),x0(2),q0,'.','color',colors(k,:),'markersize',20);

        % assign
        output{1} = hp_new;

    else
        output{1} = [];
    end

    % step direction and length
    quiver(x0(1)-p0(1),x0(2)-p0(2),p0(1),p0(2),0...
        'color',colors(k,:),'linewidth',1.5,'MaxHeadSize',0.6);

    % newton direction
    quiver(x0(1)-p0(1),x0(2)-p0(2),pN(1)/norm(pN),pN(2)/norm(pN),0.45,...
        'color',colors(k,:),'linewidth',1,'LineStyle',':');

    % negative gradient direction
    quiver(x0(1)-p0(1),x0(2)-p0(2),-g0(1)/norm(g0),-g0(2)/norm(g0),0.45,...
        'color',colors(k,:),'linewidth',1,'LineStyle','--');

    % plot new point
    plot(x0(1),x0(2),'.','color',colors(k,:),'markersize',20);

    end

end
