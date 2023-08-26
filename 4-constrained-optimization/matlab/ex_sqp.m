% ex_sqp.m
% illustration of the sequential quadratic programming (SQP) method for a
% a 2d problem with a single quadratic equality constraint
% [reference] Section 15.5 in LNO
% [course] Session 11 - Constrained Optimization (3)
close all; clear; clc

%--------------------------------------------------------------------------
% symbolic functions (Example 15.8 in LNO)
syms x1 x2 l real
x =  [x1;x2];
f = exp(3*x1) + exp(-4*x2); % you can change this function
Q = eye(2); % you can change this matrix
g = x'*Q*x - 1;

% Lagrangian
L = f - l*g;

% determine symbolic derivatives
df = gradient(f,x);
d2f = hessian(f,x);
dg = jacobian(g,x)';
dL = gradient(L,x);
d2L = hessian(L,x);

% convert to matlab functions
F = matlabFunction(f);
dF = matlabFunction(df);
d2F = matlabFunction(d2f);
G = matlabFunction(g);
dG = matlabFunction(dg);
dL = matlabFunction(dL);
d2L = matlabFunction(d2L);

%--------------------------------------------------------------------------
% setup
test = 1; % see below

switch test
    case 1
        % initial point
        x = [-1; 1]; l = -1;
    case 2
        % initial point
        x = [-0.75; 0.1]; l = 1;
    case 3
        % initial point
        x = [-1; 1]; l = 0;
end

% plot quadratic model (in 3d)
modelflag = true;

% tolerances
ConstraintTolerance = 1e-10; % constraint tolerance
OptimalityTolerance = 1e-10; % optimality tolerance
MaxIterations = 100; % maximum number of iterations

%--------------------------------------------------------------------------
%--- Sequential Quadratic Programming (SQP) method
% problem information
n = length(x); % number of variables
m = length(l); % number of constraints

% merit penalty parameter
rho = 10;

% initialize
hs = [];

% create plot
plot_helper_1(F,Q)

% go through each iteration
for k = 1:MaxIterations

    % values for Newton equations
    d2L_ = d2L(l,x(1),x(2)); % hessian matrix of Lagrangian with respect to x
    dG_ = dG(x(1),x(2)); % jacobian
    dL_ = dL(l,x(1),x(2)); % gradient of Lagrangian with respect to x
    g_ = G(x(1),x(2)); % constraint value
    F_ = F(x(1),x(2)); % function value

    % display current iteration values
    disp(" ") % spacer
    disp_helper("*iteration",k-1,[])
    disp_helper("x",x,[])
    disp_helper("L",l,[])
    disp_helper("norm(dL)",norm(dL_),100)
    disp_helper("norm(g)",norm(g_),100)
    disp_helper("M",F_ + rho*g_'*g_,6) % merit function

    % check termination conditions
    if norm(g_) <= ConstraintTolerance
        if norm(dL_) <= OptimalityTolerance
            disp("Optimal!")
            break
        end
    end

    % matrices for Newton equations
    A = [d2L_ -dG_; -dG_' zeros(m)];
    B = [-dL_; g_];

    % solve Newton equations
    sol = A\B;

    % extract parts
    p = sol(1:n);
    v = sol(n+1:end);

    % plot current iteration
    hs = plot_helper_2(x,p,dL_,d2L_,g_,dG_,modelflag,hs);

    % new estimates of the solution
    x = x + p; % variables
    l = l + v; % Lagrange multipliers

    pause

end

%--------------------------------------------------------------------------
% plot helper 1
function plot_helper_1(F,Q)

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
hf = figure; hf.Color = 'w'; hold on
N = 1e3;
x1 = linspace(-1.5,0,N);
x2 = linspace(0,1.5,N);
xlim([min(x1),max(x1)]); ylim([min(x2),max(x2)]); % axis limits
[X1,X2] = meshgrid(x1,x2);
F_ = F(X1,X2);
contour(X1,X2,F_,80)
plot_ellipse(Q)
xlabel('$x_1$'); ylabel('$x_2$'); % axis labels
ha = gca;
ha.LineWidth = 1; ha.FontSize = 18;

end

%--------------------------------------------------------------------------
% plot helper 2
function hs = plot_helper_2(x,p,dL_,d2L_,g_,dG_,modelflag,hs)

% plot current point
plot(x(1),x(2),'.','markersize',30,'color',[244, 67, 54]/255)

% compute radius to focus plot
p0 = norm(p,2);
p0 = max(0.2,p0); % minimum radius

% create grid of points around current point
p0_ = linspace(-p0,p0,1e3);
[p1,p2] = meshgrid(p0_,p0_);

% values of the quadratic model
q = zeros(size(p1));
for idx = 1:numel(p1)
    p_ = [p1(idx); p2(idx)];
    q(idx) = 0 + dL_'*p_ + 0.5*p_'*d2L_*p_;
end

% ignore values outside radius
I = p1.^2 + p2.^2 > p0^2;
q(I) = nan;

% ignore values outside plot limits
ha = gca;
I = ha.XLim(1) > p1 + x(1);
q(I) = nan;
I = ha.XLim(2) < p1 + x(1);
q(I) = nan;
I = ha.YLim(1) > p2 + x(2);
q(I) = nan;
I = ha.YLim(2) < p2 + x(2);
q(I) = nan;

% manual z limit
zlim([0 1.25])

% scale the values of the quadratic approximation to be between 0 and 1
S = min(q,[],"all");
q = q - S;
s = max(abs(q),[],"all");
q = q/s;

% delete old model
delete(hs)

% create surface plot
if modelflag

end

% create line for linear approximation of the constraint
p1_ = linspace(-1,1,500);
p2_ = -(g_ + dG_(1)*p1_)/dG_(2);

% go through each point
qc = nan(size(p1_));
for idx = 1:numel(p1_)
    p_ = [p1_(idx); p2_(idx)];
    qc(idx) = 0 + dL_'*p_ + 0.5*p_'*d2L_*p_;
end

% plot stuff
hs(1) = plot(x(1)+p1_,x(2)+p2_,'-','linewidth',2,'color',[129, 199, 132]/255);
hs(2) = plot(x(1)+p(1),x(2)+p(2),'.','markersize',30,'color',[3, 169, 244]/255); % current point

if modelflag
    hs(3) = surf(x(1)+p1,x(2)+p2,q,q,'LineStyle','none','FaceAlpha',0.75);
    hs(4) = plot3(x(1)+p1_,x(2)+p2_,(qc - S)/s,...
        '-','linewidth',2,'color',[129, 199, 132]/255); % constraint line
    hs(5) = plot3(x(1)+p(1),x(2)+p(2), (0 + dL_'*p + 0.5*p'*d2L_*p - S)/s,...
        '.','markersize',30,'color',[3, 169, 244]/255); % Lagrangian on constraint line
    view(gca,[-22.6593443093639 42.8555513022069]); % change view
end

end

%--------------------------------------------------------------------------
% plot the quadratic function x'*Q*x = 1
function plot_ellipse(Q_)

% symbolically solve for x2 given x1
syms x1 x2 a b c d real
x = [x1;x2];
Q = [a b; c d];
g = x'*Q*x;
sol = solve(g==1,x2);
SOL = matlabFunction(sol);

% compute valid x2 given x1 (might be a complex number)
x1 = linspace(-4,4,1e3)';
for k = 1:numel(x1)
    x2_(k,:) = SOL(Q_(1,1),Q_(1,2),Q_(2,1),Q_(2,2),x1(k));
end

% create vectors with only real parts
x2_ = [x2_(:,1);flipud(x2_(:,2))];
x1_ = [x1;flipud(x1)];
I = abs(imag(x2_)) > 0;
x1_(I) = [];
x2_(I) = [];
x1_ = [x1_;x1_(1)];
x2_ = [x2_;x2_(1)];

% plot ellipse
plot(x1_,x2_,'-','linewidth',2,'Color',[56, 142, 60]/255)

end

%--------------------------------------------------------------------------
% function to make it easier to display things in the command window
function disp_helper(name,number,varargin)

% default value of the number of digits
if isempty(varargin{:})
    n = 5;
else
    n = varargin{1};
end

% form string
str = strcat(string(name)," = ",mat2str(round(number,n)));

% display string
disp(str)

end