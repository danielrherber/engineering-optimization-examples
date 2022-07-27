close all; clear; clc

OPTIONS = optimoptions('fmincon');
OPTIONS.Display = 'Iter';
% OPTIONS.FiniteDifferenceType = 'central';
% OPTIONS.OptimalityTolerance = 1e-16;
% OPTIONS.StepTolerance = 0;
% OPTIONS.ConstraintTolerance = 1e-16;
% OPTIONS.FunctionTolerance = 1e-16;
% OPTIONS.ScaleProblem = true;
% OPTIONS.Algorithm = "active-set";
% OPTIONS.MaxFunctionEvaluations = 1e10;

% X0 = zeros(10,1);
X0 = [4,2,1,4,9,5,3,-2,7,0];

tic
[x,f] = fmincon(@wire_amount,X0,[],[],[],[],[],[],@building_constraints,OPTIONS);
toc

% extract
x0 = x(1); y0 = x(2);
x1 = x(3); y1 = x(4);
x2 = x(5); y2 = x(6);
x3 = x(7); y3 = x(8);
x4 = x(9); y4 = x(10);

hf = figure; hf.Color = 'w'; hold on
axis equal

% plot buildings
circle2(1,4,2) % building 1
circle2(9,5,1) % building 2
rectangle('Position',[2 -3 2 2],'FaceColor',[244, 244, 244]/255) % building 3
rectangle('Position',[6 -2 2 4],'FaceColor',[244, 244, 244]/255) % building 4

% plot wires
plot([x1,x0],[y1,y0],'-')
plot([x2,x0],[y2,y0],'-')
plot([x3,x0],[y3,y0],'-')
plot([x4,x0],[y4,y0],'-')


function h = circle2(x,y,r)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[1,1],'FaceColor',[244, 244, 244]/255);
daspect([1,1,1])
end



function w = wire_amount(x)

% extract
x0 = x(1); y0 = x(2);
x1 = x(3); y1 = x(4);
x2 = x(5); y2 = x(6);
x3 = x(7); y3 = x(8);
x4 = x(9); y4 = x(10);

% individual wire lengths
w1 = sqrt((x1-x0)^2+(y1-y0)^2);
w2 = sqrt((x2-x0)^2+(y2-y0)^2);
w3 = sqrt((x3-x0)^2+(y3-y0)^2);
w4 = sqrt((x4-x0)^2+(y4-y0)^2);

% total wire length
w = w1 + w2 + w3 + w4;

end

function [c,ceq] = building_constraints(x)

% extract
x0 = x(1); y0 = x(2);
x1 = x(3); y1 = x(4);
x2 = x(5); y2 = x(6);
x3 = x(7); y3 = x(8);
x4 = x(9); y4 = x(10);

% building 1
c(1) = (x1-1)^2 + (y1-4)^2 - 4;

% building 2
c(2) = (x2-9)^2 + (y2-5)^2 - 1;

% building 3
c(3) = 2 - x3;
c(4) = x3 - 4;
c(5) = -3 - y3;
c(6) = y3 + 1;

% building 4
c(7) = 6 - x4;
c(8) = x4 - 8;
c(9) = -2 - y4;
c(10) = y4 - 2;

% no equality constraints
ceq = []; % empty

end