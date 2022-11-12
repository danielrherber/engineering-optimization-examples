% Problem is from the Example 15.8 of Griva - Linear and Nonlinear Optimization - p.575
% minimize f(x_1,x_2) = e^(3*x_1) + e^(-4*x_2)
% subject to g(x_1,x_2)= x_1^2 + x_2^2  = 1
%
clear all; close all; clc

%% create functions and derivatives
% symbolic functions and derivatives
syms x1 x2 real
global history Q  % Global variables are declared to be used in the corresponding functions

% Keeping the log of the optimValues from OutputFcn or optimoptions
history.x = [];
history.fval = [];
history.searchdir = [];

x_ = [x1;x2]; % local symbolic variable with underscore

Q = eye(2);   % you can change this matrix
 
f = exp(3*x1) + exp(-4*x2);
% To convert symbolic function with the variables of x vector by using function handle
F = matlabFunction(f,'vars',{x_});  

constraints = @consFunc; % or directly using @consfun in associated argument in fmincon

x0 = [-1; 1];

options = optimoptions('fmincon','Algorithm','sqp','Display','iter','OutputFcn',@outfun);
[x, fval, exitflag, output, lambda] = fmincon(F, x0, [], [], [], [], [], [], constraints, options);

function [c, ceq] = consFunc(x_)
    global Q
    c = [];
    ceq = x_'*Q*x_ - 1; % Equality constraint g(x_1,x_2) = x_1^2 + x_2^2  = 1 and x = [x_1, x_2]
end

function stop = outfun(x,optimValues,state)
    global history
     stop = false;

     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval, optimValues.fval];
           history.x = [history.x, x];
         % Concatenate current search direction with searchdir
           history.searchdir = [history.searchdir, optimValues.searchdirection'];
           plot(x(1),x(2),'*');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'.
           text(x(1)+.15,x(2), num2str(optimValues.iteration));
           title('Sequence of Points Computed by fmincon');

         case 'done'
             hold off
         otherwise
     end
end

% References
% Using Symbolic Mathematics with Optimization Toolbox Solvers: 
% https://www.mathworks.com/help/optim/ug/use-symbolic-math-with-optim.html#SymbolicOptimization-2
% fmincon: 
% https://www.mathworks.com/help/optim/ug/fmincon.html
