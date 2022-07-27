close all; clear; clc;

% data
n = 100;
rng(345345);
x_data = linspace(-10,10,n);
y_data = 0.1*x_data.*cos(x_data) + 0.1*rand(1,n);

% number of activation functions
N = 6;

% starting point
X0 = rand(3*N,1);

% fminunc optimizer options
OPTIONS = optimoptions('fminunc');
OPTIONS.Display = 'iter';
OPTIONS.MaxFunctionEvaluations = 1e10; 
OPTIONS.MaxIterations = 1e10;
OPTIONS.FiniteDifferenceType = 'central';

% solve the unconstrained optimization problem
[x,f] = fminunc(@(X) objective(X,x_data,y_data,N),X0,OPTIONS);

% extract optimal parameter values
[w,a,b] = get_parameters(x,N);

% get predictions
y_pred = activation(x_data,w,a,b);  

% plot
figure; hold on
plot(x_data,y_data,'.')
plot(x_data,y_pred)

% objective function
function f = objective(x,x_data,y_data,N)

% get parameters
[w,a,b] = get_parameters(x,N);

% predict
y_pred = activation(x_data,w,a,b);

% calculate error
f = sum((y_pred - y_data).^2)/length(y_pred);

end

% get activation function parameters from optimization variables
function [w,a,b] = get_parameters(x,N)
    x = x(:); % ensure column vector
    w = x(1:N); % weights
    a = x(N+1:2*N); % amplitude
    b = x(2*N+1:end); % bias
end

% sum activation functions to compute prediction for y given x_data
function y_pred = activation(x_data,w,a,b)
    y_pred = sum(w.*1./(1+exp(a.*x_data+b)));
%     y_pred = sum(w.*(atan(a.*x_data+b)/pi + 0.5));
%     y_pred = sum(w.*0.5.*(tanh(a.*x_data+b)+1));
end