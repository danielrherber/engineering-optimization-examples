%% Managing a Session
close all % closes all figures
clear % removes all variables from the current workspace
clc % clears all the text from the Command Window

%% Basic Syntax
5+3 % addition
5-3 % subtraction
5*3 % scalar (and matrix) multiplication
3^2 % scalar (and matrix) exponentiation operator
[1 2].*[3 4] % array multiplication operator
[1 2].^[3 4] % array exponentiation operator
[1 2; 3 4]\[1; 2] % left-division, solve systems of linear equations A*x = b for x
5/3 % right-division operator
[1 2].\[3 4] % % array left-division operator
[1 2]./[3 4] % array right-division operator
[1 2]' % matrix transpose; for complex matrices, is the complex conjugate transpose
[1 2].' % array transpose
disp("Hello World") % display something to the command window
help disp % displays the help text for disp (disp could be any function/script)
% : % colon; generates regularly spaced elements and represents an entire row or column
% ( ) % parentheses; encloses function arguments and array indices; overrides precedence
% [ ] % brackets; enclosures array elements
% . % decimal point
% ... % ellipsis; line-continuation operator
% , % comma; separates statements and elements in a row
% ; % semicolon; separates columns and suppresses display
% % % percent sign; designates a comment and specifies formatting
% ' % quote sign and matrix transpose operator
% = % assignment operator

%% Special Variables and Constants
ans % most recent answer
eps % accuracy of floating-point precision
i,j % the imaginary unit sqrt(-1)
Inf % infinity
NaN % undefined numerical result (not a number)
pi % the number π, 3.1416...

%% Relational Operators
1 < 2 % less than
1 <= 2 % less than or equal to
1 > 2 % greater than
1 >= 2 % greater than or equal to
1 == 2 % equal to
1 ~= 2 % not equal to
1 & 0 % bitwise and
1 | 0 % bitwise or
~1 % logical not

%% Variables, Vectors, and Matrices
x = 3 % defining x and initializing it with a value
a = 2; b = 7; c = a*b % multiple assignments
r = [7 8 9 10 11] % row vector
c = [7;  8;  9;  10; 11] % column vector
m = [1 2 3; 4 5 6; 7 8 9] %  matrix (3x3)
whos % list variables in workspace, with sizes and types

%% Referencing the Elements of a Matrix
a = [ 1 2 3 4 5; 2 3 4 5 6; 3 4 5 6 7; 4 5 6 7 8];
a(2,5) % reference the 2nd row and 5th column of a
v = a(:,4) % extract 4th column of a
q = a(2:3,:) % extract 2nd and 3rd rows of a
a(4,:) = [] % delete 4th row of a
a(:,[1 3]) = [] % delete 1st and 3rd columns of a
a(end,:) % last row of a
a(:,end-1:end) % last and 2nd to last columns of a

%% Vector, Matrix, and Array Commands
find([1 0 1]) % finds indices of nonzero elements
length([1 2 3]) % computes number of elements
linspace(0,1,4) % creates regularly spaced vector
logspace(-1,1,5) % creates logarithmically spaced vector
max([1 4 -2]) % returns largest element
min([1 4 -2]) % returns smallest element
prod([1 4 -2]) % product of each column
size([1 2; 3 4]) % computes array size
sort([1 4 -2]) % sorts each column
sum([1 4 -2]) % sums each column
eye(3) % creates an identity matrix
ones(3,1) % creates an array of ones
zeros(3,2) % creates an array of zeros
inv([1 2; 3 4]) % computes inverse of a matrix
pinv([1 2]) % computes pseudoinverse of a matrix
rank([1 1; 1 1]) % computes rank of a matrix

%% Strings and Character Arrays
x = "Hello"; x(1) % string treats each phrase as a unit
y = 'World'; y(1) % char treats each character as a unit
z = ["Hello", "Hi"]; z(1) % string array example
w = {'World', 'Earth'}; w{1} % char array example
strcat(x,y) % concatenate strings horizontally
strsplit(x,"l") % split string at specified delimiter
strfind(y,"or") % find one string within another
strrep(x,"l","!") % find and replace substring
strcmp(x,y) % compare strings (case sensitive)
strcmpi("A","a") % compare strings (case insensitive)

%% Functions
%--- function definition example
% function y = f(x) % actual definition below
%     y = sin(x);
% end
f(1,2) % function definition usage
%--- anonymous function example
g = @(x,p) sin(x) + p; % anonymous function definition
g(1,3) % anonymous function usage
%--- pass a parameter into a function
fzero(@(x) f(x,0.5), 1)

%% Plotting Commands
figure % opens a new figure window
hold on % retains existing objects so that new objects added to figure
plot([1 2],[3 4]) % generates xy plot
title("My Title") % puts text at top of plot
xlabel("x") % adds text label to x-axis
ylabel("y") % adds text label to y-axis
legend("a") % creates a legend
text(1.4,3.8,"text") % places string in figure
exportgraphics(gcf,"figname.png") % save plot or graphics content to file

%% Decision Making and Loops
% break % terminates the loop statement and transfers execution to the statement immediately following the loop
% continue % causes the loop to skip the remainder of its body and immediately retest its condition prior to reiterating
% return % return control to invoking script or function
%--- if, elseif, else statement example
x = 2;
if x > 2
    disp("if")
elseif x == 2
    disp("elseif")
else
    disp("else")
end
%--- switch statement example
x = 2;
switch x
    case 1
        disp("1")
    case 2
        disp("2")
    otherwise
        disp("otherwise")
end
%--- while loop example
x = 3;
while x > 0
    x = x - 1;
    disp(x)
end
%--- for loop example
for k = 1:4
    x = k*k;
    disp(x)
end

%% Calculus (with the Symbolic Toolbox)
syms x y % create symbolic scalar variables and functions, and matrix variables and functions
s = sin(x); % symbolic function
diff(s,x,1) % differentiate symbolic expression or function
int(s,x) % definite and indefinite integrals
simplify(2*x + (x - 1)^2) % algebraic simplification
pretty(s^2/(s+1)) % Prettyprint symbolic expressions
subs(s,x,2) % symbolic substitution
gradient(sin(x*y) + y^2,[x y]) % gradient vector of scalar function
hessian(sin(x*y) + y^2,[x y]) % hessian matrix of scalar function
jacobian([sin(x+y); x*y],[x y]) % jacobian matrix
solve(x^2 + 2*x == 1, x) % equations and systems solver
vpasolve(x^2 + 2*x == 1, x) % solve equations numerically

%% Function Definitions
function y = f(x,p)
    y = sin(x) + p;
end

%% More Information at can be Found at:
% https://www.tutorialspoint.com/matlab/matlab_commands.htm
% https://www.mathworks.com/help/matlab/language-fundamentals.html