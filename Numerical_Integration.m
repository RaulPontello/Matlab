clear,clc
close all;

%% Finite Difference Method (FDM) %%

% Summary: This code determines de integral of the function y = f(t) using rectangle, trapezoidal and Simpson method
% Objective: Determine the integral of the function f(t) in the interval [a , b], f(t) = F'(t)
% How to use : Fill F(t) in item 1.1, fill M, C, K, x0 and v0 in item 1.2, then fill t0, tf and dt in item 1.3. Then run

%% 1) Inputs %% 

%-------------------------------------------------------------------------%
% 1.1) Definition of the function y = f(t) %
%-------------------------------------------------------------------------%

f = @(t)  t*t * exp(t^2);

%-------------------------------------------------------------------------%
% 1.2) Definition of the interval [a , b] %
%-------------------------------------------------------------------------%

a = 0;
b = 2;

%-------------------------------------------------------------------------%
% 1.3) Definition of the interval mesh t = (t1, t2,..., tn) %
%-------------------------------------------------------------------------%

dt = 1;   % dt = size of each division, dt <= b - a
n = (b - a) / dt;   % n = number of divisions of [a , b]
t = a : dt : b;   % t is the time vector


%% 2) Outputs %% 

%-------------------------------------------------------------------------%
% 2.1) Rectangle rule %
%-------------------------------------------------------------------------%

% F(t) = (b - a) * f(0.5 * (a + b))
% n divisions --> n rectangles

y = f(t);   % y = f(ti)
Ar = zeros(n , 1);   % Ar is the area using the Rectangle rule
Ar(1 , 1) = (b - a) * f(0.5 * (a + b));   % initial contidion for the area

for i = 1 : n
    Ar(i , 1) = Ar(i , 1) + dt * y(i , j);
end

%-------------------------------------------------------------------------%
% 2.2) Trapezoidal rule %
%-------------------------------------------------------------------------%

% F(t) = 0.5 * (b - a) * (f(a) + f(b)) --> F(t) = dt * y(t) --> F(ti) = dt *

At = zeros(n , 1);   % Ar is area using the Trapezoidal rule
At(1 , 1) = 0;   % initial contidion for the area

for i = 2 : n
    At(i , 1) = At(i - 1 , 1) + 0.5 * dt * (f(t(i - 1 , 1) + f(t(i , 1))));
end

%-------------------------------------------------------------------------%
% 2.3) Simpson rule %
%-------------------------------------------------------------------------%

% F(t) = 0.5 * (b - a) * (f(a) + f(b))

As = zeros(n , 1);   % Ar is area using the Simpson rule
As(1 , 1) = 0;   % initial contidion for the area

for i = 1 : 2 : n - 2
    As(i , 1) = As(i - 1 , 1) + (dt / 3) * (f(t(i , 1)) + 4 * f(t(i + 1 , 1)) + f(t(i + 2 , 1)));
end

%% 3) Results %% 

fprintf('Beginning of iteration \n\n');
 
for i = 1 : n
    fprintf('i = %1.0f ; t = %1.3f ; Ar = %1.3f ; At = %1.3f ; As = %1.3f  \n' , i , t(i , 1) , Ar(i , 1) , At(i , 1) , As(i , 1));
end

fprintf('\nIteration is over');



