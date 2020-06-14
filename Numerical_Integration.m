clear,clc
close all;

%% Numerical Integration %%

%-------------------------------------------------------------------------%
% Introduction %
%-------------------------------------------------------------------------%

% Latest update: 14/06/2020

% Summary: This code determines the integral of the function y = f(t) in the interval [a , b]
% Methods used: Rectangle, Trapezoidal and Simpson
% How to use : 
%   1) Fill f(t) in item 1.1; 
%   2) Fill a and b in item 1.2;
%   3) Fill n in item 1.3; 
%   4) Then run

%% 1) Inputs %% 

%-------------------------------------------------------------------------%
% 1.1) Definition of the function y = f(t) %
%-------------------------------------------------------------------------%

f = @(t)  t.^2 .* exp(t.^2);

%-------------------------------------------------------------------------%
% 1.2) Definition of the interval [a , b] %
%-------------------------------------------------------------------------%

a = 0;
b = 2;

%-------------------------------------------------------------------------%
% 1.3) Definition of the interval mesh t = (t1, t2,..., tN) %
%-------------------------------------------------------------------------%

n = 96;   % n = number of divisions of [a , b]

%% 2) Outputs %% 

%-------------------------------------------------------------------------%
% 2.1) Integration rules %
%-------------------------------------------------------------------------%

% n divisions --> n rectangles

Ar = zeros(n , 1);   % Ar is the area using the Rectangle rule
At = zeros(n , 1);   % Ar is area using the Trapezoidal rule
As = zeros(n , 1);   % Ar is area using the Simpson rule

for i = 1 : n
    dt = (b - a) / i;   % dt = size of each division, dt <= b - a
    t = a : dt : b;   % t is the time vector
    
    % Rectangle rule --> F(t) = dt * f(0.5 * (a + b)) %
    
    tr = a + 0.5 * dt :  dt : b - 0.5 * dt;   % tr is time vector for Rectabgle rule
    yr = f(tr);
    Ar(i , 1) = dt * sum(yr);
    
    % Trapezoidal rule --> F(t) = 0.5 * dt * (f(a) + f(b)) %
    
    tt = t;   % tt is time vector for Trapezoidal rule
    yt = zeros(1 , length(tt));
            
    for j = 1 : length(tt)
        if j == 1
            yt(j) = f(tt(j));
        elseif j == length(tt)
            yt(j) = f(tt(j));
        else
            yt(j) = 2 * f(tt(j));
        end
    end
    
    At(i , 1) = 0.5 * dt * sum(yt);
    
    % Simpson rule --> F(t) = (f(a) + 4 * f(0.5 * (a + b)) + f(b)) * dt / 3 %
    
    ts = a : 0.5 * dt : b;   % ts is time vector for Simpson rule
    ys = zeros(1 , length(ts));
    
    for j = 1 : length(ts)
        if j == 1
            ys(j) = f(ts(j));
        elseif j == length(ts)
            ys(j) = f(ts(j));
        elseif rem(j , 2) == 0
            ys(j) = 4 * f(ts(j));
        elseif rem(j , 2) ~= 0;
            ys(j) = 2 * f(ts(j));
        end
    end 
    
    As(i , 1) = (0.5 * dt / 3) * sum(ys);
end


%% 3) Results %% 

fprintf('Beginning of iteration \n\n');
 
for i = 1 : n
    fprintf('n = %1.0f ; Ar = %1.3f ; At = %1.3f ; As = %1.3f  \n' , i , Ar(i , 1) , At(i , 1) , As(i , 1));
end

fprintf('\nIteration is over');
