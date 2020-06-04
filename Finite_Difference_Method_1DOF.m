clear,clc
close all;

%% Finite Difference Method (FDM) %%

% Summary: This code solves second-order linear differential equations for damped vibrations of a single degree of freedom system using the Finite Difference Method
% Objective: Determine the solution x(t) for the equation M * x''(t) + C * x'(t) + K * x(t) = F(t)
% How to use : Fill F(t) in item 1.1, fill M, C, K, x0 and v0 in item 1.2, then fill t0, tf and dt in item 1.3. Then run

%-------------------------------------------------------------------------%
% Units used %
%-------------------------------------------------------------------------%

% Mass (M) = kg
% Damping coefficient (C) = Ns/m
% Spring coefficient (K) = N/m
% Force (F) = N
% Time (t) = s
% Displacement (x) = m
% Velocity (v) = m/s
% Acceleration (a) = m/s²
% Natural frequency (Wn) = rad/s
% Frequency (f) = Hz
% Period (T) = s

%% 1) Inputs %% 

%-------------------------------------------------------------------------%
% 1.1) Definition of the external force F(t) %
%-------------------------------------------------------------------------%

F = @(t)  0 ;   % F = f(t)

%-------------------------------------------------------------------------%
% 1.2) Definition of the (M, C, K) parameters and initial conditions %
%-------------------------------------------------------------------------%

M = 10 ; C = 35 ; K = 5;   % M = mass , C = damping coefficient , K = spring coefficient
x0 = 15 ; v0 = 0;   % initial conditions for displacement (x0) and velocity (v0)

%-------------------------------------------------------------------------%
% 1.3) Definition of the time mesh t = (t0, t1,..., tn) %
%-------------------------------------------------------------------------%

% n = number of divisions of (t0 , tf) , dt = size of each division
% dt = (tf - t0) / n

t0 = 0 ; tf = 100 ; dt = 0.01; 
n = (tf - t0) / dt;
t = zeros(n , 1);   % t is the time vector
t(1 , 1) = t0;

for i = 2 : n
    t(i , 1) = t(1 , 1) + (i - 1) * dt;
end

%% 2) Outputs %% 

%-------------------------------------------------------------------------%
% 2.1) Definition of the Wn, Tn and f %
%-------------------------------------------------------------------------%

Wn = sqrt(K / M); T = 2 * pi / Wn ; f = 1 / T;   % Wn = natural frequency , T = period , f = frequency
Cc = 2 * M * Wn;   % Cc = critical damping
z = C / Cc;   % z = damping ratio
dt_crit = T / pi;   % dt < dt_crit --->> dt < T/pi

%-------------------------------------------------------------------------%
% 2.2) Definition of the recurrence equation for x(t) = (x0, x1,...,xn) %
%-------------------------------------------------------------------------%

% OBS : xi = x(ti) = x(t = ti); i = (0, 1,..., n)
% OBS: x(t) vector must have one more element than v(t) and a(t) in order to calculate v(t) and a(t)
    
x = zeros(n , 1);   % x is the displacement vector x(t)
v = zeros(n , 1);   % v is the velocity vector v(t) = x'(t)
a = zeros(n , 1);   % a is the acceleration vector a(t) = v'(t) = x''(t)

a0 = (F(0) - C * v0 - K * x0) / M;   % a0 = x''(0)
xp = x0 - dt * v0 + (a0 * dt^2) / 2;   % xp = x(t = -1)
x(1 , 1) = x0 ; v(1 , 1) = v0 ; a(1 , 1) = a0;   % initial conditions
x(2 , 1) = (inv((M / (dt ^ 2)) + (C / (2 * dt)))) * ( (((2 * M)/(dt ^ 2)) - K) * x(1 , 1) + ((C / (2 * dt)) - (M / (dt ^ 2))) * xp + F(t(1 , 1)));

for i = 2 : n
    x(i + 1 , 1) = (inv((M / (dt ^ 2)) + (C / (2 * dt)))) * ( (((2 * M)/(dt ^ 2)) - K) * x(i , 1) + ((C / (2 * dt)) - (M / (dt ^ 2))) * x(i - 1 , 1) + F(t(i , 1)) );
    v(i , 1) = (x(i + 1 , 1) - x(i - 1 , 1)) / (2 * dt);
    a(i , 1) = (x(i + 1 , 1) - 2 * x(i , 1) + x(i - 1 , 1)) / (dt ^ 2);
end

%-------------------------------------------------------------------------%
% 2.3) Results for x(t), v(t) and a(t) %
%-------------------------------------------------------------------------%

fprintf('Start of iteration \n\n');
fprintf('Wn = %1.3f rad/s ; Tn = %1.3f s; f = %1.3f Hz\ndt = %1.3f s; dt_crit = %1.3f s; z = %1.3f\n\n' , Wn , T , f , dt , dt_crit, z);

for i = 1 : n
    fprintf('i = %1.0f ; t = %1.3f ; x = %1.3f  ; v = %1.3f  ;  a = %1.3f   \n' , i , t(i , 1) , x(i , 1) , v(i , 1) , a(i , 1));
end

fprintf('\nIteration is over');

%% 3) Ploting %% 

plot(t, x(1 : end - 1, :), 'Color', 'r', 'LineWidth', 2.5)
hold on
plot(t, v, 'Color', 'b', 'LineWidth', 2.5)
hold on
plot(t, a, 'Color', 'g', 'LineWidth', 2.5)
title('Results for x(t), v(t) and a(t)')
xlabel('Time (s)');
ylabel('[x(t)] = m,   [v(t)] = m/s,   [a(t)] = m/s²');
legend('x(t)', 'v(t)', 'a(t)')



