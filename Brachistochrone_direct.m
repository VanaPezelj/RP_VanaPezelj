
% Brachistochrone Problem using MATLAB
% This script solves the Brachistochrone problem using numerical optimization.
% It minimizes the time taken for a particle to slide under gravity between two points.

clear;
clc;

% Problem parameters
x1 = 0;  % Initial x-coordinate
y1 = 0;  % Initial y-coordinate
x2 = 2;  % Final x-coordinate
y2 = -2; % Final y-coordinate
g = 9.81; % Acceleration due to gravity (m/s^2)

% Initial guess for the parameters
N = 50; % Number of points for discretization
theta_init = linspace(pi/2, pi/4, N); % Initial guess for angles (reasonable range)

% Define the optimization problem
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1000);
theta_opt = fmincon(@(theta) objective(theta, x1, y1, x2, y2, g), ...
    theta_init, [], [], [], [], zeros(1, N), pi * ones(1, N), [], options);

% Compute the optimized path
[x_vals, y_vals] = compute_trajectory(theta_opt, x1, y1, g);

% Plot the optimized path
figure;
plot(x_vals, y_vals, 'b-', 'LineWidth', 2);
hold on;
plot(x1, y1, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(x2, y2, 'go', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Optimized Brachistochrone Path');
grid on;
legend('Brachistochrone Path', 'Start Point', 'End Point');

% Objective function: Time of travel for a given path
function T = objective(theta, x1, y1, x2, y2, g)
    [x_vals, y_vals] = compute_trajectory(theta, x1, y1, g);
    dx = diff(x_vals);
    dy = diff(y_vals);
    
    % Prevent issues with undefined or negative y-values
    valid_indices = (y_vals(1:end-1) < 0) & (dx > 0);  % Ensure y is negative and dx is positive
    if any(~valid_indices)
        T = inf; % Penalize invalid solutions
        return;
    end
    
    % Compute time integral for valid parts of the trajectory
    integrand = sqrt((1 + (dy./dx).^2) ./ (2 * g * abs(y_vals(1:end-1)) + 1e-6));
    T = sum(integrand .* dx);
end

% Compute the trajectory based on angles
function [x_vals, y_vals] = compute_trajectory(theta, x1, y1, g)
    % Compute the x and y values based on the parametric angles
    N = length(theta);
    x_vals = zeros(1, N);
    y_vals = zeros(1, N);
    x_vals(1) = x1;
    y_vals(1) = y1 - 1e-3; % Start with a small negative offset to avoid y = 0
    for i = 2:N
        dx = cos(theta(i));
        dy = -sin(theta(i));
        x_vals(i) = x_vals(i-1) + dx;
        y_vals(i) = y_vals(i-1) + dy * sqrt(2 / g * abs(y_vals(i-1)) + 1e-6); % Add small value to avoid sqrt(0)
    end
end
