% Brachistochrone Problem - Indirect Method using Euler-Lagrange equation
clear;
clc;

% Parameters
g = 9.81; % Acceleration (m/s^2)
x1 = 0;   
y1 = 0;   
x2 = 2;   
y2 = -2;  

% Initial guess for the solution
solinit = bvpinit(linspace(x1, x2, 50), @initial_guess);

% using bvp4c
sol = bvp4c(@(x, y) odefun(x, y, g), @bcfun, solinit);

% Extract solution
x_vals = sol.x;
y_vals = sol.y(1, :);

% Plot the result
figure;
plot(x_vals, y_vals, 'b-', 'LineWidth', 2);
hold on;
plot(x1, y1, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Start point
plot(x2, y2, 'go', 'MarkerSize', 8, 'LineWidth', 2); % End point
xlabel('x');
ylabel('y');
title('Brachistochrone Path (Indirect Method)');
grid on;
legend('Brachistochrone Path', 'Start Point', 'End Point');


% -----------------------
% ODE function: Euler-Lagrange equation for the Brachistochrone problem
function dydx = odefun(x, y, g)
    epsilon = 1e-6;  % Small value to avoid division by zero
    % y(1) is the position y, y(2) is the slope dy/dx
    dydx = [y(2);   % dy/dx = slope
            -g * (1 + y(2)^2) / (2 * (abs(y(1)) + epsilon))]; % Euler-Lagrange with regularization
end


% Boundary conditions function
function res = bcfun(ya, yb)
    % Boundary conditions: y(0) = 0 (start), y(x2) = -2 (end)
    res = [ya(1) - 0;   % y(0) = 0
           yb(1) - (-2)]; % y(x2) = -2
end

% Initial guess for the solution
% function yinit = initial_guess(x)
%     R = 1;  % Rough estimate for cycloid radius
%     theta = linspace(0, pi, length(x));  % Parametric angle
%     yinit = [-R * (1 - cos(theta));   % Cycloid y-coordinate
%              -R * sin(theta)];        % Cycloid slope dy/dx
% end

function yinit = initial_guess(x)
    %R = 1;  % Rough estimate for cycloid radius
    %theta = linspace(0, pi, length(x));  % Parametric angle
    yinit = [-x;   % Cycloid y-coordinate
             -1];        % Cycloid slope dy/dx
end