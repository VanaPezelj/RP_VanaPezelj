function Brachistochrone_BVP
    % Problem Parameters
    g = -9.81; % Gravity (m/s^2)
    x0 = 0; y0 = 0; % Starting position
    xf = 10; yf = -10; % Final position
    v0 = 0; % Initial velocity
    
    % Define initial guess for solution
    % p=6;
    % max_residual=100000;
    %while max_residual>6
        %p=p+1
        solinit = bvpinit(linspace(0, 1, 100), @initial_guess);

        % Set the tolerances
        options = bvpset('RelTol', 1e-5, 'AbsTol', 1e-5);

        % Solve the BVP 
        sol = bvp4c(@(x, y) odefun(x, y, g), @(ya, yb) bcfun(ya, yb, x0, y0, xf, yf, v0), solinit,options);

        
        % Evaluate the residuals at the solution points
        residuals = odefun(sol.x, sol.y, g);
        
        % Compute the maximum residual
        max_residual = sum(abs(residuals(:)));
        disp(['Maximum Residual: ', num2str(max_residual)]);    
    %end

    residuals_bc = bcfun(sol.y(:, 1), sol.y(:, end), x0, y0, xf, yf, v0);
    disp(['Boundary Condition Residuals: ', num2str(max(abs(residuals_bc)))]);

    % Extract the solution
    x_vals = sol.y(1, :); % x-position
    y_vals = sol.y(2, :); % y-position

    plot(sol.x, sol.y(3, :), 'r-', 'LineWidth', 2);
    xlabel('Time'); ylabel('Velocity');

    % Plot the solution
    figure;
    plot(x_vals, y_vals, 'b-', 'LineWidth', 2);
    hold on;
    plot(x0, y0, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(xf, yf, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Brachistochrone Trajectory');
    legend('Optimal Path', 'Start Point', 'End Point');
    grid on;
    
    % Compute total travel time
    T = trapz(sol.x, 1 ./ sol.y(3, :)); % T = integral of 1/v over time
    disp(['Total Travel Time: ', num2str(T)]);



    % R = 0.573; % Cycloid radius
    % theta = linspace(0, pi, 100);
    % x_theoretical = R * (theta - sin(theta));
    % y_theoretical = -R * (1 - cos(theta));
    % plot(x_theoretical, y_theoretical, 'k--', 'LineWidth', 2); % Overlay theoretical trajectory
    
    % error_x = max(abs(sol.y(1, :) - x_theoretical));
    % error_y = max(abs(sol.y(2, :) - y_theoretical));
    % disp(['Max Error in x: ', num2str(error_x)]);
    % disp(['Max Error in y: ', num2str(error_y)]);
end



% ------------------------------------------------------------
% Define the ODE System
function dydx = odefun(x, y, g)
    % y(1) = x, y(2) = y, y(3) = v (velocity)
    % y(4) = lambda_1, y(5) = lambda_2, y(6) = lambda_3 (costates)
    theta = atan2(y(5), y(4));            % Control angle from costates
    dydx = [y(3) * cos(theta);            % dx/dt
            y(3) * sin(theta);            % dy/dt
            g * sin(theta);               % dv/dt
            0;                            % d(lambda_1)/dt
            -y(6) * g;                    % d(lambda_2)/dt
            -y(4) * cos(theta) - y(5) * sin(theta)]; % d(lambda_3)/dt
end

% ------------------------------------------------------------
% Define Boundary Conditions
function res = bcfun(ya, yb, x0, y0, xf, yf, v0)
    % Initial conditions at t = 0
    res1 = [ya(1) - x0;  % x(0) = x0
            ya(2) - y0;  % y(0) = y0
            ya(3) - v0]; % v(0) = v0

    % Final conditions at t = tf
    res2 = [yb(1) - xf;  % x(tf) = xf
            yb(2) - yf;  % y(tf) = yf
            yb(6) - 0];      % lambda_3(tf) = 0 (example transversality condition)

    % Combine results
    res = [res1; res2];
end


function H = H_final(yb)
    % Hamiltonian at the final time (t_f)
    % Example calculation for H (modify based on your problem):
    x = yb(1); % x(t_f)
    y = yb(2); % y(t_f)
    v = yb(3); % v(t_f)
    lambda1 = yb(4); % lambda_1(t_f)
    lambda2 = yb(5); % lambda_2(t_f)
    lambda3 = yb(6); % lambda_3(t_f)
    theta = atan2(lambda2, lambda1); % Optimal control angle
    g = -9.81; % Gravity
    H = 1 + lambda1 * v * cos(theta) + lambda2 * v * sin(theta) + lambda3 * g * sin(theta);
end
% ------------------------------------------------------------
% Define Initial Guess for the Solution
function yinit = initial_guess1(x)
    % Simple linear guess for states and costates
    yinit = [x;                % x: Linear from 0 to 1
             -x;               % y: Linear from 0 to -1
             -sqrt(x);         % v: Constant initial guess
             0;                % lambda_1: Constant guess
             0;                % lambda_2: Constant guess
             1];               % lambda_3: Constant guess
end

function yinit = initial_guess(x)
    R = 0.573; % Cycloid radius approximation
    theta = pi * x; % Parametric angle
    yinit = [R * (theta - sin(theta));   % x: Cycloid x-coordinate
             -R * (1 - cos(theta));      % y: Cycloid y-coordinate
             -sqrt(x);                     % v: Quadratic growth
             1;                          % lambda_1: Constant guess
             -1;                         % lambda_2: Constant guess
             0];                         % lambda_3: Constant guess
end