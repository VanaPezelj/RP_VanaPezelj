function brachistochrone_solver
    % Problem Parameters
    g = -9.81; % Gravity (m/s^2)
    x0 = [0; 0]; % Starting point
    xf = [1; -1]; % End point
    v0 = 0; % Initial velocity
    tf_guess = 1; % Initial guess for total time
    
    % Initial guess for shooting variables
    lambda0_guess = [1; -1; 0]; % Guess for initial costates
    theta_guess = pi/4; % Guess for control (angle)

    % Solve using fsolve
    guess = [lambda0_guess; theta_guess; tf_guess];
    sol = fsolve(@(vars) shooting_residual(vars, x0, xf, v0, g), guess)

    % Extract solution
    lambda0 = sol(1:3);
    theta_opt = sol(4);
    tf_opt = sol(5);

    % Display results
    fprintf('Optimal Total Time: %.4f s\n', tf_opt);
    fprintf('Optimal Control (Theta): %.4f radians\n', theta_opt);
    fprintf('Initial Costates: %.4f, %.4f, %.4f\n', lambda0);

    % Simulate optimal trajectory
    tspan = linspace(0, tf_opt, 100);
    [~, X] = ode45(@(t, X) dynamics(t, X, theta_opt, g), tspan, [x0; v0]);

    % Plot results
    figure;
    plot(X(:, 1), X(:, 2), 'b-', 'LineWidth', 2);
    hold on;
    plot(x0(1), x0(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(xf(1), xf(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Brachistochrone Trajectory');
    legend('Optimal Path', 'Start Point', 'End Point');
    grid on;
end

function res = shooting_residual(vars, x0, xf, v0, g)
    lambda0 = vars(1:3);
    theta = vars(4);
    tf = vars(5);

    % Initial conditions
    X0 = [x0; v0];
    % Solve dynamics
    [~, X] = ode45(@(t, X) dynamics(t, X, theta, g), [0, tf], X0);

    % Residuals
    Xf = X(end, :)'; % Final state
    res = [Xf(1:2) - xf; ... % Terminal position error
           H_tf(Xf, lambda0, theta, g) - 1]; % Terminal Hamiltonian condition
end

function dxdt = dynamics(t, X, theta, g)
    x = X(1);
    y = X(2);
    v = X(3);

    dxdt = [v * cos(theta); v * sin(theta); g * sin(theta)];
end

function H = H_tf(X, lambda, theta, g)
    % Terminal Hamiltonian
    v = X(3);
    mu1=-1/(2*v)*(g*lambda(3) + v*(lambda(1) + lambda(2)));
    H = 1 + lambda(1) * v * cos(theta) + lambda(2) * v * sin(theta) + lambda(3) * g * sin(theta)+ mu1*(v*cos(theta) + v*sin(theta));
end

