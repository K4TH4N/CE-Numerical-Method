% Set parameters
Lx = 1; % Domain size in x-direction
Ly = 1; % Domain size in y-direction
Nx = 50; % Number of grid points in x-direction
Ny = 50; % Number of grid points in y-direction
dx = Lx/(Nx-1); % Grid spacing in x-direction
dy = Ly/(Ny-1); % Grid spacing in y-direction
k = 1; % Thermal conductivity
rho = 1; % Density
Re = 100;
cp = 1; % Specific heat
U = 1; % Fluid velocity
Te = 1; % Temperature of moving wall
T0 = 0; % Inlet temperature
mass_flow_rate = 1; % Mass flow rate (constant)

% Set initial temperature distribution
T = zeros(Ny,Nx); % Initialize temperature matrix
T(:,1) = T0; % Set inlet temperature
T(:,end) = Te; % Set moving wall temperature

% Set time step size
dt = 0.0001; % Use stability condition for explicit Euler method

% Initialize plot
figure;

xlim([0 Lx]);
ylim([0 Ly]);

% Solve using explicit Euler method
for n = 1:1000 % Number of time steps
    % Compute new temperature distribution
    for j = 2:Ny-1
        for i = 2:Nx-1
            T(j,i) = T(j,i) + dt/(rho*cp) * (k/dx^2 * (T(j,i+1) - 2*T(j,i) + T(j,i-1)) ...
                + k/dy^2 * (T(j+1,i) - 2*T(j,i) + T(j-1,i)) - rho*U/dx * (T(j,i) - T(j,i-1)));
        end
    end
    % Adjust inlet temperature to maintain constant mass flow rate
    T(:,1) = T(:,1) + (mass_flow_rate*dt/(rho*cp*dx))*(T(:,2) - T(:,1));
    % Plot temperature distribution
    contourf(T, 20, 'LineColor', 'none');
    title(['Time step: ', num2str(n*dt)]);
    colorbar;
    drawnow;
    
end