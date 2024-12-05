% Parameters
EL = -75e-3;             % Resting potential (mV)
V_th = -50e-3;           % Threshold potential (mV)
V_reset = -80e-3;        % Reset potential (mV)
delta_th = 2e-3;         % (mV)
Gl = 10e-9;              % (ns)
Cm = 100e-12;            % Membrane capacitance (pF)
a = 2e-9;                % (ns)
b= 0.02e-9;              % (nA)
tsra = 200e-3;           % Adaptation time constant (s)

% Simulation parameters
dt = 0.0001;
t = 0:dt:1.5;

% Input current
Iapp = zeros(size(t));
Iapp(t >= 0.5 & t <= 1.0) = 500e-12;  % 500 pA pulse from 0.5 to 1.0 s

% Initial conditions
V = EL * ones(size(t)); % Initial membrane potential
I_SRA = zeros(size(t)); % Initial adaptation conductance

% Simulation loop (Euler's method)
for i = 1:length(t)-1
    % Update dV/dt
    dVmdt = Gl * ((EL - V(i)) + delta_th * exp((V(i)-V_th) / delta_th)) - I_SRA(i) + Iapp(i);
    dVmdt = dVmdt*(1/Cm);
    V(i+1) = V(i) + dVmdt * dt;
    
    % Update dG_SRA/dt
    dISRA_dt = a * (V(i)- EL) - I_SRA;
    I_SRA(i+1) = I_SRA(i) + dISRA_dt(i) * dt;
    
    % Check for spike
    if V(i+1) > V_th
        V(i+1) = V_reset;           % Reset membrane potential
        I_SRA(i+1) = I_SRA(i+1) + b; % Increment adaptation conductance
    end
end

subplot(3,1,1);
plot(t, V * 1e3, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential');
grid on;

subplot(3,1,2);
plot(t, I_SRA * 1e9, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('G_{SRA} (nS)');
title('Adaptation Conductance');
grid on;

sgtitle('AELIF Model with Adaptation');
