% Parameters
EL = -75e-3;            % Resting potential (mV)
V_th = -50e-3;           % Threshold potential (mV)
V_reset = -80e-3;        % Reset potential (mV)
Rm = 100e6;             % Membrane resistance (Ohms)
Cm = 100e-12;           % Membrane capacitance (pF)
Ek = -80e-3;            % Adaptation reversal potential (mV)
DeltaG_SRA = 1e-9;      % Increment of adaptation conductance (nS)
tsra = 200e-3;          % Adaptation time constant (ms)

% Simulation parameters
dt = 1e-4;
t = 0:dt:1.5;

% Input current
Iapp = zeros(size(t));
Iapp(t >= 0.5 & t <= 1.0) = 500e-12;  % 500 pA pulse from 0.5 to 1.0 s

% Initial conditions
V = ones(size(t)); % Initial membrane potential
V(1) = EL;
G_SRA = zeros(size(t)); % Initial adaptation conductance

% Simulation loop (Euler's method)
for i = 1:length(t)-1
    % Update dV/dt
    dVdt = (EL - V(i)) / Rm + G_SRA(i) * (Ek - V(i)) + Iapp(i);
    dVdt = dVdt*(1/Cm);
    V(i+1) = V(i) + dVdt * dt;
    
    % Update dG_SRA/dt
    dGSRA_dt = -G_SRA(i) / tsra;
    G_SRA(i+1) = G_SRA(i) + dGSRA_dt * dt;
    
    % Check for spike
    if V(i+1) > V_th
        V(i+1) = V_reset;           % Reset membrane potential
        G_SRA(i+1) = G_SRA(i+1) + DeltaG_SRA; % Increment adaptation conductance
    end
end

% Plotting results
figure;
subplot(3,1,1);
plot(t, Iapp , 'LineWidth', 1.5);
axis([0 1.5 0 700e-12 ])
xlabel('Time (s)');
ylabel('I_{app} (pA)');
title('Input Current');
grid on;

subplot(3,1,2);
plot(t, V , 'LineWidth', 1.5);
axis([0 1.5 -0.085 -0.04])
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential');
grid on;

subplot(3,1,3);
plot(t, G_SRA , 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('G_{SRA} (nS)');
title('Adaptation Conductance');
grid on;

sgtitle('LIF Model with Adaptation');
