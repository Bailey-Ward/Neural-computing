clear

part = '6';         % Options are '2', '3', '4', '5', '6' for section 2 questions

dt = 2e-8;          % Time step
tmax=0.35;          % Maximum simulation time
t=0:dt:tmax;        % Time vector

% Set up parameters and initial conditions
ENa = 0.045;       % Reversal for sodium channels (mV)
Ek = -0.082;       % Reversal for potassium channels (mV)
El = -0.060;       % Leak reversal potential (mV)
V0 = -0.065;

Gleak = 30e-9;     % Leak conductance (nS)
GNa = 12e-6;       % Sodium conductance (uS)
Gk = 3.6e-6;       % Potassium conductance (uS)

Cm = 100e-12;       % Membrane capacitance (pF)

istart = 100e-3;    % time applied current starts
ilength= 5e-3;      % length of applied current pulse
Ibase = 0e-9;       % baseline current before/after pulse
Npulses = 1;        % number of current pulses
pulsesep = 20e-3;   % separation between current pulses
V0 = -0.0;          % initial condition for V
m0 = 0.0;           % initial condition for m
h0 = 0.0;           % initial condition for h
n0 = 0.0;           % Default initial condition for n

switch part;        
    case '2'        
        ilength= 100e-3;    % length of applied current pulse
        Ie = 0.22e-9;       % applied current     
    case '3'
        Npulses = 10;       % 10 pulses
        Ie = 0.22e-9;       % applied current for each pulse
        pulsesep = 5e-3;   % time between pulses
    case '4'
        Npulses = 10;       % 10 pulses
        Ibase = 0.6e-9;     % large baseline current 
        Ie = 0e-9;          % pulse current is below baseline
        m0= 0.05;
        h0= 0.5;
        n0=0.35;
        V0 = -0.065;
    case '5'                    
        Ibase = 0.65e-9;    % baseline current
        Ie = 1e-9;          % pulsed current
        m0= 0.05;
        h0= 0.5;
        n0=0.35;
        V0 = -0.065;
    case '6'        
        Ibase = 0.7e-9;    % baseline current
        Ie = 1e-9;          % pulsed current
        V0 = -0.065;
    otherwise
        fprintf('Variable part must be b-f')
end

% Set up the applied current vector
Iapp=Ibase*ones(size(t));       % Initialize current vector at baseline 

for pulse = 1:Npulses;          
    pulsestart = istart + (pulse-1)*pulsesep;   % Onset time of pulse
    pulsestop = pulsestart + ilength;           % Offset time of pulse
    
    % make applied current a new value for duration of current pulse
    for i=round(pulsestart/dt)+1:round(pulsestop/dt);
        Iapp(i) = Ie;
    end
end



V=zeros(size(t));       % voltage vector
V(1) = V0;              % set the inititial value of voltage

n=zeros(size(t));       % n: potassium activation gating variable
n(1) = n0;               % initialize as zero
m=zeros(size(t));       % m: sodium activation gating variable
m(1) = m0;               % initialize as zero
h=zeros(size(t));       % h: sodim inactivation gating variable
h(1) = h0;               % initialize as zero

Itot=zeros(size(t));    % to plot and look at the total current
I_Na=zeros(size(t));    % to plot and look at sodium current
I_K=zeros(size(t));     % to plot and look at potassium current
I_L=zeros(size(t));     % to plot and look at leak current

for i = 2:length(t); 
    
    Vm = V(i-1);          % membrane potential
    
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;     
    else
        alpha_m = (10^5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4 * 10^3 *exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 10^3/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      
    else;                   % potassium activation rate
        alpha_n = (10^4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.070)/0.08);     % potassium deactivation rate
      
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    
    m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    
    h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    
    n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    
    I_Na(i) = GNa*m(i)*m(i)*m(i)*h(i)*(ENa-V(i-1)); % total sodium current
    
    I_K(i) = Gk*n(i)*n(i)*n(i)*n(i)*(Ek-V(i-1)); % total potassium current
    
    I_L(i) = Gleak*(El-V(i-1));    % Leak current is straightforward
    
    Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i); % total current is sum of leak + active channels + applied current
    
    V(i) = V(i-1) + Itot(i)*dt/Cm;        % Update the membrane potential, V.
    
end

% Plot results
figure;

subplot(3,1,1);
plot(t, Iapp * 1e3, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential');
grid on;

subplot(3,1,2);
plot(t, V, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential Dynamics');
grid on;

subplot(3,1,3);
plot(t, m, 'b', 'LineWidth', 1.5); hold on;
plot(t, h, 'r', 'LineWidth', 1.5);
plot(t, n, 'g', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Gating Variables');
legend('m', 'h', 'n');
title('Gating Variable Dynamics');
grid on;