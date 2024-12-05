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
dt = 1e-4;
t = 0:dt:5;

%Input current
Iapp = 240:16:550;
Iapp = Iapp * 1e-12;

% Initial conditions
V = zeros(1, length(t));
V(1) = El;
I_SRA = zeros(size(t));
initialrate = zeros(size(Iapp)); % array to store 1/(first ISI)
finalrate = zeros(size(Iapp));   % array to store 1/(final ISI)
singlespike = zeros(size(Iapp)); % array to store "1" for only 1 spike

%Loop through applied current values
for j = 1:length(Iapp)   
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = El;
    
    I_SRA = zeros(1, length(t));
    
    spikes = zeros(size(t));
    
    %Loop through time vector
    for i = 1:length(t)-1
        if ( V(i) > V_th )         
            V(i) = V_reset;         
            I_SRA(i) = I_SRA(i) + b;   
            %neuron_fires(j) = neuron_fires(j) + 1;
            spikes(i) = 1;
        end
    
        V(i+1) = V(i)+ dt* (Gl * (El - V(i) + delta_th * exp((V(i)-V_th) / delta_th)) - I_SRA(i) + Iapp(j)) / Cm;

        dISRA_dt = a * (V(i)- El) - I_SRA/tsra;
        I_SRA(i+1) = I_SRA(i) + dt * dISRA_dt(i);
    end

    
     spiketimes = dt*find(spikes);           % extract the spike times
    
    if ( length(spiketimes) > 1 )           % if there is more than 1 spike
        ISIs = diff(spiketimes);            % ISI = interval between spikes
        initialrate(j) = 1/ISIs(1);     % inverse of first ISI
        if ( length(ISIs) > 1 )             % if there are further ISIs
            finalrate(j) = 1/ISIs(end); % inverse of final ISI
        end
        
    else
        if ( length(spiketimes) == 1 )      % if there is only one spike
            singlespike(j) = 1;         % record "1" for this trial
        end
    end
end

hold on;

plot(Iapp*1e12, finalrate,'-o', 'DisplayName', 'Initial f');

ISIindices = find(initialrate);
plot(1e12*Iapp(ISIindices),initialrate(ISIindices),'x', 'DisplayName', 'Steady-state f');

ISIindices = find(singlespike);
plot(1e12*Iapp(ISIindices),0*singlespike(ISIindices),'x', 'DisplayName', 'Initial f');

xlabel('Applied Current (nA)');
ylabel('Firing rate (Hz)');
legend('show');
title('f-I Curve for AELIF Model with Adaptation');
grid on;

hold off;
figure(1);

