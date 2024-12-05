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
t = 0:dt:5;

%Input current
Iapp = 240:16:550;
Iapp = Iapp * 1e-12;

%Initial conditions
V = zeros(1, length(t));
V(1) = El;
G_sra = zeros(1, length(t));
initialrate = zeros(size(Iapp)); % array to store 1/(first ISI)
finalrate = zeros(size(Iapp));   % array to store 1/(final ISI)
singlespike = zeros(size(Iapp)); % array to store "1" for only 1 spike

%Loop through applied current values
for j = 1:length(Iapp)   
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = El;
    
    G_sra = zeros(1, length(t));
    
    spikes = zeros(size(t));
    
    %Loop through time vector
    for i = 1:length(t)-1
        if ( V(i) > V_th )         
            V(i) = V_reset;         
            G_sra(i) = G_sra(i) + delta_G;   
            spikes(i) = 1;
        end
    
        V(i+1) = V(i) + dt * ( (El-V(i))/Rm + G_sra(i)*(Ek - V(i)) + Iapp(j))/Cm;

        G_sra(i+1) = G_sra(i) - dt*(G_sra(i)/t_sra);
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

xlabel('Applied Current (nA)');
ylabel('Firing rate (Hz)');
legend('show');
title('f-I Curve for LIF Model with Adaptation');
grid on;

hold off;
figure(1);



