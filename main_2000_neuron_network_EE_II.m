
function  [spike_E,RE,spike_time_E] = main_2000_neuron_network_EE_II(duration,step,duration_step,S_key_IE,S_key_EI,S_key_EE,S_key_II,J_E,J_EE,J_I,J_II,cross_100,Ue,Ui,stimon,stopEIWeight,sample_duration,stimon_param,plaston_param,W_IE0,W_EE0,W_EI0,W_II0) 

% =[time,W_IE,spike_E,spike_I,Synchrony,time_syn,spike_time_E,RE,W_IEN,W_IEC]
% -- Neuron Parameters --
    mew_e = 20.8; 
    sigma_e = 1; 
    mew_i = 18; 
    sigma_i = 3;
    %mew_c = 0;
    N_E = 1600;
    N_I = 400;
    Ntot = N_E + N_I;
    %Energy = zeros(duration_step,1);
    
% -- Synaptic Parameters --
    %Weight_0 = 1;
%     J_E = 100;% NOTE: J_E = J_EI
%     J_EE = 50; 
%     J_I = 260; % NOTE: J_I = J_IE
%     J_II = 50;
    %N_i = 1;
    C_E = 0.3*Ntot; 
    C_EE = Ntot;
    C_I = 0.3*Ntot;
    C_II = Ntot;
    tau_LTP = 20;
    tau_LTD = 22;   
    

    

    

    num_synapses_IE = max(max(S_key_IE));
%     num_synapses_EI = max(max(S_key_EI));
%     num_synapses_EE = max(max(S_key_EE));
%     num_synapses_II = max(max(S_key_II));
    
    W_IE = zeros(duration_step,1); 
    %W_IE_std = zeros(duration_step,1); 
    
% -- Initial Conditions --

    vE0 = 14*ones(1,N_E);
    vI0 = 14*ones(1,N_I);
    S_EI0 = zeros(1,N_E);
    S_IE0 = zeros(1,N_I);
    S_EE0 = zeros(1,N_E);
    S_II0 = zeros(1,N_I);
    X_IE0 = zeros(1,N_I);
    X_EI0 = zeros(1,N_E);
    X_II0 = zeros(1,N_I);
    X_EE0 = zeros(1,N_E);
    spt_E0 = 0;
    
    Apost0 = zeros(1,num_synapses_IE);
    Apre0 = zeros(1,num_synapses_IE);
    
    leftover_S_EI = zeros(5/step + 1 ,N_E);
    leftover_S_IE = zeros(5/step + 1 ,N_I);
    leftover_S_EE = zeros(5/step + 1 ,N_E);
    leftover_S_II = zeros(5/step + 1 ,N_I);
    ref_E = zeros(1,N_E);
    ref_I = zeros(1,N_I);    
    phif = zeros(1,N_E);
    tau_E_m = 10;
    tau_I_m = 10;
    




% -- Run --    
    
    %sample_duration = 20;%max(t_pulseE,t_pulseI);
    Synchrony = zeros(floor(duration/sample_duration),1);
    time_syn = zeros(floor(duration/sample_duration),1);

for i = 1:floor(duration/sample_duration)
    % run parameters
    comp_time = (i-1)*sample_duration;
%     if mean(W_IE0)*J_I < stopEIWeight
%         cross_100 = 0;
%     end
    ON = stimon_param*(i*sample_duration >= stimon)*cross_100; % turn on control input
    plast_on = plaston_param*(i*sample_duration >= 100); % past on/off
    
    % indexes
    a = 1 + (i>=2)*(i-1)*sample_duration/step;
    b = i*sample_duration/step;
    
    % desynchronizing input
    if stimon_param == 1
        ue = Ue(1,a:b);
        ui = Ui(1,a:b);
    else
        ue = zeros(1,b-a+1);
        ui = zeros(1,b-a+1);
    end
    percent_V_stim = 1;
    
   [timem, v_Em, v_Im, S_EIm, S_IEm, S_EEm, S_IIm, X_EIm, X_IEm, X_EEm, X_IIm, Apostm, Aprem, W_IEm, spike_Em, spike_Im, ref_Em, ref_Im, synchronym, spt_Em, phif] = ode_neuron_model(plast_on,ON,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_EE0,W_II0,mew_e,sigma_e,ue,ui,mew_i,sigma_i,J_E,J_EE,J_I,J_II,C_E,C_I,C_EE,C_II,tau_LTP,tau_LTD,step,sample_duration,N_E,N_I,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,ref_E,ref_I,tau_E_m,tau_I_m, percent_V_stim,comp_time,spt_E0,phif);
    
    
    
    % recorded variables
    time(a:b,:) = timem(1:sample_duration/step,:) + (i-1)*sample_duration;
    W_IE(a:b,1) = mean(W_IEm(1:sample_duration/step,:),2);
    spike_E(a:b,:) = spike_Em(1:sample_duration/step,:);
    spike_I(a:b,:) = spike_Im(1:sample_duration/step,:);
    Synchrony(i) = synchronym;
    time_syn(i) = sample_duration*(i);
    spike_time_E(a:b,:) = spt_Em(1:sample_duration/step,:);
    
    % generate intial condition for next run
    sample_end = sample_duration/step ;
    vE0 = v_Em(sample_end,:);
    vI0 = v_Im(sample_end,:);
    S_EI0 = S_EIm(sample_end,:);
    S_IE0 = S_IEm(sample_end,:);
    S_II0 = S_IIm(sample_end,:);
    S_EE0 = S_EEm(sample_end,:);
    X_EI0 = X_EIm(sample_end,:);
    X_IE0 = X_IEm(sample_end,:);
    X_II0 = X_IIm(sample_end,:);
    X_EE0 = X_EEm(sample_end,:);
    Apost0 = Apostm(sample_end,:);
    Apre0 = Aprem(sample_end,:);
    W_IE0 = W_IEm(sample_end,:);
    %W_EI0 = Weight_0; 
    left_sample_end = sample_end - 5/step;
    leftover_S_EI = S_EIm(left_sample_end:sample_end,:);
    leftover_S_IE = S_IEm(left_sample_end:sample_end,:);
    leftover_S_EE = S_EEm(left_sample_end:sample_end,:);
    leftover_S_II = S_IIm(left_sample_end:sample_end,:);
    spt_E0 = spt_Em(sample_end,:);
    
%     if i == floor((5*duration/sample_duration)/100)
%         disp('5% Complete')
%     end    
%     if i == floor((10*duration/sample_duration)/100)
%         disp('10% Complete')
%     end    
%     if i == floor((15*duration/sample_duration)/100)
%         disp('15% Complete')
%     end    
%     if i == floor((20*duration/sample_duration)/100)
%         disp('20% Complete')
%     end    
%     if i == floor((duration/sample_duration)/4)
%         disp('25% Complete')
%     end    
%     if i == floor(30*(duration/sample_duration)/100)
%         disp('30% Complete')
%     end    
%     if i == floor(35*(duration/sample_duration)/100)
%         disp('35% Complete')
%     end    
%     if i == floor(4*(duration/sample_duration)/10)
%         disp('40% Complete')
%     end    
%     if i == floor(45*(duration/sample_duration)/100)
%         disp('45% Complete')
%     end    
%     if i == floor((duration/sample_duration)/2)
%         disp('50% Complete')
%     end    
%     if i == floor(55*(duration/sample_duration)/100)
%         disp('55% Complete')
%     end    
%     if i == floor(6*(duration/sample_duration)/10)
%         disp('60% Complete')
%     end    
%     if i == floor(65*(duration/sample_duration)/100)
%         disp('65% Complete')
%     end    
%     if i == floor(70*(duration/sample_duration)/100)
%         disp('70% Complete')
%     end    
%     if i == floor(75*(duration/sample_duration)/100)
%         disp('75% Complete')
%     end    
%     if i == floor(8*(duration/sample_duration)/10)
%         disp('80% Complete')
%     end    
%     if i == floor(8.5*(duration/sample_duration)/10)
%         disp('85% Complete')
%     end
%     
%     if i == floor(9*(duration/sample_duration)/10)
%         disp('90% Complete')
%     end    
%     if i == floor(9.5*(duration/sample_duration)/10)
%         disp('95% Complete')
%     end

end

     
% calculated the Kuramoto Order Parameter
step = 0.1;
tf = floor(floor(duration/sample_duration)*sample_duration/step);
t = (0.1:step:tf*step);
[RE] = kuramoto_syn(spike_time_E,t,step,tf*step,N_E);
% 
% %Get E-to-I synaptic weights
% W_IEN = W_IEm(sample_duration/step,:);
% W_IEP = W_IEm(1,:);
% W_IEC = sum(abs(W_IEN-W_IEP));
end
