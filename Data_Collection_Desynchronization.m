clear -all
clearvars
clc

%Desynchronization Data
    
 %%% Parameters to Vary
    AmpV = [50]%[10,30,50,70,100]'; %10 to 100 %10, 30, 50, 70, 100
    PW1V = [0.1]%,0.5,1,1.5,2]'; %ms between 0.1ms to 2 ms  %0.1, 0.5, 1, 1.5, 2
    PW2V = [0.1]%,0.5,1,1.5,2]'; %ms between 0.1ms to 2 ms  %0.1, 0.5, 1, 1.5, 2
    PW3V = [10]%,15,20,25]'; %ms between 10 to 25 (100 Hz to 40 Hz) %10,15,20,25
    NsV = [5]%,15,25,40]'; %5 to 40 pulses 5, 15, 25, 40
    J_IV = [450]%,400,350,300,250,200,165]'; % NOTE: J_I = J_IE 450 - 165 with interval of 50
    
    aE = -1;
    aI = 1;

% -- Load Synaptic Data --

   %Load Synaptic Connection Matrix 
    load('synaptic_connection.mat');
    
    %Load Initial Synaptic Weight Matrix
    load('initial_synaptic_weights.mat');
    
    %load('W_IE0.mat')
    
    %Synaptic Strength Parameters --
    J_E = 100;% NOTE: J_E = J_EI
    J_EE = 50; 
    
    J_II = 50;
    
% -- Stimulation Input --
    step = 0.1; %ms
    sample_duration = 20; %Sample calculation (ms)
    stimon = 2000; %time at which stimulation is turned on (ms)
    
    
    cross_100 = 1; 
    stopEIWeight = 75;
   
    stimon_param = 1; %0: Stimulation OFF, 1: Stimulation ON
    plaston_param = 0; % 0 Plasticity OFF, 1: Plasticity ON
    counter = 1;
    
    %for loop for running simulation
    
    for ii = 1:size(J_IV,1)
        J_I = J_IV(ii,1);
        for jj = 1:size(NsV,1)
            Ns = NsV(jj,1);
            for kk = 1:size(PW3V,1)
                PW3 = PW3V(kk,1);
                for ll = 1:size(PW2V,1)
                    PW2 = PW2V(ll,1);
                    for mm = 1:size(PW1V,1)
                        PW1 = PW1V(mm,1);
                        for nn = 1:size(AmpV,1)
                            Amp = AmpV(nn,1);
    
                            [Ue,Ui,duration] = stim_pattern_generator(step,sample_duration,stimon,Amp,PW1,PW2,PW3,Ns);


                            duration_step = duration/step;




                            %Simulation 
    
                            [spike_E,RE,spike_time_E] = main_2000_neuron_network_EE_II(duration,step,duration_step,S_key_IE,S_key_EI,S_key_EE,S_key_II,J_E,J_EE,J_I,J_II,cross_100,Ue,Ui,stimon,stopEIWeight,sample_duration,stimon_param,plaston_param,W_IE0,W_EE0,W_EI0,W_II0); 
    
    
                            %Data Recorded - Firing Rates, Order Parameter, Spike timing   

                            %Calculating mean order parameter for 10ms before and 10ms after
	
                            idx = find(Ue ~= 0);
                            id1 = idx(1)-500:idx(1)-1;
                            nstep = PW3/step;
                            id2 = idx(end)+nstep+1:idx(end)+nstep+500;
                            Re1 = mean(RE(id1,1));
                            Re2 = mean(RE(id2,1));

                            %Calculating firing rate 

                            %Before Stimulation
                            spike_E1 = spike_E(id1,:);
                            for i = 1:size(id1,2)
                                for j = 1:1600
                                    if spike_E1(i,j) >0
                                        spike_E1(i,j) = 1;
                                    end
                                end
                            end

                        for j = 1:1600
                            lamE1(j,1) = sum(spike_E1(:,j))*1000/50; %In Hz
                        end

                        %After Stimulation
                        spike_E2 = spike_E(id2,:);
                        for i = 1:size(id2,2)
                            for j = 1:1600
                                if spike_E2(i,j) >0
                                    spike_E2(i,j) = 1;
                                end
                            end
                        end

                        for j = 1:1600
                            lamE2(j,1) = sum(spike_E2(:,j))*1000/50; %In Hz
                        end

                        %Spike timting Before and after Stimulation
                        spike_ET1 = spike_time_E(id1,:); %Before 
                        spike_ET2 = spike_time_E(id2,:); %After


                        %storing data
                        Input_data(:,counter) = [aE;aI;Amp;PW1;PW2;PW3;Ns;Re1;lamE1];
                        Label_data(:,counter) = [Re2;lamE2];
                        spike_time_input(:,:,counter) = spike_ET1;
                        spike_time_label(:,:,counter) = spike_ET2;
                        counter = counter+1;
                        clearvars spike_E RE spike_time_E idx id1 nstep id2 Re1 Re2 spike_E1 spike_E2 lamE1 lamE2
                        end
                    end
                end
            end
        end
       save(['DesyncneuraldataTrial_' num2str(ii) '.mat'],'Input_data','Label_data','spike_time_input','spike_time_label')
       clearvars Input_data Label_data spike_time_input spike_time_label
       counter = 1; 
    end


% % Plots 
% 
% tsize = size(time,1);
% figure(1)
% plot(time/1000,J_I*W_IE(1:tsize,1),'k')
% xlabel('Time (sec)')
% ylabel('Average Synpatic Weight')
% xlim([0 duration/1000])
% ylim([0 300])
% 
% figure(2)
% plot(time/1000,RE(1:tsize,1))
% ylim([0 1.1])
% 
% figure(3)
% plot(time/1000, spike_E, 'k.',time/1000, spike_I, 'r.')
% xlim([stimon/1000 time(end,1)/1000])
% ylim([0.9 2000.1])
% ylabel('Neuron Index')
% xlabel('Time (sec)')
