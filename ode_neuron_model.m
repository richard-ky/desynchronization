function [time, v_E, v_I, S_EI, S_IE, S_EE, S_II, X_EI, X_IE, X_EE, X_II, Apost, Apre, W_IE, spike_E, spike_I, ref_E, ref_I, synchrony, spt_E, phif] = ode_neuron_model(plast_on,ON1,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_EE0,W_II0,mew_e,sigma_e,ue,ui,mew_i,sigma_i,J_E,J_EE,J_I,J_II,C_E,C_I,C_EE,C_II,tau_LTP,tau_LTD,step,sample_duration,N_E,N_I,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,ref_E,ref_I,tau_E_m,tau_I_m, percent_V_stim,comp_time,spt_E0,phif)
% -- run parameters --
    duration_step = sample_duration/step;
 
% -- neuron parameters --
    tau_d = 1;
    tau_r = 1;
    vrest = 0; 
    whitenoise_E = randn(duration_step+1,N_E); % gaussian white noise
    whitenoise_I = randn(duration_step+1,N_I); % gaussian white noise
    vreset = 14;
    refractory = 2; 
     
  
     
% -- Synaptic Parameters --    
    WEI = W_EI0;
    WEE = W_EE0;
    WII = W_II0;
    syn_delay = 5;
    num_synapses_IE = max(max(S_key_IE));
    num_synapses_EI = max(max(S_key_EI));
     
% -- Plasticity Parameters --
    dApre_0 = 0.005*1;
    dApost_0 = dApre_0*1;
    Apre_i = 0;
    Apost_i = 0;
    a_LTD = -1*plast_on*1.1;
    a_LTP = 1*plast_on*1;
    eta = 0.25; 
       
% -- intialize vectors --
    % ode vectors
    dv_Edt = zeros(duration_step+1,N_E);
    dv_Idt = zeros(duration_step+1,N_I);
    dS_EIdt = zeros(duration_step+1,N_E);
    dS_IEdt = zeros(duration_step+1,N_I);
    dX_EIdt = zeros(duration_step+1,N_E);
    dX_IEdt = zeros(duration_step+1,N_I);
    dApostdt = zeros(duration_step+1,num_synapses_IE);
    dApredt = zeros(duration_step+1,num_synapses_IE);
     
    % state vectors    
    v_E = zeros(duration_step+1,N_E);
    v_I = zeros(duration_step+1,N_I);     
    S_EI = zeros(duration_step+1,N_E);
    S_EE = zeros(duration_step+1,N_E);
    S_II = zeros(duration_step+1,N_I);
    S_IE = zeros(duration_step+1,N_I);    
    X_EI = zeros(duration_step+1,N_E);
    X_EE = zeros(duration_step+1,N_E);
    X_II = zeros(duration_step+1,N_I);
    X_IE = zeros(duration_step+1,N_I);    
    Apost = zeros(duration_step+1,num_synapses_IE);
    Apre = zeros(duration_step+1,num_synapses_IE);    
    W_IE = zeros(duration_step+1,num_synapses_IE);    
    time = zeros(duration_step+1,1);
    spike_E = zeros(duration_step+1,N_E);
    spike_I = zeros(duration_step+1,N_I);
    spike_E_time = zeros(duration_step+1,N_E);
    spike_I_time = zeros(duration_step+1,N_I);
    spt_E = zeros(1,N_E);
    
     
    v_E(1,:) = vE0;
    v_I(1,:) = vI0;
     
    S_EI(1,:) = S_EI0;
    S_IE(1,:) = S_IE0;
    S_EE(1,:) = S_EE0;
    S_II(1,:) = S_II0;
     
    X_EI(1,:) = X_EI0;
    X_IE(1,:) = X_IE0;
    X_II(1,:) = X_II0;
    X_EE(1,:) = X_EE0;
     
    Apost(1,:) = Apost0;
    Apre(1,:) = Apre0;
     
    W_IE(1,:) = W_IE0;
    spt_E(1,:) = spt_E0;

% --- Stimulation Parameters ---    
    
    
    % number of E neurons stimulating
    stim_percent_E = zeros(1,N_E);
    a = 1;
    b = floor(percent_V_stim*N_E);
    stim_percent_E(1,a:b) = 1 + 0*0.1*rand(1,b);
    
    % number of I neurons stimulating
    stim_percent_I = zeros(1,N_I);
    a = 1;
    b = floor(percent_V_stim*N_I);
    stim_percent_I(1,a:b) = 1 + 0*0.1*rand(1,b);
    
    
for i = 1:duration_step
% -- time update --
    time(i+1,1) = time(i,1) + step;

% -- excitatory --
    % dxdt update
    if i > syn_delay/step
        dv_Edt(i,:) = (vrest - v_E(i,:) + J_E/C_E*S_EI(i-syn_delay/step,:) + J_EE/C_EE*S_EE(i-syn_delay/step,:) + mew_e + ON1.*ue(1,i) + sigma_e.*(tau_E_m(1,:).^0.5).*whitenoise_E(i,:))./tau_E_m(1,:);
    else
        dv_Edt(i,:) = (vrest - v_E(i,:) + J_E/C_E*leftover_S_EI(i,:) + J_EE/C_EE*leftover_S_EE(i,:) + mew_e + ON1.*ue(1,i) + sigma_e.*(tau_E_m(1,:).^0.5).*whitenoise_E(i,:))./tau_E_m(1,:);
    end
    dS_EIdt(i,:) = (-S_EI(i,:) + X_EI(i,:))/tau_d;
    dX_EIdt(i,:) = -X_EI(i,:)/tau_r;
    dS_EEdt(i,:) = (-S_EE(i,:) + X_EE(i,:))/tau_d;
    dX_EEdt(i,:) = -X_EE(i,:)/tau_r;
    
    % x update
    v_E(i+1,:) = v_E(i,:) + step*dv_Edt(i,:);
    S_EI(i+1,:) = S_EI(i,:) + step*dS_EIdt(i,:);
    X_EI(i+1,:) = X_EI(i,:) + step*dX_EIdt(i,:);
    S_EE(i+1,:) = S_EE(i,:) + step*dS_EEdt(i,:);
    X_EE(i+1,:) = X_EE(i,:) + step*dX_EEdt(i,:);
    
    % calculate kuramoto order
    spt_E(i+1,:) = spt_E(i,:);
    
    % update refractory
    ref_E(1,:) = ref_E(1,:) - step;
         
     
% -- inhibitory --
   
    % dxdt update
    if i > syn_delay/step
        dv_Idt(i,:) = (vrest - v_I(i,:) + J_I/C_I*S_IE(i-syn_delay/step,:) + J_II/C_II*S_II(i-syn_delay/step,:) + mew_i + ON1.*ui(1,i) + sigma_i.*(tau_I_m(1,:).^0.5).*whitenoise_I(i,:))./tau_I_m(1,:);
    else
        dv_Idt(i,:) = (vrest - v_I(i,:) + J_I/C_I*leftover_S_IE(i,:) + J_II/C_II*leftover_S_II(i,:) + mew_i + ON1.*ui(1,i) + sigma_i.*(tau_I_m(1,:).^0.5).*whitenoise_I(i,:) )./tau_I_m(1,:);
    end
    dS_IEdt(i,:) = (-S_IE(i,:) + X_IE(i,:))/tau_d;
    dX_IEdt(i,:) = -X_IE(i,:)/tau_r;
    dS_IIdt(i,:) = (-S_II(i,:) + X_II(i,:))/tau_d;
    dX_IIdt(i,:) = -X_II(i,:)/tau_r;
    
    % x update
    v_I(i+1,:) = v_I(i,:) + step*dv_Idt(i,:);
    S_IE(i+1,:) = S_IE(i,:) + step*dS_IEdt(i,:);
    X_IE(i+1,:) = X_IE(i,:) + step*dX_IEdt(i,:);
    S_II(i+1,:) = S_II(i,:) + step*dS_IIdt(i,:);
    X_II(i+1,:) = X_II(i,:) + step*dX_IIdt(i,:);
    
    % update refractory
    ref_I(1,:) = ref_I(1,:) - step;
     
% --- update synaptic variables ---
    % dxdt update
        dApostdt(i,:) = -Apost(i,:)/tau_LTD;
        dApredt(i,:) = -Apre(i,:)/tau_LTP;
     
    % x update
        Apost(i+1,:) = Apost(i,:) + step*dApostdt(i,:);
        Apre(i+1,:) = Apre(i,:) + step*dApredt(i,:);
        W_IE(i+1,:) = W_IE(i,:);
     
% -- check for action potentials --
 
 
    % -- excitatory spike --
    
    for k = 1:N_E
        if v_E(i+1,k) >= 20 && v_E(i,k) < 20
             
            % reset & refractory
            v_E(i+1,k) = vreset;
            ref_E(1,k) = refractory;
             
            % spike monitor
            spike_E(i+1,k) = k;
            spike_E_time(i+1,k) = time(i+1,1);
            
            % spike monitor
            spike_E(i+1,k) = k;
            spt_E(i+1,k) = time(i+1,1) + comp_time;

            % check synaptic connection E to I
            % check post-synaptic neuron
            for j = 1:N_I
                if S_key_IE(k,j) ~= 0
                    index = S_key_IE(k,j);
                    % synaptic input from E to I : _(IE)
                    X_IE(i+1,j) = X_IE(i,j) + W_IE(i,index);
                     
                    % plasticity update - "on_pre"
                    Apre(i+1,index) = Apre(i,index) + dApre_0;
                    W_IE(i+1,index) = W_IE(i,index) + eta*a_LTD*Apost(i,index);
                    
                    % max synaptic weight check                       
                    if (J_I*W_IE(i+1,index)) < 10
                        W_IE(i+1,index) = 10/J_I;
                    end
                end
            end
            
            
            % E-to-E
            % check for postsynaptic connections
            for j = 1:N_E
                if S_key_EE(k,j) ~= 0
                    index = S_key_EE(k,j);
                    X_EE(i+1,j) = X_EE(i,j) +  WEE(1,index);
                end                
            end
            
            
        elseif ref_E(1,k) >= 0
             
            % check if in refractory period
            v_E(i+1,k) = vreset;
        elseif v_E(i+1,k) < vrest
            v_E(i+1,k) = vrest;
             
        end
    end
     
    % -- inhibitory neuron spike --    
    for k = 1:N_I
        if v_I(i+1,k) >= 20 && v_I(i,k) < 20
             
            % reset & refractory
            v_I(i+1,k) = vreset;
            ref_I(1,k) = refractory;
             
            % spike monitor
            spike_I(i+1,k) = k+N_E;
            spike_I_time(i+1,k) = time(i+1,1);
             
            % check synaptic connection I to E
            for j = 1:N_E
                % synaptic input from I to E : _(EI)
                if S_key_EI(k,j) ~=0  
                    index = S_key_EI(k,j);
                    X_EI(i+1,j) = X_EI(i,j)  - WEI(1,index);
                end    
                % plasticity update - "on_post"
                if S_key_IE(j,k) ~= 0
                    index = S_key_IE(j,k);
                    Apost(i+1,index) = Apost(i,index) + dApost_0;
                    W_IE(i+1,index) = W_IE(i,index) + eta*a_LTP*Apre(i,index);
                     
                    % max synaptic weight check
                    if (J_I*W_IE(i+1,index)) > 1000
                        W_IE(i+1,index) = 1000/J_I;                   
                    end
                end
            end
            
            % check synaptic connection I to I
            for j = 1:N_I
               if S_key_II(k,j) ~= 0
                   index = S_key_II(k,j);
                   X_II(i+1,j) = X_II(i,j) - WII(1,index);
               end
            end
             
        elseif ref_I(1,k) >= 0             
            % check if in refractory period
            v_I(i+1,k) = vreset;
             
        elseif v_I(i+1,k) < vrest
            v_I(i+1,k) = vrest;
        end
    end
     
end
% --- Synchrony calculation ---
    N = N_E;
    Vcomb = zeros(sample_duration/step+1,N);
    Vcomb(:,1:N_E) = v_E;
    V1 = mean(Vcomb,2); 
    
    % variance of average volatage over whole run
    sigma_squ_v = mean(V1.^2) - (mean(V1))^2;   % sigma_v^2 = <(V(t)^2>t -[<V(t)>]^2
    
    % variance of volatage at each time step
    sigma_vi = zeros(1,N);
    sum_sig = 0;
    
    for i = 1:N
       sigma_vi(i) = mean(Vcomb(:,i).^2) - (mean(Vcomb(:,i)))^2; 
       sum_sig = sum_sig + sigma_vi(i);
    end
   
    % calculate synchrony
    syn_squ = sigma_squ_v/(sum_sig/N);
    synchrony = syn_squ^0.5;
 
end