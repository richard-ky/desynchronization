clear -all
clearvars
clc

N_E = 1600; %Number of Excitatory Neurons
N_I = 400; %Number of Inhibitory Neurons

%Network synaptic connections and initial weights
    % -- Make Random Synaptic Conncetions ---
    epsilon_E = 0.1; 
    epsilon_I = 0.1; 
    epsilon_EE = 0.1;
    epsilon_II = 0.1;
    S_key_IE = zeros(N_E,N_I);
    S_key_EI = zeros(N_I,N_E);
    S_key_EE = zeros(N_E,N_E);
    S_key_II = zeros(N_I,N_I);
    
    %--- I to E ---
    syn_count = 0;
    for pre_neuron = 1:N_I
        for post_neuron = 1:N_E
            x = rand;
            if x <= epsilon_I
                syn_count = syn_count + 1;
                S_key_EI(pre_neuron,post_neuron) = syn_count;
            else
                S_key_EI(pre_neuron,post_neuron) = 0;
            end
    
        end
    end
    
    % --- E to I ---
    syn_count = 0;
    for pre_neuron = 1:N_E
        for post_neuron = 1:N_I
            x = rand;
            if x <= epsilon_E
                syn_count = syn_count + 1;
                S_key_IE(pre_neuron,post_neuron) = syn_count;
            else
                S_key_IE(pre_neuron,post_neuron) = 0;
            end
    
        end
    end

    % --- E to E ---
    syn_count = 0;
    for pre_neuron = 1:N_E
        for post_neuron = 1:N_E
            x = rand;
            if x <= epsilon_EE && pre_neuron ~= post_neuron
                syn_count = syn_count + 1;
                S_key_EE(pre_neuron,post_neuron) = syn_count;
            else
                S_key_EE(pre_neuron,post_neuron) = 0;
            end
    
        end
    end
    

    
    % --- I to I ---
    syn_count = 0;
    for pre_neuron = 1:N_I
        for post_neuron = 1:N_I
            x = rand;
            if x <= epsilon_II && post_neuron ~= pre_neuron
                syn_count = syn_count + 1;
                S_key_II(pre_neuron,post_neuron) = syn_count;
            else
                S_key_II(pre_neuron,post_neuron) = 0;
            end
    
        end
    end
  
save('synaptic_connection.mat','S_key_EE','S_key_II','S_key_IE','S_key_EI')


%Initial Synaptic Weights

    num_synapses_IE = max(max(S_key_IE));
    num_synapses_EI = max(max(S_key_EI));
    num_synapses_EE = max(max(S_key_EE));
    num_synapses_II = max(max(S_key_II));

W_IE0 = unifrnd(0.5,1,1,num_synapses_IE); 
W_EI0 = unifrnd(0.5,1,1,num_synapses_EI); 
W_EE0 = unifrnd(0.5,1,1,num_synapses_EE);
W_II0 = unifrnd(0.5,1,1,num_synapses_II); 

save('initial_synaptic_weights.mat','W_IE0','W_EI0','W_EE0','W_II0')