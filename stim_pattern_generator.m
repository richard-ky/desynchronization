
function [Ue,Ui,duration] = stim_pattern_generator(step,sample_duration,stimon,Amp,PW1,PW2,PW3,Ns)


% Generate Stimulation Pattern
    
   
    

    V_stimE = Amp*ones(Ns,1);
    T1E = PW1*ones(Ns,1);
    T2E = PW2*ones(Ns,1);
    T3E = PW3*ones(Ns,1);
    
    V_stimI = Amp*ones(Ns,1); 
    T1I = PW1*ones(Ns,1);
    T2I = PW2*ones(Ns,1);
    T3I = PW3*ones(Ns,1);
    
    
    
    t_pulseE = T1E+T2E+T3E;
    t_pulseI = T1I+T2I+T3I;
    
    
    duration = stimon+max(sum(t_pulseI),sum(t_pulseE))+2000;% ms
   
    if mod(duration,sample_duration) ~= 0
        duration = (floor(duration/sample_duration)+1)*sample_duration+2000;
    end
    
    
    Ue = pulsatile_input_Excitatory(V_stimE,T1E,T2E,T3E,duration,step,stimon,sample_duration,sum(t_pulseI),sum(t_pulseE));
    Ui = pulsatile_input_Inhibitory(V_stimI,T1I,T2I,T3I,duration,step,stimon,sample_duration,sum(t_pulseI),sum(t_pulseE));
end