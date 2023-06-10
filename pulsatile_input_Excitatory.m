function Ue = pulsatile_input_Excitatory(V_stim,T1,T2,T3,duration,step,stimon,sample_duration,t_pulseI,t_pulseE)
Ue = zeros(1,round(duration/step));
jj = stimon/step;

%jf = duration/step;
duration1 = stimon+max(t_pulseI,t_pulseE);
jf = duration1/step;
% if mod(duration1,sample_duration) ~= 0
%     jf = duration1/step;
% else
%     jf = duration1/step;
% end
    % ue input
    counter = 0;
    counter1 = 1;
    
    for i = jj+1:jf
        counter = counter + 1;
            N1 = round(T1(counter1,1)/step);
            N2 = round((T1(counter1,1)+T2(counter1,1))/step);
            N3 = round((T1(counter1,1)+T2(counter1,1)+T3(counter1,1))/step);
        % anodic phase
        if counter > 0 && counter <= N1
            Ue(i) = -V_stim(counter1,1)*T2(counter1,1)/T1(counter1,1);
            
        end
    
        % cathodic phase
        if counter > N1 && counter <= N2
            Ue(i) = V_stim(counter1,1);
            
        end
    
        if counter > N2  + step && counter <= N3
            Ue(i) = 0;
        end
    
    
    
        if counter >= N3
            counter = 0;
            counter1 = counter1+1;
            
        end
    
    
    end


end
