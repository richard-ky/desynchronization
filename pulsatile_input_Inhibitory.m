function Ui = pulsatile_input_Inhibitory(V_stim,T1,T2,T3,duration,step,stimon,sample_duration,t_pulseI,t_pulseE)
Ui = zeros(1,round(duration/step));
jj = stimon/step;
%jf = (duration-1000)/step;
duration1 = stimon+max(t_pulseI,t_pulseE);
jf = duration1/step;
% if mod(duration1,sample_duration) ~= 0
%     jf = duration1/step;
% else
%     jf = duration1/step;
% end
    % ui input
    counter = 0;
    counter1 = 1;
    for i = jj+1:jf
        counter = counter + 1;
        N1 = round(T1(counter1,1)/step);
        N2 = round((T1(counter1,1)+T2(counter1,1))/step);
        N3 = round((T1(counter1,1)+T2(counter1,1)+T3(counter1,1))/step);
        
        % Cathodic phase
        if counter > 0 && counter <= N1
            Ui(i) = V_stim(counter1,1);
        end
    
        % Anodic phase
        if counter > N1 && counter <= N2  
            Ui(i) = -V_stim(counter1,1)*T1(counter1,1)/T2(counter1,1);
        end
    
        if counter > N2   && counter <= N3  
            Ui(i) = 0;
        end
    
        %if counter >= (2+x)*T_stim  + step && counter < (2+x+multi-1)*T_stim
        %    Ue(i) = -V_stim/multi;
        %end
    
    
        if counter >= N3
            counter = 0;
            counter1 = counter1+1;
            %Ui(i) = 0;
        end
    
    
    end


end
