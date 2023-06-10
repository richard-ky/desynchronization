function [RE] = kuramoto_syn(sptime,t,step,tf,N)
tspE = sptime;
phi = zeros(size(t,2),N);

% make vector of phases
for neuron = 1:N
        second_spike = 0;
    for i = 1:(tf)/step-1
        phi(i,neuron) = 2*pi*(t(i) - tspE(i,neuron));
        
        if tspE(i+1,neuron) ~= tspE(i,neuron)            
            if second_spike == 1
                delt = tspE(i+1,neuron) - tspE(i,neuron);
                a = floor(tspE(i,neuron)/step);
                b = floor(tspE(i+1,neuron)/step);
                unnorm = phi(a:b,neuron);
                phi(a:b,neuron) = unnorm./delt;
            end
            second_spike = 1;
        end
    end
end
RE = abs(mean((exp(1i.*phi)),2));
end