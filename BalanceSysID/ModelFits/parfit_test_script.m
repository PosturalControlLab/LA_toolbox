%% parfit test script

load test_data.mat

modelName = 'DEC_p2015.slx';

        for k=1:3; 
            stim_in{k} = mean(stim{k},1); 
            resp_in{k} = mean(resp{k},1); 
        end
        
%         [resp_sim, par_out] = pcl_parfit(stim_in, modelName, 'par_var', par_var, 'par_fix', par_fix);
        [resp_sim, par_out] = pcl_parfit(stim_in, modelName);
        
        [FD, TD] = getFRF(stim_in,resp_sim,'sr',100); 
        plotFRF(FD,TD)