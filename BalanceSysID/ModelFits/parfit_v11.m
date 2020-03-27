function [par_out, tf_sim, tf_exp, f, err] = parfit_v11(model, stim_t, resp, par_var, par_fix, do_fit)
% [par_out, tf_sim, tf_exp, f, err] = parfit_v11(model, stim, exp_com, par_var, par_fix, do_fit)
% last modified 25.06.2018
% 
% Runs parameter estimations or model simulations on any model also defined 
% in the function balanceSim.m
%
% stim, resp: stimulus and response as cells containing cycle repetitions in rows.
% par_fix and par_var are defined in the script DEC_TFpar
% do_fit: 1 executes fitting procedure
%         2 calculates transfer function of the DEC model based on input
%           parameters defined in par_var without fitting.
%
% 
% 
% 
% Lorenz Asslaender
% University of Constance
% lorenz(at)asslaender.de
% 
% future improvements: 
% - sampling rate of input data and number of frequency points considered
%   in the calculations should be defined externally.
% - Add calculation of GOF and other parameters suggested by Pasma et al.
% - function depends on function TDthreshold. Include within this function 
%   in future versions.
% 
% 

%% preallocate model parameters
        f=cell(size(resp)); tf_exp=cell(size(resp));

    % should be defined externally and input into the function in future versions
        sr = 100;

    par = par_fix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experimental transfer function and describing function for threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    nt = length(stim_t); % number of trials
    for k=1:nt
        [f{k},tf_exp{k}] = calculate_frf(stim_t{k}, resp{k}, 1000);
    end
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn = fieldnames(par_var);
    pin = cell2mat(struct2cell(par_var));
    
    
    
    if do_fit == 1
        x0 = pin(:,1); % starting value
        lb = pin(:,2); % lower bound
        ub = pin(:,3); % upper bound
    
        options=optimoptions('lsqnonlin',...
                                'DiffMinChange', 0.001,...
                                'MaxFunEvals', 1000);
                            
        [par2,~,~,~,~,~,~] = ...
                                lsqnonlin(@err_func,x0, lb, ub,options);
                            
    elseif do_fit ==0
        % get parameters to calculate transfer function without fit
        par2 = pin(:,1);
    end
    

%% calculate final function output
        for k=1:length(fn)
            eval(['par_out.' fn{k} '=' num2str(par2(k)) ';']);
        end

    tf_sim = func(par2);
    err = err_func(par2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nested subfunctions - calculation of TF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calculation of the objective function for the fit
    function err = err_func(xin)
        
        tf = func(xin);
        
        err=[];
        for m = 1:length(tf_exp)
            temp = ( abs(tf{m} - tf_exp{m})  ) ;% ./ sqrt(f{m});
            err = [err; temp'];
        end
    end
    
    % run simulation
    function tf = func(xin)
        
        for k=1:length(fn)
            eval(['par.' fn{k} '=' num2str(xin(k)) ';']);
        end

        
        nt = length(stim_t); % number of trials
        stim=[];
        
        for m=1:nt
            % format stimulus for simulation
            % two repetitions of each condition and concatenation of all stimuli
                stim = [stim; stim_t{m}'; stim_t{m}'];
                sim_time = length(stim)/sr;
                
                cl(m) = length(stim_t{m});
        end
                t = (0 : 1/sr : sim_time-1/sr)';
            % run simulink simulation
                out = sim(model, 'SrcWorkspace', 'current');

            % create output
                ct = get(out,'resp');
                st = get(out,'stim_out'); 
                
                
        for m=1:nt
                % cut simulation output into cycles, where always the first
                % cycle of each stimulus is discarded
                so_resp{m} = ct(sum(cl(1:m)) + sum(cl(1:m-1)) + 1 : 2*sum(cl(1:m)))'; 
                so_stim{m} =st(sum(cl(1:m)) + sum(cl(1:m-1)) + 1 : 2*sum(cl(1:m)))';

                [~, tf{m}] = calculate_frf(so_stim{m}, so_resp{m},sr);
        end
        % saves time domain simulation results if activated
        save('sim_outputTD','so_*')
    end
            
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

% frequency response function from time domain stimulus and response seq.
function [f_out,FRF_out] = calculate_frf(stim, resp, sr)
    
        yi = fft(stim);
        yi = yi(2:66);
        yo = fft(resp);
        yo = yo(2:66);
                
        FRF = yo./yi;
        f = (1:65)/(length(stim)/sr);
            
        f2   = f(1:2:end);
        FRF2 = FRF(1:2:end);
            
        f_out = log_smooth(f2,15);
        FRF_out = log_smooth(FRF2,15);
%         
%     % number of frequency points
%         nfp = 1:50;
%     % define frequency points being taken into account
%         FRF_out = FRF2(nfp);
%         f_out = f2(nfp);

end  



    
    