function [resp_sim, par_out, tf_sim] = parfit(stim_in, modelName, varargin)
% [par_out, tf_sim, tf_exp, f, err] = parfit(stim, resp, model, do_fit, par_var, par_fix, ...)
% last modified 5 Aug 2018
% 
% 
% Runs parameter estimations or model simulations on any model also defined 
% in the function balanceSim.m
%
% stim: stimulus as cells containing sequences in rows. Each stimulus in a
% cell is repeated twice for time domain simulations and the second
% repetition is used for analysis to avoid transients.
% 
% resp: response as cells containing cycle repetitions in rows. Average
% across cycle repetitions are used to calculate the objective function.
% 
% model
% 
% par_fix: fixed model parameters in a structure, not used for optimization
% par_var: parameters subject to optimization in a structure containing
%          [starting value, lower_bound, upper_bound].
%          par_fix and par_var are defined in the script DEC_TFpar
% 
% do_fit: 1 executes fitting procedure
%         0 calculates transfer function of the DEC model based on input
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


%% handle variable input

    p = inputParser;
    addRequired(p,'stim_in');
    addRequired(p,'modelName');
    
    % options for model simulations
    default_resp = 0;
    addOptional(p,'resp', default_resp);
    default_par_fix = 1;
    addOptional(p,'par_fix',default_par_fix);
    default_par_var = 1;
    addOptional(p,'par_var',default_par_var);
    
    % options for calculating objective function
    default_sr = 100;
    addParameter(p,'sr',default_sr,@isnumeric);
    defaultFreqPoints = [1, 2, 2];
    addParameter(p,'FreqPoints',defaultFreqPoints,@isnumeric);   
    defaultSmoothFRF = 'None';
    addParameter(p,'SmoothFRF', defaultSmoothFRF); %number of frequency points
    
    parse(p,stim_in,modelName,varargin{:})



%% handle variable input and preallocate model parameters
   f=cell(size(stim_in)); tf_exp=cell(size(stim_in));
    
    sr  = p.Results.sr;
    par = p.Results.par_fix;
    resp = p.Results.resp;
        if iscell(resp); do_fit = 1;
        else             do_fit = 0;
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experimental transfer function and describing function for threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    nt = length(stim_in); % number of trials    
    
    FreqPoints = p.Results.FreqPoints;
    if length(FreqPoints)>3 && FreqPoints(4)==1
        selFreq = FreqPoints(1):FreqPoints(2):FreqPoints(3);
    else
        for k=1:nt
            f_b = 1/length(stim_in{k})*sr;
            selFreq{k} = FreqPoints(1) : FreqPoints(2) : floor(FreqPoints(3)/f_b);
            f{k} = selFreq{k} * f_b;
        end
    end
    
    % calculate experimental FRF
    if iscell(resp)
        for k=1:nt
            [~,tf_exp{k}] = calculate_frf(stim_in{k}, resp{k}, sr, selFreq{k});
        end
    end
    
    if strcmp('tf',modelName(1:2))
        modelfun = str2func(modelName);
    else
        
        % format stimulus for simulation
        stim=[];        
        for m=1:nt
            % two repetitions of each condition and concatenation of all stimuli
            stim = [stim; stim_in{m}'; stim_in{m}'];
            cl(m) = length(stim_in{m});
        end
    end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn = fieldnames(p.Results.par_var);
    pin = cell2mat(struct2cell(p.Results.par_var));

    if do_fit
        x0 = pin(:,1); % starting value
        lb = pin(:,2); % lower bound
        ub = pin(:,3); % upper bound
    
        options=optimoptions('lsqnonlin',...
                                'DiffMinChange', 0.001,...
                                'MaxFunEvals', 1000);
                            
        [par2,~,~,~,~,~,~] = ...
                                lsqnonlin(@err_func,x0, lb, ub,options);
                            
    else
        % get parameters to calculate transfer function without fit
        par2 = pin(:,1);
    end
    

%% calculate final function output

    for k=1:length(fn)
        par_out.(fn{k}) = par2(k);
    end

    [tf_sim, resp_sim] = func(par2);
    if do_fit; par.err = err_func(par2); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nested subfunctions - calculation of TF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calculation of the objective function for the fit
    function err = err_func(xin)
        
        tf = func(xin);
        
        err=[];
        for m = 1:length(tf_exp)
            temp = ( abs(tf{m} - tf_exp{m})  ) ./ sqrt(f{m});
            err = [err; temp'];
        end
    end
    
    % run simulation
    function [tf, so_resp] = func(xin)
        
        for m=1:length(fn)
            par.(fn{m}) = xin(m);
        end
        
        if ~strcmp('slx',modelName(end-2:end))
            for n=1:nt
                tf{k} = modelfun(f{k},stim{k});
            end
            so_resp = 0;
        else
            sim_time = length(stim)/sr;         % required for Simulink settings
            t = (0 : 1/sr : sim_time-1/sr)';    % required for Simulink settings

            % run simulink simulation
            out = sim(modelName, 'SrcWorkspace', 'current');

            % create output
            ct = get(out,'resp');
            st = get(out,'stim_out'); 

            for m=1:nt
                % cut simulation output into cycles, where always the first
                % cycle of each stimulus is discarded
                so_resp{m} = ct(sum(cl(1:m)) + sum(cl(1:m-1)) + 1 : 2*sum(cl(1:m)))'; 
                so_stim{m} = st(sum(cl(1:m)) + sum(cl(1:m-1)) + 1 : 2*sum(cl(1:m)))';

                [~, tf{m}] = calculate_frf(so_stim{m}, so_resp{m},sr,selFreq{m});
            end
        end
        % saves time domain simulation results if activated
        % save('sim_outputTD','so_*')
    end
            

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

% frequency response function from time domain stimulus and response seq.
function [f,FRF] = calculate_frf(stim, resp, sr, selFreq)
    
        yi = fft(stim);
        yi = yi(selFreq+1); % +1 to exclude offset value of fft output
        yo = fft(resp);
        yo = yo(selFreq+1);
                
        FRF = yo./yi;
        f = selFreq/length(stim)*sr;
            
%         f_out = log_smooth(f2,15);
%         FRF_out = log_smooth(FRF2,15);
%         
%     % number of frequency points
%         nfp = 1:50;
%     % define frequency points being taken into account
%         FRF_out = FRF2(nfp);
%         f_out = f2(nfp);

end  



    
    