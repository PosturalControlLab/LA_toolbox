function [resp_sim, par_out, tf_sim] = pcl_TFfit(stim_in, modelName, varargin)
% [resp_sim, par_out, tf_sim] = pcl_parfit(stim_in, modelName, ...)
% last modified 6 Sept 2018
% 
% 
% Runs parameter estimations or model simulations on any model also defined 
% in the function balanceSim.m
%
% stim: stimulus as cells containing sequences in rows. Each stimulus in a
% cell is repeated twice for time domain simulations and the second
% repetition is used for analysis to avoid transients.
% 
% resp: if given, parameter estimation will be performed. Expected input is
%       the response as cells containing cycle repetitions in rows. Average
%       across cycle repetitions are used for optimization.
% 
% modelName: name of the Simulink Transferfunction model
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
% Lorenz Asslaender
% University of Constance
% lorenz(at)asslaender.de
% 
% future improvements: 
% - Add calculation of GOF and other parameters suggested by Pasma et al.


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
    resp = p.Results.resp;
    do_fit = iscell(resp);
    
    par = p.Results.par_fix;
    if ~isstruct(par); [par, ~] = getDefaultPar(modelName); end
    par_var = p.Results.par_var;
    if ~isstruct(par_var); [~, par_var] = getDefaultPar(modelName); end
    
    
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

    % format stimulus for simulation
    stim=[];        
    for m=1:nt
        % two repetitions of each condition and concatenation of all stimuli
        stim = [stim; stim_in{m}'; stim_in{m}'];
        cl(m) = length(stim_in{m});
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn = fieldnames(par_var);
    pin = cell2mat(struct2cell(par_var));

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

        for m=1:nt
%             [~, tf{m}] = calculate_frf(so_stim{m}, so_resp{m},sr,selFreq{m});
            tf{m} = simTF(f{m});
        end
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

end  



    
    