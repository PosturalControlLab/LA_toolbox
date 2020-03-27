function [par_out, tf_sim, tf_exp, f, err] = parfit_v10(model, stim, resp, par_var, par_fix, do_fit)
% [par_out, tf_sim, tf_exp, f, err] = parfit_v10(model, stim, exp_com, par_var, par_fix, do_fit)
% Version 2.0s
% last modified 14.12.2016
% 
% par_fix and par_var are defined in the script parameter_DEC_TFfit
% stim: stimulus sequences; Expected is the average across cycle
%       repetitions as a 1D cell array containing row(?) vectors
% stim: stimulus sequences; Expected is the average across cycle
%       repetitions as a 1D cell array containing row (?) vectors.
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
        sr = 50; % 
%         nfp = 1:12; % number of frequency points considered. Frequencies are base frequency and odd harmonics.

% predefine model parameters
    J=0; mgh=0; Kp=0; Kd=0; Ghs=0; ths=0; Gg=0; Flp=0; Glp=0; dt=0;
    SG=0; hTP=0; height_TP=0; Wh=0; PD=0; LT=0; Ge=0; te=0; h=0; dt_dec=0;
    Gfs=0; tfs=0; tg=0; par=0;
       
% save fixed model parameters in function workspace
    fn = fieldnames(par_fix);
    for k=1:length(fn)
        eval([fn{k} '=[' num2str(par_fix.(fn{k})) '];']);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experimental transfer function and describing function for threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    nt = length(stim); % number of trials
    for k=1:nt
        [f{k},tf_exp{k}] = calculate_frf(stim{k}, resp{k}, sr);
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
                                'MaxFunEvals', 400*9);
                            
        [par2,~,~,~,~,~,~] = ...
                                lsqnonlin(@err_func,x0, lb, ub,options);
                            
    else
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
            temp = ( abs(tf{m} - tf_exp{m})  ) ./ sqrt(f{m});
            err = [err; temp'];
        end
    end
    
    % run simulation
    function tf = func(xin)
        
        for k=1:length(fn)
            eval(['par.' fn{k} '=' num2str(xin(k)) ';']);
        end

        
        nt = length(stim); % number of trials
        for m=1:nt
            switch height_TP(m)
                case 1; hTP = 0.8;
                case 2; hTP = 1.2;
            end
            % format stimulus for simulation
            % two repetitions of each condition and concatenation of all stimuli
                stim_si = [stim{m} stim{m}]';
                sim_time = length(stim_si)/sr;
                t = (0 : 1/sr : sim_time-1/sr)';
                cl = length(stim{m});

            % run simulink simulation
                out = sim(model, 'SrcWorkspace', 'current');

            % create output
                ct = get(out,'com');
                st = get(out,'stim'); 
                so_resp{m} = ct(cl+1:end)'; 
                so_stim{m} =st(cl+1:end)';

                [~, tf{m}] = calculate_frf(so_stim{m}, so_resp{m},sr);
        end
        % saves time domain simulation results if activated
%         save('sim_outputTD','so_*')
    end
            
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

% frequency response function from time domain stimulus and response seq.
function [f_out,FRF_out] = calculate_frf(stim, resp, sr)
    
        yi = fft(stim);
        yi = yi(2:101);
        yo = fft(resp);
        yo = yo(2:101);
                
        FRF = yo./yi;
        f = (1:100)/(length(stim)/sr);
            
        f2   = decimate2(f,2);
        FRF2 = decimate2(FRF,2);
            
    % number of frequency points
        nfp = 1:10;
    % define frequency points being taken into account
        FRF_out = FRF2(nfp);
        f_out = f2(nfp);

end  



    
    