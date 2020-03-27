function [resp,t,par] = balance_sim(stim_in, modelName, varargin)
% [resp,t,par_out] = balance_sim(stim,sr,model,cond,par_in)
%
% implemented models:
%       'DEC_p2015.slx'
%           cond = 'default', 'strobe', 'eyes_open'
%       'LT_p2017.slx'
%           cond = 'default', [height of touch reference hTP in cm]
%       'IC_surfaceTilt.slx'; eyes closed condition
%           cond = 0.5, 1, 2, 4, 8 (pp tilt amplitude)
%
%
%

%% parse variable inpu
    p = inputParser;
    addRequired(p,'stim_in');
    addRequired(p,'modelName');
    
    % options for model simulations
    default_par_fix = 1;
    addOptional(p,'par_fix',default_par_fix);
    default_par_var = 1;
    addOptional(p,'par_var',default_par_var);
    default_sr = 100;
    addParameter(p,'sr',default_sr,@isnumeric);
    default_condition = 'default';
    addParameter(p,'cond',default_condition, @ischar)
    
    parse(p,stim_in,modelName, varargin{:})

    sr = p.Results.sr;
    
%% get parameters from defaults or inputs
    par = p.Results.par_fix;
    if ~isstruct(par); [par, ~] = getDefaultPar(modelName); end
    par_var = p.Results.par_var;
    if ~isstruct(par_var); [~, par_var] = getDefaultPar(modelName); end
    
    fn = fieldnames(par_var);
    for m=1:length(fn)
        par.(fn{m}) = par_var.(fn{m})(1);
    end
    
    
%% run simulation
    nocellinput=0;
    if ~iscell(stim_in); stim_in{1}=stim_in; nocellinput=1; end
    
    for k = 1:length(stim_in)
        stim = stim_in{k}';
        sim_time = length(stim)/sr;
        t = (0 : 1/sr : sim_time-1/sr)';
        out = sim(modelName, 'SrcWorkspace', 'current');
        resp{k} = out.get('resp');
    end
    
    if nocellinput; resp = resp{1}; end
end

