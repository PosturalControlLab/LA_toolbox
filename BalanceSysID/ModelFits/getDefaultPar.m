function [par_fix, par_var] = getDefaultPar(varargin)
% [par_fix, par_var] = getDefaultPar(modelName, cond, allFix, m, h, J)
% handle variable input
    p = inputParser;
    default_modelName = 'DEC_p2015.slx';
    addOptional(p,'modelName',default_modelName,@ischar);
    default_cond = 'default';
    addOptional(p,'cond',default_cond);
    default_allFix = 0;
    addOptional(p,'allFix',default_allFix);
    default_m = 69;
    addParameter(p,'m',default_m, @isnumeric);
    default_hB = 1.7;
    addParameter(p,'hB',default_hB, @isnumeric);
    default_J = 72.8;
    addParameter(p,'J',default_J, @isnumeric);
    
    parse(p,varargin{:})

    
    m = p.Results.m;
    hB = p.Results.hB;
    WT = WinterTable(m,hB);
    J = WT.J_B /180*pi;
    h = WT.h_com;

    cond = p.Results.cond;
    par_var=[];
    
switch p.Results.modelName
    case 'DEC'
        par_fix.J = J;
        par_fix.mgh = mean(m)*9.81*mean(h)  /180*pi;
        par_fix.Kp = par_fix.mgh;

        switch cond
            case 'default'
                par_var.Flp = [23 10 30];
                par_var.Glp = [0.002 0 0.1];
                par_var.PD  = [0.3 0 1];
                par_var.dt  = [0.16 0 0.3];
                par_var.Gfs = [0.75 0 1];
                par_var.tfs = [0.33 0 1];
                par_var.Gg  = [0.55 0 1];
                par_var.SG  = [1 0.5 1.5]; 
        end
    case 'DEC_pas'
        par_fix.J = J;
        par_fix.mgh = mean(m)*9.81*mean(h)  /180*pi;
        par_fix.Kp = par_fix.mgh;

        switch cond
            case 'default'
                par_var.Flp = [23 10 30];
                par_var.Glp = [0.002 0 0.1];
                par_var.PD  = [0.3 0 1];
                par_var.P  = [0.1 0 1].*par_fix.mgh;
                par_var.D  = [0.03 0 1].*par_fix.mgh;
                par_var.dt  = [0.16 0 0.3];
                par_var.Gfs = [0.75 0 1];
                par_var.tfs = [0.33 0 1];
                par_var.Gg  = [0.55 0 1];
                par_var.SG  = [1 0.5 1.5]; 
        end
        
        
    case 'DEC_LA1'
        par_fix.J = J;
        par_fix.mgh = mean(m)*9.81*mean(h)  /180*pi;
        par_fix.Kp = par_fix.mgh;

        switch cond
            case 'default'
                par_var.Flp = [15 10 30];
                par_var.Glp = [0.002 0 0.1];
                par_var.PD  = [0.3 0 1];
                par_var.dt  = [0.16 0 0.3];
                par_var.dt2 = [0.16 0 0.3];
                par_var.Gfs = [0.75 0 1];
                par_var.tfs = [0.33 0 1];
                par_var.Gg  = [0.55 0 1];
                par_var.SG  = [1 0.5 1.5]; 
                par_var.tc_leaky = [15 5 40];
        end
    case 'DEC_LA2'
        par_fix.J = J;
        par_fix.mgh = mean(m)*9.81*mean(h)  /180*pi;
        par_fix.Kp = par_fix.mgh;

        switch cond
            case 'default'
                par_var.Flp = [15 10 30];
                par_var.Glp = [0.002 0 0.1];
                par_var.PD  = [0.3 0 1];
                par_var.dt  = [0.16 0 0.3];
                par_var.dt2 = [0.16 0 0.3];
                par_var.Gfs = [0.75 0 1];
                par_var.tfs = [0.33 0 1];
                par_var.Gg  = [0.55 0 1];
                par_var.SGp  = [1 0.5 1.5]*par_fix.mgh; 
                par_var.SGd  = [1 0.5 1.5]*par_fix.mgh*0.3; 
                par_var.tc_leaky = [15 5 40];
        end
    
    case 'DEC_p2015.slx'
        par_fix.J = J;
        par_fix.mgh = mean(m)*9.81*mean(h)  /180*pi;

        par_fix.Kp = par_fix.mgh;
        
        par_fix.Kd = par_fix.mgh.*0.3;
        par_fix.SG = 0.85; 

        par_fix.K = par_fix.mgh* 0.15;     
        par_fix.D = par_fix.mgh*0.15* 0.3; 
        par_fix.Flp = 15;   
        switch cond
            case 'default'
                par_fix.K = 0;
                par_fix.D = 0;
                par_fix = rmfield(par_fix,'Flp');
                par_var.Flp = [15 10 30];
                par_fix = rmfield(par_fix,'Kd');
                par_var.Kd = par_fix.mgh.* [0.3 0 1];

                par_var.dt = [0.16 0 0.3];
                par_var.Glp = [0.13 0 10];
                par_var.Gfs = [0.75 0 1];
                par_var.tfs = [0.33 0 1];
                par_var.Gg  = [0.55 0 1];
                par_fix.tg = 0;
            case 'eyes_closed'
                par_fix.dt  = 0.16;
                par_fix.Glp = 0.13;
                par_fix.Gfs = 0.75;
                par_fix.tfs = 0.33;
                par_fix.Gg  = 0.55;
                par_fix.tg  = 0.07;
            case 'strobe'
                par_fix.dt = 0.16;
                par_fix.Glp = 0.15;
                par_fix.Gfs = 0.87;
                par_fix.tfs = 0.33;
                par_fix.Gg  = 0.46;
                par_fix.tg  = 0.02;
            case 'eyes_open'
                par_fix.dt = 0.15;
                par_fix.Glp = 0.13;
                par_fix.Gfs = 0.89;
                par_fix.tfs = 0.20;
                par_fix.Gg  = 0.46;
                par_fix.tg  = 0.03;
            case 'ec_no_tg'
                par_fix.dt = 0.16;
                par_fix.Glp = 0.13;
                par_fix.Gfs = 0.75;
                par_fix.tfs = 0.33;
                par_fix.Gg  = 0.55;
                par_fix.tg  = 0;    
        end
        
    case 'LT_p2017.slx'
        par_fix.J = J;
        mgh = mean(m)*9.81*mean(h)  /180*pi;
        par_fix.mgh = mgh;
        par_fix.h = h;

        par_fix.Flp = 20;

        par_var.Kp  = [1.1883 0.5 2];
        par_var.Kd =  [0.2967 0.1 1];
        par_var.Glp = [0.0021759 0 1];
        par_var.dt = [0.10592 0 0.5];
        par_var.Ge = [0.58129 0 1];
        par_var.te = [0.0028579 0 0.2];
        par_var.LT  = [0.28976 0 1];
        
        if ischar(cond)
            par_fix.hTP = 80;
        else
            par_fix.hTP = cond;
        end
    case 'IC_surfaceTilt'     
        par_fix.J = J;
        par_fix.mgh = mean(m)*9.81*mean(h)  /180*pi;

        par_var.Kd = [5 0 10];
        par_var.Kf = [0.1 0 1];
        par_var.Tf = [18 0 30];
        
        switch cond
            case 0.5
                par_fix.Kp = 17;
                par_fix.dt = 0.20;
                par_fix.W = 0.7;
            case 1
                par_fix.Kp = 18;
                par_fix.dt = 0.19;
                par_fix.W = 0.6;
            case 2
                par_fix.Kp = 19;
                par_fix.dt = 0.18;
                par_fix.W = 0.5;
            case {4, 'defalt'}
                par_var.Kp = [20 10 30];
                par_var.dt = [0.13 0 0.5];
                par_var.W =  [0.4 0 1];
            case 8
                par_fix.Kp = 21;
                par_fix.dt = 0.11;
                par_fix.W = 0.22;
        end

end


if p.Results.allFix
    fn = fieldnames(par_var);
    for k = 1:length(fn)
        par_fix.(fn{k}) = par_var.(fn{k})(1);
    end
    par_var=[];
end