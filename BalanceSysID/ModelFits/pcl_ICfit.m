function FD = pcl_ICfit(FD,varargin)
% [par, sim_FRF] = IC_fit(FD, SubjInfo, type)
% 
% type: different types are amplitude specific parameters and involvement
% of passive contributions (pas; NoPas)
%
%

    p = inputParser;
    addRequired(p,'FD');
    defaultSubjInfo.body_weight= 75;
    defaultSubjInfo.body_height = 175;
    addOptional(p,'SubjInfo',defaultSubjInfo,@isstruct);
    defaultPassive = 0;
    addOptional(p,'Passive', defaultPassive);

    parse(p,FD,varargin{:})
    
    Passive = p.Results.Passive;
    SubjInfo = p.Results.SubjInfo;

%% optimization    
for k=1:size(FD,1)
    for m=1:size(FD,2)
            
    % get Parameters and fit/ FRF calculation
    
    [par_fix,par_var] = IC_getPar(SubjInfo(k,m).body_weight, SubjInfo(k,m).body_height/100, Passive);
    theta_fix = [par_fix.J,par_fix.mgh];

    x=par_var;
    theta = [x.Kp(1),x.PD(1),x.dt(1),x.W(1),x.Flp(1),x.Glp(1),x.K(1),x.D(1)];
    lb = [x.Kp(2),x.PD(2),x.dt(2),x.W(2),x.Flp(2),x.Glp(2),x.K(2),x.D(2)];
    ub = [x.Kp(3),x.PD(3),x.dt(3),x.W(3),x.Flp(3),x.Glp(3),x.K(3),x.D(3)];
                
        
        
        options=optimoptions('lsqnonlin',...
                                'DiffMinChange', 0.001,...
                                'MaxFunEvals', 10000,...
                                'MaxIter', 10000);
        
        theta_out = lsqnonlin(@(theta) objectiveFunc(theta,theta_fix, FD(k,m).f, FD(k,m).FRF), theta, lb, ub, options);
    

%% output
    % transfer function
        FD(k,m).simFRF = ICmodel( FD(k,m).f,theta_out,theta_fix);
        FD(k,m).simErr = objectiveFunc(theta_out,theta_fix, FD(k,m).f, FD(k,m).FRF);
        
    % parameter structure
        theta_out=num2cell(theta_out);
        [y.Kp,y.PD,y.dt,y.W,y.Flp,y.Glp,y.K,y.D] = theta_out{:};
        FD(k,m).simPar = y;
        FD(k,m).simPar.J = par_fix.J;
        FD(k,m).simPar.mgh = par_fix.mgh;
        clear theta_out y 
    end
end
end

function tf = ICmodel(f,theta,theta_fix)

    theta = num2cell(theta);
    [Kp,PD,dt,W,Flp,Glp,K,D] = theta{:}; 

    theta_fix = num2cell(theta_fix);
    [J,mgh] = theta_fix{:};

    s = 1i*f*2*pi;

    B = 1./(J.*s.^2-mgh);
    NC = Kp.*(1 + PD.*s);
    TD = exp(-s.*dt);

    P = K+D.*s;
    F = (Glp./(Flp.*s+1) );

    tf = (P.*B + W.*NC.*B.*TD) ./ (1 - F.*NC.*TD + P.*B + NC.*B.*TD);
end
    
function err = objectiveFunc(theta,theta_fix,f,TFexp)

    tf = ICmodel(f,theta,theta_fix);
        
    err = abs( tf - TFexp ) ./sqrt(f) ;
    
end
       
function [par_fix,par_var] = IC_getPar(Bm, Bh,Passive)

% parameter_definition

    Ap = WinterTable(Bm,Bh);
    h = Ap.h_com;
    mgh = mean(Bm)*9.81*mean(h)  /180*pi;

%% fixed parameters

    par_fix.J = mean(Ap.J_B)  /180*pi;
    par_fix.mgh = mgh;

%% variable parameters

switch Passive
    case 0
        par_var.Kp = [mgh*1.3, 0.5*mgh, 3*mgh];
        par_var.PD = [0.3, 0.1, 0.8];   

        par_var.K = [0, 0, 0];
        par_var.D = [0, 0, 0];

        par_var.Flp =[15, 5, 50];
        par_var.Glp = [0.01,    0,    0.8];

        par_var.dt = [0.1, 0.05, 0.3];
        par_var.W  = [0.4,    0,    1];
    case 1
        par_var.Kp = [mgh, 0.5*mgh, 3*mgh];
        par_var.PD = [0.3, 0.1, 0.8];   

        par_var.K = [0.15*mgh, 0, mgh];
        par_var.D = [0.05*mgh, 0, mgh];

        par_var.Flp =[15, 5, 50];
        par_var.Glp = [0.01,    0,    0.8];

        par_var.dt = [0.16, 0.05, 0.3];
        par_var.W  = [0.4,    0,    1];
end



end
