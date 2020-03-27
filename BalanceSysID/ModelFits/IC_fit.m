function [par, sim_FRF, err] = IC_fit(varargin)
% [par, sim_FRF] = IC_fit(type, f, exp_FRF, Bm, Bh)
% 
% type: different types are amplitude specific parameters and involvement
% of passive contributions (pas; NoPas)
%
%

%% variable input handling

type = varargin{1};
f = varargin{2};
if nargin<3; 
    do_fit = 0; 
else do_fit = 1; 
    exp_FRF = varargin{3}; 
end

if nargin<5; 
    Bh = 1.75; 
    Bm = 75; 
else
    Bh = varargin{5}; 
    Bm = varargin{4}; 
end

%% get Parameters and fit/ FRF calculation
    
    [par_fix,par_var] = IC_getPar(Bm,Bh,type);

    theta_fix = [par_fix.J,par_fix.mgh];

        x=par_var;
        theta = [x.Kp(1),x.PD(1),x.dt(1),x.W(1),x.Flp(1),x.Glp(1),x.K(1),x.D(1)];
        lb = [x.Kp(2),x.PD(2),x.dt(2),x.W(2),x.Flp(2),x.Glp(2),x.K(2),x.D(2)];
        ub = [x.Kp(3),x.PD(3),x.dt(3),x.W(3),x.Flp(3),x.Glp(3),x.K(3),x.D(3)];
        
%% optimization
    
    if do_fit
        options=optimoptions('lsqnonlin',...
                                'DiffMinChange', 0.001,...
                                'MaxFunEvals', 10000,...
                                'MaxIter', 1000);
        
        theta_out = lsqnonlin(@(theta) objectiveFunc(theta,theta_fix,f,exp_FRF), theta, lb, ub, options);
    else
        theta_out = theta;
    end
    

%% output
    % transfer function
        sim_FRF = ICmodel(f,theta_out,theta_fix);
        err = objectiveFunc(theta_out,theta_fix,f,exp_FRF);
        
    % parameter structure
        theta_out=num2cell(theta_out);
        [y.Kp,y.PD,y.dt,y.W,y.Flp,y.Glp,y.K,y.D] = theta_out{:};
        par = y;
        par.J = par_fix.J;
        par.mgh = par_fix.mgh;
    
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
       
function [par_fix,par_var] = IC_getPar(Bm, Bh,type)

% parameter_definition

    Ap = WinterTable(Bm,Bh);
    h = Ap.h_com;
    mgh = mean(Bm)*9.81*mean(h)  /180*pi;

%% fixed parameters

    par_fix.J = mean(Ap.J_B)  /180*pi;
    par_fix.mgh = mgh;

%% variable parameters

switch type
    case 1
        par_var.Kp = [mgh*1.3, 0.5*mgh, 3*mgh];
        par_var.PD = [0.3, 0.1, 0.8];   

        par_var.K = [0, 0, 0];
        par_var.D = [0, 0, 0];

        par_var.Flp =[15, 5, 50];
        par_var.Glp = [0.01,    0,    0.8];

        par_var.dt = [0.1, 0.05, 0.3];
        par_var.W  = [0.4,    0,    1];
    case 2
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
