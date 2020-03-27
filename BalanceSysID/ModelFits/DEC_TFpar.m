function [par_fix,par_var] = DEC_TFpar(Bh, Bm)

% parameter_definition

    Ap = WinterTable(Bm,Bh);
    h = Ap.h_com;

%% fixed parameters

    par_fix.J = mean(Ap.J_B)  /180*pi;
    mgh = mean(Bm)*9.81*mean(h)  /180*pi;
    par_fix.mgh = mgh;

    par_fix.Kp = mgh;

%% variable parameters

    par_var.PD = [0.3, 0.1, 0.8];   
    par_var.SG = [0.85, 0, 2];
%     par_var.dt_bf = [0.16,    0, 0.3];

    par_var.K = [0, 0, 0]; %[0.15*mgh, 0, mgh];
    par_var.D = [0, 0, 0]; %[0.05*mgh, 0, mgh];

    par_var.Flp =[15, 0, 50];
    par_var.Glp = [0.1,    0,    0.8];

    par_var.dt_dec = [0.16,    0, 0.3];
    par_var.Gfs = [0.5,    0,  1]; %default: 0.75
    par_var.tfs = [0.9,    0,    1]; %default: 0.33
    par_var.Gg  = [0.55,    0,  1.5];
%     par_var.tg  = [0.07,    0,    0.2];
    par_var.tg = [0, 0, 0];
    
    