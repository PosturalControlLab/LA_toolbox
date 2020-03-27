function [par_out, tf_sim, tf_exp, f, fitQ] = DEC_TFfit(stim, resp, par_fix, par_var, sr, do_fit)
% [par_out, tf_sim, tf_exp, f, fitQ] = DEC_TFfit(stim, resp, par_fix, par_var, sr, do_fit)
% Performs DEC model fits to experimental data of multiple sessions.
% 
% Inputs:
% stim, resp: stimulus and response as cells containing cycle repetitions in rows.
% par_fix, par_var: fixed and variable parameters as structure. Default
% values are given in DEC_TFpar.m.
% tf_sim, tf_exp: frequency response functions of simulations and
% experiments
% sr: sampling frequency
% do_fit: boolean - if set to zero, only FRFs are calculated from default
% parameters, without any optimization.
% This function contrains a transfer function formulation of the nonlinear
% DEC model. Nonlinear elements are handled by creating lookup tables
% containing transfer function formulations of the threshold input to
% output mapping for the threshold min-max bounds.


%% preallocate model parameters
        f=cell(size(resp)); tf_exp=cell(size(resp)); Coh=cell(size(resp));
        Gth=cell(size(resp)); Gvth=cell(size(resp)); vstep=zeros(size(resp));
        J=0; mgh=0; Kp=0; K=0; D=0; Gfs=0; tfs=0; Gg=0;  
        tg=0; Flp=0; Glp=0; SG=0; dt_dec=0; dt_bf=0; PD=0;
    
%% Set fixed parameter values and bounds for fit parameters
    % save fixed model parameters in function workspace
    if ~isempty(par_fix)
        fn = fieldnames(par_fix);
        for k=1:length(fn)
            eval([fn{k} '=' num2str(par_fix.(fn{k})) ';']);
        end
    end

    % define variable model parameters and parameter bounds for fit
    if ~isempty(par_var)
        pvn = fieldnames(par_var);
        pv = cell2mat(struct2cell(par_var));
    end
       
    
%% define required parameters for TF
    for n=1:length(resp);
        [f{n}, tf_exp{n}, Coh{n}, Gth{n}, vstep(n), Gvth{n}] = getTFpar(stim{n},resp{n},pvn,pv,sr);
    end
    
%% optimization
    
    if do_fit
        options=optimset('TolX', 1e-5,'TolFun', 1e-5, 'MaxFunEvals', 400*13, 'DiffMinChange', 1e-4);
%         [par2, fval, exitflag, output] = fminsearchbnd(@err_func,pv(:,1),pv(:,2),pv(:,3),options);
        [par2,resn,res,~,~,~,jacobian] = lsqnonlin(@err_func,pv(:,1),pv(:,2),pv(:,3),options);
        
        % calculate standard error of the mean or confidence bounds
%         ci = nlparci(par2,res,'Jacobian',jacobian);
        N=length(res);
        hessian = full(jacobian'*jacobian);
        CovP = 1/N*resn*inv(hessian);
        SEM = sqrt(diag(abs(CovP)));
    else
        par2 = cell2mat(struct2cell(par_var));
        par2 = par2(:,1);
        SEM=0;
    end
    

%% calculate final function output
    for k=1:length(pvn)
        par_out.(pvn{k}) = par2(k);
    end
    tf_sim = func(par2);

    for k = 1:length(tf_sim)
        tf_sim{k} = tf_sim{k};
        tf_exp{k} = tf_exp{k};
        f{k} = f{k};
    % calculate goodness of fit
        gt = sum((abs(tf_sim{k}-tf_exp{k})).^2) ./ sum(abs(tf_exp{k}).^2);
        GOF(k) = (1-gt)*100;
    end
    
    fitQ.gof = GOF;
    if do_fit
        fitQ.resn = resn;
        fitQ.sem = SEM;
        fitQ.jacobian = jacobian;
        fitQ.hessian = hessian;
    end
    
%% nested subfunctions - calculation of error
    function err = err_func(xin)

        tf = func(xin);
        temp=zeros(length(tf_exp),length(tf_exp{1}));
        
        for m = 1:length(tf_exp)
%             temp(m,:) = ( abs(tf{m} - tf_exp{m}).^2  ) ./ (1+f{m});
            temp(m,:) = sqrt(Coh{m} ./ (1+f{m})) .* abs(log(tf_exp{m}./tf{m}));
        end
        err=temp(:);
%         err=sum(err.^2);
    end
    
%% model
    function tf = func(varargin)
        if iscell(varargin); varargin=varargin{1}; end
        % name variable parameters for model simulations according to
        % subfunction input varargin
        for m=1:length(pvn)
            eval([pvn{m} '=' num2str(varargin(m)) ';']);
        end
       
        tf=cell(size(tf_exp));
        
        for m = 1:length(tf_exp)
            s = 1i* f{m}*2*pi;

            B = 1./(J.*s.^2-mgh);
            NC = Kp + Kp*PD.*s;
            TD_dec = exp(-s.*dt_dec);
            TD_bf = exp(-s.*dt_dec);
            
%             P = K+D.*s;
            F = (Glp./(Flp.*s+1) );
%             b = Gfs* (vstep(m)-tfs)/vstep(m);
            b = Gfs* getGth(tfs, Gvth{m});
            if b<0; b=0; end
            c = Gg;%* getGth(tg,Gth{m});
            
            tf{m} = (SG.*(TD_bf - b.*TD_dec) .*NC.*B) ./ (1 -F.*NC.*TD_dec + (SG.*TD_bf + c.*TD_dec).*NC.*B);
        end
    end


end

%% other subfunctions

function [f, tf_exp, Coh, Gth, vstep, Gvth] = getTFpar(stim, resp,pvn,pv,sr)
    % define frequency points being taken into account
        nfp = 1:21;
    
    % obtain threshold bounds
        for m=1:length(pvn)
            eval([ '[' pvn{m} '_l] =' num2str(pv(m,2)) ';']);
            eval([ '[' pvn{m} '_u] =' num2str(pv(m,3)) ';']);
        end
        
	% get frequency response function
        [FD,TD] = getFRF(stim,resp,220,sr,14);
        tf_exp = FD.dec2.FRF(nfp);
        f = FD.dec2.f(nfp);
        Coh = FD.dec2.Coh(nfp);
        
    % obtain lookup tables for threshold handling
        % get Gain list for position threshold substitution
        Gth = makeGth(TD.avg_resp-mean(TD.avg_resp),tg_l,tg_u,sr);
        Gth = Gth(:,[nfp max(nfp)+1]);
        
        % get Gain list for velocity threshold substitution
        Gvth = makeGvth(TD.avg_stim,tfs_l,tfs_u,sr);
        Gvth = Gvth(:,[nfp max(nfp)+1]);
        
    % calculate stimulus velocity
        stim_vel = diff(mean(stim(:,20:end-20),1))*sr;
%         wts = [1/24; repmat(1/12,11,1); 1/24];
%         stim_vel = conv(stim_vel,wts,'valid');
    % calculate step velocity
        vstep = mean([stim_vel(stim_vel>0.1), -stim_vel(stim_vel<-0.1)]);
end

function Gth = makeGth(xin,lb,ub,sr)

    Gth=zeros(101,151);
    n=0;
    for th = lb : (ub-lb)/100 : ub
        n=n+1;
        [xth,~,~] = getSpec(TD_threshold(xin,th)',sr);

        [xnt,~,~] = getSpec(xin,sr);

        tmp = xth./xnt;
        tmp = decimate2(tmp,2);
        Gth(n,2:151) = tmp(1:150);
        Gth(n,1) = th;
    end
end

function Gth = makeGvth(xin,lb,ub,sr)

    Gth=zeros(101,151);
    n=0;
    for th = lb : (ub-lb)/100 : ub
        n=n+1;
        [xth,~,~] = getSpec(TD_vel_threshold(xin,th,sr)',sr);

        [xnt,~,~] = getSpec(xin,sr);

        tmp = xth./xnt;
        tmp = decimate2(tmp,2);
        Gth(n,2:151) = tmp(1:150);
        Gth(n,1) = th;
    end
end

function out = getGth(in,Gth)
    % get last entry that is larger then in
        ind = sum(in>=Gth(:,1));
    
    % estimate values for non-integer indices
        d = (in - Gth(ind,1)) / (Gth(ind+1,1) - Gth(ind,1));
        out =  Gth(ind,2:end) * (1-d)...
             + Gth(ind+1,2:end) * d;

end

% function out = calculate_frf(stim, com, sr)
%     
%     if size(stim,1)>size(stim,2)
%         stim=stim';
%     end
%    if size(com,1)>size(com,2)
%       com = com';
%    end
%     out.stim = stim;
%     out.com = com;
%     
%         [yi,yii,~] = getSpec(stim,sr);
%         [yo,yoo,f] = getSpec(com,sr);
%         yoi = yo.*conj(yi);
% 
%         FRF = yo./yi;
%             out.yi = yi;
%             out.yo = yo;
%             out.f = f;
%             
%             out.f2   = decimate2(f,2);
%             out.FRF2 = decimate2(FRF,2);
% 
%         yoi = mean(decimate2(yoi,2));
%         yii = mean(decimate2(yii,2));
%         yoo = mean(decimate2(yoo,2));
%         
%         out.Coh2=(abs(yoi).^2)./(yii.*yoo);
% end  

