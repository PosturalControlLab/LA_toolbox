clear all; close all

addpath([map(2) 'DEC_two_PRTS'])
% load subject specific parameters
    info = info_DTP;
% get experimental data
    load('DTP_comData.mat')

tic
for sn = 1:11;
    sp = WinterTable(info(sn+3).subj.weight,info(sn+3).subj.height);
    mgh(sn) = sp.m_B*sp.h_com*9.81/180*pi;
   
    [par_fix,par_var] = ICpar_v1(info(sn+3).subj.height, info(sn+3).subj.weight);
       
    k=0;
    for n=1:6
        k=k+1;
        [par_out(sn,k), tfso{sn,k}, tfeo{sn,k}, f{sn,k}] = ICfit_v1(TD(sn,n).stim, TD(sn,n).resp, par_fix, par_var, 1);
        vstim = abs(cdiff(TD(1,n).avg_stim)./0.01);
        vstep(k) = mean(vstim(vstim>0.2));
    end

end
toc

%%
for k=1:6
    struct2table(par_out(:,k))
    p = cell2mat(struct2cell(par_out(:,k)));
    Kp(k,:)  = p(1,:);
    PD(k,:)  = p(2,:);
    dt(k,:)  = p(3,:);
    W(k,:)   = p(4,:);
    Flp(k,:) = p(5,:);
    Glp(k,:) = p(6,:);
end


%%
for k=1:6
    p = cell2mat(struct2cell(par_out(:,k)));
    ICf(k) = mean(p(1,:).*p(4,:)./mgh);
    ICf_std(k) = std(p(1,:).*p(4,:)./mgh);
end

vstep_ideal = [0.3704    0.7327    1.4655    1.2442    0.7162    0.3703];

errorbar(1./vstep_ideal,ICf,ICf_std,'k.')
axis([0 5 0 1])

parA = [[1,1,1]; 1./vstep(1:3)]'\ICf(1:3)';
parB = [[1,1,1]; 1./vstep(4:6)]'\ICf(4:6)';

return
%% plot changes of PD

for k=1:6
    p = cell2mat(struct2cell(par_out(:,k)));

%     PDrel(k) = mean(p(2,:)./p(1,:)./mgh);
%     PDrel_std(k) = std(p(2,:)./p(1,:)./mgh);
    PD(k,:) = p(2,:);
    PDstd(k) = std(p(2,:));
end
%     errorbar(PD,PDstd,'k.')

    
%% 
figure
for k = 1: length(tfso)
    semilogx(eo(k).f2s,eo(k).Mag2s,'r'); hold on
    semilogx(so(k).f2s,so(k).Mag2s,'k-');
    semilogx(tff{k},abs(tfso{k}),'k--');
    
    if k==1
        legend('exp. TF', 'TF fit',  'SL fit')
        title(['TF and SL fit to experimental data; ' ind{cond}])
    end
end

keep cond

