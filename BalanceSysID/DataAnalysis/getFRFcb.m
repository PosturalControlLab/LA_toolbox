function [FD, TD] = getFRFcb(FD,TD,varargin)
% FD = getFRFcb(TD,FD,nbs,ci)
% TD,FD: output of getFRF script
% nbs: number of bootstraps; default is 400
% ci: confidence interval; default is 0.95

    if nargin>2
        nbs = varargin{1};
    else
        nbs = 400;
    end
    if nargin>3
        ci = varargin{2};
    else
        ci = 0.95;
    end

    for k = 1:length(FD)
        % confidence bounds for data
            out = bootstrp(nbs,@mfrf,FD(k).yi,FD(k).yo);
            bstrG = out(:,1:round(end/3));
            bstrP = out(:,round(end/3)+1:round(end/3*2));
            bstrC = out(:,round(end/3*2)+1:end);

            FD(k).cbMag = bs_getCL(bstrG, ci);
            FD(k).cbPha = bs_getCL(bstrP, ci);
            FD(k).cbCoh = bs_getCL(bstrC, ci);
            
        % confidence bounds for TD data
            out = bootstrp(nbs,@mean,TD(k).stim);
            TD(k).cbStim = bs_getCL(out, ci);
            out = bootstrp(nbs,@mean,TD(k).resp);
            TD(k).cbResp = bs_getCL(out, ci);
    end
    
end


function out = mfrf(yi,yo)

    Myo = mean(yo,1);
    Myi = mean(yi,1);
    FRF = Myo./Myi;

    G = abs(FRF);
    P = smooth_phase(phase(FRF).*180./pi);
    
    yoi = yo.*conj(yi);
    yii = yi.*conj(yi);
    yoo = yo.*conj(yo);
    Coh=(abs(mean(yoi,1)).^2)./(mean(yii,1).*mean(yoo,1));
    
    out = [G, P, Coh];
    
end

function cb = bs_getCL(list,cl)
    % function cb = bs_getCL(list,cl)
    % outputs the confidence bounds based on the bootstrap input table 'list'.
    %
    % bootstraps are in rows (Dim=1)
    %

    ilb = round(size(list,1)*(1-cl)/2);
    if ilb==0; ilb=1; end
    iub = round(size(list,1)*((1-cl)/2+cl));

        list = sort(list,1);
        cb(1,:) = list(ilb,:);
        cb(2,:) = list(iub,:);
end