function [FD,TD] = getFRF(stim,resp,varargin)
% function [FD,TD] = getFRF(stim,resp,sr)
% 
%   Inputs:
%       stim:  Input time series; matrix with each cycle in one row
%       resp:  Output time series; matrix with each cycle in one row
%       
%   Optional:
%       sr:    sampling rate [numeric]
%   
%   Parameters:
%       'FreqPoints' [first selected entry, interval, maxFreq (Hz), type]
%           Define which frequencies are selected in output structure
%           [FirstFreqPoint, SelectionInterval, MaximalFreq(Hz), type]
%           default: [1,2,2], selects base frequency, followed by every
%           second frequency until last frequency below 2 Hz.
%           type: if given, other selection methods can be chosen.
%               default (as above); 
%               type=1: Highest selected fequency point is given as index
%                   Example: [1,2,40,1] selects every second frequency
%                   point from 1 to 40 (i.e. highest frequency is 40 times
%                   the base frequency)
%
%       'SmoothFRF', nfpts
%           nfpts: integer giving the numer of frequencies for output
%           Performs a logarithmic smoothing across frequencies
%       
%       'SmoothPhase', true/false
%           Unwraps phase. Default is true.
%
%       'getOtnesCL', true/false
%           Calculate FRF confidence limits. Default is true.
%           Output: cl_otnes: [1x1 struct]      FRF confidence limits (Otnes & Enochson 1972)
%           
%        'getBstrpCL', nbstrp
%           Calculates FRF and TD confidence limits. Calculation takes a
%           while. nbstr is the number of bstr. samples.
%           nbstr=0: no calculation
%           nbstr=1: calculatin with 400 bootstrap samples
%           nbstr>1: calculation with nbstr bootstrap samples
%           Output: FD. cbMag; cbPha; cbCoh;
%                   TD. cbStim; cbResp
%
%   Frequency domain output
%       FD.
%                  f: [1 x nf double]   frequency vector
%                 yi: [nc x nf double]  input amplitude spectra of individual cycles
%                 yo: [nc x nf double]  output amplitude spectra of individual cycles
%            yi_mean: [1 x nf double]   input amplitude spectrum averaged across cycles
%            yo_mean: [1 x nf double]   output amplitude spectrum averaged across cycles
%                yii: [25x16 double]    input power spectra of individual cycles
%                yoo: [25x16 double]    output power spectra of individual cycles
%           yii_mean: [1 x nf double]   input power spectrum averaged across cycles
%           yoo_mean: [1 x nf double]   output power spectrum averaged across cycles
%           yoi_mean: [1 x nf double]   cross power spectrum averaged across cycles
%                FRF: [1 x nf double]   frequency response function
%                Mag: [1 x nf double]   Magnitude of FRF
%                Pha: [1 x nf double]   Phase of FRF
%                Coh: [1 x nf double]   Coherence
%         TotalPower: [double]          Power of total spectrum
%    RemnantsAllFreq: [1 x nf* double]  Remnants (van der Kooij and Peterka 2005)
%       RemnantPower: [double]          Power of Remnants
%           Remnants: [1 x nf double]   Remnants of selected frequencies
%     RemnantsNotSel: [1 x nf* double]  Remnants of not selected frequencies
%      f_NotSelected: [1 x nf* double]  not selected frequencies
%
%   Time domain output (nc: number of cycles; ns: number of samples):
%       TD.
%       PP_resp_cyc: [nc x 1 double]    peak-to-peak of resp angle for each cycle
%              cl95: [1 x ns double]    95% conficence limits from standard deviation
%          RMS_resp: [double]           rms of average response (deg)
%          RMS_stim: [double]           rms of average stimulus (deg)
%           PP_resp: [double]           peak-to-peak value of average response (deg)
%           PP_stim: [double]           peak-to-peak value of average stimulus (deg)
%              resp: [nc x ns double]   response all cycles
%              stim: [nc x ns double]   stimulus all cycles
%          variance: [1 x ns double]    variance of response across cycles
%          avg_stim: [1 x ns double]    stimulus averaged across cycles
%          avg_resp: [1 x ns double]    response averaged across cycles
%                 t: [1 x ns double]    time vector
%        PowerTotal: [double]           signal power of response
%     PowerPeriodic: [double]           signal power of response, periodic component
%     PowerRemnants: [double]           signal power of response, non-periodic component
% 
%
%
%
% (c) Lorenz Assländer, lorenz@asslaender.de,   23-Sept-2015
%                                               updated Aug-2018

%% handle multiple trial inputs in the form of cells
% this section calls the function k times with each trial as non-cell input
    if iscell(stim)
        for k = 1:length(stim)
            [FD(k),TD(k)] = getFRF(stim{k},resp{k},varargin{:});
        end
        return
    end

%% input and default values

    p = inputParser;
    addRequired(p,'stim');
    addRequired(p,'resp');
    default_sr = 1000;
    addOptional(p,'sr',default_sr,@isnumeric);
    defaultFreqPoints = [1, 2, 2];
    addParameter(p,'FreqPoints',defaultFreqPoints,@isnumeric);   
    defaultSmoothFRF = 'None';
    addParameter(p,'SmoothFRF', defaultSmoothFRF); %number of frequency points
    defaultSmoothPhase = 1;
    addParameter(p,'SmoothPhase', defaultSmoothPhase);
    defaultBstrpCL = 0;
    addParameter(p,'getBstrpCL', defaultBstrpCL, @isnumeric);
    defaultOtnesCL = 'false';
    addParameter(p,'getOtnesCL', defaultOtnesCL, @isboolean);
    
    parse(p,stim,resp,varargin{:})

%% calculate temporal parameters of input
        noc = size(stim,1); % number of cycles      
        cycleLength = size(stim,2);
        samprate = p.Results.sr;
        cycleTime = cycleLength/samprate;
        
        fbase = 1/cycleTime;
        f = fbase:fbase:fbase*cycleLength/2;
        
        FreqPoints = p.Results.FreqPoints;
        if length(FreqPoints)>3 && FreqPoints(4)==1
            selFreq = FreqPoints(1):FreqPoints(2):FreqPoints(3);
        else
            fmax = find(f>FreqPoints(3),1) -1;
            selFreq = FreqPoints(1):FreqPoints(2):fmax;
        end
        
%% time domain describing variables (rms, pp,...)
		RMS_cyc=zeros(1,noc); % preallocate
        for cyc=1:noc
			respc=resp(cyc,:);
			RMS_cyc(cyc)=sqrt(mean((respc-mean(respc)).^2));	     % calculate rms value of resp angle for each cycle
        end % cyc

		TD.PP_resp_cyc=max(resp,[],2)-min(resp,[],2);				     % calc peak-to-peak of resp angle for each cycle
        
            
		avg_stim=mean(stim,1);	% average stimulus
		avg_resp=mean(resp,1);	% average response
		resp_offset=mean(avg_resp);					% mean angle of subject resp
        
        
		TD.cl95=1.96*(std(resp,1))/sqrt(cyc);	% 95% confidence limits about mean response pos
		TD.RMS_resp=rms(avg_resp-resp_offset); % rms of average response (deg)
		TD.RMS_stim=rms(avg_stim-mean(avg_stim)); % rms of average stimulus (deg)
		TD.PP_resp=max(avg_resp)-min(avg_resp);	% peak-to-peak value of average response (deg)
		TD.PP_stim=max(avg_stim)-min(avg_stim);	% peak-to-peak value of average stimulus (deg)
        TD.resp=resp;
        TD.stim=stim;
        TD.variance=var(resp,1);
        TD.avg_stim=avg_stim;
        TD.avg_resp=avg_resp;
        TD.t = (1:length(avg_resp))/samprate;
                
        
        
%% calculation of spectra 
    [yi,yii] = getSpec(stim,samprate);
    [yo,yoo] = getSpec(resp,samprate);
    yoi = yo.*conj(yi);
    yoi = 1./samprate./2.*size(stim,2) .* yoi; % scale cross spectrum by same factor as power spectra are scaled in getSpec
    
    % reduce to selected frequencies
        f  = f(selFreq);
        yi = yi(:,selFreq);
        yo = yo(:,selFreq);
        yii = yii(:,selFreq);
        yoo = yoo(:,selFreq);
        yoi = yoi(:,selFreq);

	% smooth across frequencies
%         if isnumeric(p.Results.SmoothFRF)
%             f   = log_smooth(f,p.Results.SmoothFRF);                % frequency vector
%             yi  = log_smooth(yi,p.Results.SmoothFRF);               % input amplitude spectra of individual cycles
%             yo  = log_smooth(yo,p.Results.SmoothFRF);               % output amplitude spectra of individual cycles
%             yii = yi.*conj(yi) .* 1./samprate./2.*size(stim,2); % recalculate power spectra and scale
%             yoo = yo.*conj(yo) .* 1./samprate./2.*size(stim,2); % recalculate power spectra and scale
%             yoi = yi.*conj(yo) .* 1./samprate./2.*size(stim,2); % recalculate cross power spectra and scale
%         end
%         
    % mean spectra
        yi_mean=mean(yi,1);
		yo_mean=mean(yo,1);
        
        yoi_mean=mean(yoi,1);
        yii_mean=mean(yii,1);
        yoo_mean=mean(yoo,1);
        
        
        
%% Calculate FRF, Magnitude and Phase of FRF, as well as Coherence
    
        % FRF from position data - Pintelon & Schoukens eq 2-17
        FRF=yo_mean./yi_mean;
        Mag=abs(FRF);
        Pha=phase(FRF).*180./pi;
        
        if p.Results.SmoothPhase
            Pha = smooth_phase(Pha);
        end
        
        Coh=(abs(yoi_mean).^2)./(yii_mean.*yoo_mean);
        
    % FD output for all frequency points
        FD.f           = f;        % frequency vector
        FD.yi          = yi;       % input amplitude spectra of individual cycles
        FD.yo          = yo;       % output amplitude spectra of individual cycles
        FD.yi_mean     = yi_mean;  % input amplitude spectrum averaged across cycles
        FD.yo_mean     = yo_mean;  % output amplitude spectrum averaged across cycles
        FD.yii         = yii;      % input power spectra of individual cycles
        FD.yoo         = yoo;      % output power spectra of individual cycles
        FD.yii_mean    = yii_mean; % input power spectrum averaged across cycles
        FD.yoo_mean    = yoo_mean; % output power spectrum averaged across cycles
        FD.yoi_mean    = yoi_mean; % cross power spectrum averaged across cycles
        FD.FRF         = FRF;      % frequency response function
        FD.Mag         = Mag;      % Magnitude of FRF
        FD.Pha         = Pha;      % Phase of FRF
        FD.Coh         = Coh;      % Coherence
        
        if isnumeric(p.Results.SmoothFRF)
            fn = fieldnames(FD);
            for k=1:length(fn)
                FD.(fn{k}) = log_smooth(FD.(fn{k}), p.Results.SmoothFRF);
            end
        end
        
        [FD, TD] = getRemnants(FD, TD, samprate, selFreq);

        if p.Results.getOtnesCL
            FD.cl_otnes = FRF_cl_otnes(Coh, Mag, yoo_mean, yii_mean, ones(1,length(Coh)), noc);
        end
        if p.Results.getBstrpCL~=0
            if p.Results.getBstrpCL==1
                [FD, TD] = getFRFcb(FD,TD,200);
            else
                [FD, TD] = getFRFcb(FD,TD,p.Results.getBstrpCL);
            end
        end
        
end

function cl_otnes = FRF_cl_otnes(Coh, Mag, yoo_mean, yii_mean, smooth_pts, cyc)
% [cl_otnes] = FRF_otnes(Coh, Mag, dof)
% 
% Function gives the 95% confidence limits for Coherence, Gain and Phase
% as calculated in Otnes & Enochson 1972 "Digital Time Series Analysis"
%
% basic assumption is a confidence "circle" around the estimated FRF point
% in the complex plane, which is then transferred into gain and phase
% confidence limits. the circle radius is estimated based on the
% f-distribution.
% 
% Inputs:
%       Coh: Coherence
%       Mag: Gain values
%       yoo_mean: mean Output power spectrum
%       yii_mean: mean Input power spectrum
%       dof: degrees of freedom, calculated as the number of averaged points. 
%               i.e. 2 times the number of cycles times the number of 
%               averaged frequencies.
% 
% Outputs:
%       cl_otnes.Coh_cl95_low
%       cl_otnes.Coh_cl95_high
%       cl_otnes.Pha_cl95
%       cl_otnes.Mag_cl95

% get F-distribution table
        dof = smooth_pts.*cyc.*2;
        Fdist = finv(.95,2,dof-2);
%   without statistics toolbox, use presaved matfile:
%         load('f_distr_095_2_x.mat')
%         Fdist = f_distr_095_2_x(dof-2);

% calculate error Radius
        r = sqrt((((2*(1-Coh)).*Fdist)./(dof-2)).*(yoo_mean./yii_mean));

        cl_otnes.Pha_cl95 = 180/pi*asin(r./Mag);
        cl_otnes.Mag_cl95 = r;
        
% Coherence confidence interval calculations
    z = atanh(sqrt(Coh));
    Za = .95; % 95% confidence level

    b = 1./(dof./2-2);
    sigma1 = sqrt(1./(dof./2-2));

    % correction for Coh<.3
    sig_corr = ones(size(Coh));
    Esigma = 1 - .004.^(1.6.*Coh+.22);
    sig_corr(Coh<.3) = Esigma(Coh<.3);
    sigma2 = sigma1 .*sig_corr;

    cl_otnes.Coh_cl95_low = abs(tanh(z-b-sigma2.*Za).^2-Coh);
    cl_otnes.Coh_cl95_high = abs(tanh(z-b+sigma2.*Za).^2-Coh);

end

function [FD, TD] = getRemnants(FD, TD, sf, selFreq)
    % Expected input: signal (data, sf), where data should be cut into cycles
    % of the according repetitions and f is the record sampling frequency
    % sev: stimulus evoked (periodic component of sway).

    % frequency resolution
        N = size(TD.resp,2);
        df = sf / N;
        dt = 1/sf;
        
        x=zeros(size(TD.resp));
        y=zeros(size(TD.resp,1),N/2);
    % detrend of signal, removes offset
        for i=1:size(TD.resp,1)
            x(i,:)=detrend(TD.resp(i,:),'constant');
        %     x(i,:)=detrend(x(i,:));
            y1 = fft(x(i,:));
            y(i,:) = 2.* y1(1,2:N/2+1);
        end

    % x=data; % without detrending the signal
        Sxx = dt./2./N .* abs(y).^2;
        FD.TotalPower = mean(sum(Sxx,2))*df;

    % mean spectral density; mean across cycles. For calculation of the power of
    % the periodic part, same arguments as above.
        Msev = dt./2./N .* abs(mean(y,1)).^2;
        FD.PeriodicPower = sum(Msev) *df;

  % calculation of remnants
    % gives variation of the signal in the frequency domain in terms of
    % power per frequency squared
        for i=1:size(x,1)
            rx(i,:)= dt./2./N .* abs(y(i,:)-mean(y,1)).^2;
        end

    % mean power of the variation
        FD.RemnantsAllFreq =  sum(rx,1) ./ (size (rx,1) - 1);

    % calculation of the mean remnant power.
        FD.RemnantPower = sum (FD.RemnantsAllFreq) .* df;

    % stimulus and non stimulus evoked frequencies
        FD.Remnants = FD.RemnantsAllFreq(selFreq);
        ind = ones(size(FD.RemnantsAllFreq)); ind(selFreq) = 0; %not selected freq.
        FD.RemnantsNotSelected = FD.RemnantsAllFreq(logical(ind));
        
        f=df:df:length(FD.RemnantsAllFreq)*df;
        FD.Remnants_f_NotSelected = f(logical(ind));

  % calculations in time domain
    % total sway power
        TD.PowerTotal = mean(mean(x.^2,2));

    % Mean signal in time domain
        Mx = mean(x,1);
    % power of periodic component
        TD.PowerPeriodic = mean(Mx.^2);

    % calculation of remnants and mean
        for i=1:size(x,1)
            rT(i,:) = abs(x(i,:) - Mx).^2;
        end
    % power of remnants
        TD.PowerRemnants = sum(mean(rT,2)) / (size(rT,1));


end





