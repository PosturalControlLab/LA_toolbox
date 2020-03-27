function plotFRF(FD,TD,varargin)
% plotFRF(FD,TD,leg,fname)
% plots results analysed using the script 'getFRF'
% if multiple structures in TD and FD are given as input, all will be
% plotted in the output figure
%
% TD/FD: output structures from function getFRF
% leg: legend entries expected input is a cell array
% fname: filename; if input is provided, plot is stored as pdf in current
% folder

%% handle variable input
    p = inputParser;
    addRequired(p,'FD');
    addRequired(p,'TD');
    defaultType = 0;
    addOptional(p,'type',defaultType,@isnumeric);
    
    defaultColors={'k';'b';'c';'r';'g';'y';'k';'b';'c';'r';'g';'y'};
    addParameter(p,'Colors',defaultColors,@iscell);
    defaultLegend = 0;
    addParameter(p,'Legend',defaultLegend,@iscell);
    defaultFilename = 0;
    addParameter(p,'Filename',defaultFilename,@ischar);
    defaultAxis = 0;
    addParameter(p,'Axis',defaultAxis,@isnumeric);

    parse(p,FD,TD,varargin{:})

tmax=0; amax=0; pmin=0; pmax=0;

% figure
pc = p.Results.Colors;
type = p.Results.type;
for k=1:length(TD)
    subplot(3,2,1)
        switch type
            case 0; plot(TD(k).t, TD(k).avg_stim, pc{k}); 
            case 1; 
                if ~isfield(FD,'cbMag')
                    [FD,TD] = getFRFcb(FD,TD); % calculate FD confidence bounds
                end
                plot(TD(k).t, TD(k).avg_stim, pc{k},'LineWidth',2); hold on
                plot(TD(k).t, TD(k).cbStim(1,:), pc{k});
                plot(TD(k).t, TD(k).cbStim(2,:), pc{k});
            case 2
                if ~isfield(FD,'cbMag')
                    [FD,TD] = getFRFcb(FD,TD); % calculate FD confidence bounds
                end
                shadedErrorBar(TD(k).t, TD(k).avg_stim, abs(TD(k).cbStim-repmat(TD(k).avg_stim,2,1)), 'lineprops', pc{k}); hold on
            case 3
                plot(TD(k).t, TD(k).avg_stim, pc{k}, 'LineWidth', 3); hold on
                plot(TD(k).t, TD(k).stim, pc{k}, 'LineWidth', 0.1);
        end
        hold on
        xlabel('time (s)')
        ylabel('stim')
        xlim([0 TD(k).t(end)])
        
    subplot(3,2,2)
        switch type
            case 0; plot(TD(k).t, TD(k).avg_resp, pc{k});
            case 1 
                plot(TD(k).t, TD(k).avg_resp, pc{k},'LineWidth',2); hold on
                plot(TD(k).t, TD(k).cbResp(1,:), pc{k});
                plot(TD(k).t, TD(k).cbResp(2,:), pc{k});
            case 2
                shadedErrorBar(TD(k).t, TD(k).avg_resp, abs(TD(k).cbResp-repmat(TD(k).avg_resp,2,1)), 'lineprops', pc{k}); hold on
            case 3
                plot(TD(k).t, TD(k).avg_resp, pc{k}, 'LineWidth', 2); hold on
                plot(TD(k).t, TD(k).resp, pc{k}, 'LineWidth', 0.1);
        end
        hold on
        xlabel('time (s)')
        ylabel('resp')
        xlim([0 TD(k).t(end)])
        
    subplot(3,2,3)
        switch type
            case 0; semilogx(FD(k).f, FD(k).Mag, pc{k}, 'Marker','.');
            case 1 
                semilogx(FD(k).f, FD(k).Mag, pc{k},'LineWidth',2); hold on
                semilogx(FD(k).f, FD(k).cbMag(1,:), pc{k});
                semilogx(FD(k).f, FD(k).cbMag(2,:), pc{k});
            case 2
                shadedErrorBar(FD(k).f, FD(k).Mag, abs(FD(k).cbMag([2 1],:)-repmat(FD(k).Mag,2,1)), 'lineprops', pc{k}); hold on
                set(gca,'XScale','log'); 
            case 3
                semilogx(FD(k).f, FD(k).Mag, pc{k}); hold on
        end
        hold on
        xlabel('frequency (Hz)')
        ylabel('gain (deg/deg)')
    subplot(3,2,4)
        switch type
            case 0; semilogx(FD(k).f, FD(k).Coh, pc{k}, 'Marker','.');
            case 1 
                semilogx(FD(k).f, FD(k).Coh, pc{k},'LineWidth',2); hold on
                semilogx(FD(k).f, FD(k).cbCoh(1,:), pc{k});
                semilogx(FD(k).f, FD(k).cbCoh(2,:), pc{k});
            case 2
                shadedErrorBar(FD(k).f, FD(k).Coh, abs(FD(k).cbCoh([2 1],:)-repmat(FD(k).Coh,2,1)), 'lineprops', pc{k}); hold on
                set(gca,'XScale','log'); 
            case 3
                semilogx(FD(k).f, FD(k).Coh, pc{k}); hold on
        end
        hold on
        xlabel('frequency (Hz)')
        ylabel('coherence')
    subplot(3,2,5)
        switch type
            case 0; semilogx(FD(k).f, FD(k).Pha, pc{k}, 'Marker','.');
            case 1 
                semilogx(FD(k).f, FD(k).Pha, pc{k},'LineWidth',2); hold on
                semilogx(FD(k).f, FD(k).cbPha(1,:), pc{k});
                semilogx(FD(k).f, FD(k).cbPha(2,:), pc{k});
            case 2
                shadedErrorBar(FD(k).f, FD(k).Pha, abs(FD(k).cbPha([2 1],:)-repmat(FD(k).Pha,2,1)), 'lineprops', pc{k}); hold on
                set(gca,'XScale','log'); 
            case 3
                semilogx(FD(k).f, FD(k).Pha, pc{k}); hold on
        end
        hold on
        xlabel('frequency (Hz)')
        ylabel('phase (deg)')

    tmax = max([tmax, TD(k).t(end)]);
    amax = max([amax, max(FD(k).Mag)]);
    pmin = min([pmin, min(FD(k).Pha)]);
%     pmax = max([pmax, max(FD(k).Pha(1:4))]);
    pmax = max([pmax, max(FD(k).Pha)]);
end

if iscell(p.Results.Legend)
    subplot(3,2,4)
%     plot(10,10); hold on
%     axis([0 10 0 1])
    legend(p.Results.Legend,'Position',[0.8 0.2 0.05 0.05])
end
        

    
if length(p.Results.Axis)==1 && p.Results.Axis==0
    fmin = 1/tmax;
    fmax = 2;
    amin = 0;
    amax = ceil(amax*2)/2;
elseif length(p.Results.Axis)==1 && p.Results.Axis==1
    fmin = 0.0165;
    fmax = 2;
    amin = 0;
    amax = 5;
else
    fmin = p.Results.Axis(1);
    fmax = p.Results.Axis(2);
    amin = p.Results.Axis(3);
    amax = p.Results.Axis(4);
end

    pmin = -300;
    pmax = 100;
%     pmin = floor(pmin/50)*50;
%     pmax = ceil(pmax/50)*50;


    subplot(3,2,1); axis([0 tmax -2 4])
    subplot(3,2,2); axis([0 tmax -4 4])
    subplot(3,2,3); axis([fmin fmax amin amax])
    subplot(3,2,4); axis([fmin fmax 0 1])
    subplot(3,2,5); axis([fmin fmax pmin pmax])

if ischar(p.Results.Filename)
%     fname = varargin{4};
%     saveas(gcf, [fname '.eps'], 'psc2')
%     saveas(gcf, [fname '.png'])
    make_pdf(gcf, p.Results.Filename)
end

% 
% for k=1:length(TD)
%     subplot(3,2,2)
%         errBar = abs(TD(k).cb - [TD(k).avg_resp; TD(k).avg_resp]);
%         shadedErrorBar(TD(k).t, TD(k).avg_resp, errBar, pc{k}); hold on
%         
%     subplot(3,2,3)
%         errBar = abs(FD(k).cbMag - [FD(k).Mag;FD(k).Mag]);
%         shadedErrorBar(FD(k).f, FD(k).Mag, errBar, pc{k}); hold on
%     
%     subplot(3,2,4)
%         errBar = [FD(k).cl_otnes.Coh_cl95_high; FD(k).cl_otnes.Coh_cl95_high];
%         shadedErrorBar(FD(k).f, FD(k).Coh, errBar, pc{k}); hold on
%     
%     subplot(3,2,5)
%         errBar = abs(FD(k).cbPha - [FD(k).Pha;FD(k).Pha]);
%         shadedErrorBar(FD(k).f, smooth_phase(FD(k).Pha), errBar, pc{k}); hold on
%     
%     subplot(3,2,6)
%         plot(10,10,pc{k}); hold on
% 
% end
    