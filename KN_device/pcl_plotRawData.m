function pcl_plotRawData(varargin)
% pcl_plotRawData(filename,save_boolean)
%
% plots raw and processed data of file 'filename'
% save_boolean can be 
%       0 (not saved)
%       1 (saved as pdf)
%       2 (saved as png)
%
%

    if nargin<2
        sp = 0;
    else
        sp = varargin{2};
    end
    
    filename = varargin{1};

    data = pcl_getData('data',filename);
    com  = pcl_getData('com',filename);
    xcom = pcl_getData('xcom',filename);
    cop  = pcl_getData('cop',filename);
    pf   = pcl_getData('pf',filename);
    t    = pcl_getData('t',filename);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(511); plot(t, pf);   axis([0 max(t) -4 6]);  ylabel('platform (deg)'); title(filename,'Interpreter','None');
    subplot(512); plot(t, com);  axis([0 max(t) 0 10]);  ylabel('com (deg)')
    subplot(513); plot(t, cop);  axis([0 max(t) -5 15]); ylabel('cop (cm)')
    subplot(514); plot(t, xcom); axis([0 max(t) -5 15]); ylabel('com (cm)')
    subplot(515); plot(t, data(4:7,:)); axis([0 max(t) -5 5]); ylabel('SR angle (deg)'); legend('SR1','SR2','SR3','SR4')
    
    switch sp
        case 0
        case 1
            [~,name,~] = fileparts(filename);
            if ~exist([pwd filesep 'plots' filesep 'rawData'],'dir'); mkdir plots rawData; end
            make_pdf(gcf,[pwd filesep 'plots' filesep 'rawData' filesep name],'h')
        case 2
            [~,name,~] = fileparts(filename);
            if ~exist([pwd filesep 'plots' filesep 'rawData'],'dir'); mkdir plots rawData; end
            saveas(gcf,[pwd filesep 'plots' filesep 'rawData' filesep name '.png'])
    end
        
  
    
    
    
    