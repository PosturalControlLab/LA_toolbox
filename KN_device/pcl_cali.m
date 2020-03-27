function pcl_cali(varargin)
% pcl_cali(FileName)
% 
% No input: script asks for calibration trial and trials to save
% calibration values for later com analysis.
% 
% FileName as input: either cell array or string of files to save
% calibration values for later com analysis. Calibration trial is again
% input via user interface.
% 
% Lorenz Asslaender Aug 2018
% lorenz(a)asslaender.de

[calFile,calPathName] = uigetfile('*.mat','Select Calibration Trial','/home/PosturoLab/Projekte/','MultiSelect','off');
if calFile == 0; return; end

if nargin==1
    FileName = varargin{1};
    PathName = '';
elseif nargin == 0
    [FileName,PathName] = uigetfile('*.mat','Select Trials associated with Calibration Trial',calPathName,'MultiSelect','on');
end
    
load([calPathName calFile],'t')

x_hip = pcl_getData('x_hip', [calPathName calFile]);
x_sho = pcl_getData('x_sho', [calPathName calFile]);
x_cop = pcl_getData('cop', [calPathName calFile]);

%% cut out parts, where sway rods exceeded measurement bounds
valid = x_sho < 14.15;

cut_points = length(x_sho) - sum(valid);
if cut_points > 0; disp(['Sway-rod hit limit. Number of points ignored: ' num2str(cut_points)]); end

x_hip = x_hip(valid);
x_sho = x_sho(valid);
x_cop = x_cop(valid);
t = t(valid);

%% run calibration routine

[OFF, A1, A2] = COM_calibration(t, x_hip',x_sho',x_cop',1);

save([calPathName calFile],'OFF','A1','A2','-append')

if ~iscell(FileName) && ~ischar(FileName) && FileName == 0 
    return
elseif iscell(FileName)
    for k=1:length(FileName)
        save([PathName FileName{k}],'OFF','A1','A2','-append')
    end
else
        save([PathName FileName],'OFF','A1','A2','-append')
end




