function [out, out_cyc] = pcl_getData(varargin)
% out = pcl_loadData(resp,fname);
%
% ToDo: implement cut to cycle
%       create file to add CycleLenth post hoc to trials

addpath(genpath('/home/lorenz/Matlab/LA_toolbox'))
%% handle variable input
    if nargin<2
        [FileName,PathName] = uigetfile('*.mat','Select Trial(s)','/home/PosturoLab/Experimente/','MultiSelect','on');
        if ~iscell(FileName) && ~ischar(FileName) && FileName == 0; return; end
        if iscell(FileName)
            for k=1:length(FileName)
            fname{k} = [PathName FileName{k}];
            end
        else
            fname = [PathName FileName];
        end
    else
        fname=varargin{2};
    end

    if nargin<1
        resp='com';
    else
        resp=varargin{1};
    end


%% loop to handle multiple files
if iscell(fname); k_max=length(fname); else k_max=1; end
fname_l = fname;
out=[];

for k=1:k_max
    if iscell(fname_l)
        fname=fname_l{k};
    end
    
    load(fname,'data');
    load(fname,'SubjInfo');


%% select appropriate output
switch resp
    case 'data'
        out_t = data;
    case 't'
        load(fname,'t')
        out_t = t';
    case 'xcom'
        [~,out_t] = getCOM(fname);
    case 'com'
        [out_t,~] = getCOM(fname);
    case 'x_hip'
        out_t = SubjInfo.hip_hook_distance * tand(data(:,4));
    case 'x_sho'
        out_t = SubjInfo.shoulder_hook_distance * tand(data(:,5));
    case 'hip'
        msgbox('not implemented yet'); return
    case 'ls'
        msgbox('not implemented yet'); return
    case 'ts'
        msgbox('not implemented yet'); return
    case 'torque'
        out_t = data(:,8) * 36.9;
    case 'cop'
        out_t = data(:,8) * 36.9 / (SubjInfo.body_weight*9.81) * 100;
    case 'sr1'
        out_t = data(:,4);
    case 'sr2'
        out_t = data(:,5);
    case 'sr3'
        out_t = data(:,6);
    case 'sr4'
        out_t = data(:,7);
    case 'stim'
        load(fname,'stim');
        out_t = stim;
        if isrow(out_t); out_t = out_t'; end
    case 'pf'
        out_t = data(:,3);
end
    dl(k) = length(out_t);
    out{k} = out_t';    

    try
        load(fname,'CycleLength')
        cl(k) = CycleLength;
    catch
        cl(k)=0;
    end
end

%% format output and cut to cycles
    out_cyc = []; make_cyc = 1;
    if ~all(cl==cl(1)); disp('Selected files do not have the same cycle length.'); make_cyc = 0; end
    if ~all(cl~=0); disp('At least one selected file has no defined cycle length.'); make_cyc = 0; end
    if strcmp(resp,'data'); make_cyc = 0; end    
    
    if make_cyc
    for k=1:length(out)
        out_cyc_t = reshape(out{k}, CycleLength, size(out{k},2)/CycleLength);
        out_cyc = [out_cyc; out_cyc_t'];
    end
    end
    
    % change out from cell to matrix if all trials have the same length
    if all(dl==dl(1)); out = cell2mat(out'); end

    
end


function [a_com,x_com] = getCOM(fname)
    
        load(fname,'OFF','A1','A2')
    if ~exist('OFF','var')
        answer = questdlg('No calibration values found. Select calibration file?');
        if strcmp(answer,'Yes')
            [OFF, A1, A2] = pcl_cali(fname);
        else
            x_com = [];
            a_com = [];
            return
        end
    end
    
    x_hip = pcl_getData('x_hip',fname);
    x_sho = pcl_getData('x_sho',fname);
    load(fname,'sr');
    load(fname,'SubjInfo');

    x_com = OFF + A1*x_hip + A2*x_sho;
    ap = WinterTable(SubjInfo.body_weight,SubjInfo.body_height/100);
    a_com = asind(x_com / (ap.h_com*100)) ;
    
    x_com = x_com'; a_com = a_com';
end