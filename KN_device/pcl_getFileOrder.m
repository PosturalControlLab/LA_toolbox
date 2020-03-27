function fo = pcl_getFileOrder(varargin)

if nargin==1
    FileName = varargin{1};
    PathName = '';
elseif nargin == 0
    [FileName,PathName] = uigetfile('*.xlsx','Select xlsx file containing file order.','/home/PosturoLab/Projekte/','MultiSelect','on');
end
    
out = readcell([PathName filesep FileName]);

for k=1:size(out,2)
    for subj = 2:size(out,1)
        if k==1
            fo(subj-1).SubjNum = out{subj,1};
        else
            fo(subj-1).(out{1,k}) = out{subj,k};
        end
    end
end
