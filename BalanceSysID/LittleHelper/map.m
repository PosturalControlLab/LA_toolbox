function map = map(varargin)

switch computer
    case 'GLNXA64'
        map ='/media/lorenz/E/Cloud/';
        [~,hostname]= system('hostname');
        if strcmp(strcat(hostname),'t460s')
            map='/home/lorenz/';
        end
    case 'PCWIN64'
        map ='E:\Cloud\';
end

if nargin>0
    switch varargin{1}
        case 1; map = [map 'Matlab' filesep];
        case 2; map = [map 'Experimente' filesep];
        case 3; map = [map 'Dokumente' filesep];
        case 'dtp'; map = [map 'Experimente' filesep 'DEC_two_PRTS' filesep];
        case 'dcf'; map = [map 'Experimente' filesep 'strobdcf' filesep];
    end
end