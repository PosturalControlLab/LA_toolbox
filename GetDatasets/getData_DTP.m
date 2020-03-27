function [output1, output2] = getData_DTP(subject,trial,resp)
% info=struct;
n=0; m=0;

addpath(genpath([map 'Experimente' filesep 'DEC_two_PRTS']))
info = info_DTP;

for subj=subject
    m=m+1;
switch resp
    case 'xcom'
        [x,~] = getCOM;
    case 'com'
        [~,x] = getCOM;
    case 'ts'
        resp_num  = 46; %channelnumber of response
    case 'hip'
        resp_num = 49;
    case 'ls'
        resp_num = 47;
    case 'cop'
        resp_num = 35;
    case 'm1'
        resp_num = 3;
    case 'm2'
        resp_num = 9;
    case 'm3'
        resp_num = 15;
    case 'm4'
        resp_num = 21;    
    case 'stim'
        resp_num=38;
    case 'trig'
        resp_num=44;
    case 'pf'
        resp_num=51;
end

        
    x = data(:,resp_num);
    
%     if sum(resp_num==[46 49 3 9 21])
%         x=smooth_spikes(x);
%         x=smooth_opto(x);
%     end

%% select cycles 2-10 and resample to ideal cycle length
    trig = getTrigger_DTP(subj,trial);
    
    x = x(trig(2):trig(end));
    
    % define ideal cycle length
    if strcmp(trial,'warmup')
        lcyc = length(makeStim_DTP(2,1));
    else
        lcyc = length(makeStim_DTP(trial,1));
    end
    ncyc = round(length(x)/lcyc);
	% resample
        output1(m,:) = resample(x,ncyc*lcyc,length(x)); %time drift between reizgen and recording PC
    % cut to cycles
        for i = 0:ncyc-1
            n=n+1;
            output2(n,:) = output1(m,(1:lcyc)+i*lcyc);
        end
    keep n m output1 output2 subj subject info trial resp
end


function [xcom,com] = getCOM
if strcmp(trial,'warmup')
    calib = info(subj).warmup.calib(1);
else
    calib = info(subj).file.calib(trial);
end

OFF = info(subj).calfile.values(calib,1);
A1 = info(subj).calfile.values(calib,2);
A2 = info(subj).calfile.values(calib,3);

if strcmp(trial,'warmup')
    [data,~] = load_record_file_v2(info(subj).warmup.name{1},'','',0);
else
    [data,~] = load_record_file_v2(info(subj).file.name{trial},'','',0);
end
    
x_hip = smooth_spikes(data(:,15));
x_sho = smooth_spikes(data(:,9)); %x_sho=butterworth_filter(x_sho,100,5,[],2,'low');
% x_sho = smooth_opto(x_sho);
% x_pf = data(:,21);

xcom = (OFF + A1 * x_hip + A2 * x_sho)/100;

Spar =  WinterTable(info(subj).subj.weight, info(subj).subj.height);

com = asin(xcom./Spar.h_com).*180./pi;

end
end


function output = getTrigger_DTP(subj,tr)

    info = info_DTP;
    addpath(genpath([map 'Experimente' filesep 'DEC_two_PRTS']))

    if strcmp(tr,'warmup')
        [data,~] = load_record_file_v2(info(subj).warmup(1).name{1},'','',0);
    else
        [data,~] = load_record_file_v2(info(subj).file.name{tr},'','',0);
    end

        trig = data(:,44);
        dtrig = diff(trig(101:end-100));

        % trigger threshold - no output if trigger difference value is below 50
        trig_pos = find(abs(dtrig)>50)+100;

        % silence period - no output if distance between triggers is below 200 samples
        sp = diff([trig_pos; length(data)])>200;

        output = trig_pos(sp);
end

function stim = makeStim_DTP(tr,n)

switch tr
    case 1
        [~,xi]=pseudogen3(4,[0 0 -1 1],[2 0 1 1],30);	% 80 state sequence
        xi=xi/(max(xi)-min(xi))*1;
    case 2
        [~,xi]=pseudogen3(4,[0 0 -1 1],[2 0 1 1],30);	% 80 state sequence
        xi=xi/(max(xi)-min(xi))*2;
        % vel = 0.7407
    case 3
        [~,xi]=pseudogen3(4,[0 0 -1 1],[2 0 1 1],30);	% 80 state sequence
        xi=xi/(max(xi)-min(xi))*4;
    case 4
        [~,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 2 0 2],10);	% 242 state sequence %10
        xi=xi/(max(xi)-min(xi));
        xi=xi*0.3675/0.336*2;
        % vel = 1.287
    case 5
        [~,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 2 0 2],10);	% 242 state sequence %10
        xi=xi/(max(xi)-min(xi));
        xi=xi*0.3675/0.336*1.1512;
        % vel = 0.7407
    case 6
        [~,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 2 0 2],10);	% 242 state sequence %10
        xi=xi/(max(xi)-min(xi));
        xi=xi*0.3675/0.336*1.1512*.5;
    case 7
        [~,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 2 0 2],30);	% 242 state sequence %10
        xi=xi/(max(xi)-min(xi));
        xi=xi/0.1961*0.7407;
end


stim = [];
for k=1:n
        stim = [stim; xi];
end


end


