function [data_all, data_cycles, stim_all, stim_cycles] = getData_v3_strobdcf(subjects,trials,resp)

%% get calibrated data of subjects
addpath(genpath([map 'Experimente' filesep 'mat_files']));
addpath(genpath([map 'Experimente' filesep 'subfunctions']));


run([map 'Experimente' filesep 'strobdcf' filesep 'filelists_strobdcf_120910.m']);

switch resp
    case 'com'
        resp_num  = 52; %channelnumber of response
    case 'ts'
        resp_num  = 46;
    case 'tl'
        resp_num = 49;
    case 'ls'
        resp_num = 47;
    case 'cop'
        resp_num = 35;
    case 'head'
        resp_num = 3;
    case 'shoulder'
        resp_num = 9;
    case 'm4'
        resp_num = 21;
    case 'b4'
        resp_num=51;
end


data_all=zeros(36300:2*length(subjects));
stim_all=zeros(36300:2*length(subjects));
data_cycles = zeros(10*length(subjects),6050);
stim_cycles = zeros(10*length(subjects),6050);
%%
count = 0;cy_count=0;
for m = subjects %subject number
  for n = 1:2 %test/retest number

    for i= trials %number of trials
      if ~strcmp(filelist{m,n}{i,1}, 'xxxxx'); %datei fehlt in Liste-> ueberspringen
      if ~exist([filelist{m,n}{i,1} '.mat'], 'file'); filelist{m,n}{i,1}
      else
        count = count+1; %Kontrolle wie viele dateien ausgewertet wurden

          x=load([filelist{m,n}{i,1} '.mat']); % ermoeglicht zugriff auf .data

            data=x.(filelist{m,n}{i,1}).data;
            
            if size(data,2)==51; data = com_korrektur(x, data, filelist{m,n}{i,1}); end 
            com_h(count) = get_com_h(x, filelist{m,n}{i,1});
            
            if resp_num == 9
                h2=x.(filelist{m,n}{i,1}).info.h2;
                data_all(:,count) = asin(data(:,9)./(h2*100)).*180./pi;
            else
                data_all(:,count) = data(:,resp_num);
            end
            stim_all(:,count) = data(:,38);
            data_all_files{count,1} = [filelist{m,n}{i,1} ' ' filelist{m,n}{i,2}];
            
            
            for cy = 1:5 % cycles 2-6 are selected!
                
                cy_count=cy_count+1;
%                 data_cycles(cy_count,:) = data(cy*6050+1:(cy+1)*6050,resp_num);
                data_cycles(cy_count,:) = data_all(cy*6050+1:(cy+1)*6050,count);
                stim_cycles(cy_count,:) = data(cy*6050+1:(cy+1)*6050,38);
            end
            
           keep filelist stim_all stim_cycles data_all data_all_files resp_num m n i count subjects trials cy_count data_cycles com_h
      
      end
      end
    end
  end
end

data_all=data_all';
end

function out = get_com_h(Matrix, filename_b)
    sh = getfield(Matrix, filename_b, 'info', 'h2') *100;
    hh = getfield(Matrix, filename_b, 'info', 'h3') *100;
    diff23 = sh-hh;
    
    h_COM_L = 0.553*hh;
    h_COM_T = hh+0.626*(diff23+5); %Glenohumeral Gelenk liegt 5cm �ber Schultermarker

    h_COM_B = h_COM_L*0.32+h_COM_T*0.678;
    out=h_COM_B/100;
end


% count