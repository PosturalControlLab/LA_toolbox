%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     % Filter %
% 
% function data = butterworth_filter (data,sampling_frequency,low_limit,high_limit,order,type)
% type: 'low'; 'high'; 'band_pass'; 'band_stop'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% als Filter wird ein Butterworth Filter verwendet. Um die Zeitliche
% Verschiebung zu eliminieren wird das Signal noch einmal umgedreht, von
% der anderen Seite gefiltert und zur�ckgedreht.

%Filterfunktion: (Ordnung,Grenzfrequenz/Nyquist-frequenz,type)

function data = butterworth_filter (data,sampling_frequency,low_limit,high_limit,order,type)

switch type
case 'low'        %order is 2* selected order
    [b,a] = butter(order,low_limit/(sampling_frequency/2));
%     Hd = dfilt.df2t(b,a);        
case 'high'       %order is 2* selected order
    [b,a] = butter(order,high_limit/(sampling_frequency/2),'high');
case 'band_pass'  %order is 4* selected order
    [b,a] = butter(order,[low_limit high_limit]/(sampling_frequency/2));
case 'band_stop'  %order is 4* selected order
    [b,a] = butter(order,[low_limit high_limit]/(sampling_frequency/2),'stop');    
end

data = filtfilt(b,a,data);