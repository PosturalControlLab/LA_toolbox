function [Sx,Sxx,f]=getSpec(varargin)
% function [Sx,Sxx,f]=getSpec(x,Fs)
% 
% Inputs:
%       x:      input time series (time along 2nd Dimension)
%       Fs:     sampling frequency of time series (pts/sec)
%
% Outputs:
%       Sx:     amplitude spectrum; complex values
%       Sxx:    power spectrum
%       f:      frequencies (Hz)
%
% Sx is scaled such that the amplitude of a sine input is given by abs(Sx)
% Sxx is scaled such that the integrated power density Sxx is equal to the 
% mean power of the time domain input. sum(Sxx*df) = mean(data^2). 
% df is the frequency bandwidth accounted for by each frequency point.
%
% (c) Lorenz Assländer, lorenz@asslaender.de, 23-Sept-2015

x = varargin{1};
Fs = varargin{2};

if nargin > 2
    pb = varargin{3};
else
    pb = 0;
end

N=size(x,2);		% # time series sample points

f=(1:N/2)./N.*Fs;   % frequency points for the output

fk=fft(x,[],2);
y=fk(:,2:N/2+1).*2; % half sided spectrum

% scaling to yield Sx, such that abs(Sx) = A
Sx=1/N.*y;

% scaling to yield Sxx, such that Sxx = mean(data^2)
Sxx = 1./Fs./2./N .* abs(y).^2;

if pb
%     figure
    subplot(2,1,1)
    semilogx(f,abs(Sx)); hold on
    ylabel('amplitude')
    subplot(2,1,2)
    loglog(f,Sxx); hold on
    xlabel('frequency (Hz)')
    ylabel('signal power')
end
