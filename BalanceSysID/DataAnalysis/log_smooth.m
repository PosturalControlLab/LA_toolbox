function [xout,npoints] = log_smooth(xin,n)
% function [xout,npoints] = log_smooth(xin,n)
%
% performs an averaging across frequencies with approximately logarithmic
% bin widths. Averaging is performed along Dim2.
%
% xin: linearly scaled values
% n: number of desired output points - if length(xin) is too short, the
%    number of output points will be smaller then n
%
% xout: vector of output values. Each output point is calculated as the 
%       mean across all points in a bin.
% npoints: number of averaged points of each output value xout.

n2=n-1; npoints=0;
while length(npoints)<n
    n2=n2+1;
    ls = logspace(log10(1),log10(size(xin,2)),n2);
    npoints = [1 diff(round(ls))]; % get bin width
    npoints=npoints(npoints~=0); % kick out frequencies which are double
end



% find first point where bins contain more the one data point, i.e. next
% bin contains averaging
% n2 = find(npoints~=1, 1)-1;
% ls = logspace(log10(n2),log10(size(xin,2)), n-n2);
% npoints = [ones(1,n2) diff(round(ls))]; % get bin width


bins = cumsum(npoints);

xout = zeros(size(xin,1),length(bins));
for k = 1:length(bins)        
    if npoints(k)==1
        xout(:,k)=xin(:,bins(k));
    else
        xout(:,k) = mean(xin(:,bins(k-1)+1:bins(k)),2);
    end
    
end
