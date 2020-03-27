function x = cdiff(x)%CDIFF	Central Difference function.  If X is a vector % [x(1) x(2) ... x(n)], then CDIFF(X) returns a vector % of central differences between every second element% [x(2)-x(1)  (x(3)-x(1))/2  (x(4)-x(2))/2  ... %  (x(n)-x(n-2))/2  x(n)-x(n-1)].  %% For time series, divide result x by deltat for proper% velocity scaling where deltat is time interval between% adjacent points x(i) and x(i+1)%% If X is a matrix, the differences are calculated down% each column.%% The first and last elements in CDIFF are the simple% differences between adjacent elements, and the returned% vector or matrix has the same dimensions as the original.%	J.N. Little 8-30-85	[m,n] = size(x);	if m == 1		y = x(3:n) - x(1:n-2);  x = [x(2)-x(1) y./2 x(n)-x(n-1)];	else		y = x(3:m,:) - x(1:m-2,:);  x = [x(2,:)-x(1,:); y./2; x(m,:)-x(m-1,:)];	end