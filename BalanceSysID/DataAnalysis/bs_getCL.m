function cb = bs_getCL(varargin)
% function cb = bs_getCL(list,cl)
%
% outputs the confidence bounds based on the bootstrap input table 'list'.
%
% confidence limit default is 95%, i.e. cl=0.95, therefore 0.025 and 0.975
% values of the sorted list.
% bootstraps are in rows (Dim=1)
%

if nargin<2
    cl=0.95;
else
    cl=varargin{2};
end
list = varargin{1};

ilb = round(size(list,1)*(1-cl)/2);
if ilb==0; ilb=1; end
iub = round(size(list,1)*((1-cl)/2+cl));

    list = sort(list,1);
    cb(1,:) = list(ilb,:);
    cb(2,:) = list(iub,:);





