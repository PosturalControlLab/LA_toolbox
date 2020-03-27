function output = TD_threshold(input, threshold)

if size(input,1) < size(input,2)
    input = input';
end

output=zeros(size(input));
nzp = input >  abs(threshold); % positive values above threshold
nzn = input < -abs(threshold); % negative values above threshold

    output(nzp)=input(nzp)-abs(threshold);
    output(nzn)=input(nzn)+abs(threshold);
    
if size(input,1) < size(input,2)
    output = output';
end