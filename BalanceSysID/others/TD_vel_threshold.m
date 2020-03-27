function output = TD_vel_threshold(input, threshold,sr)

if size(input,1) < size(input,2)
    input = input';
end

input=cdiff(input).*sr;

output=zeros(size(input));
nzp = input >  abs(threshold); % positive values above threshold
nzn = input < -abs(threshold); % negative values above threshold

    output(nzp)=input(nzp)-abs(threshold);
    output(nzn)=input(nzn)+abs(threshold);
    
output=cumsum(output)./sr;
    
if size(input,1) < size(input,2)
    output = output';
end