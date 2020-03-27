function boxPlot(data0, lineWidth, width)
% boxPlot(data0) - plot box-whiskers diagram, accept multiple columns
% Arguments: data0 -  unsorted data, mxn, m samples, n columns
%            lineWidth -  line thickness in the plot default = 1;
%            width -  the width of the box, default = 1;
% Returns:	 
% Notes: each column is considered as a single set	


    if(nargin < 3)
        width = 1;
    end;
    if(nargin < 2)
        lineWidth = 1;
    end;


    [m n] = size(data0);

    data = sort(data0, 1); % ascend
    
    
    q2 = median(data, 1);
    
    if(rem(m,2) == 0)
        
        upperA = data(1:m/2,:);
        lowA =  data(m/2+1:end,:);
        
    else
        
        upperA = data(1:round(m/2), :);
        lowA =  data(round(m/2):end, :);  
        
    end;
    
    q1 = median(upperA, 1);
    q3 = median(lowA, 1);
    
    min_v = data(1,:);
    max_v = data(end,:);
    
    draw_data = [max_v; q3; q2; q1; min_v];
    
    % adjust the width
    drawBox(draw_data, lineWidth, width);


return;


function drawBox(draw_data, lineWidth, width)

    n = size(draw_data, 2);

    unit = (1-1/(1+n))/(1+9/(width+3));
    
        
    hold on;       
    for i = 1:n
        if i==1 || i==7; lineWidth=2; else lineWidth=1.5; end
        v = draw_data(:,i);
        
        % draw the min line
%         plot([i-unit, i+unit], [v(5), v(5)], 'LineWidth', lineWidth, 'Color', 'k');
        % draw the max line
%         plot([i-unit, i+unit], [v(1), v(1)], 'LineWidth', lineWidth, 'Color', 'k');
        % draw middle line
        plot([i-unit, i+unit], [v(3), v(3)], 'LineWidth', 4, 'Color', 'k');
        % draw vertical line
        plot([i, i], [v(5), v(4)], 'LineWidth', .5, 'Color', 'k');
        plot([i, i], [v(2), v(1)], 'LineWidth', .5, 'Color', 'k');
        % draw box
        plot([i-unit, i+unit, i+unit, i-unit, i-unit], [v(2), v(2), v(4), v(4), v(2)], 'LineWidth', lineWidth, 'Color', 'k');
        
    end;

return;