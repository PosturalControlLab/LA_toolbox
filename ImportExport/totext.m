function totext(xpoints,ypoints,filename)
%save datapoints to .txt-file e.g. to make plots with tikZ

    M = [xpoints; ypoints]';
    dlmwrite([filename '.txt'], M, 'delimiter', ' ', 'precision', '%e');
    
end

