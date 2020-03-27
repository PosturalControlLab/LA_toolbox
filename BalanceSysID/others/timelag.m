function dt = timelag(data_a, data_b, sf)
% data_a earlier as data_b -> positive sign

d_xc = xcorr(data_a, data_b);

if max(d_xc) > abs(min(d_xc))
    pm_l = (d_xc == max(d_xc));
    pm = find(pm_l);
else
    pm_l = (d_xc == min(d_xc));
    pm = find(pm_l);
end

dt = (length(data_a)-pm)/sf;