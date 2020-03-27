function [OFF, A1, A2] = COM_calibration(t, x_hip,x_sho,x_cop,pl)
% function [OFF, A1, A2] = calibrate_com_v1(t, x_hip,x_sho,x_cop,pl)
% 
% Derives calibratin parameters for COM from a two segment model.
% t is the time vector, x_hip and x_sho are hip and shoulder markers.
% x_cop is the corresponding center of pressure shift.
% pl is a boolen for plot ON (1) or OFF (2)
%
% To derive COM data from marker displacements during trials use:
% com = ones(size(x_hip))*OFF + x_hip*A1 + x_sho*A2


par=[ones(size(x_hip)),x_hip,x_sho] \ x_cop;
OFF=par(1); A1=par(2); A2=par(3);

if pl
    figure; 
    subplot(2,1,1); plot(t,[x_sho,x_hip, x_cop]); legend('x shoulder','x hip','x cop')
    subplot(2,1,2); plot(t,[x_cop, [ones(size(x_hip)),x_hip,x_sho]*par]); hold on
    text(5,12,{['A1: ' num2str(A1)]; ['A2: ' num2str(A2)]; ['OFF: ' num2str(OFF)]})
end