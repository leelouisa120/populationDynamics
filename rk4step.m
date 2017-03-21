function [time, y] = rk4step (rhs, time_in, y_in, opt)
% runs through numerical method
% Input:
%   rhs = calls the rhs function
%   time_in = current time
%   y_in = current value of dependent variable
%   opt = struct containing the time step

h = opt.dt;

k1 = rhs(time_in, y_in);
k2 = rhs(time_in+(h/2), y_in + h/2*k1);
k3 = rhs(time_in+(h/2), y_in + h/2*k2);
k4 = rhs(time_in + h, y_in + h*k3);
y = y_in + h*((k1/6) + (k2/3) + (k3/3) + (k4/6));
time = time_in + h;