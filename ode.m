function [time_out, y_out] = ode (rhs, time_interval, y_init, opt)
% collect data as program runs
% time_interval = array = [t_init, t_final]

y_out = zeros(1, length(time_interval));
time_out = time_interval(1) : opt.dt : time_interval(2);
y_out(1) = y_init;


for ii = 1: length(time_out)-1
    [~, y_out(ii+1)] = rk4step(rhs, time_out(ii), y_out(ii), opt);
end

if time_out(end)< time_interval(2)
    opt.dt = time_interval(2)-time_out(end);
    time_out(end+1) = time_interval(2);
    [~, y_out(end+1)] = rk4step(rhs, time_out(end), y_out(end), opt);
end
end