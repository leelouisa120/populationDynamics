%% Problem Function (rhs.m)
function yp = rhs (time,y)
% Computes the right-hand side of equation given t and P(t)

% Simplified Equation
% yp = -y + 2*y^2;

% Full Equation
a = -1;
B = 2;
g_0 = -0.9;
g_1= -0.0001;

yp = y*(a+(B*y)+(g_0+g_1*time)*y^2);
end


%% Method Function (rk4step.m)
function [time, y] = rk4step (rhs, time_in, y_in, opt)
% runs through numerical method
% Input:
%   rhs = calls the rhs function
%   time_in = current time
%   y_in = current value of dependent variable
%   opt = struct containing the time step

% Initialize Variable
h = opt.dt;

% Runge-Kutta Method
k1 = rhs(time_in, y_in);
k2 = rhs(time_in+(h/2), y_in + h/2*k1);
k3 = rhs(time_in+(h/2), y_in + h/2*k2);
k4 = rhs(time_in + h, y_in + h*k3);
y = y_in + h*((k1/6) + (k2/3) + (k3/3) + (k4/6));
time = time_in + h;
end


%% Driver Function (ode.m)
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

%% Validation of Code (script file)

% Numerical Solution
time_step = [0.0001,0.0002, 0.0005, 0.001, 0.0025, 0.005, 0.01];
error = zeros (length(time_step));
for jj = 1: length(time_step)
    
    opt = struct ('dt', time_step(jj), 'method', @rk4step);
    y_init = 0.6;
    time_interval = [0, 1.79];
    
    
    [time_out, y_out] = ode (@rhs, time_interval, y_init, opt);
    
    p(jj) = plot(time_out,y_out);
    hold on;
    
    % Error
    error(jj) = abs(y_out(end) - P(end));
end

% Analytical Solution
P = zeros(1,length(time_out));

for ii = 1: length(time_out)
    A = ((2*y_init) - 1)/y_init;
    P(ii) = 1/(2-A*exp(1)^(time_out(ii)));
end

AnaGraph = plot(time_out,P);
hold on;

% Graphing Runge-Kutta and Analytical Solutions
legend([p(1), p(2), p(3), p(4), p(5), p(6), p(7), AnaGraph], '0.0001','0.0002', '0.0005', '0.001', '0.0025', '0.005', '0.01', 'Analytical');
title('Runge-Kutta and Analytical Solutions')
xlabel('Time')
ylabel('Population')

% Graphing Error
loglog(time_step, error);
xlabel('Time Step');
ylabel('Error');
title('Error vs. Time Step');
