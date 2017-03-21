clear all;
 %% Numerical Solution
    time_step = [0.0001,0.0002, 0.0005, 0.001, 0.0025, 0.005, 0.01];
    error = zeros (length(time_step));
for jj = 1: length(time_step)
   
    opt = struct ('dt', time_step(jj), 'method', @rk4step);
    y_init = 0.6;
    time_interval = [0, 1.79];
    
    
    [time_out, y_out] = ode (@rhs, time_interval, y_init, opt);
    
    p(jj) = plot(time_out,y_out);
    hold on;


    %% Error
        %error(jj) = abs(y_out(end) - P(end));
end
%% Analytical Solution
    P = zeros(1,length(time_out));
    
    for ii = 1: length(time_out)
        A = ((2*y_init) - 1)/y_init;
        P(ii) = 1/(2-A*exp(1)^(time_out(ii)));
    end
    
    AnaGraph = plot(time_out,P);
    hold on;

    legend([p(1), p(2), p(3), p(4), p(5), p(6), p(7), AnaGraph], '0.0001','0.0002', '0.0005', '0.001', '0.0025', '0.005', '0.01', 'Analytical');
    title('Runge-Kutta and Analytical Solutions')
    xlabel('Time')
    ylabel('Population')
%  loglog(time_step, error);
%  xlabel('Time Step');
%  ylabel('Error');
%  title('Error vs. Time Step');
