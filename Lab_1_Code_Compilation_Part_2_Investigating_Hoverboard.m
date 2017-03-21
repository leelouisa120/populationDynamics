%% Finding Initial Condition P_0 using Stable Critical Point for t = 0
clear all;

a = -1;
B = 2;
g_0 = -0.9;
% g_0 = -1.1;     Understanding Behavior of Solution
g_1= -0.0001;

% Finding the Roots
P_0 = (-B + (B^2 - (4*g_0*a))^(1/2))/(2*g_0);
P_02 = (-B - (B^2 - (4*g_0*a))^(1/2))/(2*g_0)

% Checking whether critical points are stable
critpoint = B + (2*g_0*P_0);
stablecritpoint2 = B + (2*g_0*P_02);

%% 4.1-Compute for t = [0,950] using time_step = 0.8
time_step = 0.8;

    opt = struct ('dt', time_step, 'method', @rk4step);
    y_init = P_02;
    time_interval = [0, 950];

    
    [time_out, y_out] = ode (@rhs, time_interval, y_init, opt);

% Plotting    
    figure;
    subplot(1,2,1);
    figure1 = plot(time_out,y_out);
    title('RK4 for t = [0,950]');
    xlabel('Time');
    ylabel('Population');
    axis([0 1500 0 1.5]);

  
%% 4.2-Compute for t = [0,1300] using time_step = 0.8
time_step = 0.8;

    opt = struct ('dt', time_step, 'method', @rk4step);
    y_init = P_02;
    time_interval = [0, 1300];
    
    
    [time_out, y_out] = ode (@rhs, time_interval, y_init, opt);
    
% Plotting    
    subplot(1,2,2)
    figure2 = plot(time_out,y_out);
    title('RK4 for t = [0,1300]');
    xlabel('Time');
    ylabel('Population');

%% 4.3- Measure delta_P2 in P(t = 1057) when time_step decreases by factor of 2

% Value for P(t = 1057) with time_step of 0.8
time_step=0.8;
    opt=struct('dt',time_step, 'method', @rk4step);
    y_init=P_02;
    time_interval=[0,1057];
    [~, y_out]=ode(@rhs,time_interval,y_init, opt);
    p_1057=y_out(end);

% Initialization
p_old = p_1057;
p_new = p_old+.01;
time_step = 0.8;
delta_p = 0;
time_step_v=[];
rel_change=[];

% While loop to check for when error is > 0.0001
 while ((abs(p_new-p_old))/p_old) > 0.0001
    p_old=p_new;
    opt=struct('dt',time_step, 'method', @rk4step);
    y_init=P_02;
    time_interval=[0,1057];
    [~, y_out]=ode(@rhs,time_interval,y_init, opt);
    p_new = y_out(end);
    loglog(time_step,(abs(p_new-p_old))/p_old,'bo');
    time_step_v(end+1) = time_step;
    rel_change(end+1) = (abs(p_new-p_old))/p_old;
    xlabel('time step');
    ylabel('relative change in population');
    title('Change in Quigonian Population vs. Time Step');
    hold on;
    % New values for time_step
    time_step=time_step/2;
    delta_p=(abs(p_new-p_old))/p_old;
 end
 
 % Plotting the timestep that produces a relative change that is < 0.0001   
 p_old=p_new;
 opt=struct('dt',time_step, 'method', @rk4step);
 y_init=P_o;
 time_interval=[0,1057];
 [time_out, y_out]=ode(@rhs,time_interval,y_init, opt);
 p_new=y_out(end);
 loglog(time_step,(abs(p_new-p_old))/p_old,'bo');
 
 fprintf('To get a relative change in population of less than 0.0001, time step is: %4.4f\n', time_step);


%% 4.5- Final Solution and Graph
 time_step = [0.1,0.5,2.0,3.0,4.0];

for jj = 1: length(time_step)
   
    opt = struct ('dt', time_step(jj), 'method', @rk4step);
    y_init = P_02;
    time_interval = [0, 1300];
    
    
    [time_out, y_out] = ode (@rhs, time_interval, y_init, opt);
    
    td(ii)=plot(time_out, y_out);
    hold on;
end
 
legend([td(1),td(2),td(3),td(4),td(5)],'0.1','0.5','2.0','3.0','4.0');
title('Time-Dependent Quigonian Population');
xlabel('time');
ylabel('population');

%% Bonus
time_step=[0.05];
opt=struct('dt',time_step, 'method', @rk4step);
y_init=P_02;
time_interval=[0,1058];
[time_out, y_out]=ode(@rhs,time_interval,y_init, opt);
td(1)=plot(time_out, y_out,'-b');
hold on

time_step=[0.05];
opt=struct('dt',time_step, 'method', @rk4step);
y_init=P_02;
time_interval=[0,1000];
[time_out, y_out]=ode(@rhs,time_interval,y_init, opt);
td(2)=plot(time_out, y_out,'-r');
hold on

legend([td(1),td(2)],'t=1057','t=1000');
title('Time-Dependent Quigonian Population');
xlabel('time');
ylabel('population');
hold on