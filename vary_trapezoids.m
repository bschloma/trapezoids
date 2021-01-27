% vary_trapezoids.m

% Script for testing how different trapezoids lead to different accumulated
% mRNA signals. Define linear arrays for r, t_on, and t_off values and loop through them
%. Create a 2 tile plot with trapezoid signals on the left, mRNA signals on the
% right.

%% params
% time in hours. fix the period and 2 of the trapezoid params. 
period = 0.5;
r = [5,8,20];
inverse_r = 1./r;
t_on = 0.1.*ones(size(r));
t_off = period - t_on - 2.*inverse_r;
if sum(t_off < 0) > 0
    disp(['t_off is negative!'])
    return
end

Tmax = 4;               % total time to simulation
dt = 0.001;             % time step. Note: instabilities appear to abound. check timestep robustness.
decay_rate = 10;        % decay rate (1/hr)

tvec = 0:dt:Tmax;       % time vector


%% loop over params, compute, and plot
figure; hold on;

% arrays for color scheme. fade magenta to cyan.
reds = linspace(1,0,numel(r));
blues = ones(1,numel(r));
greens = linspace(0,1,numel(r));

% trapezoids
subplot(1,2,1); hold on;

for i = 1:numel(r)
    
    [trapezoid_signal] = make_trapezoid_signal(r(i),t_on(i),t_off(i),Tmax,dt);
    
    this_color = [reds(i) greens(i) blues(i)];
    plot(tvec,trapezoid_signal,'-','linewidth',3,'color',this_color)
    
    
    
end

set(gca,'fontsize',24,'linewidth',4)
xlabel('time (hours)','fontsize',24)
ylabel('MS2','fontsize',24)
axis([0,Tmax,0,1])

% mrna
subplot(1,2,2); hold on;

for i = 1:numel(r)
    
    [trapezoid_signal] = make_trapezoid_signal(r(i),t_on(i),t_off(i),Tmax,dt);
    [mrna] = integrate_trapezoid_signal(trapezoid_signal,decay_rate,Tmax,dt);
    
    this_color = [reds(i) greens(i) blues(i)];
    
    plot(tvec,mrna,'-','linewidth',3,'color',this_color)
  
end

set(gca,'fontsize',24,'linewidth',4)
xlabel('time (hours)','fontsize',24)
ylabel('accumulated mRNA','fontsize',24)

    
