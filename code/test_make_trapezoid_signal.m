% test_make_trapezoid_signal.m

% time in hours
r = 4;
t_off = 0.4;
t_on = 0.4;
Tmax = 4;
dt = 0.0001;
tvec = 0:dt:Tmax;
decay_rate = 10;

[trapezoid_signal] = make_trapezoid_signal(r,t_on,t_off,Tmax,dt);
[mrna] = integrate_trapezoid_signal(trapezoid_signal,decay_rate,Tmax,dt);

figure; hold on;
plot(tvec,trapezoid_signal,'b-','linewidth',3)
plot(tvec,mrna./max(mrna),'r','linewidth',3)
set(gca,'fontsize',24,'linewidth',4)
xlabel('time (hours)','fontsize',24)
ylabel('signal','fontsize',24)