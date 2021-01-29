% test_predicted_protein.m

period = 0.5;
r = 5;
t_on = 0.2;
t_off = period - t_on;
translation_rate = 4.2*60;
mrna_decay_rate = 10;        % decay rate (1/hr)
protein_decay_rate = 12;%0.1*20*mrna_decay_rate;
protein_maturation_delay = 0*5/60;
maturation_rate = (1/10)*60;
slope = 30;

Tmax = 6;               % total time to simulation
dt = 0.001;             % time step. Note: instabilities appear to abound. check timestep robustness.

tvec = 0:dt:Tmax;       % time vector

[trapezoid_signal] = make_trapezoid_signal(r,t_on,t_off,Tmax,dt,slope);
[mrna] = integrate_trapezoid_signal(trapezoid_signal,mrna_decay_rate,Tmax,dt);
protein = compute_protein_signal_from_mrna(mrna,translation_rate,protein_decay_rate,Tmax,dt);
mature_protein = compute_mature_protein_signal_from_total_protein(protein,maturation_rate,protein_decay_rate,Tmax,dt);

figure; hold on;
plot(tvec,mature_protein + 8,'c','linewidth',3)
plot(tvec,protein - 0,'m','linewidth',3)
plot(tvec,trapezoid_signal./max(trapezoid_signal).*2,'k-','linewidth',3)
%plot(tvec,mrna./max(mrna),'m','linewidth',3)
set(gca,'fontsize',24,'linewidth',4)
xlabel('time (hours)','fontsize',24)
ylabel('signal (AU)','fontsize',24)
legendcell = {'mature fluorescent protein','total protein','MS2'};
legend(legendcell,'location','nw','fontsize',16)
axis([0.9,2.5,0,18])