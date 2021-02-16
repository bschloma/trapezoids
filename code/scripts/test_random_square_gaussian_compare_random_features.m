% test_random_square_gaussian_compare_random_features.m

l_assemble = 1;
l_statistics = 1;
period = 0.5;
amplitude = 5;
duration = 0.3;
translation_rate = 4.2*60;
mrna_decay_rate = 10;        % decay rate (1/hr)
protein_decay_rate = 12; 
protein_maturation_delay = 0*5/60;
maturation_rate = (1/10)*60;
Tmax = 6;
dt = 0.00001;             % time step. Note: instabilities appear to abound. check timestep robustness.

tvec = 0:dt:Tmax;       % time vector

std_period = [0.01,0.2,0.4].*period;
std_amplitude = [0.01,0.2,0.4].*amplitude;
std_duration = [0.01,0.2,0.4].*duration;




%% time domain, just compare 3

if l_assemble
    figure; hold on;
    for i = 1:numel(std_amplitude)
        
        
        [square_signal,~,~,~] = make_random_square_signal_gaussian(period,amplitude,duration,Tmax,dt,std_period(i),std_amplitude(i),std_duration(i));
        [mrna] = integrate_trapezoid_signal(square_signal,mrna_decay_rate,Tmax,dt);
        protein = compute_protein_signal_from_mrna(mrna,translation_rate,protein_decay_rate,Tmax,dt);
        mature_protein = compute_mature_protein_signal_from_total_protein(protein,maturation_rate,protein_decay_rate,Tmax,dt);
        
        subplot(1,3,i); hold on;
        plot(tvec,mature_protein + 11,'c','linewidth',3)
        plot(tvec,protein + 3,'m','linewidth',3)
        plot(tvec,square_signal,'k-','linewidth',3)
        %plot(tvec,mrna,'g','linewidth',3)
        set(gca,'fontsize',24,'linewidth',4)
        xlabel('time (hours)','fontsize',24)
        ylabel('signal (AU)','fontsize',24)
        if i==1
            legendcell = {'mature fluorescent protein','total protein','MS2'};
            legend(legendcell,'location','nw','fontsize',16)
            title('deterministic bursting','fontsize',24)
        elseif i==2
            title('weakly random duration and amplitude','fontsize',24)
        elseif i==3
            title('strongly random duration and amplitude','fontsize',24)
        end
        axis([0.9,4.5,0,25])
        
        
        
    end
end

%% statistics
if l_statistics
    % assumes square wave
    Tmax = 100;
    tvec = 0:dt:Tmax;
    period_cell = cell(1,numel(std_amplitude));
    duration_cell = cell(1,numel(std_amplitude));
    amplitude_cell = cell(1,numel(std_amplitude));
    markersize = 80;
    for i = 1:numel(std_amplitude)
        
        [square_signal,period_cell{i},duration_cell{i},amplitude_cell{i}] = make_random_square_signal_gaussian(period,amplitude,duration,Tmax,dt,std_period(i),std_amplitude(i),std_duration(i));

% older code for computing statistics from signal. Now I just output them from the main function.        
%         diff_sig = diff(square_signal);
%         on_times = tvec(diff_sig > 0);
%         off_times = tvec(diff_sig < 0);
%         on_ids = find(diff_sig>0);
%         on_ids = on_ids + 1;
%         amplitudes = square_signal(on_ids);
%         
%         amplitude_cell{i} = amplitudes;
%         period_cell{i} = diff(on_times);
%         
%         durations = zeros(1,numel(on_times));
%         for k = 1:numel(durations)
%             this_off_time = off_times(find(off_times>on_times(k),1));
%             durations(k) = this_off_time - on_times(k);
%         end
%         duration_cell{i} = durations;
        
    end
    
    all_statistics_cell = {period_cell,duration_cell,amplitude_cell};
    
    % plot dot plot
    figure; hold on;
    counter = 0;
    color_cell = {[0,0,1],[1,0,0],[0.5,0.5,0.5]};
    sig_x = 0.1;
    bar_width = 0.5;
    for j = 1:numel(all_statistics_cell)
        
        this_statistic = all_statistics_cell{j};
        for i = 1:numel(this_statistic)
            if i==1
                this_mean = mean(this_statistic{i});
            end
            counter = counter+1;
            these_data_points = this_statistic{i};
            these_data_points = these_data_points./this_mean;
            
            this_x = counter.*ones(size(these_data_points)) + sig_x.*randn(size(these_data_points));
            thiscolor = (i./numel(this_statistic)).*color_cell{j};
            h = scatter(this_x,these_data_points,markersize,thiscolor,'filled');
            alpha(h,0.05);
            
            errorbar(counter,mean(these_data_points),std(these_data_points),'k-','linewidth',3)
            mean_bar_x = linspace(counter - 0.5*bar_width,counter+0.5*bar_width,5);
            mean_bar_y = mean(these_data_points).*ones(size(mean_bar_x));
            plot(mean_bar_x,mean_bar_y,'k-','linewidth',3)
            
        end
        
        counter = counter + 1;
        
        
    end
    
    set(gca,'fontsize',24,'linewidth',4,'xtick',[2,6,10],'xticklabel',{'period','duration','amplitude'})
    axis([0,12,0,3])
    ylabel('normalized burst parameter','fontsize',24)
end
    

    
    

