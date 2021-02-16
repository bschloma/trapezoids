% make_random_square_signal_gaussian.m
%
% Make a periodic, symmetric trapezoid signal of height=1 to model MS2 traces. 
% Uses Euler integration of a piecewise derivative. 

% MODIFIED with new parameterization. fixed slope param steep, and r
% becomes the height.

function [trapezoid_signal,all_periods,all_durations,all_amplitudes] = make_random_square_signal_gaussian(mean_period,mean_amplitude,mean_duration,Tmax,dt,std_period,std_amplitude,std_duration)


% time and signal arrays
tvec = 0:dt:Tmax;
trapezoid_signal = zeros(1,numel(tvec));

num_cycles = 0;
this_period = mean_period;


this_amplitude = mean_amplitude;
this_duration = mean_duration;
this_period = mean_period;

all_durations = this_duration;
all_periods = this_period;
all_amplitudes = this_amplitude;

time_of_last_cycle = 0;

% loop over time
for s = 2:numel(tvec)
    
    if tvec(s) > time_of_last_cycle + this_period-dt
        
        this_period = mean_period + std_period*randn(1);
        if this_period < dt
            this_period = dt;
        end
         

        this_duration = mean_duration + std_duration*randn(1);
        if this_duration > this_period
            this_duration = this_period - 2*dt;
        elseif this_duration < dt
            this_duration = 2*dt;
        end
        
        
        this_amplitude = mean_amplitude + std_amplitude*randn(1);
        if this_amplitude < 0
            this_amplitude = eps;
        end
        
        
        
        time_of_last_cycle = tvec(s);
        
        all_durations = [all_durations,this_duration];
        all_periods = [all_periods,this_period];
        all_amplitudes = [all_amplitudes,this_amplitude];
    end
       
    % call update function (below)
    trapezoid_signal(s) = update_random_square_signal_gaussian(this_period,this_amplitude,this_duration,tvec(s-1));
    

    
end




end

function [trapezoid_signal_out] = update_random_square_signal_gaussian(this_period,this_amplitude,this_duration,t)
% piecewise-defined square wave
% 

% For simplicity, define derivative using a shifted time variable "one_cycle_t" such that trapezoid
% starts at one_cycle_t = 0. Do this using t modulo the period.
%
one_cycle_t = mod(t,this_period);

% piecewise definition
if one_cycle_t < this_duration 
    trapezoid_signal_out = this_amplitude;
elseif one_cycle_t < this_period
    trapezoid_signal_out = 0;
else
    disp('error in update_trapezoid_signal: one_cycle_t should be less than period');
    trapezoid_signal_out = NaN;     % need to pass something to function output
    return
end
    



end