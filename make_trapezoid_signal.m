% make_trapezoid_signal.m
%
% Make a periodic, symmetric trapezoid signal of height=1 to model MS2 traces. 
% Uses Euler integration of a piecewise derivative. See accompanying graphic for
% definition of t_on, t_off, and r. Oscillation period = 2/r + t_on +
% t_off. Signal starts at zero.
%
% 

function [trapezoid_signal] = make_trapezoid_signal(r,t_on,t_off,Tmax,dt)
% time and signal arrays
tvec = 0:dt:Tmax;
trapezoid_signal = zeros(1,numel(tvec));

% loop over time
for s = 2:numel(tvec)
   
    trapezoid_signal_prior = trapezoid_signal(s-1);
    
    % call update function (below)
    trapezoid_signal(s) = update_trapezoid_signal(trapezoid_signal_prior,r,t_on,t_off,dt,tvec(s-1));
    
end



end

function [trapezoid_signal_out] = update_trapezoid_signal(trapezoid_signal_prior,r,t_on,t_off,dt,t)
% updates trapezoid_signal with Euler integration of a piecewise derivative
% (dfdt).
% 
% 
% period = 2/r + t_on + t_off;
% 
% num_cycles_occured = floor(t/period);
% 
% % piecewise definition
% if t < (1/r) + num_cycles_occured*period
%     dfdt = r;
% elseif t < 1/r + t_off + num_cycles_occured*period
%     dfdt = 0;
% elseif t < t_off + 2/r + num_cycles_occured*period
%     dfdt = -r;
% elseif t < period+ num_cycles_occured*period
%     dfdt = 0;
% else
%     disp('error in update_trapezoid_signal: one_cycle_t should be less than period');
%     trapezoid_signal_out = NaN;     % need to pass something to function output
%     return
% end
% 


% For simplicity, define derivative using a shifted time variable "one_cycle_t" such that trapezoid
% starts at one_cycle_t = 0. Do this using t modulo the period.
%
period = 2/r + t_on + t_off;
one_cycle_t = mod(t,period);

% piecewise definition
if one_cycle_t < (1/r) 
    dfdt = r;
elseif one_cycle_t < 1/r + t_off 
    dfdt = 0;
elseif one_cycle_t < 2/r + t_off
    dfdt = -r;
elseif one_cycle_t < period
    dfdt = 0;
else
    disp('error in update_trapezoid_signal: one_cycle_t should be less than period');
    trapezoid_signal_out = NaN;     % need to pass something to function output
    return
end
    
% Euler intergration
trapezoid_signal_out = trapezoid_signal_prior + dt*dfdt;

% manually force trapezoid to be bounded by 1. Otherwise numerical errors
% lead to a drift.
if trapezoid_signal_out > 1
    trapezoid_signal_out = 1;
end




end