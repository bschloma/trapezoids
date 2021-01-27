% integrate_trapezoid_signal.m
%
% Euler integrate a given trapezoid (MS2) signal, representing mRNA
% production, together with a decay term, to get accumulated mRNA.
% If f(t) = trapezoid signal,
%       d(mRNA)/dt = f(t) - (decay_rate)*mRNA
%

function [mrna] = integrate_trapezoid_signal(trapezoid_signal,decay_rate,Tmax,dt)

% time and signal arrays
tvec = 0:dt:Tmax;
mrna = zeros(1,numel(tvec));

% inherit initial condition from trapezoid (usually 0).
mrna(1) = trapezoid_signal(1);

% loop over time
for s = 2:numel(tvec)
   
    % Euler integrate
    mrna(s) = mrna(s-1) + dt*(trapezoid_signal(s-1) - decay_rate*mrna(s-1));
    
end





end