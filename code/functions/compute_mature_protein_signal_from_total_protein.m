% compute_mature_protein_signal_from_total_protein.m
%
% Euler integrate a given mRNA signal, together with a decay term, to get protein signal.
% If f(t) = trapezoid signal,
%       d(protein)/dt = (translation_rate)*(mRNA) - (decay_rate)*protein
%

function [mature_protein] = compute_mature_protein_signal_from_total_protein(protein,maturation_rate,decay_rate,Tmax,dt)

% time and signal arrays
tvec = 0:dt:Tmax;
mature_protein = zeros(1,numel(tvec));

% loop over time
for s = 2:numel(tvec)
   
    % Euler integrate
    mature_protein(s) = mature_protein(s-1) + dt*(maturation_rate*protein(s-1) - decay_rate*mature_protein(s-1));
    
end





end