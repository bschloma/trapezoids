% compute_protein_signal_from_mrna.m
%
% Euler integrate a given mRNA signal, together with a decay term, to get protein signal.
% If f(t) = trapezoid signal,
%       d(protein)/dt = (translation_rate)*(mRNA) - (decay_rate)*protein
%

function [protein] = compute_protein_signal_from_mrna(mRNA,translation_rate,decay_rate,Tmax,dt)

% time and signal arrays
tvec = 0:dt:Tmax;
protein = zeros(1,numel(tvec));

% loop over time
for s = 2:numel(tvec)
   
    % Euler integrate
    protein(s) = protein(s-1) + dt*(translation_rate*mRNA(s-1) - decay_rate*protein(s-1));
    
end





end