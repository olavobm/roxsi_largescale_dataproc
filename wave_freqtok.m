function k = wave_freqtok(frequency, H, fbounds, Hbound)
%% k = WAVE_FREQTOK(frequency, H, fbounds, Hbound)
%
%   inputs
%       -
%       -
%       -
%       -
%
%   outputs
%       -
%
%
%
%
% WAVE_FREQTOK.m does NOT deal with useful calculation of
% having 
%
%
%
%


%% Define parameters

% Acceleration of gravity
g = 9.8;

% From Hz to radians per second
if ~exist('fbounds', 'var')
    fbounds = [(1/5000), (1/0.5)];
end
% % fbounds = 2*pi.*fbounds;   % recent change -- should remove, right?!

% Bound at very shallow water
if ~exist('Hbound', 'var')
    Hbound = 0;
end


%% Dispersion relationship
% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, freq, H) g*k*tanh(k*H) - (2*pi*freq)^2;


%%

%
Npts_H = length(H);
Npts_freq = length(frequency);


%% Check which frequencies are good
% for the calculation of wavenumber

%
lgoodfreq = ~isnan(frequency) & ...
            (frequency >= fbounds(1)) & ...
            (frequency <= fbounds(2));


%% Convert frequency to wavenumber

%
k = NaN(Npts_H, Npts_freq);

%
for i1 = 1:Npts_H

    %
    H_aux = H(i1);

    %
    if (H_aux < Hbound)
        continue
    end

    %
    for i2 = 1:Npts_freq

        %
        if lgoodfreq(i2)

            %
            disp_rel_eval = @(k) disp_rel(k, frequency(i2), H_aux);
    
            %
            try
            k(i1, i2) = fzero(disp_rel_eval, fbounds);
            catch
                keyboard
            end
        end
    end

end



