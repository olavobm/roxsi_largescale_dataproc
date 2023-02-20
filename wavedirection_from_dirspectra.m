function [meandir_f, meandirspread_f] = wavedirection_from_dirspectra(theta, See_theta, D_theta)
%% [meandir_f, meandirspread_f] = WAVEDIRECTION_FROM_DIRSPECTRA(theta, See_theta, D_theta)
%
%   inputs
%       -
%       -
%       -
%
%   outputs
%       -
%       -
%
%
%
%
%
%


%%

%
Nfrequency= size(See_theta, 1);
Ntimestamps = size(See_theta, 3);


%%
% Matrices where row dimension should be frequency
% and column dimension should be time

% %
% meandir_f = NaN(Nfrequency, Ntimestamps);
% bulkstats.meandirspread_f = bulkstats.meandir_f;

% They correspond theta0 and sigma0 in usual notation
% and I'll use this to make the code clearer
theta0 = NaN(Nfrequency, Ntimestamps);
sigma0 = theta0;


%%

% % % What are these???
% % theta0 = bulkstats.meandir_f;
% % sigma0 = theta0;

% the mean direction (and spread) that is (are) used to compute a more refined estimate(s)???



%%

% Loop over time
for i1 = 1:Ntimestamps

    % Loop over frequencies
    for i2 = 1:Nfrequency

        % Create a circular version of the directional spread (???)
        D = D_theta(i2, :, i1);
        D = D(:);
        
        % Get the maximum of the directional distribution function
        [MAXD, ind_max] = max(D);
        %
        guess = theta(ind_max);

        %
        if ~isnan(MAXD)    % WHEN IS IND_MAX NAN ???
            
            % Use fzero to solve for theta0 (per timestep, per frequency)
            theta0(i2, i1) = fzero(@(theta0) trapz(theta, sin(theta - theta0).*D), guess);
            
            % Compute directional spread -- why this expression???
            sigma0(i2, i1) = rad2deg(sqrt(4.*trapz(theta, ...
                                                   (sin((theta - theta0(i2, i1))./2).^2).*D)));
% % %         else % isn't this unnecessayr??
% % %             theta0(i1, i2) = NaN;
% % %             sigma0(i1, i2) = NaN;
        end
        
    end
end

%
theta0 = mod(theta0, 2*pi);    % ???


%%

%
meandir_f = theta0;
meandirspread_f = sigma0;




