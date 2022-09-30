function [bottomdepth] = waveprestoH(presbottom, dtsampling, dtavg, freqcutoff)
%% function [eta, Hmo, Tp] = WAVEPRESTOH(presbottom, dtsampling, dtavg, freqcutoff)
%
%   inputs
%       - presbottom: timeseries (vector) of bottom pressure in dbar.
%       - dtavg: in seconds.
%       - dtavg: in seconds.
%       - freqcutoff: 
%
%   outputs
%       - bottomdepth: timeseries of bottom depth, computed from
%                      pressure data and linear wave theory.
%
%
% eta is sea surface elevation, Hmo=significant wave hieght, Tp=peak period
% computes wave height from pressure measurements
% p=pressure(decibars or meters); hi =instrument height (meters)
% dt= sampling period (seconds)
%
% With pressure data at 2 Hz, this code takes about 7 seconds per week
% of data to run (in my 2015 laptop, 3.1 GHz dual-core Intel Core i7,
% 8GB of memory). This running time is linear (i.e. 4 weeks of data at
% 2 Hz should take 28 seconds to run).
%
% Olavo Marques (Sep/2022), adapted (and highly modified) from
% Loraine K. Chial/ Jamie MacMahan 10/30/01.


%%

%
rho0 = 1030;
g = 9.8;


%%

%
Npts = length(presbottom);

% % % Frequency resolution
% % df = 1/(N.*dtsampling)


%%

% Should the data be detrended for eta??? --
% the difference should be tiny


%%

%
if ~exist('dtavg', 'var')
    % Will just compute of the full timeseries
    dtavg = dtsampling*Npts;
end

%
Nptsavg = dtavg/dtsampling;

% Frequency resolution per averaging segment
df = 1/(Nptsavg*dtsampling);


%%

%
ind_begin_seg = 1 : (Nptsavg) : Npts;
%
ind_end_seg = [(ind_begin_seg(2:end) - 1), Npts];

% If the last segment is too short, then merge the last two segments
% Or I could just get rid of the last segment because it's
% annoying and not helpful to deal with that exception


% Get number of segments
Navgsegments = length(ind_begin_seg);

% % %
% % if ((ind_end_seg(end) - ind_begin_seg(end) + 1) == Nptsavg)
% %     %
% %     Navgwholesegments = Navgsegments;
% % % %     ind_
% % %
% % else
% %     Navgwholesegments = Navgsegments - 1;
% % end


%% Frequency vector

% Frequency vector, in Hz (does NOT start at 0)
frequency = (1:floor(Nptsavg/2)) .* df;

%
lbelowcutoff = (frequency <= freqcutoff);

%
Nbeforecutoff = length(find(lbelowcutoff));


%% Compute low-frequency pressure and SSH
% from hydrostatic pressure

% % %
% % presbottom_butlastsegment = reshape(presbottom(ind_begin_seg(1):ind_end_seg(end-1)), Nptsavg, Navgsegments);

% All but the last segment, which likely
% has different number of points
pressure_reshaped = reshape(presbottom(ind_begin_seg(1):ind_end_seg(end-1)), Nptsavg, (Navgsegments-1));

% Compute mean pressure (and include the last segment)
mean_pressure = [mean(pressure_reshaped, 1), ...
                 mean(presbottom(ind_begin_seg(end):end))];

% May or may not be integer
centerdummy = (ind_begin_seg + ind_end_seg)./2;
% Now interpolate to the same grid as pressure
pressure_lowfrequency = interp1(centerdummy, mean_pressure, 1:Npts);

% Compute depth with the approximation that the pressure
% is hydrostatic (at time scales equal or longer than dtavg)
depth_lowfrequency_perseg = mean_pressure*1e4 ./ (rho0*g);
depth_lowfrequency = pressure_lowfrequency*1e4 ./ (rho0*g);


%% Function handle to compute wavenumbers

% (PS: frequency in Hz and k in radians per meter)
disp_rel = @(k, H, freq) g*k*tanh(k*H) - (2*pi*freq)^2;


%% Compute wavenumbers (in radians per meter) from linear wave theory

%
k_matrix = NaN(length(frequency), (Navgsegments-1));
% % k_matrix = NaN(Nbeforecutoff, Navgsegments);  % or this smaller one?

% Loop over segments
tic
for i1 = 2:(Navgsegments-1)

    % Loop over frequencies (until the cutoff)
    for i2 = 1:Nbeforecutoff

        %
        disp_rel_solver = @(k) disp_rel(k, depth_lowfrequency_perseg(i1), frequency(i2));

        % In radians per meter
        try
        k_matrix(i2, i1) = fzero(disp_rel_solver, [(2*pi/(5000)), (2*pi/(1))]);
        catch
            keyboard
        end

    end
end
toc


%% Try ffting without a loop

%
pressure_reshaped = reshape(presbottom(ind_begin_seg(1):ind_end_seg(end-1)), Nptsavg, (Navgsegments-1));
%
pressure_meanreshaped = reshape(pressure_lowfrequency(ind_begin_seg(1):ind_end_seg(end-1)), Nptsavg, (Navgsegments-1));
%
pressure_reshaped_demeaned = pressure_reshaped - pressure_meanreshaped;

%
pressure_fft_coefs = fft(pressure_reshaped_demeaned, [], 1);


% Zero out zero frequency high frequencies
pressure_fft_coefs(1, :) = 0;

%
matrix_k_forcoefs = [NaN(1, (Navgsegments-1)); k_matrix(1:(floor((Nptsavg-1)/2)), :); flipud(k_matrix)];
matrix_waterdepth = repmat(depth_lowfrequency_perseg(1:end-1), Nptsavg, 1);

%
pres2ssh = cosh(matrix_k_forcoefs .* matrix_waterdepth);
pres2ssh(isnan(pres2ssh)) = 0;

%
ssh_fft_coefs = pres2ssh .* pressure_fft_coefs;


% % % % Compute SSH fft coefficients
% % % ssh_fft_coefs = pressure_fft_coefs;
% % % ssh_fft_coefs() = pres2ssh .* ssh_fft_coefs();
% % % 
% % % 
% % %     Kp=[1 cosh(k.*hi)./cosh(k.*h)];
% % %     Kp(lo(end)+1:nf+1)=0;
% % %     
% % %     Kp=[Kp(1:(N/2+1)) fliplr(Kp(2:(N/2)))]; 


% Do inverse fft to compute timeseries
% of elevation
SSH_var_timeseries = real(ifft(ssh_fft_coefs, [], 1));

%%

% Reshape to a timeseries
SSH_timeseries = SSH_var_timeseries(:);

% Add the low frequency part
bottomdepth = SSH_timeseries + depth_lowfrequency(1:length(SSH_timeseries)).';


%%

return

%%

%
for i = 1:Navgsegments

    %
    inds_sub_aux = ind_begin_seg(i):ind_end_seg(i);




    % Compute fft
    pressure_fft = fft(presbottom());

    % Normalization factors?

    % zero out high frequencies
    

    % Use transfer function to go from
    % Fourier coefficients of pressure
    % coefficients of SSH
     
    % Do inverse fft to compute timeseries
    % of elevation

    

end





%%







%%
return

%%

%
if rem(N, 2)==0
    pz=detrend(p);
    P=fft(pz)/N;
    f = ([1:1:N/2+1]' - 1) * df;
    lo=find(f<=0.15);
    nf = length(f)-1;
    f_lo=f(lo);
    
    keyboard
    k=wavek_test(1./f_lo(2:end),abs(h));
    Kp=[1 cosh(k.*hi)./cosh(k.*h)];
    Kp(lo(end)+1:nf+1)=0;
    
    Kp=[Kp(1:(N/2+1)) fliplr(Kp(2:(N/2)))];
    
    Q=P(1:end)./Kp';
    Q=denanize(Q,0);
    eta=real(ifft(Q)).*N;
    
% %     S=2.*Q(1:N/2+1).*conj(Q(1:N/2+1))./df;
% %     sum_spec=sum(S).*df;
% %     variance=var(eta,1);
% %     Hmo=4.01.*sqrt(sum_spec);
% %     [u]=find(S(3:end)==max(S(3:end)));
% %     peak_f=f(u);
% %     Tp=1./peak_f;

else
    %
    pz=detrend(p);
    P=fft(pz)/N;
    f = ([1:1:N/2+1]'-1) * df;
    lo=find(f<=0.15);
    nf = length(f);
    f_lo=f(lo);
    k=wavek_test(1./f_lo(2:end),abs(h));

    Kp=[1 cosh(k.*hi)./cosh(k.*h)];
    Kp(lo(end)+1:nf+1)=0;
    % keyboard
    Kp=[Kp(1:(N/2+.5)) fliplr(Kp(2:(N/2+.5)))];
    %keyboard

    Q=P(1:end)./Kp';
    Q=denanize(Q,0);
    eta=real(ifft(Q)).*N;
    
% %     S=2.*Q(1:N/2+.5).*conj(Q(1:N/2+.5))./df;
% %     sum_spec=sum(S).*df;
% %     variance=var(eta,1);
% %     Hmo=4.01.*sqrt(sum_spec);
% %     [u]=find(S(3:end)==max(S(3:end)));
% %     peak_f=f(u);
% %     Tp=1./peak_f;
end    

Hmo = 1;
Tp = 1;

