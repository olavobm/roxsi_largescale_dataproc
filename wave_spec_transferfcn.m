function Tfcn = wave_spec_transferfcn(varstransfer, zhabobs, k, H)
%% Tfcn = WAVE_SPEC_TRANSFERFCN(varstransfer, zhabobs, k, H)
%
%   inputs
%       - varstransfer:
%       - k:
%       - H:
%       - zhabobs:
%
%   outputs
%       - Tfcn
%
%
%
%
%
%
%
% Olavo Badaro Marques.

%%

list_available_transfers = {'ptoeta', 'etatop'};

%
if ~any(strcmp(list_available_transfers, varstransfer))
    error(['Input specifying transfer function not ' ...
           'accepted by ' mfilename() '.'])
end


%% Check dimensions

%
if ~isequal(size(k), size(H))
    %
    if ~isvector(H) && ~isvector(k)
        %
        error('H and k have inconsistent dimensions.')
    end
end


%% Set correct dimensins for k and H if necessary

%
if isvector(H)
    %
    H = H(:);    % make sure it's a column vector

    %
    if size(k, 1) ~= length(H)
        k = k.';
        if size(k, 1) ~= length(H)
            %
            error('H and k have inconsistent dimensions.')
        end
    end

    % Turn H into an array of the same dimension as k
    H = repmat(H, 1, size(k, 2));
end


%% Compute transfer function

%
if strcmp(varstransfer, 'ptoeta')
    %% From pressure to vertical elevation
    Tfcn = cosh(k .* H) ./ cosh(k * zhabobs);
%
elseif strcmp(varstransfer, 'etatop')
    %% From vertical elevation to pressure
    Tfcn = cosh(k * zhabobs) ./ cosh(k .* H);
end

%% Square to do the conversion between spectra
Tfcn = Tfcn.^2;


