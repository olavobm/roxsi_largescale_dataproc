function dirfilesout = Signature_orderfiles(dirmatdata, mooringID)
%% dirfilesout = SIGNATURE_ORDERFILES(dirmatdata, mooringID)
% 
%   inputs
%       - dirmatdata: directory with *.mat files that have Signature1000
%                     data from the large-scale array of ROXSI 2022.
%       - mooringID: three-character mooring ID.
%
%   outputs
%       - dirfilesout: output of dir function.
%
%
% Function SIGNATURE_ORDERFILES.m outputs the same as the dir.m
% function, but in the correct order of Signature1000 files.
% These files have a few different filename formats, and bad
% formatting such that dir.m does not output in the correct order.



%%

%
if (strcmp(mooringID, 'A01')) || (strcmp(mooringID, 'B15')) || ...
   (strcmp(mooringID, 'C01')) || strcmp(mooringID, 'X11')

    %
    dir_1_aux = dir(fullfile(dirmatdata, ['*_ROXSIsig' mooringID '_00*.mat']));
    %
    dir_2_aux = dir(fullfile(dirmatdata, ['*_ROXSIsig' mooringID '_1.mat']));
    dir_3_aux = dir(fullfile(dirmatdata, ['*_ROXSIsig' mooringID '_2.mat']));

    %
    dirfilesout = [dir_1_aux; dir_2_aux; dir_3_aux];
    

%
elseif strcmp(mooringID, 'B10')
    %
    dir_1_aux = dir(fullfile(dirmatdata, '*_ROXSI_B10_*.mat'));
    dir_2_aux = dir(fullfile(dirmatdata, '*_ROXSI_B10.mat'));

    %
    dirfilesout = [dir_1_aux; dir_2_aux];

%
elseif strcmp(mooringID, 'B13') || strcmp(mooringID, 'X05')


    %
    dir_all_aux = dir(fullfile(dirmatdata, ['*_ROXSIsig' mooringID '_*.mat']));
    %
    dirfilesout = dir_all_aux;
    %
    for i = 1:length(dir_all_aux)
        %
        dirfilesout(i) = dir(fullfile(dirmatdata, ['*_ROXSIsig' mooringID '_' num2str(i) '.mat']));
    end


%
elseif strcmp(mooringID, 'B17')

    %
    dir_1_aux = dir(fullfile(dirmatdata, '*_ROXSIsigB17_*.mat'));
    dir_2_aux = dir(fullfile(dirmatdata, '*_ROXSIsigB17.mat'));

    %
    dirfilesout = [dir_1_aux; dir_2_aux];


    %
else
    error('!!!!')
end

