function dirpath = repo_dirpath('')
%%
%
%
%

%%

%
if nargin==0
    %
    full_filename = mfilename('fullpath');
    %
    ind_sep = find(full_filename, filesep);

    %
    dirpath = full_filename(1:(ind_sep(end)-1));
end
