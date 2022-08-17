function dirpath = repo_dirpath()
%% dirpath = repo_dirpath()
%
% REPO_DIRPATH.m is a function with no arguments
% that returns the directory of where this function
% is located. The idea is that this function is
% in outermost directory of this repository.

%
full_filename = mfilename('fullpath');

%
ind_sep = strfind(full_filename, filesep);

%
dirpath = full_filename(1:(ind_sep(end)-1));

