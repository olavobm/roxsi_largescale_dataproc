%% Script to extract some log info from RBR sensors (SoloDs and SoloTs).
%
% Extract info from programming and recovery log files.
%
% Here is the general timeline:
%   - 2022-06-07: programmming of SoloD's on Jamie's computer (file ruskin_jamie_programming.log)
%   - 2022-06-08: programmming of SoloT's on Jamie's computer (file ruskin_jamie_programming.log)
%
%   - 2022-06-14: programmming of 6 SoloD's on Falk's computer (file ruskin_26July2022.log)
%   - 2022-07-22: recovery of both??? on Falk's computer (file ruskin_26July2022.log)
%   - 2022-07-26: recovery of remaining SoloD's and SoloT's (on B10at) on
%                 Falk's computer (file ruskin_26July2022.log).
%
% Something wrong happened with SoloD's 077802 (A08p) and 077822 (C04p),
% such that we do not have data from them.





%% Jamie's SoloD's and day when they were programmed

% All, but the two we don't have data
list_Jamies_SoloD = {'077810'; '077818'; '077814'; '077807'; '077815'; ...
                     '077816'; '077589'; '077590'; '077817'; '077587'; ...
                     '077848'; '077591'; '077845'; '077813'; ...
                     '077808'; '077586'; '077272'; '077594'; ...
                     '077595'; '077809'; '077819'; '077804'; '077274'};

% % % All SoloD's
% % list_Jamies_SoloD = {'077810'; '077818'; '077814'; '077807'; '077815'; ...
% %                      '077816'; '077589'; '077590'; '077817'; '077587'; ...
% %                      '077848'; '077802'; '077591'; '077845'; '077813'; ...
% %                      '077808'; '077586'; '077272'; '077594'; '077822'; ...
% %                      '077595'; '077809'; '077819'; '077804'; '077274'};

% % % The two we don't have data
% % list_Jamies_SoloD = {'077802'; '077822'};

%
day_program_Jamies_SoloD = '2022-06-07';

% only 077274 in Falk's log on 07/22


%% Falk's SoloD's

% Falk's SoloDs (on a different programming/log file)
list_Falks_SoloD = {'124020'; '124022'; '124073'; '124039'; '124015'; '124017';};

%
day_program_Falks_SoloD = '2022-06-14';


%% Jamie's SoloT's and day when they were programmed
% (all SoloT's that were deployed belong to Jamie)

%
list_SoloT = {'077090'; '077151'; '077296'; '076976'; '076990'; '077051'; '077187'; ...    % X02t
              '077083'; '077157'; '077184'; '076997'; '077082'; '077299'; ...    % X05at
              '077002'; '077164'; '077081'; '077004'; '077155'; '077064'; '076988'; '077073'; '077242'; '077297'; '077269'; ...    % B02at
              '077243'; '077241'; '076979'; '076958'; '077089'; '077031'; '077069'; '077003'; '077226'; ...    % B04at
              '077032'; '076897'; '077025'; '076972'; '077057'; '076980'; '077079'; '077172'; '076893'; '077279'; '077274'; ...    % B06at
              '077068'; '077292'; '076954'; '076925'; '076939'; '077168'; '077156'; '077050'; '077056'; ...    % B08at
              '077277'; '076905'; '076932'; '076964'; '077283'; '077195'; '077190'};    % B10at

%
day_program_SoloT = '2022-06-08';

%
day_recovery_SoloT_1 = '2022-07-22';    % for all, except for those on B10at, which is probably on the 26th
day_recovery_SoloT_2 = '2022-07-26';    % for those on B10at

%
alldays_recovery_SoloT = [repmat({day_recovery_SoloT_1}, 53, 1); ...
                          repmat({day_recovery_SoloT_2}, 7, 1)];

% % % Checking just those on B10at
% % list_SoloT = list_SoloT(54:end);
% % alldays_recovery_SoloT = alldays_recovery_SoloT(54:end);


%%

addpath('/Users/olavobm/Documents/ROXSI_Postdoc/MyResearch/ROXSI/Common_Code/LargeScale_Data_2022/code_proc/')


%%
% ----------------------------------------------------
% -------------- PROGRAMMING/DEPLOYMENT --------------
% ----------------------------------------------------


%% Get programming log info from Jamie's SoloD's

%
for i = 1:length(list_Jamies_SoloD)
    bla = RBR_extract_loginfo('ruskin_jamie_programming.log', day_program_Jamies_SoloD, list_Jamies_SoloD{i});
end


%% Get programming log info from Falk's SoloD's

%
for i = 1:length(list_Falks_SoloD)
    bla = RBR_extract_loginfo('ruskin_26July2022.log', day_program_Falks_SoloD, list_Falks_SoloD{i});
end


%% Get programming log info from SoloT's

%
for i = 1:length(list_SoloT)
    bla = RBR_extract_loginfo('ruskin_jamie_programming.log', day_program_SoloT, list_SoloT{i});
end


%%
% ----------------------------------------------------
% --------------------- RECOVERY ---------------------
% ----------------------------------------------------

%% SoloD's

% % % All of Jamie's SoloD's in this log file (except for those
% % % we don't have data, where this log file does not have any
% % % line whatsoever with the SN of these 2 instruments)
% % for i = 1:length(list_Jamies_SoloD)
% %     bla = RBR_extract_loginfo('ruskin.log', '2022-07-22', list_Jamies_SoloD{i});
% % end


% % % One of Falk's SoloD's
% % bla = RBR_extract_loginfo('ruskin_124039.log', '2022-07-26', '124039');


% 5 of Falk's SoloD's
for i = 1:length(list_Falks_SoloD)
    bla = RBR_extract_loginfo('ruskin.log', '2022-07-22', list_Falks_SoloD{i});
end


%% SoloT's

%
for i = 1:length(list_SoloT)
    bla = RBR_extract_loginfo('ruskin_26July2022.log', alldays_recovery_SoloT{i}, list_SoloT{i});
end




