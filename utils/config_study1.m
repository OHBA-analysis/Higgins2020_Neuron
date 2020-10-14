
%% SPECIFY SUBJECTS
SpecifySubjectDetails

%% GENERAL STUFF
is.fs = filesep;
disk_address='/Volumes/CamsHD2/YunzheData/Replaydata4Cam';
%cd(disk_address);
%is.rootBehav = [disk_address,'/BehavData/']; %'C:\Users\yliu\Dropbox\Projects\UCL\Experiment03_MEG_Replay\BehavData\';
is.rootBehav = '/Users/chiggins/data/YunzheData/Replaydata4Cam/BehavData/';
is.rootMEG = [disk_address,'/MEGData/'];%'C:\Users\yliu\Dropbox\Projects\UCL\Experiment03_MEG_Replay\MEGData\';
is.networkPath = [disk_address,'/MEGData/RawData/'];%'C:\Users\yliu\Dropbox\Projects\UCL\Experiment03_MEG_Replay\MEGData\RawData\';   % where to pull the unpreprocessed MEG data from
is.smthdPath =[disk_address,'/MEGData/DF3Data/'];%'C:\Users\yliu\Dropbox\Projects\UCL\Experiment03_MEG_Replay\MEGData\DF3data\';
is.OPTPath = [disk_address,'/MEGData/OPTData/'];%C:\Users\yliu\Dropbox\Projects\UCL\Experiment03_MEG_Replay\MEGData\OPTdata\';
is.AnalysisPath = [disk_address,'/Analysis/'];
is.pNames = {'epsilon', 'beta', 'gammaG', 'gammaS', 'psi4', 'psiA'};


is.smoothFact = 6;
% note that we wish to train on a single timepoint's data, then test
% generalisation across time; this amounts to a training window of is.tss
% offsets from the training times specified in is.whichTimes
is.tss = -15:50; 
is.whichTimes = [31,32,33,34,35,36,37,38];
%note 31 denotes is.tss(31), or t=140msec, as training point
is.msPerSample = 1000 * is.smoothFact / 600; 
is.lgncy = 1:60;  % the latencies to consider cross-correlation at, in units of samples
is.Ltimes = is.lgncy*is.msPerSample/1000;
is.nShuf = 20; %20;%20;%29;
is.ENL1 = 0.001:0.001:0.01;   % L1=lambda*alpha   L2=(lambda-lambda*alpha)/2
is.whichSubj = 1:is.nSubj;%11:30;  % subject11 is the first non-pilot. 21 and 26 are excluded for huge artifacts.
is.highpass = 1; % decide which hiph pass spec is used for the further analysis
is.bandpass = [];
is.eyeBlinkSens = 0;  % change to 1 if we want to only focus on the sensory channels.
is.Zebcorrection = 0; % change to 1 if we want to also do Zeb's correction on the preprocessed data.
is.delete = 0; % if set to 1, will directly delete the sample from the data set, set to 0 is the same as Zeb.
is.prestimnormalisation=0; % if set to one, trials are onyl normalised using the mean and std from points pre stimulus presentation to avoide bias
is.scalefunc=0;
is.balancedata=1;
% is.exclu = [21; 26]; % sessions to exclude. 21 and 26 have large artifacts. i'm not sure if we can exclude 11. (indexed in the 1:30 scheme)
% is.dupes = [12 13; 15 18; 14 19; 22 28; 25 27; 23 30];   % two sessions from the same subject to merge (indexed in the 1:30 scheme)

is.behavdatafilenames={'26-Jun-2017\02_Replay_0001.mat','05-Sep-2017\04_Replay_0001.mat',...
    '12-Sep-2017\10_Replay_0001.mat','14-Sep-2017\11_Replay_0001.mat',...
    '14-Sep-2017\12_Replay_0001.mat','19-Sep-2017\14_Replay_0001.mat',...
    '21-Sep-2017\16_Replay_0001.mat','27-Sep-2017\19_Replay_0001.mat',...
    '28-Sep-2017\20_Replay_0001.mat','29-Sep-2017\21_Replay_0001.mat',...
    '05-Oct-2017\22_Replay_0001.mat','06-Oct-2017\23_Replay_0001.mat',...
    '10-Oct-2017\24_Replay_0001.mat','24-Oct-2017\27_Replay_0001.mat',...
    '24-Oct-2017\28_Replay_0001.mat','24-Oct-2017\29_Replay_0001.mat',...
    '25-Oct-2017\30_Replay_0001.mat','26-Oct-2017\31_Replay_0001.mat',...
    '26-Oct-2017\32_Replay_0001.mat','28-Oct-2017\36_Replay_0001.mat',...
    '29-Oct-2017\38_Replay_0001.mat'};