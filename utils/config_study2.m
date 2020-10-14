
%% SPECIFY SUBJECTS
SpecifySubjectDetails_StudyII;
%% GENERAL STUFF
%data_dir = get_homedir();
data_dir = strrep(rawdatadir,'/MEGData/','');
is.fs = filesep;
is.rootBehav = [data_dir,'\BehavData\'];
is.rootMEG = rawdatadir;
is.networkPath = [data_dir,'\MEGData\'];   % where to pull the unpreprocessed MEG data from
%is.smthdPath = ' ';
is.OPTPath = [data_dir,'\OPTdata\'];
%is.BFPath = 'C:\Users\yliu\Dropbox\Projects\UCL\Experiment04_Replay_StrLearn\Data\Beamforming\';
is.AnalysisPath = [data_dir,'\Analysis\'];

is.pNames = {'epsilon', 'beta', 'gammaG', 'gammaS', 'psi4', 'psiA'};
is.smoothFact = 6; 
is.tss = -15:100;

is.whichTimes = [31:50]; % this is for stim code, from 150ms to 350ms 31:50
is.whichTimesPos = [42:51, 67:76]; % this is for pos code, from 250ms to 350ms 42:51; from 500ms to 600ms 67:76; 
is.whichTimesSeq = [28:37, 67:76]; % this is for seq code, from 100ms to 200ms 28:37; from 500ms to 600ms 67:76; 

is.msPerSample = 1000 * is.smoothFact / 600; 
is.lgncy = 1:60;  % the latencies to consider cross-correlation at, in units of samples
is.Ltimes = is.lgncy*is.msPerSample/1000;
is.nShuf = 20; %20;%20;%29;
is.ENL1 = 0.001:0.001:0.01;   % L1=lambda*alpha   L2=(lambda-lambda*alpha)/2
is.whichSubj = 1:is.nSubj;%11:30;  % subject11 is the first non-pilot. 21 and 26 are excluded for huge artifacts.
is.highpass = 0.5; % decide which hiph pass spec is used for the further analysis
is.bandpass = [];
is.eyeBlinkSens = 0;  % change to 1 if we want to only focus on the sensory channels.
is.Zebcorrection = 1; % change to 1 if we want to also do Zeb's correction on the preprocessed data.
is.delete = 0; % if set to 1, will directly delete the sample from the data set, set to 0 is the same as Zeb.
is.prestimnormalisation=0; % if set to one, trials are onyl normalised using the mean and std from points pre stimulus presentation to avoide bias
is.scalefunc=1;
is.balancedata=1;
is.goodsub =[1,3,5:10,11:17,20:22,24:27]; % new analysis

pd =pwd;
if strcmp(pd(1),'/')
    % mac: change all address identifiers:
    for i=1:length(is.rootBehav)
        if strcmp(is.rootBehav(i),'\')
            is.rootBehav(i)='/';
        end
    end
    for i=1:length(is.rootMEG)
        if strcmp(is.rootMEG(i),'\')
            is.rootMEG(i)='/';
        end
    end
    for i=1:length(is.networkPath)
        if strcmp(is.networkPath(i),'\')
            is.networkPath(i)='/';
        end
    end
    for i=1:length(is.OPTPath)
        if strcmp(is.OPTPath(i),'\')
            is.OPTPath(i)='/';
        end
    end
    for i=1:length(is.AnalysisPath)
        if strcmp(is.AnalysisPath(i),'\')
            is.AnalysisPath(i)='/';
        end
    end
end

% and omit bad subjects:
is.fnDate = is.fnDate(is.goodsub);
is.fnSID = is.fnSID(is.goodsub);
is.fnBehav = is.fnBehav(is.goodsub);
is.fnTrainClass = is.fnTrainClass(is.goodsub);
is.fnMEG = is.fnMEG(is.goodsub);
is.MEGruns = is.MEGruns(is.goodsub);
is.taskversion = is.taskversion(is.goodsub);
is.subjects_to_do = 1:length(is.goodsub);


