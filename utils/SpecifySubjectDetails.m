%% MEGruns: 'sss' = simplyShowStimuli, 'lci' = lexicallyCuedImagination, 'str' = Structure Learning, 'rst' = resting state
%%          'rwd' = reward learning,  'tst' = testing phase, 'tsg' = testing phase for generalization,

is.nSubj = 1;
is.fnDate{is.nSubj} = '2017-06-26'; is.fnSID{is.nSubj} = '02';
is.fnBehav{is.nSubj} = {'0002'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05019_Yunzhe_20170626'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=5; % 5.0 (3 pairs, with fixed order), rest in the beginning: first session 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [2]
is.fnDate{is.nSubj} = '2017-09-05'; is.fnSID{is.nSubj} = '04';
is.fnBehav{is.nSubj} = {'0004'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05095_Yunzhe_20170905'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=5; % 5.0 (3 pairs, with fixed order), rest in the beginning: first session 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [3]  
 
is.fnDate{is.nSubj} = '2017-09-12'; is.fnSID{is.nSubj} = '10';
is.fnBehav{is.nSubj} = {'0010'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05110_Yunzhe_20170912'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [4] 
 
is.fnDate{is.nSubj} = '2017-09-14'; is.fnSID{is.nSubj} = '11';
is.fnBehav{is.nSubj} = {'0011'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05116_Yunzhe_20170914'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [5] 
 
is.fnDate{is.nSubj} = '2017-09-14'; is.fnSID{is.nSubj} = '12';
is.fnBehav{is.nSubj} = {'0012'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05117_Yunzhe_20170914'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [6] 

is.fnDate{is.nSubj} = '2017-09-19'; is.fnSID{is.nSubj} = '14';
is.fnBehav{is.nSubj} = {'0014'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05121_Yunzhe_20170919'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [7] 

is.fnDate{is.nSubj} = '2017-09-21'; is.fnSID{is.nSubj} = '16';
is.fnBehav{is.nSubj} = {'0016'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05128_Yunzhe_20170921'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [8] 

is.fnDate{is.nSubj} = '2017-09-27'; is.fnSID{is.nSubj} = '19';
is.fnBehav{is.nSubj} = {'0019'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05132_Yunzhe_20170927'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [9] 

is.fnDate{is.nSubj} = '2017-09-28'; is.fnSID{is.nSubj} = '20';
is.fnBehav{is.nSubj} = {'0020'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05133_Yunzhe_20170928'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [10] 

is.fnDate{is.nSubj} = '2017-09-29'; is.fnSID{is.nSubj} = '21';
is.fnBehav{is.nSubj} = {'0021'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05134_Yunzhe_20170929'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [11]

is.fnDate{is.nSubj} = '2017-10-05'; is.fnSID{is.nSubj} = '22';
is.fnBehav{is.nSubj} = {'0022'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05135_Yunzhe_20171005'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [12]

is.fnDate{is.nSubj} = '2017-10-06'; is.fnSID{is.nSubj} = '23';
is.fnBehav{is.nSubj} = {'0023'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05136_Yunzhe_20171006'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [13]

is.fnDate{is.nSubj} = '2017-10-10'; is.fnSID{is.nSubj} = '24';
is.fnBehav{is.nSubj} = {'0024'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05138_Yunzhe_20171010'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [14] 

is.fnDate{is.nSubj} = '2017-10-24'; is.fnSID{is.nSubj} = '27';
is.fnBehav{is.nSubj} = {'0027'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05151_Yunzhe_20171024'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [15] 
is.fnDate{is.nSubj} = '2017-10-24'; is.fnSID{is.nSubj} = '28';
is.fnBehav{is.nSubj} = {'0028'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05152_Yunzhe_20171024'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [16] 

is.fnDate{is.nSubj} = '2017-10-24'; is.fnSID{is.nSubj} = '29';
is.fnBehav{is.nSubj} = {'0029'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05153_Yunzhe_20171024'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [17]

is.fnDate{is.nSubj} = '2017-10-25'; is.fnSID{is.nSubj} = '30';
is.fnBehav{is.nSubj} = {'0030'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05155_Yunzhe_20171025'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [18] 

is.fnDate{is.nSubj} = '2017-10-26'; is.fnSID{is.nSubj} = '31';
is.fnBehav{is.nSubj} = {'0031'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05159_Yunzhe_20171026'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [19]

is.fnDate{is.nSubj} = '2017-10-26'; is.fnSID{is.nSubj} = '32';
is.fnBehav{is.nSubj} = {'0032'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05160_Yunzhe_20171026'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;


%% [20] 
is.fnDate{is.nSubj} = '2017-10-28'; is.fnSID{is.nSubj} = '36';
is.fnBehav{is.nSubj} = {'0036'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05165_Yunzhe_20171028'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [21] 

is.fnDate{is.nSubj} = '2017-10-29'; is.fnSID{is.nSubj} = '38';
is.fnBehav{is.nSubj} = {'0038'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'Replay_0001'};
is.fnMEG{is.nSubj} = 'mg05167_Yunzhe_20171029'; 
is.MEGruns{is.nSubj} = {'sss' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'rwd' 'rst' 'tst' 'tsg'}; 
is.taskversion{is.nSubj}=6; % 6.0 (random order across participants), rest after FL 
is.nSubj = is.nSubj + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


is.nSubj = is.nSubj - 1;