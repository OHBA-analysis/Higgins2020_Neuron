%% need to insert the preproc files before this!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path
%     
basedir='/Users/chiggins/';
cd([basedir 'Documents/MATLAB/osl/osl-core']);
osl_startup;
wd = [basedir,'data/Neuron2020/'];
% Add netlab and fmt
addpath( fullfile(osldir,'ohba-external','netlab3.3','netlab') );
addpath( fullfile(osldir,'ohba-external','fmt') );

download_path = [basedir,'Documents/MATLAB/Neuron2020/'];
addpath(genpath(download_path));

TESTRUN = true; % when this is on, just run full script for one subject only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generic filename / parameter details:
session_name{1} = 'CanonicalRS/';
nscans{1} = 55;
Fs_to_run{1} = 250;
nSj{1} = 55;
rawdatapath{1} = '/Volumes/CamsHD2/NottinghamRS/raw_data';

session_name{2} = 'Study1/';
nscans{2} = 42;
Fs_to_run{2} = [250,600];
nSj{2} = 21;
rawdatapath{2} = '/Volumes/CamsHD2/YunzheData/Replaydata4Cam/RawData';

session_name{3} ='Study2/';
nscans{3} = 44;
Fs_to_run{3} = [250,600];
nSj{3} = 22;
rawdatapath{3} = '/Volumes/CamsHD2/YunzheData/StrLearn_MEGexp/MEGData';

session_name{4} = 'Study1_FLI/';
nscans{4} = 2*21;
Fs_to_run{4} = 250;
nSj{4} = 21;

session_name{5} = 'Study2_FLI/';
nscans{5} = 3*22;
Fs_to_run{5} = 250;
nSj{5} = 22;
   

%% Filter, beamform, parcellate and signflip:

for whichstudy=1:5

    spmfilesdir=[wd,session_name{whichstudy} ];
    
    if TESTRUN;nsub=1;else, nsub=nscans{whichstudy};end

    for iFreq=1:length(Fs_to_run{whichstudy})
        % general settings:
        S=[];
        %downsampling frequency
        Fs=Fs_to_run{whichstudy}(iFreq);
        if Fs==600
            preproc_name='600Hz/';
            S.freq_range=[0 0];
        else
            preproc_name='250Hz/';
            S.freq_range=[1 45];
        end
        
        S.do_prepare=1; 
        S.do_hmm=0; 
        S.do_spectral_estimation=0;
        S.preproc_name=preproc_name;
        S.session_name=session_name{whichstudy}; 
        S.spmfilesdir=spmfilesdir;
        S.prep_sessions_to_do=1:nsub;
        S.hmm_sessions_to_do=1:nsub;
        S.num_embeddings=14;
        
        S.parcellations_dir = [osldir, '/parcellations'];

        S.K = 12;

        if whichstudy==1
            disp(S);
            [hmm hmmfname hmmoptions settings_prepare] = run_full_hmm_pipeline_neuron2020edit(S);
            disp(hmmfname);
            signfliptemplatesubj = [spmfilesdir,preproc_name,'bfnew_1to45hz/sfold_giles_symmetric_f_session',int2str(settings_prepare.templatesubj)];
        else
            S.signfliptemplatesubj =signfliptemplatesubj;
            [hmm hmmfname hmmoptions settings_prepare] = run_full_hmm_pipeline_neuron2020edit(S);
        end
    end
end

%% Now fit HMM:

% keep same S settings as above, just change actions performed:
S.do_prepare=0; 
S.do_hmm=1; 
S.do_spectral_estimation=0;

for whichstudy = 1:3
    
    spmfilesdir=[wd,session_name{whichstudy} ];
    preproc_name='250Hz/';
    S.preproc_name=preproc_name;
    S.session_name=session_name{whichstudy}; 
    S.spmfilesdir=spmfilesdir;
        
    if whichstudy==1
        if isfield(S,'templateHMM')
            S = rmfield(S,'templateHMM');
        end
        for i=1:Nruns
            S.hmm_name=int2str(i); 
            disp(S);
            [hmm hmmfname hmmoptions settings_prepare] = run_full_hmm_pipeline_neuron2020edit(S);
            disp(hmmfname);
        end
    else
        for i=1:Nruns
            templatestring = 'usingtemplate';
            S.hmm_name = [int2str(i),templatestring];
            S.templateHMM = [wd,session_name{1},preproc_name,'/hmm_1to45hz/'...
                    'hmm',int2str(i),'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(S.K),'_big1_dyn_modelhmm.mat'];
            [hmm hmmfname hmmoptions settings_prepare] = run_full_hmm_pipeline_neuron2020edit(S);
        end
    end
end

%% determine model with lowest free energy to keep:

for whichstudy=1:3
    if whichstudy>1
        templatestring = 'usingtemplate';
    else
        templatestring = '';
    end
    for i=1:Nruns
        hmmfname = [wd,session_name{whichstudy},preproc_name,'/hmm_1to45hz/'...
                    'hmm',int2str(i),templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(S.K),'_big1_dyn_modelhmm.mat'];
        temp = load(hmmfname);
        FEcompare(i,whichstudy) = temp.hmm.fehist(end);
    end
end
[~,bestmodel] = min(sum(FEcompare,2));

save([wd,'bestmodel.mat'],'bestmodel','FEcompare');



%% Analyses for each figure looped over different studies:
%K = S.K;

bestmodel=5;
K=12;

for whichstudy=3
    
    %NeuronFig2Analyses;
    clearvars -except K whichstudy bestmodel wd session_name nscans Fs_to_run nSj
    NeuronFig3Analyses
    clearvars -except K whichstudy bestmodel wd session_name nscans Fs_to_run nSj
    NeuronFig4Analyses
    clearvars -except K whichstudy bestmodel wd session_name nscans Fs_to_run nSj
%     if whichstudy>1
%         NeuronFig5Analyses
%         clearvars -except K whichstudy bestmodel wd session_name nscans Fs_to_run nSj
%     end
    close all;
end