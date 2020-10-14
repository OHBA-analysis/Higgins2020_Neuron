%% main_preproc
%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path
    
basedir='/Users/chiggins/';
cd([basedir 'Documents/MATLAB/osl/osl-core']);
download_path = [basedir,'Documents/MATLAB/Neuron2020/'];
addpath(genpath(download_path));
wd = [basedir,'/data/Neuron2020/'];
TESTRUN = false; % when this is on, just run full script for one subject only


%% generic filename / parameter details:
session_name{1} = 'CanonicalRS/';
nsubjects{1} = 55;
Fs_to_run{1} = 250;
rawdatapath{1} = '/Volumes/CamsHD2/NottinghamRS/raw_data/';

session_name{2} = 'Study1/';
nscans{2} = 42;
Fs_to_run{2} = [250,600];
nSj{2} = 21;
rawdatapath{2} = '/Volumes/CamsHD2/YunzheData/Replaydata4Cam/RawData/';

session_name{3} ='Study2/';
nscans{3} = 44;
Fs_to_run{3} = [250,600];
nSj{3} = 22;
rawdatapath{3} = '/Volumes/CamsHD2/YunzheData/StrLearn_MEGexp/MEGData/';

session_name{4} = 'Study1_FLI/';
nsubjects{4} = 2*21;
Fs_to_run{4} = 250;

session_name{5} = 'Study2_FLI/';
nsubjects{5} = 3*22;
Fs_to_run{5} = 250;

%%%%%%%%%%%%%%%%%%
%% Settings

do_convert=true;

for whichstudy = [3] %this denotes the study being analysed - either 1 (Nottingham data) or 2 or 3 (Yunzhe's analysis)


    for iFreq=2:length(Fs_to_run{whichstudy})
        %downsampling frequency
        Fs=Fs_to_run{whichstudy}(iFreq);
        if Fs==600
            preproc_name='600Hz/';
        else
            preproc_name='250Hz/';
        end
        
        % Initialise OSL
        osl_startup;

        % Add netlab and fmt
        addpath( fullfile(osldir,'ohba-external','netlab3.3','netlab') );
        addpath( fullfile(osldir,'ohba-external','fmt') );

        % MRI conversions need freesurfer:
        addpath('/Applications/freesurfer/matlab/');

        %%%%%%%%%%%%%%%%%%
        %% SET UP FILE NAMES

        studydir = [wd,session_name{whichstudy}]; % this is the directory the new files will all be saved into
        spmfilesdir = [studydir,preproc_name];
        mkdir(spmfilesdir);
        africadir   = [spmfilesdir 'africa/'];

        if whichstudy==2 || whichstudy==4
            subjects = 1:21;
            config_study1;
            subjects_to_do = subjects;
            rawdatadir = rawdatapath{2};
        elseif whichstudy==3 || whichstudy==5
            subjects = 1:22;
            rawdatadir = rawdatapath{3};
            config_study2;
            subjects_to_do = subjects;
        elseif whichstudy==1
            config_study0;
            subjects = is.subjects_to_do;
            rawdatadir = rawdatapath{1};
        end
        
        if TESTRUN;subjects=1;end

        if whichstudy~=1
            ctf_files            = cell(2*length(subjects),1);
            spm_files            = cell(2*length(subjects),1);
            structural_files        = cell(2*length(subjects),1);
            pos_files               = cell(2*length(subjects),1);
        else
            fields_to_workspace(is);
        end
        
        if TESTRUN;subjects_to_do=subjects_to_do(1);end
        
        if whichstudy>1
            for j = 1:length(subjects)

                inputsubjDir = [strrep(is.fnDate{j},'-','') '/' is.fnMEG{j} '_'];    
                if whichstudy<4
                    % resting state scans:
                    sessions = find(ismember(is.MEGruns{j},'rst')); % find resting state scans
                else
                    % functional localiser scans:
                    sessions = find(ismember(is.MEGruns{j},'lci'));
                end

                for k = 1:length(sessions)

                    session = num2str(sessions(k), '%02d');
                    s = length(sessions)*(j-1) + k;
                    dsfile = [rawdatadir inputsubjDir session '.ds'];
                    if ~isempty(dsfile)
                        ctf_files{s} = dsfile;
                        % set up a list of SPM MEEG object file names (we only have one here)
                        spm_files{s}    = [spmfilesdir '/RSDataRun' int2str(s) '.mat'];
                    end
                    
                    % no structural files:
                    structural_files{s} = [];
                    
                    % list of head position files
                    pfname=[rawdatadir 'pos_files/' inputsubjDir '.pos'];
                    pf = dir(pfname); 
                    if ~isempty(pf)
                        pos_files{s}=pfname;
                    end

                end
            end
            subjects_to_do = 1:(length(sessions)*length(subjects));
        end
        
        %% CONVERT FROM .ds TO AN SPM MEEG OBJECT:

        if do_convert

            S2=[];
            mkdir( [spmfilesdir, '/results'])
            for ss=1:length(subjects_to_do), % iterates over subjects    
                s=subjects_to_do(ss);

                S2.outfile=spm_files{s};

                % The conversion to SPM will show a histogram of the event codes
                % and correspond to those listed below in the epoching section
                D = osl_import(ctf_files{s},S2);

                % crop:
                if whichstudy~=1
                        S3 = struct;
                        %S.D = fullfile(localPath,'spmeeg.mat');
                        S3.D = spm_files{s};
                        S3.prefix='';            
                        event = ft_read_event(ctf_files{s});
                        sample = [event(find(strcmp('UPPT001', {event.type}))).sample]';
                        S3.timewin = [0,round(sample(end)/600*1000)+5000];
                        try
                            D=spm_eeg_crop(S3);
                        catch
                            D=spm_eeg_crop(S3);
                        end
                        D = D.chantype(find(strcmp(D.chanlabels,'UADC001')),'EOG');
                        D = D.chanlabels(find(strcmp(D.chanlabels,'UADC001')),'EOG1');

                        D = D.chantype(find(strcmp(D.chanlabels,'UADC002')),'EOG');
                        D = D.chanlabels(find(strcmp(D.chanlabels,'UADC002')),'EOG2');

                        D = D.chantype(find(strcmp(D.chanlabels,'UADC003')),'EOG');
                        D = D.chanlabels(find(strcmp(D.chanlabels,'UADC003')),'EOG3');

                        D.save();
                end
            end

            %%%%%%%%%%%%%%%%%%%
            % Sort out fiducials:
            if whichstudy==1
                for ss=1:length(subjects_to_do) % iterates over subjects    
                    f=subjects_to_do(ss);

                    spm_file = prefix(spm_files{f},'');
                    D = spm_eeg_load(spm_file);

                    fID = fopen(pos_files{f});
                    fid_data = textscan(fID, '%s %f %f %f');
                    fclose(fID);

                    fid_new = [];

                    % Fiducials:
                    fid_inds = [find(strcmpi(fid_data{1},'nasion'))
                    find(strcmpi(fid_data{1},'left'))
                    find(strcmpi(fid_data{1},'right'))];

                    fid_new.fid.pnt = [fid_data{2}(fid_inds) fid_data{3}(fid_inds) fid_data{4}(fid_inds)] * 10;
                    fid_new.fid.label = {'nas';'lpa';'rpa'};

                    % Headshape:
                    hs_inds = setdiff(2:length(fid_data{1}),fid_inds);
                    fid_new.pnt = [fid_data{2}(hs_inds) fid_data{3}(hs_inds) fid_data{4}(hs_inds)] * 10;
                    %nose = fid_new.pnt(:,3) < -10;
                    %fid_new.pnt(nose,:) = [];

                    fid_new.unit = 'mm';

                    % Labels:
                    fid_labels = fid_data{1};

                    D = fiducials(D,fid_new);
                    D.save;

                end
            end
        end

        
        %% OPT

        opt=[];
        opt.dirname = [studydir preproc_name '_opt/'];
        runcmd(['rm -rf ' opt.dirname]);

        clear spm_files_in structural_files_in;

        for ss=1:length(subjects_to_do), % iterates over subjects  

            f=subjects_to_do(ss);

            spm_files_in{ss} = prefix(spm_files{f},'');
            structural_files_in{ss} = structural_files{f};     

        end

        opt.datatype='ctf';
        opt.spm_files=spm_files_in;

        opt.downsample.do=1;
        opt.downsample.freq=Fs;
        
        % We will bandpass filter later, not here
        opt.highpass.do = 0;
        opt.highpass.cutoff = 0.1;
        opt.mains.do = 0;

        opt.africa.todo.ica=0;
        opt.africa.todo.ident=0;
        opt.africa.todo.remove=0;

        opt.africa.ident.func = @identify_artefactual_components_auto;
        opt.africa.ident.kurtosis_wthresh=0.2;
        opt.africa.ident.max_num_artefact_comps=2;
        opt.africa.precompute_topos=1;

        opt.bad_segments.do=1;
        opt.bad_segments.event_significance=0.05;

        opt.coreg.do=true;
        if whichstudy==1  
            opt.coreg.mri=structural_files_in;
        else
            opt.coreg.mri=cell(length(subjects_to_do),1);
        end
        opt.coreg.use_rhino=0;
        opt.coreg.useheadshape=0;
        opt.coreg.forward_meg='MEG Local Spheres';
        opt.epoch.do=0;
        opt.outliers.do = 0;
    
        % run:
        opt=osl_run_opt(opt);

        sessions_to_do = 1:(length(subjects_to_do)*length(subjects));
        if TESTRUN;sessions_to_do=1;end
        spm_links=cell(length(sessions_to_do),2);
        for ss=1:length(spm_files),

            S=[];
            S.D=opt.results.spm_files{sessions_to_do(ss)};
            S.outfile=[studydir preproc_name 'PreprocessedData/' '_session' num2str(ss)];
            Dnew=spm_eeg_copy(S);

            % write log of how original spm files are linked to copies
            spm_links{ss,1}=S.D;
            spm_links{ss,2}=S.outfile;

        end
        
        % remove the remaining build files:
        for ss=1:length(spm_files)
            delete(spm_files{ss});
            delete(strrep(spm_files{ss},'.mat','.dat'));
        end
        for ss = 1:length(opt.results.spm_files)
            delete([opt.results.spm_files{ss},'.mat']);
            delete([opt.results.spm_files{ss},'.dat']);
        end
    end
end