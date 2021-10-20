function [hmm, hmmfname, hmmoptions, settings_prepare, hmm_input_spm_files, hmm_input_epoched_spm_files] = run_full_hmm_pipeline_notts(S)

% [hmm, hmmfname, hmmoptions, settings_prepare] = run_full_hmm_pipeline(S)

try preproc_name=S.preproc_name; catch, error('S.preproc_name needed'); end
try session_name=S.session_name; catch, error('S.session_name needed'); end
try prep_sessions_to_do=S.prep_sessions_to_do; catch, error('S.prep_sessions_to_do needed'); end
try hmm_sessions_to_do=S.hmm_sessions_to_do; catch, hmm_sessions_to_do=prep_sessions_to_do; end
try spmfilesdir=S.spmfilesdir; catch, error('S.spmfilesdir needed'); end

try do_prepare=S.do_prepare; catch do_prepare=1; end
try do_hmm=S.do_hmm; catch do_hmm=1; end
try do_spectral_estimation=S.do_spectral_estimation; catch do_spectral_estimation=1; end

try parcellations_dir=S.parcellations_dir; catch parcellations_dir='/Users/woolrich/Dropbox/vols_scripts/hmm_misc_funcs/parcellations'; end
try num_embeddings=S.num_embeddings; catch num_embeddings  = 12; end
try hmm_name=S.hmm_name; catch hmm_name  = ''; end
try freq_range=S.freq_range; catch, freq_range=[1 45]; end

% build file names
clear spm_files epoched_spm_files;

session_name2 = 'PreprocessedData/';%(1:end-1);

for ss=1:length(prep_sessions_to_do)
    session=prep_sessions_to_do(ss);
    spm_files{ss}=[spmfilesdir preproc_name session_name2 '_session' num2str(session)];
    epoched_spm_files{ss}=[spmfilesdir preproc_name session_name2 '_session' num2str(session) '_epoched'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

settings_prepare=[];

settings_prepare.spm_files=spm_files;
settings_prepare.sessions_to_do=1:length(spm_files);

settings_prepare.freq_range=freq_range;

% setup parcellation
settings_prepare.parc_dir=parcellations_dir;

%settings_prepare.dirname=[spmfilesdir session_name ];%'_' session_name '/'];
settings_prepare.dirname=[strrep(spmfilesdir,session_name,''), session_name , preproc_name ];%'_' session_name '/'];

%settings_prepare.parcellation.parcellation_to_use='test';
settings_prepare.parcellation.parcellation_to_use='giles';

%settings_prepare.parcellation.orthogonalisation='innovations_mar';
%settings_prepare.parcellation.innovations_mar_order=14;
settings_prepare.parcellation.orthogonalisation='symmetric';

settings_prepare.signflip.num_iters=1500;
settings_prepare.signflip.num_embeddings=num_embeddings;
if isfield(S,'signfliptemplatesubj')
    settings_prepare.templatesubj=S.signfliptemplatesubj;
end
    
settings_prepare.do_beamform=1;
settings_prepare.do_parcellation=1;
D = spm_eeg_load(spm_files{1});
%D.fsample = 250;
if D.fsample < 600
    settings_prepare.do_signflip=1;
    settings_prepare.do_signflip_diagnostics=1;
else
    settings_prepare.do_signflip=0;
    settings_prepare.do_signflip_diagnostics=0;
    settings_prepare.freq_range = [0 0]; % this a hack to avoid any filtering
end
settings_prepare.do_hilbert=0;

if do_prepare
    [~,templatesubj] = prep_parcellated_data( settings_prepare );
    settings_prepare.templatesubj = templatesubj;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup HMM

%hmm_mode='envelope';
hmm_mode='raw';

% build file names for signflip_parcellated_files for HMM input:
parc_spmprefix=settings_prepare.parcellation.parcellation_to_use;
switch settings_prepare.parcellation.orthogonalisation
    case 'innovations_mar'
        parc_prefix   = [parc_spmprefix '_' settings_prepare.parcellation.orthogonalisation num2str(settings_prepare.parcellation.innovations_mar_order) '_'];
    otherwise
        parc_prefix   = [parc_spmprefix '_' settings_prepare.parcellation.orthogonalisation '_'];
end

signflipped_files={};
bf_files={};
enveloped_files={};

bfdir = [settings_prepare.dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];

bf_files=[];
for ss = 1:length(hmm_sessions_to_do)        
   session = hmm_sessions_to_do(ss);
   [~, fname]=fileparts(spm_files{session});
   bf_files{ss}=[bfdir fname];
   bf_files{ss}=prefix(bf_files{ss},'f');
   signflipped_files{ss}=prefix(bf_files{ss},['sfold_' parc_prefix]);
   unsignflipped_files{ss}=prefix(bf_files{ss},[parc_prefix]);
   
   enveloped_files{ss}=prefix(bf_files{ss},['h' parc_prefix]);
   
   hmm_input_epoched_spm_files{ss}=epoched_spm_files{session};
end    

% check signflipped_files exist
signflipped_files_ext = cellfun(@strcat,signflipped_files,repmat({'.mat'},1,length(signflipped_files)),'uniformoutput',0);
isfile = cell2mat(cellfun(@exist,signflipped_files_ext,repmat({'file'},1,length(signflipped_files)),'uniformoutput',0));

% check enveloped_files exist
%enveloped_files_ext = cellfun(@strcat,enveloped_files,repmat({'.mat'},1,length(enveloped_files)),'uniformoutput',0);
%isfile = cell2mat(cellfun(@exist,enveloped_files_ext,repmat({'file'},1,length(enveloped_files)),'uniformoutput',0));

disp(['Number of sessions is ' num2str(length(isfile))]);

if (any(~isfile))
    warning('Invalid signflipped_files');
end


switch hmm_mode   
    case 'raw'    
        if D.fsample==600
            hmm_input_spm_files=unsignflipped_files;
        else
            hmm_input_spm_files=signflipped_files;
        end

    case 'envelope'    
        hmm_input_spm_files=enveloped_files;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HMM

hmmoptions=[];

switch hmm_mode
    
    case 'raw'
        % use tideh (time delay embedded HMM) on raw data        

        % dimensionality reduction
        D=spm_eeg_load(hmm_input_spm_files{1});
        %D = spm_eeg_load
        nparcels=size(D.parcellation.weights,2);
        hmmoptions.prepare.pcadim           = 80;

        % embedding 
        hmmoptions.prepare.embed.do         = 1;
        hmmoptions.prepare.embed.num_embeddings = num_embeddings;
           
    case 'envelope'        
        
        % no dimensionality reduction
        hmmoptions.prepare.pcadim           = 0;
       
        % no embedding this time
        hmmoptions.prepare.embed.do         = 0;
        hmmoptions.prepare.embed.num_embeddings = 0;
end

% K defines the number of states
if isfield(S,'K')
    hmmoptions.hmm.K                    = S.K;
else
    hmmoptions.hmm.K                    = 12;
end

if isfield(S,'templateHMM') 
    temp = load(S.templateHMM);
    hmmoptions.hmm.hmm              = temp.hmm;
    hmmoptions.hmm.updateObs = 0;
    hmmoptions.hmm.BIGcyc = 1; % this the critical param: it means there is just one call to hmmdecode
end

hmmoptions.prepare.normalisation    = 'voxelwise';
hmmoptions.prepare.whiten           = 1; 
hmmoptions.prepare.savePCmaps       = 0;
%hmmoptions.prepare.max_ntpts        = 40000;

hmmoptions.hmm.dynamic_model_type   = 'hmm';
if isfield(S,'dynamic_model_type')
    hmmoptions.hmm.dynamic_model_type   = S.dynamic_model_type;
end
%hmmoptions.hmm.dynamic_model_type   = 'vbrt';

hmmoptions.hmm.initcyc              = 60;
hmmoptions.hmm.initrep              = 4;

hmmoptions.hmm.big                  = 1;
hmmoptions.hmm.BIGNbatch            = 10;
hmmoptions.hmm.name                 = hmm_name;

hmmoptions.output.method       = 'pcorr';
hmmoptions.output.use_parcel_weights = 0;
hmmoptions.output.assignment   = 'hard';
 
% setup filenames
 
hmmoptions.prepare.filename     = ['hmm' hmmoptions.hmm.name ...
                                '_parc_' parc_prefix ...                                
                                '_pcdim' num2str(hmmoptions.prepare.pcadim) ...
                                '_' num2str(hmmoptions.prepare.normalisation) ...
                                '_embed' num2str(hmmoptions.prepare.embed.num_embeddings)];
hmmoptions.hmm.filename        = [hmmoptions.prepare.filename ...
                                '_K' num2str(hmmoptions.hmm.K)...
                                '_big' num2str(hmmoptions.hmm.big)...
                                '_dyn_model' hmmoptions.hmm.dynamic_model_type];
hmmoptions.output.filename     = [hmmoptions.hmm.filename '_output'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run HMM
hmm = [];hmmfname = [];
if do_hmm==1

    hmmoptions.todo.prepare  = 1;
    hmmoptions.todo.hmm      = 1;
    hmmoptions.todo.output   = 0;                                                                       
    
    hmmoptions.hmmdir = [settings_prepare.dirname 'hmm_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];
    [HMMresults_raw_flips] = teh_groupinference_parcels(hmm_input_spm_files,hmmoptions);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load previously run HMM

    hmmdir=[settings_prepare.dirname 'hmm_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];
    load([hmmdir hmmoptions.hmm.filename]);

    hmmfname=[hmmdir hmmoptions.hmm.filename];

    % add some settings to hmm for use later

    D=spm_eeg_load(hmm.data_files{1});
    hmm.parcellation= D.parcellation;
    hmm.parcellation.file = strrep(hmm.parcellation.S.parcellation,'/Users/woolrich/Dropbox/vols_scripts/hmm_misc_funcs',...
        osldir);

    sres=nii.get_spatial_res(hmm.parcellation.file);
    gridstep=sres(1);

    hmm.parcellation.mask=[osldir '/std_masks/MNI152_T1_' num2str(gridstep) 'mm_brain'];

    save([hmmdir hmmoptions.hmm.filename],'hmm');
end
%%

if do_spectral_estimation
    
    hmmdir=[settings_prepare.dirname 'hmm_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];
    load([hmmdir hmmoptions.hmm.filename]);

    hmmfname=[hmmdir hmmoptions.hmm.filename];
    
    storage_dir=[hmmfname '_store/'];

    mkdir(storage_dir);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Estimate spectra and cross spectra using multitaper on concatenated data

    do_run=1;

    S=[];
    S.parcellated_filenames=hmm.data_files;
    S.normalisation='voxelwise';
    S.assignment='hard';
    S.global_only=false;
    S.embed.do=0;
    S.embed.rectify=false;

    S.netmat_method=@netmat_spectramt;
    S.netmat_method_options.fsample=hmm.fsample;
    S.netmat_method_options.fband=freq_range;
    S.netmat_method_options.type='coh';
    S.netmat_method_options.full_type='full';
    S.netmat_method_options.var_normalise=false;
    S.netmat_method_options.reg=2; % higher is less reg
    S.netmat_method_options.order=0;

    if do_run

        [ state_netmats_mt ] = hmm_state_netmats_teh_concat( hmm, S );

        save([storage_dir '/state_netmats_mt' num2str(floor(S.netmat_method_options.reg)) ...
            '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment '_' ...
            'global' num2str(S.global_only)], '-v7.3', 'state_netmats_mt');
    else
        load([storage_dir '/state_netmats_mt_' num2str(floor(S.netmat_method_options.reg)) ...
        '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment  '_' ...
        'global' num2str(S.global_only)],'state_netmats_mt');

    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Estimate spectra and cross spectra using multitaper on each session separately

    do_run=1;

    S=[];
    S.parcellated_filenames=hmm.data_files;
    S.normalisation='voxelwise';
    S.assignment='hard';
    S.global_only=false;
    S.embed.do=0;
    S.embed.rectify=false;

    S.netmat_method=@netmat_spectramt;
    S.netmat_method_options.fsample=hmm.fsample;
    S.netmat_method_options.fband=freq_range;
    S.netmat_method_options.type='coh';
    S.netmat_method_options.full_type='full';
    S.netmat_method_options.var_normalise=false;
    S.netmat_method_options.reg=2;
    S.netmat_method_options.order=0;

    if do_run

        [ state_netmats_mtsess ] = hmm_state_netmats_teh( hmm, S );

        save([storage_dir '/state_netmats_mtsess_' num2str(floor(S.netmat_method_options.reg)) ...
            '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment '_' ...
            'global' num2str(S.global_only)], '-v7.3', 'state_netmats_mtsess');

    else

        load([storage_dir '/state_netmats_mtsess_' num2str(floor(S.netmat_method_options.reg)) ...
        '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment  '_' ...
        'global' num2str(S.global_only)],'state_netmats_mtsess');

    end
end

end

