function [ spm_files_preproc,template_subj ] = prep_parcellated_data( settings )

try freq_range=settings.freq_range; catch, error('freq_range'); end
try parcellation_to_use=settings.parcellation.parcellation_to_use; catch, error('parcellation_to_use'); end

%settings.parcellation.parcellation_to_use='giles';
%settings.parcellation.parcellation_to_use='giles_3pcc';
%settings.parcellation.parcellation_to_use='giles_3pcc_ips';
%settings.parcellation.parcellation_to_use='aal_cort';

try sessions_to_do=settings.sessions_to_do; catch, error('sessions_to_do'); end
try dirname=settings.dirname; catch, error('dirname'); end
try spm_files=settings.spm_files; catch, error('spm_files'); end
try parc_dir=settings.parc_dir; catch, error('parc_dir'); end
try tmp=settings.parcellation.orthogonalisation; catch, settings.parcellation.orthogonalisation='innovations_mar'; end
try innovations_mar_order=settings.parcellation.innovations_mar_order; catch, innovations_mar_order=14; end
try num_iters=settings.signflip.num_iters; catch, num_iters=500; end
try num_embeddings=settings.signflip.num_embeddings; catch, num_embeddings=10; end
try modalities=settings.modalities; catch, modalities={'MEGGRAD'}; end

try do_signflip_diagnostics=settings.do_signflip_diagnostics; catch, do_signflip_diagnostics=settings.do_signflip; end

spm_files_preproc=[];

% setup some names
switch settings.parcellation.parcellation_to_use
    case 'test'
        parc_file= [parc_dir '/test'];
    case 'giles'
        parc_file= [parc_dir '/fmri_d100_parcellation_with_PCC_reduced_2mm'];
        parc_file = [parc_file '_ss5mm_ds8mm.nii.gz']; 
    case 'giles_3pcc'
        parc_file = [parc_dir '/fmri_d100_parcellation_with_3PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
    case 'giles_3pcc_ips'
        parc_file = [parc_dir '/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm.nii.gz'];
    case 'aal_cort'
        parc_file = [parc_dir '/aal2mni_cortical_4d_8mm.nii.gz'];
    case 'harvox_cort'
        parc_file = [parc_dir '/HarvOx_hemis_prob_2mm_8mm.nii.gz'];
    case 'sam49'
        parc_file = [parc_dir '/SParcellation_49_8mm_avw.nii.gz'];
    case 'dkc'
        parc_file = [parc_dir '/dk_cortical.nii.gz'];
end
parc_spmprefix=settings.parcellation.parcellation_to_use;

switch settings.parcellation.orthogonalisation
    case 'innovations_mar'
        parc_prefix   = [parc_spmprefix '_' settings.parcellation.orthogonalisation num2str(innovations_mar_order) '_'];
    otherwise
        parc_prefix   = [parc_spmprefix '_' settings.parcellation.orthogonalisation '_'];
end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filter and beamform

if settings.do_beamform

    bfdir = [settings.dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];

    % Copy into new dir, bandpass filter & beamform
    mkdir(bfdir);

    for ss = 1:length(sessions_to_do)

        session = sessions_to_do(ss);

        % Band pass filter:
        S      = [];
        S.D    = spm_files{session};
        S.band = 'bandpass';
        S.freq = freq_range;
        S.prefix = 'f';
        if all(S.freq >0)
            D=spm_eeg_filter(S);
        else
            D = spm_eeg_load(S.D)
            D = D.copy([dirname,'PreprocessedData/',S.prefix,D.fname]);
        end

        % move into bfdir
        Dnew=D.move(bfdir);

        bf_files=fullfile(Dnew);

        % Beamform:
        S                   = [];
        S.modalities        = modalities;
        S.timespan          = [0 Inf];
        S.pca_order         = 120;
        S.type              = 'Scalar';
        S.inverse_method    = 'beamform';
        S.prefix            = '';

        mni_coords        = osl_mnimask2mnicoords(fullfile(osldir,'/std_masks/MNI152_T1_8mm_brain.nii.gz'));
        osl_inverse_model(bf_files,mni_coords,S);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parcellate 

if settings.do_parcellation

   % setup bf file list
    bfdir = [dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];

    bf_files=[];
    for ss = 1:length(sessions_to_do)        
       session = sessions_to_do(ss);
       [~, fname]=fileparts(spm_files{session});
       bf_files{session}=[bfdir fname];
       bf_files{session}=prefix(bf_files{session},'f');
    end

    %%%%%%%%%%%%
    % parcellate the data

    parcellated_Ds=[];

    for ss = 1:length(sessions_to_do)
        session = sessions_to_do(ss);
        S                   = [];
        S.D                 = bf_files{session};
        S.parcellation      = parc_file;
        S.orthogonalisation = settings.parcellation.orthogonalisation;
        S.innovations_mar_order = innovations_mar_order;
        S.method            = 'spatialBasis';
        S.normalise_voxeldata = 0;
        S.prefix=parc_prefix;
        [parcellated_Ds{ss},parcelWeights,parcelAssignments] = osl_apply_parcellation(S);

        parcellated_Ds{ss}.parcellation.weights=parcelWeights;
        parcellated_Ds{ss}.parcellation.assignments=parcelAssignments;

        parcellated_Ds{ss}.save;
    end

    % D=spm_eeg_load(options.Ds{1}); tmp=std(D(:,:,1),[],2);fslview(nii_quicksave(tmp,'davve'))
    % tmp=std(newDs{1}(:,:,1),[],2);fslview(ROInets.nii_parcel_quicksave(tmp,parcelAssignments,'dave'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% enveloping

if settings.do_hilbert % needed for HMM on envelopes

    %%%%%%%%%
    % setup bf file list
    bfdir = [dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];

    bf_files=[];
    for ss = 1:length(sessions_to_do)        
       session = sessions_to_do(ss);
       [~, fname]=fileparts(spm_files{session});
       bf_files{session}=[bfdir fname];
       bf_files{session}=prefix(bf_files{session},'f');
    end

    %%%%%%%%%
    % setup parcellated_files list
    parcellated_files=[];
    for ss = 1:length(sessions_to_do)
        session = sessions_to_do(ss);
        parcellated_files{ss}=prefix(bf_files{session},parc_prefix);
    end

    for ss=1:length(sessions_to_do)

        subnum=sessions_to_do(ss);
        S=[];
        S.D = parcellated_files{ss};
        S.winsize = 1/40; %secs
        %S.winsize = 0; %secs
        S.downsample=0;
        S.remove_edge_effects=1;
        S.prefix  ='h';

        Dnew = osl_hilbenv(S);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sign flip stuff
parcellated_files = cell(length(spm_files),1);
if settings.do_signflip || do_signflip_diagnostics

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find template subject
    clear state_netmats_cov_preflipped bf_files parcellated_files;


    % setup bf file list
    bfdir = [dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];

    bf_files=[];
    for ss = 1:length(sessions_to_do)        
       session = sessions_to_do(ss);
       [~, fname]=fileparts(spm_files{session});
       bf_files{session}=[bfdir fname];
       bf_files{session}=prefix(bf_files{session},'f');
    end

    
    %%%%%%%%%
    % plot pre-flipped results

    parcellated_files=[];
    for ss = 1:length(sessions_to_do)
        session = sessions_to_do(ss);
        parcellated_files{ss}=prefix(bf_files{session},parc_prefix);
    end
    
    % input template, if one has been nominated:
    if isfield(settings,'templatesubj')
        parcellated_files{session+1} = settings.templatesubj;
    end
    
    %%%%%%%%%%
    % establish a good template subject
    S=[];
    S.concat = [];
    S.concat.protocol='none';
    S.concat.embed.do=1;
    S.concat.embed.num_embeddings=num_embeddings;
    S.concat.embed.rectify=false;
    S.concat.whiten=1;
    S.concat.normalisation='voxelwise';
    S.concat.pcadim=-1;
    S.netmat_method=@netmat_cov;

    state_netmats_cov_preflipped = hmm_full_global_cov( parcellated_files, S );
else
    state_netmats_cov_preflipped = [];
end

clear template_subj;

% assess which subject is the best template:
state_netmats=state_netmats_cov_preflipped;

modes={'none','abs'};
diag_offset=15;
if ~isfield(settings,'templatesubj');
    metric_global=zeros(length(state_netmats),length(state_netmats),length(modes));

    for mm=1:length(modes)
        for subj=1:length(state_netmats)
           for subj2=1:length(state_netmats)
                if subj2~=subj
                    metric_global(subj,subj2,mm)=matrix_distance_metric(state_netmats{subj}.global.netmat_full, state_netmats{subj2}.global.netmat_full,diag_offset,modes{mm},[]);
                end
           end
        end
    end

    tmp=sum(metric_global(:,:,2),2);

    [~, template_subj]=max(tmp);
else
    template_subj = length(parcellated_files);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  do actual sign flip

if settings.do_signflip

    bfdir = [dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];

    if 1

        % sign flipping settings
        S=[];
        S.roinets_protocol=settings.parcellation.orthogonalisation;
        S.innovations_mar_order = innovations_mar_order;            
        S.Ds=parcellated_files;
        S.num_iters=num_iters;
        S.prefix='sfold_';
        S.num_embeddings=num_embeddings;
        S.subj_template=template_subj;
        [ signflipped_files_out, sign_flip_results ] = find_sign_flips( S );

        sign_flip_results.signflipped_files=signflipped_files_out;
        sign_flip_results.energies_group=mean(sign_flip_results.energies,2);
        sign_flip_results.energies=sign_flip_results.energies(1:20:end,:);
        save([bfdir 'sign_flip_results_sfold_' parc_prefix],'-struct','sign_flip_results','-v7.3');

    else

        S=[];
        S.Ds = parcellated_files;
        S.prefix = 'sf_';
        S.options=[];
        S.options.maxlag=4; % max lag to consider, considering that we are including lagged autocovariance matrices in the calculation (default to 4).
        S.options.noruns=50; % how many random initialisations will be carried out (default to 50).
        S.options.maxcyc=200; % for each initialization, maximum number of cycles of the greedy algorithm (default to 100 * N * no. of channels).
        [Dsnew, flips, scorepath] = osl_resolve_sign_ambiguity(S);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot some sign flip diagnostics

if do_signflip_diagnostics

    bfdir = [dirname 'bfnew_' num2str(freq_range(1)) 'to' num2str(freq_range(2)) 'hz/'];
    sign_flip_results=load([bfdir 'sign_flip_results_sfold_' parc_prefix]);

    signflip_parcellated_files=[];
    for ss = 1:length(sessions_to_do)
        subnum = sessions_to_do(ss);
        signflip_parcellated_files{ss}=prefix(bf_files{subnum},['sfold_' parc_prefix]);
    end
    
    % input template, if one has been nominated:
    if isfield(settings,'templatesubj')
        signflip_parcellated_files{ss+1} = settings.templatesubj;
    end
    
    S=[];
    S.concat = [];
    S.concat.protocol=settings.parcellation.orthogonalisation;
    S.innovations_mar_order = innovations_mar_order;            
    S.concat.embed.do=1;
    S.concat.embed.num_embeddings=num_embeddings;
    S.concat.embed.rectify=false;
    S.concat.whiten=1;
    S.concat.normalisation='voxelwise';
    S.concat.pcadim=-1;
    S.netmat_method=@netmat_cov;

    [ state_netmats_cov_signflipped ] = hmm_full_global_cov( signflip_parcellated_files, S );

    subj_template_no=template_subj;

    print_fname=[bfdir '/sign_flip_plot'];

    plot_sign_flip_results(state_netmats_cov_preflipped,state_netmats_cov_signflipped, subj_template_no, freq_range, sign_flip_results, print_fname);

    disp('Sign flip diagonostic plot saved to:');
    disp(print_fname);

    spm_files_preproc=signflip_parcellated_files;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

