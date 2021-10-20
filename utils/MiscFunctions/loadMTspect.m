function [psd,coh,f] = loadMTspect(basedir,K,templatestring,permuteSTCs);
if nargin<4
    permuteSTCs = false; % note this is a flag for testing whether the estimation is biased; if true, it scrambles the STCs
end
if permuteSTCs
    PMstring = '_permuted';
else
    PMstring = '';
end
    
hard = 0; % 1 for hard state assignment, 0 for soft
if hard
     fulldir = [basedir, '/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm_store/',...
         'state_netmats_mtsess_2_vn0_hard_global0.mat'];

    load(fulldir);

    state_netmats=state_netmats_mtsess;
    NK=length(state_netmats{1}.state);
    num_nodes=size(state_netmats{1}.state{1}.netmat,1);
    num_freqs=length(state_netmats{1}.state{1}.spectramt.f);
    nsubjects=length(state_netmats);

    hmmfile = [basedir,'/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
    load(hmmfile,'new_state_ordering');
    
    psd=zeros(nsubjects,NK,num_freqs,num_nodes,num_nodes);
    coh = zeros(nsubjects,NK,num_freqs,num_nodes,num_nodes);
    for ss=1:length(state_netmats)
        for kk=1:length(state_netmats{1}.state)           
            psd(ss,kk,:,:,:)=state_netmats{ss}.state{new_state_ordering(kk)}.spectramt.psd; 
            coh(ss,kk,:,:,:)=state_netmats{ss}.state{new_state_ordering(kk)}.spectramt.coh; 
        end
        % add global on the end
        %psds(ss,kk+1,:,:,:)=state_netmats{ss}.global.spectramt.psd;
    end

    f = state_netmats{1}.state{1}.spectramt.f;
else
    fulldir = [basedir, '/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm_store/',...
        PMstring,'state_netmats_mtsess_2_vn0_soft_global0.mat'];
    if exist(fulldir)
        load(fulldir);
    else
        disp('Soft state spectra not found; recomputing for all data now:')
        
        % load hmm:
        hmmfile = [basedir,'/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
        load(hmmfile);
        
        mtfilename = [basedir, '/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm_store/',...
            PMstring,'state_netmats_mtsess_2_vn0_soft_global0.mat'];
        
        S=[];
        S.parcellated_filenames=hmm.data_files;
        S.normalisation='voxelwise';
        S.assignment='soft';
        S.global_only=false;
        S.embed.do=0;
        S.embed.rectify=false;

        % S.netmat_method=@netmat_spectramt;
        % S.netmat_method_options.fsample=hmm.fsample;
        % S.netmat_method_options.fband=freq_range;
        % S.netmat_method_options.type='coh';
        % S.netmat_method_options.full_type='full';
        % S.netmat_method_options.var_normalise=false;
        % S.netmat_method_options.reg=2;
        % S.netmat_method_options.order=0;

        logtrans=0;
        assignment=S.assignment;

        normalisation = S.normalisation;
        embed=S.embed;

        % and Diego's code setup:
        fitmt_subj = cell(length(hmm.data_files),1);
        %d = length(options_mt.embeddedlags) - 1; 
        
        Hz=250;
        options_mt = struct('Fs',Hz); % Sampling rate - for the 25subj it is 300
        options_mt.fpass = [1 45];  % band of frequency you're interested in
        options_mt.tapers = [4 7]; % taper specification - leave it with default values
        options_mt.p = 0; %0.01; % interval of confidence  
        options_mt.win = 2 * Hz; % multitaper window
        options_mt.to_do = [1 0]; % turn off pdc
        options_mt.order = 0;
        options_mt.embeddedlags = [0];%-7:7;

        for subnum=1:length(hmm.data_files)
            disp(['Computing for subj ' num2str(subnum)]);
            Dfname = hmm.data_files{subnum};
            if contains(Dfname,'WooliePipeline')
                Dfname = strrep(Dfname,'/Users/chiggins/data/YunzheData/Replaydata4Cam/WooliePipeline/spm/sept2019_eo',basedir);
            end
            Dp = spm_eeg_load(Dfname);
            embed.tres=1/Dp.fsample;

            datap = osl_teh_prepare_data(Dp,normalisation,logtrans,[],embed);
            datap=datap';
            nparcels=Dp.nchannels;
            hmm_sub = hmm; 
            hmm_sub.statepath = hmm.statepath(hmm.subj_inds==subnum); 
            hmm_sub.gamma = hmm.gamma(hmm.subj_inds==subnum,:);
            if permuteSTCs
                hmm_sub.gamma = [flipud(hmm_sub.gamma),ones(length(hmm_sub.gamma),1)];
            end
            fit = hmmspectramt(datap,length(datap),hmm_sub.gamma,options_mt);
            if subnum==1
                f = fit.state(1).f;
                psd=zeros(length(hmm.data_files),hmm_sub.K,length(f),nparcels,nparcels);
                coh=zeros(length(hmm.data_files),hmm_sub.K,length(f),nparcels,nparcels);
            end
            for jj = 1:length(fit.state)
                psd(subnum,jj,:,:,:) = fit.state(jj).psd;
                coh(subnum,jj,:,:,:) = fit.state(jj).coh;
            end
            clear fit

        
        end
        dir1 = [basedir, '/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm_store/'];
        if ~isfolder(dir1)
            mkdir(dir1);
        end
        save( mtfilename ,'psd','coh','f')
    end
    hmmfile = [basedir,'/hmm_1to45hz/hmm',templatestring,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
    load(hmmfile,'new_state_ordering');
    psd_temp = psd;
    coh_temp = coh;
    for k=1:length(new_state_ordering)
        psd(:,k,:,:,:) = psd_temp(:,new_state_ordering(k),:,:,:);
        coh(:,k,:,:,:) = coh_temp(:,new_state_ordering(k),:,:,:);
    end
end
end