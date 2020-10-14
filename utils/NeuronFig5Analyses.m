%% Replay Paper Figure 5: High Frequency Oscillations

if whichstudy==1
    CanonicalRSN = false;
else
    CanonicalRSN = true; % canonical refers to RSNs trained on noittingham data and refitted
end

% add osl_braingraph:
addpath(genpath([osldir,'/osl-core/old']));

% Define colours to use in state plots
color_scheme = set1_cols();

% Define sample rate
sample_rate = 250;

if CanonicalRSN
    template_string = [int2str(bestmodel),'usingtemplate'];
else
    template_string = [int2str(bestmodel),''];
end

parc_name='Giles';%'NatComms_wHippocampus';DiegoNatComms, or Giles


studydir = [wd,session_name{whichstudy},'250Hz/'];


% Load in run indices:
load([studydir,'hmm_1to45hz/hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'],'hmmT','subj_inds');
hmmdir = [studydir,'hmm_1to45hz/'];
load([hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat']);
hmm = hmm_permutestates(hmm,new_state_ordering);
Gamma = hmm.gamma;
scan_T = cell2mat(hmmT);

savebase = fullfile( [wd,session_name{whichstudy},'Fig5_',template_string,'/' ])
if ~exist( savebase )
    mkdir(savebase);
end

parc_file = ['fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
parc = parcellation(parc_file);
statelabels={'1','2','3','4','5','6','7','8','9','10','11','12'};

%enter directories of high frequency raw data:
if whichstudy==2
    datapath600Hz = [wd,session_name{whichstudy},'600Hz/bfnew_0to0hz/'];
elseif whichstudy==3
    datapath600Hz = [wd,session_name{whichstudy},'600Hz/bfnew_0to0hz/'];
elseif whichstudy==1
    
end


%% First need to process 600Hz data:
outfile600Hz = fullfile( datapath600Hz, 'prepared_data.mat' )
if ~exist(outfile600Hz)
    
    % We'll need to collect the data, T and epoch information
    data = [];               % HMM ready dataset
    T = [];                  % Length of continuous good segments
    R = [];                  % Indices of single run within data
    B = cell(nscans{whichstudy},1);      % Indices of bad samples per session
    trl = cell(nscans{whichstudy},1);    % epoch info per segment
    runlen = zeros(nscans{whichstudy},1);          % Length of run per good segment

    % Define HMM folder 
    hmm_folder600Hz = datapath600Hz;
    
    if ~exist( hmm_folder600Hz )
        mkdir( hmm_folder600Hz );
    end
    outfile600Hz = fullfile( datapath600Hz, 'prepared_data' );
    
    % need to re-order names:
    fnames = 'giles_symmetric_f_session';
    
    for ind = 1:nscans{whichstudy}

        fprintf(['\nProcessing ',fnames, int2str(ind)]);
        %-------------------------------
        % continuous file
        D_600 = spm_eeg_load([hmm_folder600Hz, fnames, int2str(ind)]);
        
        runlen(ind) = size(D_600,2);

        %-------------------------------
        % get data and orthogonalise
        dat = D_600(:,:,1);
        %dat = ROInets.remove_source_leakage(dat, 'symmetric');

        %-------------------------------
        % Get badsamples - Note these have already been synchronised to match low freq bad segments:
        runlen(ind) = size(dat,2);
        bs = ~good_samples( D_600 );

        % find single good samples - bug when we have consecutive bad segments
        xx = find(diff(diff(bs)) == 2)+1;
        if ~isempty(xx)
            bs(xx) = 1;
        end

        % store bad samples
        B{ind}=find(bs);

        % indices of good samples
        good_inds=setdiff(1:runlen(ind),B{ind});

        % remove bad samples,
        dat = dat(:,good_inds);

        if any(bs)
            t_good = ~bs;
            db = find(diff([0; t_good(:); 0]));
            onset = db(1:2:end);
            offset = db(2:2:end);
            t = offset-onset;

            % sanity check
            if size(dat,2) ~= sum(t)
                disp('Mismatch between Data and T!!');
            end
        else
            t = size(dat,2);
        end

        %--------------------------------
        % Store info

        offset = sum(T);

        R = cat(1,R,[offset+1 offset+size(dat,2)]);

        T = cat(1,T,t);

        data = cat(2,data,dat);

    end
    save( outfile600Hz, 'data', 'R', 'T', 'B', 'runlen', '-v7.3' );
else
    load( outfile600Hz, 'data', 'R', 'T', 'B', 'runlen');
end
data=data';
T600 = T;

%% now load HMM model and realign:
datadir = [studydir,'bfnew_1to45hz/'];


if whichstudy==3
    if 1;~isfile([wd,'GenericReplayData/STUDYII_ReplayOnset/600HzMasks.mat'])
        [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datapath600Hz);
        save([wd,'GenericReplayData/STUDYII_ReplayOnset/600HzMasks.mat'],'maskA','maskB','triggerpoints','goodsamples')
    end
    [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);
    masks_100Hz = load([wd,'GenericReplayData/STUDYII_ReplayOnset/100HzMasks.mat'],'maskA');
    %[maskA,maskB,triggerpoints,goodsamples] = getSubjectRestingMasks_quinn(s600);
    masks_600Hz = load([wd,'GenericReplayData/STUDYII_ReplayOnset/600HzMasks.mat'],'maskA','maskB','triggerpoints','goodsamples');
    badsub = zeros(nscans{whichstudy},1);
elseif whichstudy==2
    if 1;~isfile([wd,'GenericReplayData/STUDYI_ReplayOnset/600HzMasks.mat'])
        [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datapath600Hz);
        save([wd,'GenericReplayData/STUDYI_ReplayOnset/600HzMasks.mat'],'maskA','maskB','triggerpoints','goodsamples')
    end
    [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);
    masks_600Hz = load([wd,'GenericReplayData/STUDYI_ReplayOnset/600HzMasks.mat'],'maskA','maskB','triggerpoints','goodsamples');
    masks_100Hz = load([wd,'GenericReplayData/STUDYI_ReplayOnset/100HzMasks.mat'],'maskA');
    badsub = zeros(nscans{whichstudy},1);
end


%% Align data:

newdata600Hz = [];
Gamma_600 = [];
R_new = zeros(size(R));
for iSj=1:nscans{whichstudy}
    if iSj==1
        R_new(1,1)=1;
    else
        %R_new(iSj,1) = max(R_new(:,2))+1;
        R_new(iSj,1) = R_new(iSj-1,2)+1;
    end
    fprintf(['Processing data from subj ',int2str(iSj),'\n']);
    dattemp = zeros(length(masks_600Hz.goodsamples{iSj}),size(data,2));
    dattemp(masks_600Hz.goodsamples{iSj},:) = data(R(iSj,1):R(iSj,2),:);
    if ~ismember(iSj,find(badsub))
        D = spm_eeg_load(hmm.data_files{iSj-sum(badsub(1:iSj))});
        GS = good_samples(D);
        GamTemp = zeros(length(GS),size(Gamma,2));
        GamTemp(GS,:) = Gamma(hmm.subj_inds==iSj-sum(badsub(1:iSj)),:);
        GamNew = resample(GamTemp,600,250);
        % correct odd edge effects:
        GamNew(GamNew>1)=1;
        GamNew(GamNew<0)=0;
        % also correct rounding errors in size:
        p = min([length(dattemp),length(GamNew)]);
        debugrec1(iSj,[1:2]) = [length(dattemp),length(GamNew)];
        if length(dattemp)-length(GamNew)>5
            error('Something wrong in alignment code');
            badsubs(iSj)=1;
            R_new(iSj,2)=0;
        else
            Gamma_600 = [Gamma_600;GamNew(1:p,:)];
            newdata600Hz = [newdata600Hz;dattemp(1:p,:)];
            R_new(iSj,2) = R_new(iSj,1)-1+p;
        end
    else
        R_new(iSj,2) = R_new(iSj-1,2);
    end
end

R = R_new(~badsub,:);

if size(Gamma_600,1) ~= size(newdata600Hz,1)
	warning('The size of data and Gamma do not match');
end


%% Learn high frequency spectra associated with each RSN-state:

% These options specify the how the spectra will be computed. Full details can
% be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#spectra

spec_options = struct();
spec_options.fpass = [1 160];
spec_options.p = 0; % no confidence intervals
spec_options.to_do = [0 0]; % no coherence - set to [1 0] to do coherence but not pdc 
spec_options.embeddedlags = [0];
spec_options.Fs = 600;%D.Fsample;
spec_options.win = 2*spec_options.Fs;
spec_options.Nf = 256;
nparcels = parc.n_parcels;nstates=hmm.K;

mt_outfile = fullfile( savebase, sprintf('embedded_HMM_K%d_spectra',hmm.K));
    
if ~isfile(mt_outfile)
    for iSubj=1:length(R)
        fprintf(['Computing state spectra for subj ',int2str(iSubj),'\n']);
        subj_data = newdata600Hz(R(iSubj,1):R(iSubj,2),:);

        fit = hmmspectramt(subj_data,R(iSubj,2)-R(iSubj,1),Gamma_600(R(iSubj,1):R(iSubj,2),:),spec_options);
        if iSubj==1
            nres = size(fit.state(1).psd,1);
            psd = zeros(length(R),nstates,nres,nparcels,nparcels);
        end

        for jj = 1:nstates
            psd(iSubj,jj,:,:,:) = fit.state(jj).psd;
        end
        clear fit
    end
    % Save the MT outputs
    save( mt_outfile ,'psd','-v7.3')
else
    load( mt_outfile ,'psd')
end

%% and plot:
nparcels = 38;
I = find(eye(nparcels));

f_band = linspace(spec_options.fpass(1),spec_options.fpass(2),size(psd,3));
Subj_psd = mean(psd(:,:,:,I),4);
if whichstudy>1
    % account for each subject's two runs:
    Subj_psd = squeeze(mean(reshape(Subj_psd,[2,nSj{whichstudy},size(Subj_psd,2),size(Subj_psd,3)]),1));
end
PSD_to_plot = squeeze(mean(Subj_psd,1));
PSD_ste = squeeze(std(Subj_psd,[],1)./sqrt(size(Subj_psd,1)));

figure('Position',[440 430 535 362]);
for k=hmm.K:-1:1
    if mod(k,2)==0,ls='--',else,ls='-',end
    semilogy(f_band,PSD_to_plot(k,:),'LineWidth',2,'Color',color_scheme{k},'LineStyle',ls);hold on;
    shadedErrorBar(f_band,PSD_to_plot(k,:),PSD_ste(k,:),{'LineWidth',2,'Color',color_scheme{k},'LineStyle',ls},1);hold on;
    h(k) = plot(NaN,NaN,'Color',color_scheme{k},'LineWidth',2,'LineStyle',ls);
end
for k=1:K,h(k).DisplayName=['RSN-State ',int2str(k)];end
leg=legend(h,'Location','EastOutside');
plot4paper('Frequency (Hz)','PSD');
YL = ylim();
grid on;grid minor;grid minor;

xlim([45,160]);
ylim([1e-3,3e-3]);
print([savebase, 'PSD_per_state','xaxislimited'],'-depsc');


%% work out statistical significance:

f_analyse = find(f_band>52 & f_band<148);
for ifreq=1:length(f_analyse)
    pval_anova(ifreq) = anova1(Subj_psd(:,:,f_analyse(ifreq)));
    pval_anova_noDMN(ifreq) = anova1(Subj_psd(:,[1,3:12],f_analyse(ifreq)));
    close all;
end

%% and spatially localising:

Subj_to_do = size(psd,1);%[1:(2*nSj)];
Subj_psd = squeeze(mean(psd(Subj_to_do,:,:,I),1));
for i=1:3
    HFOpower(:,:,1) = squeeze(mean(Subj_psd(:,f_band>52 & f_band<98,:),2));
    HFOpower(:,:,2) = squeeze(mean(Subj_psd(:,f_band>102 & f_band<148,:),2));
end

for k=1:hmm.K
    toplot=HFOpower(k,:,2);
    %toplot(toplot<prctile(toplot,75))=NaN;
    plot_surf_summary_amended(parc,toplot)
    str = (['State ',int2str(k), ', power 100-150Hz']);
    a = annotation('textbox',[0.5,0.5,0.5,0.5],'String',str,'FitBoxToText','on');
    set(a,'FontSize',15);
    set(a,'Position',[0.3720 0.8555 0.1875 0.1094]);
    print([savebase 'PSD_K',int2str(k),'_102to148Hz'],'-depsc')
end


%% And investigate spectra around epoched replay events:
% 
% outfile600Hz = fullfile(hmm_folder600Hz, 'embedded_hmm_data.mat' );
% load( outfile600Hz, 'data', 'R', 'T', 'B', 'runlen');
% data = data';
% T600 = T;

outfile600Hz = fullfile( datapath600Hz, 'prepared_data.mat' )
load( outfile600Hz, 'data', 'R', 'T', 'B', 'runlen');
data = data';
T600 = T;
%% Need to convert replay times into a state timecourse, much like Gamma

%addpath(genpath('/Users/chiggins/Google Drive/MATLAB/3.0 Yunzhe Exp'));
%s=study([config.analysisdir,'spm_sss_600Hz_processed'],parc_name);
%[~,~,triggerpoints,goodsamples] = getSubjectRestingMasks_quinn(s);

% Load replay scores
if whichstudy==2
    load([wd,'GenericReplayData/STUDYI_ReplayOnset/StrReplayOnset'],'ToRall'); % the replay scores for first resting state
    replayScores(1:2:(2*21),:) = ToRall;
    load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
    replayScores(2:2:(2*21),:) = ToRall;
else
    load([wd,'GenericReplayData/STUDYII_ReplayOnset/STUDYII_PreplayOnset'],'ToRall'); % the replay scores for first resting state
    replayScores(1:2:(2*22),:) = ToRall;
    load([wd,'GenericReplayData/STUDYII_ReplayOnset/STUDYII_ReplayOnset'],'ToRall');
    replayScores(2:2:(2*22),:) = ToRall;
end
%
t_window=10;

DataAll=[];T_cropped=[];GammaAll=[];
for iSj=1:nscans{whichstudy}    
    fprintf(['\nProcessing scan ',int2str(iSj)]);
    dattemp = zeros(length(masks_600Hz.goodsamples{iSj}),size(data,2));
    dattemp(masks_600Hz.goodsamples{iSj},:) = data(R(iSj,1):R(iSj,2),:);
    DataSjCropped = dattemp(masks_600Hz.triggerpoints{iSj},:);
    
    % epoch by replay onset:
    Gamma_sj = zeros(size(DataSjCropped,1),2);
    for i_gam=1:2
        if i_gam==1
             replaytimes = replayScores(iSj,:) > prctile(replayScores(iSj,:),99);
        else
            Sj_perm = randi(nscans{whichstudy});
            replaytimes = replayScores(Sj_perm,:) > prctile(replayScores(Sj_perm,:),99);
        end
        t_i = find(replaytimes);
        Q = length(DataSjCropped) / length(replaytimes);
        t_g = round(t_i*Q);

        if any(t_g>size(DataSjCropped,1)-t_window)
            t_g(t_g>size(DataSjCropped,1)-t_window) = [];
            t_i(t_i>size(DataSjCropped,1)-t_window) = [];
        end
        
        if any(t_g-2*t_window<0)
            t_g(t_g-2*t_window<0) = [];
            t_i(t_g-2*t_window<0) = [];
        end
        % omit any replay events for which data is NaN (bad segments in
        % different pipeline)
        to_omit=[];
        for i=1:length(t_g)
            if any(isnan(DataSjCropped((t_g(i)-2*t_window):(t_g(i)+2*t_window))))
                to_omit=[to_omit,i];
            end
        end
        t_g(to_omit) = [];

        % create faux state timecourse around replay events
        for i=1:length(t_g)
              Gamma_sj((t_g(i)-t_window):(t_g(i)+t_window),i_gam)=1;
        end
    end

    % also control for bad samples:
    Gamma_sj(~masks_600Hz.goodsamples{iSj}(masks_600Hz.triggerpoints{iSj}),:)=0;
    
    % omit bad samples:
    BS = isnan(DataSjCropped(:,1));
    DataSjCropped(BS,:) = [];
    Gamma_sj(BS,:) = [];
    T_cropped = [T_cropped,size(DataSjCropped,1)];
    
    DataAll = [DataAll;DataSjCropped];
    GammaAll = [GammaAll;Gamma_sj];
end

%% and run multitaper:
t_window=10;
mt_outfile = fullfile( savebase, ['Replay_evoked_spectra_twindow_',int2str(t_window)]);
    
spec_options = struct();
spec_options.fpass = [1 160];
spec_options.p = 0; % no confidence intervals
spec_options.to_do = [0 0]; % no coherence - set to [1 0] to do coherence but not pdc 
spec_options.embeddedlags = 0;
spec_options.Fs = 600;%D.Fsample;
spec_options.win = 2*spec_options.Fs;
spec_options.Nf = 256;

nparcels = size(DataAll,2);
nstates = size(GammaAll,2);
if ~isfile(mt_outfile)
    for iSubj=1:length(T_cropped)
        fprintf(['Computing replay spectra for Sj ',int2str(iSubj),'\n']);
        
        t_in = (sum(T_cropped(1:(iSubj-1))) +1):(sum(T_cropped(1:iSubj)));
        subj_data = DataAll(t_in,:);
        subj_Gamma = GammaAll(t_in,:);
        
        fit = hmmspectramt(subj_data,T_cropped(iSubj),subj_Gamma,spec_options);
        I = logical(eye(nparcels));
        if iSubj==1
            nres = size(fit.state(1).psd,1);
            psd = zeros(length(R),nstates,nres,nparcels);
        end
        for jj = 1:nstates
            psd(iSubj,jj,:,:,:) = fit.state(jj).psd(:,I);
        end
        if any(isnan(psd(:)))
            error('Problem in psd spectra');
        end
        clear fit
    end
    % Save the MT outputs
    save( mt_outfile ,'psd','-v7.3')
else
    load( mt_outfile ,'psd')
end

%% and plot replay-evoked PSD:
% psd is subj x states x frequencies x parcels

I = find(eye(nparcels));
f_band = linspace(spec_options.fpass(1),spec_options.fpass(2),size(psd,3));
Subj_psd = mean(psd(:,:,:,:),4);
Subj_psd = squeeze(mean(reshape(Subj_psd,[2,nSj{whichstudy},size(Subj_psd,2),size(Subj_psd,3)]),1));
PSD_to_plot = squeeze(mean(Subj_psd,1));

PSD_ste = squeeze(std(Subj_psd,[],1)./sqrt(nSj{whichstudy}));
figure('Position',[440 430 686 362]);
semilogy(f_band,(PSD_to_plot),'LineWidth',2);
colors_replay = {0.05*[1 1 1],0.55*[1 1 1]};
clear h;
for k=1:2
    semilogy(f_band,PSD_to_plot(k,:),'LineWidth',2,'Color',colors_replay{k});hold on;
    shadedErrorBar(f_band,PSD_to_plot(k,:),PSD_ste(k,:),{'LineWidth',2,'Color',colors_replay{k}},1);hold on;
    h(k) = plot(NaN,NaN,'Color',colors_replay{k},'LineWidth',2);
end
h(1).DisplayName='Replay evoked PSD';
h(2).DisplayName='Baseline';
leg=legend(h,'Location','EastOutside');

plot4paper('Frequency (Hz)','PSD');
title('Replay Evoked Power');
grid on;grid minor;grid minor;
axis square;
xlim([45,160]);

ylim([min(PSD_to_plot(:)),max(squash(PSD_to_plot(:,140:end)))])
print([savebase, '\Replay_evoked_psd'],'-depsc');
