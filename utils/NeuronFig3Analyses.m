
%% Replay Paper Figure 3:

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

parc_name='Giles';

studydir = [wd,session_name{whichstudy},'250Hz/'];

% Load in HMM results:
load([studydir,'hmm_1to45hz/hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'],'hmmT','subj_inds');
hmmdir = [studydir,'hmm_1to45hz/'];
load([hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat']);
hmm = hmm_permutestates(hmm,new_state_ordering);
Gamma = hmm.gamma;
scan_T = cell2mat(hmmT);

savebase = fullfile( [wd,session_name{whichstudy},'Fig3_',template_string,'/' ])
if ~exist( savebase )
    mkdir(savebase);
end


%% SECTION 1: PLOT EXAMPLE SUBJECT SEGMENT
% load alignment info:
if whichstudy>1
    if whichstudy==3
        datadir = [studydir,'bfnew_1to45hz/'];
        [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);
        masks_100Hz = load([wd,'GenericReplayData/STUDYII_ReplayOnset/100HzMasks.mat'],'maskA');
        load([wd,'GenericReplayData/RedAll_1StudyII.mat'])
        % Select subject, run and time period to plot:
        ex_sj = 12;
        ex_t = [101:12100] +6000;%90 seconds
    else
        datadir = [studydir,'bfnew_1to45hz/'];
        [maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);
        masks_100Hz = load([wd,'GenericReplayData/100HzMasks.mat'],'maskA');
        % Select subject, run and time period to plot:
        ex_sj = 30;
        ex_t = 8000+[1:12000];
        load([wd,'GenericReplayData/RedAll_1.mat'])
    end
    % load classifier scores:
    scores = permute(RedAll_1,[3,4,2,1]);
    scores = permute(scores(:,:,:),[3,1,2]);

    % for each combination, compute replay probability:
    t_lag = 4; % peak lag found by Liu et al
    for k=1:K;leglabels{k} = ['RSN-State ',int2str(k)];end


    clear replay_pair replay_pair_logical
    for i=[1,2,3,5,6,7]
        replay_pair(i,:) = scores(ex_sj,ex_t,i) .* scores(ex_sj,ex_t+t_lag,i+1);
        replay_pair_logical(i,:) = replay_pair(i,:) > prctile(replay_pair(i,:),99.9);
    end

    %remove redundant index:
    replay_pair(4,:)=[];
    replay_pair_logical(4,:) =[];

    icolormap = 3;
    figure('Position',[143 360 1273 438]);
    LineFormat = struct()
    LineFormat.Color = [0.3 0.3 0.3];
    LineFormat.LineWidth = 2;
    LineFormat.LineStyle = '-';
    subplot(2,1,1);
    plotSpikeRaster(replay_pair_logical,'PlotType','vertline','LineFormat',LineFormat);
    xlim([1,length(ex_t)]);
    set(gca,'XTick',[3000:3000:length(ex_t)]);
    set(gca,'XTickLabel',[30:30:(length(ex_t)/100)]);
    set(gca,'YTick',[1:6]);
    set(gca,'YTickLabel',{'A\rightarrowB','B\rightarrowC','C\rightarrowD', ...
        'A''\rightarrowB''','B''\rightarrowC''','C''\rightarrowD'''});
    plot4paper('Time (sec)','Replayed Sequence')
    legend({'Replay'},'Location','EastOutside')

    % Align subject's state timecourse:
    R = [[1;1+cumsum(scan_T(1:end-1))'],cumsum(scan_T(1:end))'];
    GamToPlot = NaN(length(goodsamples{ex_sj}),K);
    GamToPlot(goodsamples{ex_sj},:) = Gamma(R(ex_sj-1,2)+1:R(ex_sj,2),:);
    GamToPlot = GamToPlot(triggerpoints{ex_sj},:);

    % and plot:
    clear State_windowed;
    for ik=1:K
        State_windowed(:,ik) = conv(GamToPlot(:,ik),ones(1,sample_rate),'same');
    end

    subplot(2,1,2);
    if sample_rate==250
        ex_t_toplot = round(ex_t(1)/100*250):round(ex_t(end)/100*250);
    end
    toplot = (State_windowed(ex_t_toplot,:));
    toplot = toplot ./ repmat(sum(toplot,2),1,size(toplot,2));
    plotStateSequence(fliplr(toplot),1,icolormap);
    xlim([1,length(ex_t_toplot)]);
    t_step = (30*250)-2;%floor(length(toplot)/3);
    set(gca,'XTick',[t_step:t_step:length(toplot)]);
    set(gca,'XTickLabel',[30:30:(length(ex_t)/100)]);
    plot4paper('Time (sec)','Active RSN-State')
    legend({'State 1','State 2'},'Location','EastOutside')
    print([savebase,'/ReplayBursting',int2str(icolormap)],'-depsc')
    set(gcf,'Position',[143 360 861 438]);
    print([savebase,'/ReplayBursting',int2str(icolormap),'squash'],'-depsc')


    %% separately save an appropriate colorbar:
    temp = 1:12;
    figure('Position', [906 567 94 231]);
    yyaxis right;imagesc(temp');
    set(gca,'XTick',[]);
    imagesc(temp');
    yyaxis left;,set(gca,'YTick',[]);
    yyaxis right;
    set(gca,'YTick',[1:12]);
    set(gca,'YTickLabel',leglabels)
    set(gca,'YColor',[0 0 0]);
    for cmapind=1:4
        if cmapind==1
            colors = (parula(K));
        elseif cmapind==2
            colors = (winter(K));
        elseif cmapind==3
            colors = (hot(K));
        elseif cmapind==4
            colors = (copper(K));
        end
        colormap(colors);
        print([savebase,'ExCmap',int2str(cmapind)],'-depsc');
    end
end

%% SECTION 2: HMM Temporal statistics
%
% Here we compute the global temporal statistics for the HMM:
    
% Load in HMM results
load([hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat']);
hmm = hmm_permutestates(hmm,new_state_ordering);

% Create basepath for saving results
savebase = fullfile( studydir, ['figures',template_string ],[sprintf('HMM_K%i',K)],'Figure3/'); 
if ~exist( savebase )
    mkdir(savebase);
end


% account for delay embedding in state gammas
Gamma = hmm.gamma;

options.Fs=sample_rate;
% Compute temporal stats

% Fractional Occupancy is the proportion of time spent in each state
FO = getFractionalOccupancy( Gamma, scan_T, options,2);
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma, scan_T, [],[],[],false);
ITmerged = cellfun(@mean,IT);clear IT
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma, scan_T, [],[],[],false);
LTmerged = cellfun(@mean,LT); clear LT

% Make summary figures
fontsize = 18;

figure;subplot(111);
distributionPlot(FO,'showMM',2,'color',{color_scheme{1:size(FO,2)}});
set(gca,'YLim',[0 1.1*max(FO(:))],'FontSize',fontsize)
title('Fractional Occupancy');plot4paper('RSN-State','Proportion');grid on;
print([savebase 'temporalstats_FO'],'-depsc')

figure;subplot(111);
distributionPlot(LTmerged ./ sample_rate * 1000,'showMM',2,'color',{color_scheme{1:size(FO,2)}})
title('Life Times');plot4paper('RSN-State','Time (ms)');grid on;
YL = 1.1*max(LTmerged(:))./ sample_rate * 1000;
set(gca,'YLim',[0 YL],'FontSize',fontsize,'FontSize',fontsize);
print([savebase 'temporalstats_LT'],'-depsc')

figure;subplot(111);
distributionPlot(log10(ITmerged ./ sample_rate),'showMM',2,'color',{color_scheme{1:size(FO,2)}})
title('Interval Times');plot4paper('RSN-State','Time (secs)');grid on
YL(2) =1.5* max(mean(log10(ITmerged ./ sample_rate)));
YL(1) = min(squash(log10(ITmerged ./ sample_rate)));
set(gca,'YLim',YL,'FontSize',fontsize)
set(gca,'YTick',log10([0.05,0.1,0.5,1,5,10]))
y_labels = get(gca,'YTickLabel');
for i=1:length(y_labels)
    y_labels{i}=num2str(10.^(str2num(y_labels{i})),1);
end
set(gca,'YTickLabels',y_labels);
print([savebase 'temporalstats_IT_logscale'],'-depsc')

figure;subplot(111);
distributionPlot(ITmerged ./ sample_rate,'showMM',2,'color',{color_scheme{1:size(FO,2)}})
title('Interval Times');plot4paper('RSN-State','Time (secs)');grid on
YL(2) =1.5* max(mean((ITmerged ./ sample_rate)));
YL(1) = 0;
set(gca,'YLim',YL,'FontSize',fontsize)
print([savebase 'temporalstats_IT'],'-depsc')

close all;


%% SECTION 2: Replay Temporal statistics
%
% Here we analyse the Fano Factor of the replay timecourse:
if whichstudy>1
    % Fano Factor plots:

    % load replay timecourse:
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

    thresholdlevel = 99; %percentile cutoff for replay
    T_session = size(replayScores,2);
    replayscores_processed = zeros(nSj{whichstudy},2*T_session);

    % setup window sizes to compute FF from:
    windows = round(exp(linspace(log(2),log(10000),100)));
    scores_windowed = nan(2*T_session,length(windows),nSj{whichstudy});

    %%
    clear FF_sj_perm;

    if whichstudy==2
        permfile = [wd,'GenericReplayData/STUDYI_ReplayOnset/FanoFactor_Perms.mat']; 
    else
        permfile = [wd,'GenericReplayData/STUDYII_ReplayOnset/FanoFactor_Perms.mat'];
    end
    if ~exist(permfile)
        % generate surrogate data by permuting the intervals between events:
        nperms = 1001;
        clear replayscores_processed;
            for iSj=1:nSj{whichstudy}
            fprintf(['Processing subj: ',int2str(iSj),'\n']);
            nsess = (iSj-1)*2+1;
            replayscores_processed(iSj,1:T_session) = double(replayScores(nsess,:)>prctile(replayScores(nsess,logical(masks_100Hz.maskA{nsess})),thresholdlevel));
            maskfull = masks_100Hz.maskA{nsess};
            nsess = (iSj-1)*2+2;
            replayscores_processed(iSj,T_session+1:2*T_session) = double(replayScores(nsess,:)>prctile(replayScores(nsess,logical(masks_100Hz.maskA{nsess})),thresholdlevel));
            maskfull = [maskfull;masks_100Hz.maskA{nsess}];
            maskfull(T_session+1) = 0;
            % extract intervals:
            intervals = find(replayscores_processed(iSj,:));
            intervals = [intervals(1),intervals(2:end) - intervals(1:end-1)];
            %remove segments marked as bad:
            badintervals = find(maskA{iSj}==0);
            to_remove = false(length(intervals),1);
            for i=2:length(intervals)
                if any(sum(intervals(1:i-1))<=badintervals & badintervals<=sum(intervals(1:i)))
                    to_remove(i)=true;
                end
            end
            to_remove(1)=true; %offset that was kept for alignment with badsamples
            intervals(to_remove)=[];
            for iperm = 1:nperms
                % we will save in the first row the non-permuted intervals:
                if iperm>1 
                    int_permuted = intervals(randperm(length(intervals)));
                else
                    int_permuted = intervals;
                end
                replayscores_permuted = [1,zeros(1,sum(int_permuted))];
                replayscores_permuted(1+cumsum(int_permuted))=1;
                scores_windowed = zeros(1+sum(int_permuted),length(windows));
                for i=1:length(windows)
                    windowlength=windows(i);
                    % non-overlapping windows:
                    scores_windowed = conv(replayscores_permuted,ones(windowlength,1),'same');
                    temp = scores_windowed(round(windowlength/2):windowlength:end);
                    mu = nanmean(temp);
                    sigsq = nanvar(temp);
                    FF_sj_perm(iperm,i,iSj) = sigsq/mu;
                end
            end
        end
        %save for quicker access later
        save(permfile,'FF_sj_perm'); 
    else
        load(permfile,'FF_sj_perm'); 
    end

    FF_sj_run = squeeze(FF_sj_perm(1,:,:));
    FF_sj_perm = FF_sj_perm(2:end,:,:);

    % Fano Factor Plot
    colors = set1_cols();
    remove_outliers = false;
    if remove_outliers
        outliers = nanmean(FF_sj_run) - nanmean(nanmean(FF_sj_run)) >  1.5*std(nanmean(FF_sj_run));
        FF_sj_run(:,outliers) = [];
        FF_sj_perm(:,:,outliers)=[];
    else
        outliers = false(nSj{whichstudy},1);
    end
    displaypoints = [1:100];
    figure('Position',[293 225 736 421]);
    clear h;
    semilogx(windows/sample_rate,mean(FF_sj_run(:,:),2));hold on;
    plot4paper('Window length (seconds)','Fano factor');
    axis square;
    title('Replay: Fano factor vs window size');
    AX = gca;
    AX.XMinorGrid = 'off';
    grid on;hold on;
    shadedErrorBar(windows/sample_rate,mean(FF_sj_run,2),std(FF_sj_run(displaypoints,:),[],2)./sum(outliers==0));
    semilogx(windows/sample_rate,mean(FF_sj_run(:,:),2),'LineWidth',2,'Color','black');hold on;
    h(1) = line(NaN,NaN,'Color','black','LineWidth',2);

    xlim([10e-2,30]);
    if whichstudy==2
        ylim([0,8]);
    else
        ylim([0,13]);
    end
    set(gca,'XTick',[0.1,1,10,30]);
    % find mean permutation FF:
    mu = squeeze(mean(FF_sj_perm,1));
    mu = mean(mu,2);
    hold on;
    mu_ste = squeeze(mean(FF_sj_perm,3));
    mu_ste_upper = prctile(mu_ste,97.5); 
    semilogx(windows/sample_rate,mu_ste_upper,'LineWidth',1,'Color','k','LineStyle','--');
    semilogx(windows/sample_rate,ones(1,length(windows/100)));
    h(2) = line(NaN,NaN,'Color',[0 0 0],'LineStyle','--','LineWidth',1);
    h(1).DisplayName = ['Replay Fano Factor'];
    h(2).DisplayName = ['IID Significance level'];
    leg=legend(h,'Location','EastOutside');

    print([savebase,'/FanoFactorVsWindowSize_Replay_intervals'],'-depsc')
end
%% SECTION 3: Fano Factor for state sequences:
% recompute windows given new sample rate:
if whichstudy>1
    windows = round(exp(linspace(log(2),log(10000),100)));
    clear FF_sj_perm_states;
    if ~exist([savebase,'FanoFactor_Perms_state.mat'])
        % generate surrogate data by permuting the intervals between events:
        nperms = 2;
        for iSj=1:nSj{whichstudy}
            fprintf(['Processing subj: ',int2str(iSj),'\n']);
            nsess = (iSj-1)*2+1;
            temp=zeros(length(goodsamples{nsess}),1);
            temp(goodsamples{nsess}) = hmm.statepath(hmm.subj_inds==nsess);
            nsess = (iSj-1)*2+2;
            temp2 = zeros(length(goodsamples{nsess}),1);
            temp2(goodsamples{nsess}) = hmm.statepath(hmm.subj_inds==nsess);
            vpath_sj = [0;temp;0;temp2;0];
            mask = [0,goodsamples{nsess-1},0,goodsamples{nsess},0];
            for ik=1:K
                ind_off = find(diff(vpath_sj==ik)==-1);ind_on = find(diff(vpath_sj==ik)==1);
                badintervals = find(mask==0);
                to_remove = false(length(ind_off),1);
                for i=2:length(ind_on)
                    if any(ind_off(i-1)<=badintervals & badintervals<=ind_off(i))
                        to_remove(i)=true;
                    end
                end
                ind_on(to_remove)=[];ind_off(to_remove)=[];
                intervals_all{ik,iSj} = ind_on(2:end)-ind_off(1:end-1);
            end

            for ik=1:K
                intervals = intervals_all{ik,iSj};
                for iperm = 1:nperms
                    if iperm>1
                        int_permuted = intervals(randperm(length(intervals)));
                    else
                        int_permuted = intervals;
                    end
                    replayscores_permuted = [1,zeros(1,sum(int_permuted))];
                    replayscores_permuted(1+cumsum(int_permuted))=1;
                    for i=1:length(windows)
                        windowlength = windows(i);
                        % non-overlapping windows:
                        scores_windowed = conv(replayscores_permuted,ones(windowlength,1),'same');
                        temp = scores_windowed(round(windowlength/2):windowlength:end);
                        mu = nanmean(temp);
                        sigsq = nanvar(temp);
                        FF_sj_perm_states(iperm,i,ik,iSj) = sigsq/mu;
                    end
                end
            end
        end
        %save for quicker access later
        save([savebase,'FanoFactor_Perms_state'],'FF_sj_perm_states','intervals_all'); 
    else
        load([savebase,'FanoFactor_Perms_state'],'FF_sj_perm_states'); 
    end

    FF_sj_run_k = squeeze(FF_sj_perm_states(1,:,:,:));
    FF_sj_perm_k = FF_sj_perm_states(2:end,:,:,:);

    %% Fano Factor Plot

    remove_outliers = false;

    K_to_plot = [1:K];
    K_to_emphasise = [1:2];%new_order(1:2);
    if remove_outliers
        outliers = [1:nSj{whichstudy}]==6;
        FF_sj_run_k(:,:,outliers) = [];
        FF_sj_perm(:,:,outliers)=[];
    else
        outliers = false(nSj{whichstudy},1);
    end
    displaypoints = [1:100];
    figure();
    for ik=1:K
        if mod(ik,2)==1
             ls = '-';
         else
             ls='--';
        end
        if ~isempty(intersect(K_to_emphasise,ik))
            lw = 4;tp=1;
        else
            lw=1;tp=0.2;
        end
        semilogx(windows/sample_rate,mean(FF_sj_run_k(:,ik,:),3),'LineWidth',0.5);hold on;
        mu=mean(FF_sj_run_k(:,ik,:),3);
        ste=std(FF_sj_run_k(:,ik,:),[],3)./sqrt(size(FF_sj_run_k,3));
        shadedErrorBar(windows/sample_rate,mu,ste,{'Color',colors{ik},'LineWidth',lw,'LineStyle',ls},0.2);
        hold on;
        h(ik) = line(NaN,NaN,'Color',colors{ik},'LineWidth',lw,'LineStyle',ls);
    %errorbar(windows(displaypoints)/100,mean(FF_sj_run(displaypoints,:),2),std(FF_sj_run(displaypoints,:),[],2)./sum(outliers==0));hold on;
    %semilogx(windows/100,mean(FF_sj_run(:,:),2),'LineWidth',2);
    end
    xlim([10e-2,30]);%ylim([0,7]);
    set(gca,'XTick',[0.1,1,10,30]);
    for k=1:K,h(k).DisplayName=['RSN-State ',int2str(k)];end
    leg=legend(h,'Location','EastOutside');
    plot4paper('Window length (seconds)','Fano factor');
    axis square;
    %title('HMM States: Fano factor vs window size');
    AX = gca;
    AX.XMinorGrid = 'off';
    grid on;hold on;
    %shadedErrorBar(windows/100,mean(FF_sj_run_k,3),std(FF_sj_run(displaypoints,:),[],2)./sum(outliers==0));
    %semilogx(windows/100,mean(FF_sj_run(:,:),2),'LineWidth',2);
    print([savebase,'/FanoFactorVsWindowSize_States'],'-depsc')

    %% parametric tests:

    for W_totest = 1:100
        groupstats = (squeeze(FF_sj_run_k(W_totest,:,:)));
        groupinfo = repmat([1:12]',1,nSj{whichstudy});
        [p_anova(W_totest),tbl,stats] = anova1(groupstats(:),groupinfo(:));
        close all;
        % now run one way t-tests:
        for istate=1:12
            [~,pval_ttest(W_totest,istate)] = ttest2(groupstats(groupinfo==istate),groupstats(groupinfo~=istate));
        end
    end

    %%
    figure('Position',[1 378 1440 420]);subplot(1,3,1);
    semilogy(windows/sample_rate,p_anova);title('OneWay ANOVA p value for differences between states');hold on;
    semilogy(windows([1,end])/sample_rate,[0.05,0.05],'LineStyle','--');
    semilogy(windows([1,end])/sample_rate,[0.05,0.05]/K,'LineStyle','--');
    semilogy(windows([1,end])/sample_rate,[0.05,0.05]/(K*length(windows)),'LineStyle','--');
    plot4paper('WindowSize','Pvalue');
    xlim([0,30]);
    legend({'pval','p=0.05','Bonf corrected for states','Bonfcorrected for states+windows'},'Location','SouthEast');
    subplot(1,3,2);
    semilogy(windows/sample_rate,pval_ttest(:,1));title('OneWay ttest for state 1 larger than all others');hold on;
    semilogy(windows([1,end])/sample_rate,[0.05,0.05],'LineStyle','--');
    semilogy(windows([1,end])/sample_rate,[0.05,0.05]/K,'LineStyle','--');
    semilogy(windows([1,end])/sample_rate,[0.05,0.05]/(K*length(windows)),'LineStyle','--');
    plot4paper('WindowSize','Pvalue');
    xlim([0,30]);
    subplot(1,3,3);
    semilogy(windows/sample_rate,pval_ttest(:,2));title('OneWay ttest for state 2 larger than all others');hold on;
    semilogy(windows([1,end])/sample_rate,[0.05,0.05],'LineStyle','--');
    semilogy(windows([1,end])/sample_rate,[0.05,0.05]/K,'LineStyle','--');
    semilogy(windows([1,end])/sample_rate,[0.05,0.05]/(K*length(windows)),'LineStyle','--');
    plot4paper('WindowSize','Pvalue');
    legend({'pval','p=0.05','Bonf corrected for states','Bonfcorrected for states+windows'},'Location','SouthEast');
    print([savebase,'/FanoFactorVsWindowSize_States_parametrictests'],'-depsc')
    xlim([0,30]);

    pvalmaxanova = max(p_anova((windows/sample_rate)<30))
    pvalmaxttests= max(pval_ttest((windows/sample_rate)<30,:),[],1)
end
%% SECTION 4: REPLAY INTERVALS CONDITIONED ON ACTIVE STATE
if whichstudy>1
    % compute replay intervals given active state:
    IT_given_state=cell(nSj{whichstudy},K);
    IT_subj=cell(nSj{whichstudy},1);
    for iSj=1:nSj{whichstudy}
        fprintf(['Processing subj: ',int2str(iSj),'\n']);
        nsess = (iSj-1)*2+1;
        replayscores_processed(iSj,1:T_session) = double(replayScores(nsess,:)>prctile(replayScores(nsess,logical(masks_100Hz.maskA{nsess})),thresholdlevel));
        maskfull = masks_100Hz.maskA{nsess};
        if nsess==25 && whichstudy==2 % one short session:
            maskfull = [maskfull;zeros(18000,1)];
        end
        nsess = (iSj-1)*2+2;
        replayscores_processed(iSj,T_session+1:2*T_session) = double(replayScores(nsess,:)>prctile(replayScores(nsess,logical(masks_100Hz.maskA{nsess})),thresholdlevel));
        maskfull = [maskfull;masks_100Hz.maskA{nsess}];
        maskfull(T_session+1) = 0;
        intervals = find(replayscores_processed(iSj,:));

        replaytimes = round(sample_rate/100 * intervals);
        %align state timecourses:
        nsess = (iSj-1)*2+1;
        temp1 = zeros(length(goodsamples{nsess}),K);
        temp1(goodsamples{nsess},:) = Gamma(R(nsess,1):R(nsess,2),:);
        temp1 = temp1(triggerpoints{nsess},:);
        if nsess==25 && whichstudy==2
            temp1 = [temp1;zeros(45000,K)];
        end
        nsess = (iSj-1)*2+2;
        temp2 = zeros(length(goodsamples{nsess}),K);
        temp2(goodsamples{nsess},:) = Gamma(R(nsess,1):R(nsess,2),:);
        temp2 = temp2(triggerpoints{nsess},:);
        Gam_long = [temp1;temp2];

        for i_int=1:length(replaytimes)-1
            if any(all(Gam_long(replaytimes(i_int):replaytimes(i_int+1),:)==0,2) )
                % bad segments between replay times; ignore
            else

                [~,activestate] = max(Gam_long(replaytimes(i_int),:));
                IT_given_state{iSj,activestate} = cat(2,IT_given_state{iSj,activestate},(replaytimes(i_int+1)-replaytimes(i_int))./ sample_rate);

    %             % also track intervals in consecutive order with
    %             % corresponding gamma info:
    %             if i_int~=length(replaytimes)-1 && ~any(all(Gam_long(replaytimes(i_int+1):replaytimes(i_int+2),:)==0,2) )
    %                 IT_subj{iSj} = cat(1,IT_subj{iSj},(replaytimes(i_int+2)-replaytimes(i_int+1))./ sample_rate);
    %                 ARterms = [];
    %                 for t=1:ARorder
    %                     if i_int+1-t>0
    %                         ARterms = [ARterms,(replaytimes(i_int+2-t)-replaytimes(i_int+1-t))./ sample_rate];
    %                     else
    %                         ARterms = [ARterms,0];
    %                     end
    %                 end
    %                 DM_subj{iSj} = cat(1,DM_subj{iSj},[Gam_long(replaytimes(i_int+1),1:K-1),ARterms]);
    %             end

                %replay_given_state(replaytimes(i),activestate)=1;
                %replay_given_state(replaytimes(i+1),activestate)=1;
            end
        end
    end

    ITmerged = cellfun(@mean,IT_given_state);
    IT_mean_all = mean(ITmerged,2);

    % remove entries with very low sampling:
    samplelimit=10;
    for i1=1:nSj{whichstudy}
        for i2=1:hmm.K
            Lmat(i1,i2) = length(IT_given_state{i1,i2});
        end
    end
    labels = repmat([1:hmm.K],nSj{whichstudy},1);
    sufficientlysampled = Lmat>=samplelimit;
    goodsamps = ITmerged(sufficientlysampled);
    labels = labels(sufficientlysampled);
    [pval_anova_states,~,stats] = anova1(goodsamps,labels);

    for ik=1:K
        labels2 = labels==ik;
        [~,pval_ind(ik)] = ttest2(goodsamps(labels==ik), goodsamps(labels~=ik),'Tail','left');
    end
    ITmerged = [IT_mean_all,ITmerged];


        % Make summary figures
    %     fontsize = 18;
    %     figure;subplot(111);
    %     distributionPlot(log10(ITmerged ./ sample_rate),'showMM',2,'color',{[0.8,0.8,0.8],set1_cols{1:K}})
    %     title('Interval Times');plot4paper('State','Time (secs)');grid on
    %     %YL(2) =1.5* max(mean(log10(ITmerged ./ sample_rate)));
    %     %YL(1) = min(squash(log10(ITmerged ./ sample_rate)));
    %     %set(gca,'YLim',YL,'FontSize',fontsize)
    %     %set(gca,'YTick',log10([0.05,0.1,0.5,1,5,10]))
    %     set(gca,'XTick',1:K+1);
    %     set(gca,'XTickLabels',{'Mean',1:K});
    %     y_labels = get(gca,'YTickLabel');
    %     for i=1:length(y_labels)
    %         y_labels{i}=num2str(10.^(str2num(y_labels{i})),1);
    %     end
    %     set(gca,'YTickLabels',y_labels);
    %     box on;
        %print([savebase '_temporalstats_IT_logscale'],'-depsc')

    figure;subplot(111);
    distributionPlot(ITmerged ,'showMM',2,'color',{[0.8,0.8,0.8],color_scheme{1:K}});

    ylim([0,2.5]);
    for ik=1:K;
        YL = ylim();
        if pval_ind(ik)<0.05
            hold on;
            plot(ik+1,YL(2)*0.95,'k*','MarkerSize',10,'linewidth',2)
        end
    end
    title('Replay Interval Times given Active State');
    plot4paper('State','Time (secs)');grid on
    set(gca,'XTick',1:K+1);
    set(gca,'XTickLabels',{'Mean',1:K});

    %YL(2) =1.5* max(mean((ITmerged ./ sample_rate)));
    %YL(1) = 0;
    %set(gca,'YLim',YL,'FontSize',fontsize)
    box on;
    print([savebase '/StateConditionalReplayIntervals_mean'],'-depsc')

    figure;subplot(111);
    cpl=[0 0 0 ];for k=1:K;cpl = cat(1,cpl,color_scheme{k});end
    B = boxplot(ITmerged ,'color',cpl,'PlotStyle','compact','Widths',1);
    h=findobj(gca,'tag','Box'); % Get handles for outlier lines.
       set(h,'LineWidth',15); % Change symbols for all the groups.
       %set(h(1),'MarkerEdgeColor','b');
    %boxplot(ITmerged ,'color',cpl,'BoxStyle','filled','Widths',1);
    h=findobj(gca,'tag','Whisker');
    set(h,'LineWidth',4);
    ylim([0,2.5]);
    for ik=1:K;
        YL = ylim();
        if pval_ind(ik)<0.05
            hold on;
            plot(ik+1,YL(2)*0.95,'k*','MarkerSize',10,'linewidth',2)
        end
    end
    title('Replay Interval Times given Active State');
    plot4paper('State','Time (secs)');grid on
    set(gca,'XTick',1:K+1);
    set(gca,'XTickLabels',{'Mean',1:K});
    line([0,K+2],[mean(IT_mean_all),mean(IT_mean_all)],'Color','black','LineStyle','--');
    xlim([0.5,K+1.5]);
    %YL(2) =1.5* max(mean((ITmerged ./ sample_rate)));
    %YL(1) = 0;
    %set(gca,'YLim',YL,'FontSize',fontsize)
    box on;
    print([savebase '/StateConditionalReplayIntervals_boxplot'],'-depsc')

    % or plot all events:
    for ik=1:K
        ITmerged_all{ik} = cat(2,[],IT_given_state{:,ik});
    end
    ITmerged2 = ITmerged_all;
    figure;subplot(111);
    distributionPlot(ITmerged2 ,'showMM',2,'color',{[0.8,0.8,0.8],color_scheme{1:K}},'globalNorm',2)
    title('Replay Interval Times given Active State');
    plot4paper('State','Time (secs)');grid on
    %set(gca,'XTick',1:K+1);
    %set(gca,'XTickLabels',{'Mean',1:K});
    box on;
    ylim([0,5]);
    print([savebase '/StateConditionalReplayIntervals_all'],'-depsc')
end

%% SECTION 5: STATE EVOKED REPLAY INCIDENCE
if whichstudy>1
    replayrate = zeros(nSj{whichstudy},K);
    for iSj=1:nSj{whichstudy}
        fprintf(['Processing subj: ',int2str(iSj),'\n']);

        % load replay scores:
        nsess = (iSj-1)*2+1;
        replayscores_processed(iSj,1:T_session) = double(replayScores(nsess,:)>prctile(replayScores(nsess,logical(masks_100Hz.maskA{nsess})),thresholdlevel));
        maskfull = masks_100Hz.maskA{nsess};
        nsess = (iSj-1)*2+2;
        replayscores_processed(iSj,T_session+1:2*T_session) = double(replayScores(nsess,:)>prctile(replayScores(nsess,logical(masks_100Hz.maskA{nsess})),thresholdlevel));
        maskfull = [maskfull;masks_100Hz.maskA{nsess}];
        maskfull(T_session+1) = 0;

        repscore = replayscores_processed(iSj,:);

        % load viterbi paths:
        temp1 = zeros(length(goodsamples{nsess}),1);
        temp1(goodsamples{nsess}) = hmm.statepath(R(nsess,1):R(nsess,2));
        temp1 = temp1(triggerpoints{nsess},:);
        if nsess==25 && whichstudy==1
            temp1 = [temp1;zeros(45000,K)];
        end
        nsess = (iSj-1)*2+2;
        temp2 = zeros(length(goodsamples{nsess}),1);
        temp2(goodsamples{nsess}) = hmm.statepath(R(nsess,1):R(nsess,2));
        temp2 = temp2(triggerpoints{nsess});
        vpath = [temp1;temp2];

        % now compute Poisson rate of replay for each state:
        for ik=1:K
            t_on = find(diff([0;vpath==ik;0])==1)-1;
            t_off = find(diff([0;vpath==ik;0])==-1)-1;
            % adjust for sampling rate difference:
            t_on = t_on * 100/250;
            t_off = t_off * 100/250;
            t_in = false(size(repscore,2),1);
            for t=1:length(t_on)
                t_in = t_in | [[1:60000]'>t_on(t) & [1:60000]'<t_off(t)];
            end
            replayrate(iSj,ik) = mean(repscore(t_in));
        end

     end

    replayrate = [mean(replayrate,2),replayrate];

    for ik=1:K
        notK = setdiff(1:K,ik);
        [~,pval_ind(ik)] = ttest2(replayrate(:,ik+1),squash(replayrate(:,notK+1)),'Tail','right');
    end

    %%
    figure('Position',[440 536 560 262]);subplot(111);
    cpl=[0 0 0 ];for k=1:K;cpl = cat(1,cpl,color_scheme{k});end
    B = boxplot(replayrate ,'color',cpl,'PlotStyle','compact','Widths',1);
    h=findobj(gca,'tag','Box'); % Get handles for outlier lines.
       set(h,'LineWidth',15); % Change symbols for all the groups.
       %set(h(1),'MarkerEdgeColor','b');
    %boxplot(ITmerged ,'color',cpl,'BoxStyle','filled','Widths',1);
    h=findobj(gca,'tag','Whisker');
    set(h,'LineWidth',4);
    O = 0.15;
    for ik=1:K;
        YL = ylim();
        if pval_ind(ik)<0.0005
            hold on;
            plot((ik+1)- O ,YL(2)*0.95,'k*','MarkerSize',10,'linewidth',2)
            plot((ik+1)+ O ,YL(2)*0.95,'k*','MarkerSize',10,'linewidth',2)
        end
    end

    title('Replay Rate given Active State');
    plot4paper('State','Time (secs)');grid on
    set(gca,'XTick',1:K+1);
    set(gca,'XTickLabels',{'Mean',1:K});
    line([0,K+2],[mean(replayrate(:)),mean(replayrate(:))],'Color','black','LineStyle','--');
    xlim([0.5,K+1.5]);
    %YL(2) =1.5* max(mean((ITmerged ./ sample_rate)));
    %YL(1) = 0;
    %set(gca,'YLim',YL,'FontSize',fontsize)
    box on;
    set(gca,'YTick',[0:0.01:0.05]);
    print([savebase '/Replayrate_givenRSN_boxplot'],'-depsc')
    set(gcf,'Position',[440 379 560 419]);
    print([savebase '/Replayrate_givenRSN_boxplot_taller'],'-depsc')

    %% histogram replay scores and plot a sample replay timecourse:
    for iSj=1:2:size(replayScores)
        RS = replayScores(iSj,:);
        RS = RS./max(RS);
        xs = 0:0.02:1;
        [xblock,xs] = histcounts(RS,xs);
        if iSj==1
            xblock_all = xblock;
        else
            xblock_all = xblock_all + xblock;
        end
        m(iSj) = max(RS);
    end
    figure('Position',[440 515 503 283]);
    bar(xs(2:end),xblock,1);
    ylim([0,2000]);
    set(gca,'YTick',[0:500:1500,1750,2000]);
    set(gca,'YTickLabel',{'0','500','1000','1500','//','2.7e4'});
    plot4paper('Replay Score','Frequency');
    %title('Histogram of replay scores');
    print([savebase '/ReplayScoreHistogram'],'-depsc')

    RS = replayScores(1,:);
    RS = RS./max(RS);
    t = 0:0.01:15;
    figure('Position',[440 515 503 283]);plot(t,RS(9000:10500),'LineWidth',2,'Color','black');
    plot4paper('Time (sec)','Replay Score');
    print([savebase '/ReplayScoreSample'],'-depsc')
end