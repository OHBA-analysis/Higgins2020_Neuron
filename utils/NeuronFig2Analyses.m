%% script to run analyses related to figure 2 of Neuron2020 Paper

%% general setup for this study data:
if whichstudy==1
    CanonicalRSN = false;
else
    CanonicalRSN = true; % canonical refers to RSNs trained on nottingham data and refitted
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

studydir = [wd,session_name{whichstudy},'250Hz/'];

savebase = fullfile( [wd,session_name{whichstudy},'Fig2_',template_string,'/' ])
if ~exist( savebase )
    mkdir(savebase);
end


hmmdir = [studydir,'hmm_1to45hz/'];
hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
load(hmmfile);    % no need to permute here as will do below 
prepdatafile = [hmmdir,'hmm_parc_giles_symmetric__pcdim80_voxelwise_embed14.mat'];
load(prepdatafile,'hmmT','subj_inds');
scan_T = cell2mat(hmmT);

parc_file = ['fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
parc = parcellation(parc_file);
for k=1:hmm.K;statelabels{k}={'RSN-State ',int2str(k)};end

%% SECTION 1: Analyse Transition matrix and order states:
    
% Plot MDS map of state network:
disttoplot = plotMDS_states(hmm);
[~,new_state_ordering] = sort(disttoplot(:,1));

if any(new_state_ordering(2:end) ~= new_state_ordering(1:end-1)+1)
    hmm = hmm_permutestates(hmm,new_state_ordering);
    hmmfile = [hmmdir,'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
    save(hmmfile,'new_state_ordering','-append');
    disttoplot = disttoplot(new_state_ordering,:);
end
Gamma = hmm.gamma;

figure('Position', [440 519 391 279]);
for ik=1:hmm.K
    for jk=1:hmm.K
        if ik~=jk
            line([disttoplot(ik,1),disttoplot(jk,1)],[disttoplot(ik,2),disttoplot(jk,2)],...
                'color',0.5*[1,1,1]);hold on;
        end
    end
end
for ik=1:hmm.K
    scatter1 = scatter(disttoplot(ik,1),disttoplot(ik,2),400,...
        'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik}); 
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    if ik==10
        scatter1 = scatter(disttoplot(ik-1,1),disttoplot(ik-1,2),400,...
        'MarkerFaceColor',color_scheme{ik-1},'MarkerEdgeColor',color_scheme{ik-1}); 
        hold on
        % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
        scatter1.MarkerFaceAlpha = 0.5;
        text(disttoplot(ik-1,1)-0.03,disttoplot(ik-1,2),int2str(ik-1),'FontSize',12,'FontWeight','bold');hold on;
    end
    if ik<10
        text(disttoplot(ik,1)-0.03,disttoplot(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    else
        text(disttoplot(ik,1)-0.05,disttoplot(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    end
end
axis square
axis off
print([savebase '/network_transmat_allconnections'],'-depsc')

% and plot transmat as matrix:
figure('Position',[440 527 294 271]);
imagesc(hmm.P);caxis([0,1.1*max(squash(hmm.P-diag(diag(hmm.P))))]);
colormap(hot)
plot4paper('State at time t+1','State at time t');
title('HMM Transition Matrix');
axis square;
set(gca,'YTick',[1:2:11]);
print([savebase '/network_transmat_matrix'],'-depsc');
set(gcf,'Position',[440 527 458 271]);colorbar;
print([savebase '/network_transmat_matrix_colorbar'],'-depsc');

%% SECTION 2: STATE REPLAY ALIGNMENT
if whichstudy > 1
    datadir = [studydir,'bfnew_1to45hz/'];
    
    % load in alignment data:
    [~,~,triggerpoints,goodsamples] = getSubjectMasks(datadir);
    
    % load replay times:
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
   
    pthreshold = 0.05/hmm.K; %display threshold
    emphasis_states = 1:6;
    emphasis_string= 'SigEmphasis';
    
    t_window = sample_rate/2;
    [~,betas_replay,betas_replay_norm] = plotReplayFig1(Gamma,replayScores,t_window,goodsamples,triggerpoints,scan_T,pthreshold,1,color_scheme,emphasis_states);
    YL = ylim();
    if whichstudy==2;YL(2) = YL(2)+0.02;ylim(YL);end
    line([0,0],YL,'LineWidth',2,'Color','red','LineStyle','--');
    
    gcf;
    set(gcf,'Position',[440 411 967 387]); % this one for legend = 'RSN-state'; above just for 'state'
    pause(0.5);
    print([savebase '_replayAlignment',emphasis_string],'-depsc')
    set(gcf,'Position',[440 236 588 563]); % this one for legend = 'RSN-state'; above just for 'state'
    print([savebase '_replayAlignment_tall',emphasis_string],'-depsc')

    % Supplementary figure: correlation analysis with transmat row:
    xs = setdiff([1:12],2);
    transprobs = log10(hmm.P(2,:)');
    b0 = mean(betas_replay_norm(126,:,:),3)';
    [rho,pval,rho_upper,rho_downer] = corrcoef(transprobs(xs),b0(xs),'alpha',0.05);
    figure('Position', [440 432 398 366]);
    for i=1:12
        scatter(transprobs(i),b0(i),50,color_scheme{i},'filled');hold on;
    end
    plot4paper('Transition probability from RSN-state 2','Replay evoked RSN-state probability');
    box on;
    axis square;
    ylim(YL);
    set(gca,'YTick',-0.05:0.05:0.1);
    xlim([-6,-1]);
    set(gca,'XTick',[-6:1:-1]);
    for i=1:6; labelstemp{i} = ['10^{-',int2str(7-i),'}'];end
    set(gca,'XTickLabel',labelstemp);
    % add correlation line with confidence interval:
    x = [-6:-1];
    y = mean(b0(xs)) + rho(2,1)*std(b0(xs))*(x-mean((transprobs(xs))))./std((transprobs(xs)));
    plot(x,y,'Color','black')
    TB1 = annotation('textbox','Position',[0.3555 0.4891 0.3 0.0628]);
    TB1.String = {['\rho = ',num2str(rho(2,1),2)]};
    TB1.FontSize = 18;
    TB1.LineStyle = 'none';
    print([savebase 'Correlation_st2Transprob',emphasis_string],'-depsc')
    
    YL_replayalignment = ylim;
        
    % run direct analysis:
    if whichstudy==2
        % find cluster times:
        for k=1:K
            dat = permute(betas_replay_norm(:,k,:),[3,2,1]);
            thresh = 3;
            nP=5000;
            [corrp(:,k),tstats_cluster(:,k)] = osl_clustertf(dat,thresh,nP);
        end
        clustertimes = corrp>(1-1e-3);
        save([savebase,'Study1_clustertimes.mat'],'corrp','clustertimes','tstats_cluster');
    elseif whichstudy==3
        load([strrep(savebase,'Study2','Study1'),'Study1_clustertimes.mat'],'corrp','clustertimes','tstats_cluster');
        corrp2 = zeros(K,1);
        for k=1:K
            dat = permute(betas_replay_norm(:,k,:),[3,2,1]);
            thresh = 3;
            nP=5000;
            [~,tstats_cluster2(:,k)] = osl_clustertf(dat,thresh,nP);
            threshKDirect = min(tstats_cluster2(clustertimes(:,k),k));
            if ~isempty(threshKDirect)
                [~,~,nulldist] = osl_clustertf_direct(dat,threshKDirect,nP);
                % nulldist is vector of largest cluster with t stat given by
                % threshKdirect:
                corrp2(k) = 1-sum(sum(clustertimes(:,k))>nulldist) ./ nP;
            end
            
        end
    end
        
    % and FLI evoked response:
    if whichstudy==2
        datadirFLI = [strrep(studydir,'Study1','Study1_FLI'),'bfnew_1to45hz/'];
        hmm_fli_file = [strrep(hmmdir,'Study1','Study1_FLI'),'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
    elseif whichstudy==3
        datadirFLI = [strrep(studydir,'Study2','Study2_FLI'),'bfnew_1to45hz/'];
        hmm_fli_file = [strrep(hmmdir,'Study2','Study2_FLI'),'hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm.mat'];
    end
    clear mg;
    for k=1:nscans{whichstudy};mg(k,:) = mean(hmm.gamma(hmm.subj_inds==k,:));end
    
    YL=[];ylimstring = '';
    YL = YL_replayalignment;
    ylimstring = 'scaled';

    emphasis_states = 1:6;
    emphasis_string= 'SigEmphasis';

    [~,betas_FLI,betas_FLI_norm] = plotFLIEvokedResponse(datadirFLI,hmm_fli_file,mg,color_scheme,YL,emphasis_states);
    set(gcf,'Position',[440 411 931 387]);
    gcf;xlim([-0.5,0.5]);
    yl = ylim();
    line([0.2,0.2],yl,'LineWidth',2,'Color','red','LineStyle','--');
    print([savebase '_FLIAlignment',ylimstring,emphasis_string],'-depsc')
    gcf;xlim([-0.3,0.7]);
    print([savebase '_FLIAlignment_late',ylimstring,emphasis_string],'-depsc')
    gcf;xlim([-0.5,0.5]);
    set(gcf,'Position',[440 235 549 563]);pause(0.1)
    print([savebase '_FLIAlignment_tall',ylimstring,emphasis_string],'-depsc');pause(0.1);
    gcf;xlim([-0.3,0.7]);
    print([savebase '_FLIAlignment_tall_late',ylimstring,emphasis_string],'-depsc')
        
    
    
    % and do a paired t test on inst distribution:
    save([savebase,'betas_replayFLI.mat'],'betas_replay','betas_FLI','betas_replay_norm','betas_FLI_norm');
    figure();
    t_FLI = 800; % timepoint corresponding to training data
    t_replay = t_window+1;%timepoint corresponding to replay
    markerpoint(1) = max(squash(nanmean(betas_FLI_norm(t_FLI,:,:),3) + nanstd(betas_FLI_norm(t_FLI,:,:),[],3)./sqrt(size(betas_FLI_norm,3))));
    markerpoint(2) = max(squash(nanmean(betas_replay_norm(t_replay,:,:),3) + nanstd(betas_replay_norm(t_replay,:,:),[],3)./sqrt(size(betas_replay_norm,3))));
    markerpoint = 1.15*max(markerpoint);
    for k=1:hmm.K
        errorbar(k-0.175,nanmean(betas_FLI_norm(t_FLI,k,:),3),std(betas_FLI_norm(t_FLI,k,:),[],3)./sqrt(size(betas_FLI_norm,3)),'+','LineWidth',1,'Color',color_scheme{k});hold on;
        errorbar(k+0.175,nanmean(betas_replay_norm(t_replay,k,:),3),std(betas_replay_norm(t_replay,k,:),[],3)./sqrt(size(betas_FLI_norm,3)),'*','LineWidth',2,'Color',color_scheme{k});
        [~,pval_paired(k)] = ttest(betas_FLI_norm(t_FLI,k,:),betas_replay_norm(t_replay,k,:),'tail','left');
        if pval_paired(k)<0.001 && k<6 % only plot first six states
            plot(k-0.15,markerpoint,'*k','LineWidth',1.5,'MarkerSize',10)
            plot(k+0.15,markerpoint,'*k','LineWidth',1.5,'MarkerSize',10)
        elseif pval_paired(k)<0.05 && k<6
            plot(k,markerpoint,'*k','LineWidth',1.5,'MarkerSize',10)
        end
    end
    clear h;
    h(1) = errorbar(k+6,nanmean(betas_replay_norm(t_replay,k,:),3),std(betas_replay_norm(t_replay,k,:),[],3)./sqrt(size(betas_replay_norm,3)),'*','LineWidth',2,'Color','black');
    h(2) = errorbar(k+5,nanmean(betas_FLI_norm(t_FLI,k,:),3),std(betas_FLI_norm(t_FLI,k,:),[],3)./sqrt(size(betas_FLI_norm,3)),'+','LineWidth',1,'Color','black');hold on;
    plot4paper('RSN-State','RSN-State probability');
    legend(h,{'Replay evoked','Stimulus evoked'},'Location','SouthWest')
    %title('Comparison of state distribution during encoding and replay');
    set(gca,'XTick',[1:hmm.K]);
    if whichstudy==1;
        set(gca,'YTick',[-0.2:0.1:0.2]);
    else
        set(gca,'YTick',[-0.2:0.04:0.2]);
    end
    grid on;
    xlim([0,hmm.K+1]);
    %ylim([-0.155,0.195]);
    print([savebase '_FLIReplayComparison_ste'],'-depsc')
    set(gcf,'Position',[440 475 560 323]);
    print([savebase '_FLIReplayComparison_ste_Long'],'-depsc')
    set(gcf,'Position',[440 235 445 563]);pause(0.1)
    print([savebase '_FLIReplayComparison_ste_tall'],'-depsc')
    
    
end

%% SECTION 3: WIDEBAND POWER AND COHERENCE MAPS

% load spectral info:
[psd,coh,f] = loadMTspect(studydir,K,template_string);
    
% check for bad subjects:
BadSubj = any(isnan(psd(:,:) + coh(:,:)),2);
psd(BadSubj,:,:,:,:) = [];
coh(BadSubj,:,:,:,:) = [];

nparcels = parc.n_parcels;
net_mean = zeros(nparcels,size(psd,2));
f_max = ceil(30*size(psd,3)/50);
f_band = f>1 & f<30; % corresponds to 1-30 Hz
for kk = 1:size(psd,2)
    tmp = squeeze(nanmean(nanmean(abs(psd(:,kk,f_band,:,:)),3),1) ); %note this includes all cross spectra; just averages over subjects and frequencies
    net_mean(:,kk) = (diag(tmp));
end

% unclear the below are any use:
for kk=1:K
       toplot = net_mean(:,kk);
       psdthresh = prctile(abs(toplot),50);
       CL = max(abs(toplot(:))) * [0 1];
       toplot(abs(toplot)<psdthresh) = NaN;
       f2 = plot_surf_summary_neuron(parc,toplot,0,false,'enclosing',[],[],CL);
       print([savebase '/St',int2str(kk),'_power'],'-depsc');
       close(f2);
end
%% remove noise with 2 mode NNMF:
nnmf_outfileWB = fullfile( hmmdir, ['embedded_HMM_K',int2str(K),template_string,'_nnmfWB']);
if ~isfile(nnmf_outfileWB)
    S=[];
    S.psds=psd;

    % for wideband
    S.maxP=2; S.maxPcoh=2; 

    nnmfWB_res = teh_spectral_nnmf(S);
    nnmfWB_res = rmfield(nnmfWB_res,'coh');

    save(nnmf_outfileWB,'nnmfWB_res')
    
    % Visualise the mode shapes
    for ii = 1:2
        figure('Position',[100 100*ii 256 256])
        h = area(nnmfWB_res.nnmf_coh_specs(ii,:));
        h.FaceAlpha = .5;
        h.FaceColor = [.5 .5 .5];
        grid on;axis('tight');
        set(gca,'YTickLabel',[],'FontSize',14);
        set(gca,'XTick',20:20:100);set(gca,'XTickLabel',[10:10:50]);
        print([savebase '/NNMFMode',int2str(ii)],'-depsc')
    end
else
    load( nnmf_outfileWB );
end
    
%%
% wideband power plots:

net_mean = zeros(nparcels,size(Gamma,2));
thresh_mean = zeros(size(Gamma,2),1);
for k = 1:size(Gamma,2)
    net_mean(:,k) = squeeze(nnmfWB_res.nnmf_psd_maps(k,1,:))';
end
    
% Plot 4 way power plots:
for kk=1:K
   toplot = net_mean(:,kk);%-mean(net_mean,2);
   psdthresh = prctile(abs(toplot),50);
   %CL = psdthresh*[-1.25,1.25];
   CL = max(abs(squash(net_mean(:,:)))) * [0 1];
   %CL = [min(toplot(:)), max(toplot(:))];
   toplot(abs(toplot)<psdthresh) = NaN;
   %f2 = parc.plot_surface(psds,0,false,'enclosing',CL);
   f2 = plot_surface_4way(parc,toplot,0,false,'enclosing',[],[],CL);
   print([savebase '/St',int2str(kk),'_power_WB_4wayplot'],'-depsc');
   close(f2);

end
    
% and plot colorbar separately:
figure();
cm = cat(2,linspace(.5, 1 ,63)',linspace(0, 1 ,63)',linspace(0, 0 ,63)');
c=flipud(cm);
cm = cat(1,[.6 .6 .6],cm);
colormap(cm);
caxis(CL);
colorbar;axis off

print([savebase '/St_colorbar_power_WB_4wayplot'],'-depsc');
    

%% Plot glass brain networks
    

for kk = 1:K
    figure();
    graph = abs(squeeze(nnmfWB_res.nnmf_coh_maps(kk,1,:,:)));
    tmp=squash(triu(graph));
    inds2=find(tmp>1e-10);
    data=tmp(inds2);

    S2=[];
    S2.data=squash(data);
    S2.do_fischer_xform=0;
    S2.pvalue_th=0.05/(38.^2); 
    S2.do_plots=0;
    graph_ggm=teh_graph_gmm_fit(S2);

    th=graph_ggm.normalised_th;
    graph=graph_ggm.data';

    if th<1.96 % less than 2 stds from mean
        graph(graph<th)=NaN;
        graphmat=nan(nparcels, nparcels);
        graphmat(inds2)=graph;
        graph=graphmat;
    else
        % few sparse connections, do not plot:
        graph = nan(nparcels);
    end
                
    if all(isnan(graph(:)))
        graph(1,1)=1;
    end
    parcelFile = fullfile(osldir,'parcellations',parc_file);
    spatialRes = 8; 
    spatialMap = nii.quickread(parcelFile, spatialRes); 
    mni_coords = find_ROI_centres(spatialMap, spatialRes, 0); 

    % plot 
    nROIs = size(graph,2);
    colorLims = [th th+1]; 
    sphereCols = repmat([30 144 255]/255, nROIs, 1); 
    edgeLims = [4 8];

    osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 

    view([0 90]);
    zoom(1);
    print([savebase '/WB_Coherence',int2str(kk),'_power'],'-depsc');

    viewZ = {[270,0],[-270,0],[0,90]};
    figure('Color', 'w','Position',[547 100 577 453]);
    ax(1) = axes('Position',[0 0.5 0.5 0.5]);
    ax(2) = axes('Position',[0.55 0.5 0.5 0.5]);
    ax(3) = axes('Position',[0.27 0.1 0.5 0.5]);
    
    % and plot 3 way brain graphs:
    for iplot=1:3
        axes(ax(iplot))
        osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 
        view(ax(iplot),viewZ{iplot})
        colorbar('hide') 
    end
    print([savebase '/WB_Coherence',int2str(kk),'_3way'],'-depsc');
    close all;
end
        
    

