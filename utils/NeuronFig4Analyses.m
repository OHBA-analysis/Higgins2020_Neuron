%% Replay Paper Figure 4: State spectral profiles

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

savebase = fullfile( [wd,session_name{whichstudy},'Fig4_',template_string,'/' ])
if ~exist( savebase )
    mkdir(savebase);
end

parc_file = ['fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
parc = parcellation(parc_file);
statelabels={'1','2','3','4','5','6','7','8','9','10','11','12'};

%% SECTION 1: STATE SPECIFIC FREQUENCY BREAKDOWN:

% first we must check that the soft computed frequency spectra have been 
% computed, and if not, compute them:
mtfilename = [studydir, '/hmm_1to45hz/hmm',template_string,'_parc_giles_symmetric__pcdim80_voxelwise_embed14_K',int2str(K),'_big1_dyn_modelhmm_store/',...
    'state_netmats_mtsess_2_vn0_soft_global0.mat'];
if ~exist(mtfilename)
    error('Soft state timecourses not found - rerun this analysis!');
end
    
[psd,coh] = loadMTspect(studydir,K,template_string);

% check for bad subjects:
BadSubj = any(isnan(psd(:,:) + coh(:,:)),2);
psd(BadSubj,:,:,:,:) = [];
coh(BadSubj,:,:,:,:) = [];

nparcels = parc.n_parcels;

%% plot psd and coh as function of freq averaged over parcels:

f_band=1:size(psd,3);
offdiags = eye(nparcels)==0;
psd_all = zeros(K,length(f_band));
coh_all= zeros(K,length(f_band));
psd_cross= zeros(K,length(f_band));
for kk=1:K
    G = squeeze( mean(abs(coh(:,kk,f_band,offdiags)),4));
    P = squeeze(mean(abs(psd(:,kk,f_band,~offdiags)),4));
    P_off=squeeze(mean(abs(psd(:,kk,f_band,offdiags)),4));
    coh_all(kk,:) = mean(G,1);
    coh_ste(kk,:) = std(G,[],1)./sqrt(length(G));
    psd_all(kk,:) = mean(P,1);
    psd_ste(kk,:) = std(P,[],1)./sqrt(length(P));
    psd_cross(kk,:) = mean(P_off,1);
    psd_c_ste(kk,:) = std(P_off,[],1)./sqrt(length(P_off));
end
colorscheme = set1_cols();

figure('Position',[440 508 708 290]);
for k=1:K
    if mod(k,2)==1
         ls{k} = '-';
     else
         ls{k}='--';
    end
    shadedErrorBar(0.5*f_band,psd_all(k,:),psd_ste(k,:),{'LineWidth',2,'LineStyle',ls{k},'Color',colorscheme{k}},1);hold on;
    statelabels{k} = ['State ',int2str(k)];
    h(k) = plot(NaN,NaN,'Color',colorscheme{k},'LineWidth',2,'LineStyle',ls{k});
end
grid on;
title('PSD per state');
plot4paper('Frequency');
for k=1:K,h(k).DisplayName=['State ',int2str(k)];end
leg=legend(h,'Location','EastOutside');
print([savebase '/CrossStatePSD_long'],'-depsc')

figure('Position',[440 508 708 290]);
for k=1:K;
    if mod(k,2)==1
         ls{k} = '-';
     else
         ls{k}='--';
    end
    shadedErrorBar(0.5*f_band,coh_all(k,:),coh_ste(k,:),{'LineWidth',2,'LineStyle',ls{k},'Color',colorscheme{k}},1);hold on;
    h(k) = plot(NaN,NaN,'Color',colorscheme{k},'LineWidth',2,'LineStyle',ls{k});
end
grid on;title('Coherence per state');
plot4paper('Frequency');
for k=1:K,h(k).DisplayName=['State ',int2str(k)];end
leg=legend(h,'Location','EastOutside');
print([savebase '/CrossStateCOH_long'],'-depsc')

%% Infer NNMF spectral decomposition:
    
% Spectral Mode NNMF
nnmf_outfile = fullfile( savebase, ['embedded_HMM_K',int2str(K),template_string,'_nnmf']);
    
% Compute the NNMF
if ~isfile([nnmf_outfile,'.mat'])

    S = [];
    S.psds = psd(:,:,:,:,:);
    S.maxP=4;
    S.maxPcoh=4;
    S.do_plots = true;
    nnmf_res  = run_nnmf( S, 10 );
    nnmf_res = rmfield(nnmf_res,'coh');

    save(nnmf_outfile,'nnmf_res')

    % Visualise the mode shapes
    for ii = 1:4
        figure('Position',[100 100*ii 256 256])
        h = area(nnmf_res.nnmf_coh_specs(ii,:));
        h.FaceAlpha = .5;
        h.FaceColor = [.5 .5 .5];
        grid on;axis('tight');
        set(gca,'YTickLabel',[],'FontSize',14);
        set(gca,'XTick',20:20:100);set(gca,'XTickLabel',[10:10:50]);
        print([savebase '/NNMFMode',int2str(ii)],'-depsc')
    end
else
    load( nnmf_outfile );
end

%% Spectral Mode Power Plots

net_mean = zeros(nparcels,4,size(Gamma,2));
for k = 1:size(Gamma,2)
    diagselectormatrix = find(eye(nparcels));
    for jj=1:3
        temp = squeeze(mean(psd(:,k,:,diagselectormatrix),1)) ...
            .*repmat(nnmf_res.nnmf_coh_specs(jj,:)',[1,nparcels]);
        net_mean(:,jj,k) = sum(temp,1);
    end
end
    
% Plot 4 way power plots:
for kk=1:K
    for jj=1:3
       
       
       toplot = net_mean(:,jj,kk);
       psdthresh = prctile(abs(toplot),90);
       CL = max(abs(squash(net_mean(:,jj,:)))) * [0 1];
       toplot(abs(toplot)<psdthresh) = NaN;
       f2 = plot_surface_4way(parc,toplot,0,false,'enclosing',[],[],CL);
       print([savebase '/St',int2str(kk),'_power_nnmfmode_',int2str(jj),'_10prctile'],'-depsc');
       subplot(1,2,2);
       c1 = colorbar;caxis(CL);
       set(c1,'Ticks',c1.Ticks([1,end]));
       set(gcf,'Position',[547 360 577 193])
       print([savebase '/St',int2str(kk),'_power_nnmfmode_',int2str(jj),'_10prctile_cbar'],'-depsc');
       close(f2);
    end
end
    
%% Summary plots of network coherences
num_nodes=parc.n_parcels;
S = [];
S.psds = psd(:,:,:,:,:);
S.maxP=4;
S.maxPcoh=4;
S.do_plots = false;
for kk=[1:K]        
    viewZ = {[270,0],[-270,0],[0,90]};
    for pp=1:S.maxPcoh-1
        figure('Color', 'w','Position',[547 100 577 453]);
        ax(1) = axes('Position',[0 0.5 0.5 0.5]);
        ax(2) = axes('Position',[0.55 0.5 0.5 0.5]);
        ax(3) = axes('Position',[0.27 0.1 0.5 0.5]);
        graph=abs(squeeze(nnmf_res.nnmf_coh_maps(kk,pp,:,:)));
        tmp=squash(triu(graph));
        inds2=find(tmp>1e-10);
        data=tmp(inds2);


        S2=[];
        S2.data=squash(data);
        S2.do_fischer_xform=false;
        S2.pvalue_th=0.05./(nparcels.^2); 
        S2.do_plots=0;
        graph_ggm=teh_graph_gmm_fit(S2);

        th=graph_ggm.normalised_th;
        graph=graph_ggm.data';

        graphrec{kk,pp} = graph_ggm;
        
        if th<1.96 % less than 2 stds from mean
            graph(graph<th)=NaN;
            graphmat=nan(nparcels, nparcels);
            graphmat(inds2)=graph;
            graph=graphmat;
        else
            % few sparse noisy connections, do not plot:
            graph = nan(nparcels);
            graph(1,1)=1;
        end

        % get node co-ordinates parcellation
        parcelFile = fullfile(osldir,'parcellations',parc_file);
        spatialRes = 8; 
        spatialMap = nii.quickread(parcelFile, spatialRes); 
        mni_coords = find_ROI_centres(spatialMap, spatialRes, 0); 

        % plot 
        nROIs = size(graph,2);
        colorLims = [th-1 th+1]; 
        sphereCols = repmat([30 144 255]/255, nROIs, 1); 
        edgeLims = [4 8];
        for iplot=1:3
            axes(ax(iplot))
            osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 
            view(ax(iplot),viewZ{iplot})
            colorbar('hide') 
        end
        print([savebase '/CoherencePlot_K',int2str(kk),'_mode',int2str(pp)],'-depsc')
    close all;  
    end

end   
    
    
%% PSD and Coherence scatter plots:
I = logical(eye(nparcels));
for pp=1:3
    figure('Position',[702 470 453 320]);
    for kk=1:K
        scatter_coh(:,kk) = abs(squeeze(sum(nnmf_res.nnmf_coh_maps(kk,pp,:,:),4)))/nparcels;
        scatter_psd(:,kk) = abs(squeeze(nnmf_res.nnmf_psd_maps(kk,pp,:)));
        scatter(scatter_psd(:,kk),scatter_coh(:,kk),20,'MarkerFaceColor',colorscheme{kk});hold on;
        scatterplotlabels{kk} = ['RSN-State ',int2str(kk)];
    end
    xlim([0,0.1]);ylim([0,0.01]);
    plot4paper('PSD','Coherence');
    box on;axis square;
    legend(scatterplotlabels,'Location','EastOutside');
    Ylabel = get(gca,'YTick');
    for i=1:length(Ylabel);Ylabelstr{i} = num2str(Ylabel(i));end
    set(gca,'YTickLabel',Ylabelstr);
    print([savebase,'ScatterPlot_mode',int2str(pp)],'-depsc');
    %close all;

    pval_anova_nb_coh(pp) = anova1(scatter_coh)
    pval_anova_nb_psd(pp) = anova1(scatter_psd)
    for kk=1:K
        [~,pval_ttest_nb_coh(pp,kk)] = ttest(scatter_coh(:,kk),mean(scatter_coh(:,setdiff([1:K],kk)),2),'Tail','right');
        [~,pval_ttest_nb_psd(pp,kk)] = ttest(scatter_psd(:,kk),mean(scatter_psd(:,setdiff([1:K],kk)),2),'Tail','right');
    end
end

% and same for wideband:
nnmf_outfileWB = fullfile( hmmdir, ['embedded_HMM_K',int2str(K),template_string,'_nnmfWB']);
load( nnmf_outfileWB );
%figure('Position',[440 519 408 279]);
 figure('Position',[702 470 453 320]);
for kk=1:K
    scatter_coh(:,kk) = abs(squeeze(sum(nnmfWB_res.nnmf_coh_maps(kk,1,:,:),4)))/nparcels;
    scatter_psd(:,kk) = abs(squeeze(nnmfWB_res.nnmf_psd_maps(kk,1,:)));
    scatter(scatter_psd(:,kk),scatter_coh(:,kk),20,'MarkerFaceColor',colorscheme{kk});hold on;
    scatterplotlabels{kk} = ['RSN-State ',int2str(kk)];
end
plot4paper('PSD','Coherence');
xlim([0,0.1]);ylim([0,0.01]);
box on;
legend(scatterplotlabels,'Location','EastOutside');
Ylabel = get(gca,'YTick');
for i=1:length(Ylabel);Ylabelstr{i} = num2str(Ylabel(i));end
set(gca,'YTickLabel',Ylabelstr);
axis square;
print([savebase,'ScatterPlot_wideband'],'-depsc');
%close all;

% compute stats:

pval_anova_wb_coh = anova1(scatter_coh)
pval_anova_wb_psd = anova1(scatter_psd)
for kk=1:K
    [~,pval_ttest_wb_coh(k)] = ttest2(scatter_coh(:,kk),squash(scatter_coh(:,setdiff([1:K],kk))));
    [~,pval_ttest_wb_psd(k)] = ttest2(scatter_psd(:,kk),squash(scatter_psd(:,setdiff([1:K],kk))));
end

%% PART 2: HMM-REGULARISED REPLAY SPECTRAL INFO

% extract instantaneous and time-aligned Fractional Occupancies, per
% subject:
if whichstudy >1
    
    betadistfile = fullfile( strrep(studydir,'250Hz',''), ['Fig2_',template_string ],'/betas_replayFLI.mat'); 
    load(betadistfile);
    nROIs = 38;

    for i=1:2
        if i==1
            % Find replay evoked stats using this study's learned spectra
            spect_dir = [wd,session_name{whichstudy},'250Hz/'];
            savebase_rep = [savebase,'ReplayEvoked_ReplayRSSpectra/'];
            mkdir(savebase_rep);
            template_stringtemp = template_string;
            nnmf_dir = [wd,session_name{whichstudy},'Fig4_',template_string];
        else
            % Find replay evoked stats using the canonical spectra
            spect_dir = [wd,session_name{1},'250Hz/'];
            savebase_rep = [savebase,'ReplayEvoked_CanonicalSpectra/'];
            template_stringtemp = int2str(bestmodel);
            mkdir(savebase_rep);
            nnmf_dir = [wd,session_name{1},'Fig4_',int2str(bestmodel)];
        end
        
        % set colorbar limits - limits determined by basis set used
        %if whichstudy==1
        CA1 = [-0.0005    0.0033;...
            -0.0005    0.0031 ];
        CA2 = [-0.0042    0.0034;...
            -0.0024    0.0066];
        %end
        if whichstudy==3 && i==1
            CA2(1,:) = [-0.0024    0.002];
        end
        % Plot replay evoked PSD and coherence, over all channels:
        [psd,coh,f] = loadMTspect(spect_dir,K,template_stringtemp);  
        
        if i==1
            % account for multiple runs per subject:
            psd = squeeze(mean(reshape(abs(psd),[2,nSj{whichstudy},K,size(psd,3),nROIs,nROIs]),1));
            coh = squeeze(mean(reshape(abs(coh),[2,nSj{whichstudy},K,size(psd,3),nROIs,nROIs]),1));
        else
            psd = mean(abs(psd),1);
            coh = mean(coh,1);
        end
        
        psd_zscore = demean(psd,2);
        coh_zscore = demean(coh,2);
        
        I=logical(eye(nparcels));
        f_n = size(psd_zscore,3);
        t_window=125;
        
        for iFLI=1:2
            if iFLI==1
                FO_replay = betas_replay;
                t_selected = 1:size(betas_replay,1);
                idstring='_Replay';
            else
                FO_replay = betas_FLI(801-t_window:801+t_window,:,:);
                t_selected = 1:size(betas_replay,1);
                idstring='_FLI';
            end
        
            clear psd_reg coh_reg CL
            for it=t_selected
                temp = mean(psd_zscore(:,:,:,I),4).*repmat(permute(FO_replay(it,:,:),[3,2,1]),[1,1,f_n]);
                psd_reg(it,:,:,:) = squeeze(sum(temp,2));
                temp = mean(coh_zscore(:,:,:,~I),4).*repmat(permute(FO_replay(it,:,:),[3,2,1]),[1,1,f_n]);
                coh_reg(it,:,:,:) = squeeze(sum(temp,2));
            end
            psd_reg = mean(psd_reg,2);
            coh_reg = mean(coh_reg,2);
            %
            figure('Position',[440 439 942 359]);
            subplot(1,2,1);
            imagesc(flipud(squeeze(psd_reg)'));colorbar;
            colormap(hot);
            if iFLI==1
                set(gca,'XTick',[1:(t_window/2):(2*t_window)+1]);
                set(gca,'XTickLabel',[-t_window:(t_window/2):(2*t_window)]/sample_rate);
                caxis(CA1(i,:));
            else
                caxis(CA1(i,:));
                set(gca,'XTick',[1:(t_window/2):(2*t_window)+1]);
                set(gca,'XTickLabel',[-0.3:0.25:0.8]);
            end
            set(gca,'YTick',[1,9:10:89]);
            set(gca,'YTickLabel',fliplr([0:5:45]));
            axis square;
            plot4paper('Time from Replay Event','Frequency');
            title('PSD');
            subplot(1,2,2);
            imagesc(flipud(squeeze(coh_reg)'));colorbar;
            colormap(hot);
            if iFLI==1
                set(gca,'XTick',[1:(t_window/2):(2*t_window)+1]);
                set(gca,'XTickLabel',[-t_window:(t_window/2):(2*t_window)]/sample_rate);
                caxis(CA2(i,:));
            else
                caxis(CA2(i,:));
                set(gca,'XTick',[1:(t_window/2):(2*t_window)+1]);
                set(gca,'XTickLabel',[-0.3:0.25:0.8]);
            end
            set(gca,'YTick',[1,9:10:89]);
            set(gca,'YTickLabel',fliplr([0:5:45]));
            axis square;
            plot4paper('Time from Replay Event','Frequency');
            title('Coherence');
            print([savebase_rep '/ReplayEvokedPSD_COH_nozscore',idstring],'-depsc');
            
            % Now plot replay regularised coherence and power:
            % use normalised (relative) power and coherence:
            if iFLI==1
                FO_replay = betas_replay_norm;
            else
                FO_replay = betas_FLI_norm(801-t_window:801+t_window,:,:);
            end

            t_replay = t_window+1;

            temp = psd(:,:,:,I).*repmat(permute(FO_replay(t_replay,:,:),[3,2,1]),[1,1,f_n]);
            temp = squeeze(sum(temp,2));
            psd_replay = squeeze(mean(temp,1));


            temp = coh(:,:,:,~I).*repmat(permute(FO_replay(t_replay,:,:),[3,2,1]),[1,1,f_n]);
            temp = squeeze(sum(temp,2));
            coh_replay = squeeze(mean(temp,1));

            % load nnmf results:
            nnmf_outfile = fullfile(nnmf_dir, ['embedded_HMM_K',int2str(K),template_stringtemp,'_nnmf']);
            load( nnmf_outfile );
            NN = nnmf_res;clear nnmf_res;
            
            psd_replay_wb = pinv(NN.nnmf_coh_specs')*psd_replay;
            CL{1} = [0    0.0009158; 0    0.0019; 0    0.0001535];
            CL{2} = [0    0.0010; 0    0.0013; 0    0.0002141];
            if whichstudy==3 & i==1;
                CL{1} = [0    0.0005352;0    0.0011;0    0.00006280];
            end
            psd_replay_wb = pinv(NN.nnmf_coh_specs')*psd_replay;
            for ii=1:size(NN.nnmf_coh_specs,1)-1
                toplot = psd_replay_wb(ii,:);
                toplot(toplot<prctile(toplot,90))=NaN;
                toplot(toplot<0)=realmin;
                %CL = max(abs(squash(psd_replay))) * [0 1];
                %f2 = plot_surface_4way(parc,psd_replay(25,:),0,false,'enclosing',[],[],[]);
                f2 = plot_surf_summary_neuron(parc,toplot,0,false,'enclosing',[],[],CL{i}(ii,:));
                print([savebase_rep '/PSDPlot_Replay_mode',int2str(ii),'_10prctile',idstring],'-depsc')
            end

            
            coh_replay_wb = pinv(NN.nnmf_coh_specs')*coh_replay;
            viewZ = {[270,0],[-270,0],[0,90]};
            for pp=1:size(NN.nnmf_coh_specs,1)-1

                toplot=coh_replay_wb(pp,:);

                graph = zeros(nparcels);graph(~I)=toplot;

                figure('Color', 'w','Position',[547 100 577 453]);
                ax(1) = axes('Position',[0 0.5 0.5 0.5]);
                ax(2) = axes('Position',[0.55 0.5 0.5 0.5]);
                ax(3) = axes('Position',[0.27 0.1 0.5 0.5]);
                tmp=squash(triu(graph));
                inds2=find(tmp>1e-10);
                data=tmp(inds2);
                
                if length(data)>0
                    S2=[];
                    S2.data=squash(data);
                    S2.do_fischer_xform=false;
                    S2.pvalue_th=0.05./(nparcels.^2); 
                    S2.do_plots=0;
                    graph_ggm=teh_graph_gmm_fit(S2);

                    if iFLI==1
                        th(pp)=graph_ggm.orig_th;
                    end
                    graph=S2.data';
                else
                    graph = nan(nparcels);
                    inds2 = find(ones(nparcels));
                end

                if th(pp)<1.96 % less than 2 stds from mean
                    graph(graph<th(pp))=NaN;
                    graphmat=nan(nparcels, nparcels);
                    graphmat(inds2)=graph;
                    graph=graphmat;
                else
                    % few sparse connections, do not plot:
                    graph = nan(nparcels);
                end

                % get node co-ordinates parcellation
                parcelFile = fullfile(osldir,'parcellations',parc_file);
                spatialRes = 8; 
                spatialMap = nii.quickread(parcelFile, spatialRes); 
                mni_coords = find_ROI_centres(spatialMap, spatialRes, 0); 

                % plot 
                nROIs = size(graph,2);
                colorLims = [th(pp) th(pp)*1.2]; 
                sphereCols = repmat([30 144 255]/255, nROIs, 1); 
                edgeLims = [4 8];
                for iplot=1:3
                    axes(ax(iplot))
                    osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 
                    view(ax(iplot),viewZ{iplot})
                    colorbar('hide') 
                end
                print([savebase_rep '/CoherencePlot_Replay_mode',int2str(pp),idstring],'-depsc')
                close all;  
            end
        end
    end
end
