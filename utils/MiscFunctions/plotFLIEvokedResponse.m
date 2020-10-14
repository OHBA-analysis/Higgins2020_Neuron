function [fig1,Regmodel_params_unnormal,Regmodel_params] = plotFLIEvokedResponse(datadirFLI,hmm_fli_file,mg,set1_cols,YL,states_to_emphasise);

if nargin < 5
    YL=[];
end
if nargin < 6
    states_to_emphasise = [];
end
% load in alignment data:
[~,~,triggerpoints,goodsamples] = getSubjectMasks(datadirFLI);
load(hmm_fli_file);
if contains(hmm_fli_file,'Study1')  
    main_hmmfile = strrep(hmm_fli_file,'Study1_FLI','Study1');
    stimmarker = 'UPPT001_down';
    numsess = 2;
    nSj=21;
else
    main_hmmfile = strrep(hmm_fli_file,'Study2_FLI','Study2');
    stimmarker = 'FIL_UPPT001_down';
    numsess = 2;
    nSj=22;
end
load(main_hmmfile,'new_state_ordering');
hmm = hmm_permutestates(hmm,new_state_ordering);
t_ind=0;
for isub = 1:(nSj*numsess)
    D = spm_eeg_load(hmm.data_files{isub});
    ev = D.events;
    %Gamma_subj = hmm.gamma(subj_inds==isub,:);
    % warning: subject inds are not correct!
    GS = good_samples(D);
    Gamma_subj = zeros(size(GS,1),hmm.K);
    Gamma_subj(GS,:) = hmm.gamma(t_ind+[1:sum(GS)],:);
    t_ind = t_ind+sum(GS);
    mask = false(1,length(ev));
    for i=1:length(ev)
        mask(i) = strcmp(ev(i).type,stimmarker) && isnumeric(ev(i).value) && mod(ev(i).value,10)==1;
    end
    ntr = sum(mask); % number of trials
    tr_inds = find(mask);
    Gam_subj_evoked = nan(4*hmm.fsample+1,hmm.K,ntr);
    for i=1:ntr
        t = round(ev(tr_inds(i)).time*hmm.fsample);
        if t-3*hmm.fsample>1 & t+hmm.fsample<length(Gamma_subj)
            Gam_subj_evoked(:,:,i) = Gamma_subj(t-3*hmm.fsample : t+hmm.fsample,:);
        end
    end
    if nargin<3 || isempty(mg)
        Regmodel_params(:,:,isub) = nanmean(Gam_subj_evoked,3) - nanmean(Gamma_subj);
        %Regmodel_params(:,:,isub) = nanmean(Gam_subj_evoked,3) - nanmean(Gam_subj_evoked(751,:,:),3);
        %Regmodel_params(:,:,isub) = nanmean(Gam_subj_evoked,3) - nanmean(nanmean(Gam_subj_evoked,3),1);
    else
        Regmodel_params(:,:,isub) = nanmean(Gam_subj_evoked,3) - mg(isub,:);
    end
    Regmodel_params_unnormal(:,:,isub) = nanmean(Gam_subj_evoked,3);
end
%%
%figure();plot(mean(Regmodel_params,3),'LineWidth',2)
Regmodel_params = squeeze(nanmean(reshape(Regmodel_params,[4*hmm.fsample+1,hmm.K,numsess,isub/numsess]),3));
t_stats = nanmean(Regmodel_params,3)./(nanstd(Regmodel_params,[],3)./sqrt(size(Regmodel_params,3)));
Regmodel_params_unnormal= squeeze(nanmean(reshape(Regmodel_params_unnormal,[4*hmm.fsample+1,hmm.K,numsess,isub/numsess]),3));

fig1 = figure('Position',[440 235 549 563])
% colors_alt=jet(hmm.K);
% for k=1:hmm.K
%     colors{k}=colors_alt(k,:);
% end
t = linspace(-3,1,size(Regmodel_params,1));
for k=1:hmm.K
     if mod(k,2)==1
         ls = '-';
     else
         ls='--';
     end
    if ~isempty(states_to_emphasise) && any(ismember(states_to_emphasise,k))
        lw = 2.5;
    else
        lw = 1.5;
    end
    %plot(nanmean(GamBlock(:,k,:),3),'LineWidth',2,'LineStyle',ls);hold on;
    mu=nanmean(Regmodel_params(:,k,:),3);
    ste=nanstd(Regmodel_params(:,k,:),[],3)./sqrt(size(Regmodel_params,3));
    shadedErrorBar(t,mu,ste,{'Color',set1_cols{k},'LineWidth',lw,'LineStyle',ls},1);
    hold on;
    h(k) = plot(NaN,NaN,'Color',set1_cols{k},'LineWidth',lw,'LineStyle',ls);
end
 
n=size(Regmodel_params,3);
pthreshold = 0.05/hmm.K;
hyp = tcdf((t_stats),(n)-1) > (1-pthreshold);

% OR do cluster sign flipping test:
t_included = 751-125 + [0:250];
for k=1:hmm.K
    dat = permute(Regmodel_params(t_included,k,:),[3,2,1]);
    thresh = 3;
    nP=5000;
    [corrp(t_included,k),tstats_cluster(t_included,k)] = osl_clustertf(dat,thresh,nP);
end
hyp = corrp>(1-1e-3);

if ~isempty(YL)
    ylim(YL);
end
ylow = ylim;
if isempty(YL)
    ylow(1) = min(squash(nanmean(Regmodel_params,3)))*1.1;
    yvals = ylow(1) * [1.05:0.05:2];%,1.1,1.15,1.15,1.15,1.20,1.25,1.3,1.35];
else
    yvals = ylow(1) * [0.64:0.03:1];
end
t_window = (size(Regmodel_params,1)-1)/2;
lengththresh = 2;
for iSt=1:hmm.K
    if 1%iSt~=7
        sig_line = find(diff([0;hyp(:,iSt);0]));
        S = reshape(sig_line,[2,length(sig_line)/2]);
        S(:,diff(S)<lengththresh) = [];
        sig_line = S(:);
        sig_line(sig_line>length(hyp))=length(hyp);
        
        %yval = [ylow-0.05+iSt*0.005];
        %yval = [ylow(1)*(1+0.05*iSt)];
        yval = yvals(iSt);
        for k=1:length(sig_line)/2
            %if [sig_line(k*2)-sig_line(k*2-1)]>8
                line([t(sig_line(k*2-1)),t(sig_line(k*2))],[yval,yval],'Color',set1_cols{iSt},'LineWidth',2);
            %end
        end
    end
end
if isempty(YL)
    ylim([yval(1)*1.05,ylow(2)]);
end
for k=1:hmm.K,h(k).DisplayName=['State ',int2str(k)];end
leg = legend(h,'Location','EastOutside');
plot4paper('Time (sec) from stimulus presentation','State probability')
grid on

end