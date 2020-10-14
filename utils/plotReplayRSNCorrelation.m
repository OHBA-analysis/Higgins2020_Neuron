function fig1 = plotReplayRSNCorrelation(Gamma,replayScores,t_window,s,goodsamples,triggerpoints,R,pthreshold,method,classifiersregression)
% note last argument, method, denotes whether correlation is over ALL valid
% replay points (method==1); the FIRST replay point in a burst (method==2);
% or the LAST replay point in a burst (method==3).
if nargin<9
    method=1;
end
if nargin<10
    classifiersregression=0;
end
% also load raw scores:
load('/Volumes/CamsHD2/YunzheData/Replaydata4Cam/Analysis/classifiers/YZ_sigmoidtransformscores/RedAll_1.mat')
scores = permute(RedAll_1,[3,4,2,1]);
scores = permute(scores(:,:,:),[3,1,2]);
K=size(Gamma,2);

GamBlock = NaN(2*t_window+1,K,s.n);
if length(triggerpoints{1})>40000
    Fs = 250;
else
    Fs = 100;
end
t=[-t_window:t_window]./Fs;
colors = utils.set1_cols;
colors_alt=jet(K);
for k=1:K
    colors{k}=colors_alt(k,:);
end

if ~all(R(2:end,1)==R(1:end-1,2)+1)
    % hack to correct some uncorrected offsets in R:
    R(2:end,1) = R(1:end-1,2)+1;
end
for iSj=1:s.n
    % convert Gamma to padded timecourse
    Gam_long = NaN(length(goodsamples{iSj}),K);
    Gam_long(goodsamples{iSj},:) = Gamma(R(iSj,1):R(iSj,2),:);
    Gam_long = Gam_long(triggerpoints{iSj},:);
%     % epoch by replay onset:
     replaytimes = replayScores(iSj,:) > prctile(replayScores(iSj,:),99);
    %eliminate adjacent points
    replaytimes= [0,diff(replaytimes)]==1;
    %eliminate border points:
    replaytimes(1:t_window)=0;replaytimes([end-t_window]:end)=0;
%     t_i = find(replaytimes);
    %extract epochs:
    if method==1
        t_i = findBurst(replaytimes,'all');
    elseif method==2
        t_i = findBurst(replaytimes,'first',200);
        %Smooth gamma for longer timescales:
        for k=1:size(Gam_long,2),Gam_long(:,k)=conv(Gam_long(:,k),ones(100,1),'same');end
        for k=1:size(Gam_long,2),Gam_long(:,k)=fliplr(conv(fliplr(Gam_long(:,k)),ones(100,1),'same'));end
        Gam_long = (Gam_long - repmat(nanmean(Gam_long),length(Gam_long),1)) ./ repmat(nanstd(Gam_long),length(Gam_long),1); 
    elseif method==3
        t_i = findBurst(replaytimes,'last',200);
        for k=1:size(Gam_long,2),Gam_long(:,k)=conv(Gam_long(:,k),ones(100,1),'same');end
        for k=1:size(Gam_long,2),Gam_long(:,k)=fliplr(conv(fliplr(Gam_long(:,k)),ones(100,1),'same'));end
        Gam_long = (Gam_long - repmat(nanmean(Gam_long),length(Gam_long),1)) ./ repmat(nanstd(Gam_long),length(Gam_long),1); 
    end
    
    %upsample for state timecourses at higher sampling rate than replay
    %scores:
    if length(Gam_long)>length(replaytimes) 
        Q = length(Gam_long) / length(replaytimes);
        t_g = round(t_i*Q);
    else
        Q = 1;
        t_g = t_i;
    end
    
    if any(t_g>size(Gam_long,1)-t_window)
        t_g(t_g>size(Gam_long,1)-t_window) = [];
        t_i(t_i>size(Gam_long,1)-t_window) = [];
    end
    
    if classifiersregression
        %regress out effect of raw classifiers across whole state timecourse:
        Gam_long(isnan(Gam_long))=0;
        DM = [squeeze(scores(iSj,1:length(Gam_long),:)),ones(length(Gam_long),1)];
        betas = pinv(DM)*Gam_long;
        Gam_long = Gam_long - DM*betas;
    end
    
    X_subj = Gam_long-repmat(nanmean(Gam_long),length(Gam_long),1);
    Gammat = NaN(2*t_window+1,K,length(t_i));
    scoremat = NaN(2*t_window+1,8,length(t_i));
    for iEv=1:length(t_i)
        t_sample_g = [t_g(iEv)-t_window : t_g(iEv)+t_window];
        Gammat(:,:,iEv) = X_subj(t_sample_g,:);
        if Q==1 && classifiersregression
            scoremat(:,:,iEv) = squeeze(scores(iSj,t_sample_g,:));
        else
            % do not do for now:
            scoremat(:,:,iEv) = zeros(2*t_window+1,8);
        end
    end
    if 1
        GamBlock(:,:,iSj) = nanmean(Gammat,3);
    else
        % regress out raw classifier score affect at each timepoint:
        Gammat(isnan(Gammat))=0;
        for it=1:2*t_window+1
            %DM = [ones(length(t_i),1),(squeeze(scoremat(it,:,:))')];
            DM = [ones(length(t_i),1),max(squeeze(scoremat(it,:,:))',[],2)];
            betas = pinv(DM)*(squeeze(Gammat(it,:,:)))';
            GamBlock(it,:,iSj) = betas(1,:);
        end
    end
end
t_stats = nanmean(GamBlock,3)./(nanstd(GamBlock,[],3)./sqrt(size(GamBlock,3)));

fig1 = figure('Position',[440 235 549 563])
for k=1:K
     if k<=K/2
         ls = '-';
     else
         ls='--';
    end
    %plot(nanmean(GamBlock(:,k,:),3),'LineWidth',2,'LineStyle',ls);hold on;
    mu=nanmean(GamBlock(:,k,:),3);
    ste=nanstd(GamBlock(:,k,:),[],3)./sqrt(size(GamBlock,3));
    shadedErrorBar(t,mu,ste,{'Color',colors{k},'LineWidth',1.5,'LineStyle',ls},1);
    hold on;
    h(k) = plot(NaN,NaN,'Color',colors{k},'LineWidth',1.5,'LineStyle',ls);
 end
%leg=[];for k=1:K;leg{k}=int2str(k);end


% compute and plot significance bars:
%subplot(2,length(K_options),length(K_options)+iK);
%significant_points=t_stats;

% for k=1:K; for j=1:length(t)
%     [hyp(j,k),pvalstest(j,k)] = ttest(GamBlock(j,k,:));
% end;end

% alternatively manual 2 sided ttest:
n=size(GamBlock,3);
hyp = tcdf(abs(t_stats),n-1) > (1-pthreshold/2);

ylow = ylim;
ylow(1) = min(squash(nanmean(GamBlock,3)))*1.05;
yvals = ylow(1) * [1.05:0.05:2];%,1.1,1.15,1.15,1.15,1.20,1.25,1.3,1.35];
for iSt=1:K
    if 1%iSt~=7
        sig_line = find(diff([0;hyp(:,iSt);0]));
        sig_line(sig_line==(2*t_window+2))=2*t_window;
        
        %yval = [ylow-0.05+iSt*0.005];
        %yval = [ylow(1)*(1+0.05*iSt)];
        yval = yvals(iSt);
        for k=1:length(sig_line)/2
            %if [sig_line(k*2)-sig_line(k*2-1)]>8
                line([t(sig_line(k*2-1)),t(sig_line(k*2))],[yval,yval],'Color',colors{iSt},'LineWidth',2);
            %end
        end
    end
end
ylim([yval*1.05,ylow(2)]);
for k=1:K,h(k).DisplayName=['State ',int2str(k)];end
leg=legend(h,'Location','EastOutside');
plot4paper('Time (sec) from replay event','State probability')
grid on
end