function [fig1,FO_replay,GamBlock] = plotReplayFig1(Gamma,replayScores,t_window,goodsamples,triggerpoints,scan_T,pthreshold,method,set1_cols,states_to_emphasise)
% note last argument, method, denotes whether correlation is over ALL valid
% replay points (method==1); the FIRST replay point in a burst (method==2);
% or the LAST replay point in a burst (method==3).
if nargin<8
    method=1;
end
if nargin<10
    states_to_emphasise = [];
end

if nargin<9
    download_path = '/Users/chiggins/Documents/MATLAB/QuinnFrontiers/Quinn2018_TaskHMM/';
    addpath(genpath(download_path));
    set1_cols = utils.set1_cols;
end

nSj = length(goodsamples);

K=size(Gamma,2);

GamBlock = NaN(2*t_window+1,K,nSj);
if length(triggerpoints{1})>40000
    Fs = 250;
else
    Fs = 100;
end

t=[-t_window:t_window]./Fs;

% colors_alt=jet(K);
% for k=1:K
%     colors{k}=colors_alt(k,:);
% end
R = [[1;1+cumsum(scan_T(1:end-1)')],[cumsum(scan_T(1:end)')]];
if ~all(R(2:end,1)==R(1:end-1,2)+1)
    % hack to correct some uncorrected offsets in R:
    R(2:end,1) = R(1:end-1,2)+1;
end
for iSj=1:nSj
    % convert Gamma to padded timecourse
    Gam_long = NaN(length(goodsamples{iSj}),K);
    Gam_long(goodsamples{iSj},:) = Gamma(R(iSj,1):R(iSj,2),:);
    Gam_long = Gam_long(triggerpoints{iSj},:);
%     % epoch by replay onset:
    if method==1
         replaytimes = replayScores(iSj,:) > prctile(replayScores(iSj,:),99);
    elseif method==2
        replaytimes = replayScores(iSj,:) < prctile(replayScores(iSj,:),1);
    end
    % replaytimes = replayScores(iSj,:) < prctile(replayScores(iSj,:),1);
    %eliminate adjacent points
    replaytimes= [0,diff(replaytimes)]==1;
    %eliminate border points:
    replaytimes(1:t_window)=0;replaytimes([end-t_window]:end)=0;
     t_i = find(replaytimes);
   
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
    
    Gamma_mu = nanmean(Gam_long);
    X_subj = Gam_long-repmat(Gamma_mu,length(Gam_long),1);
    Gammat = NaN(2*t_window+1,K,length(t_i));
    for iEv=1:length(t_i)
        t_sample_g = [t_g(iEv)-t_window : t_g(iEv)+t_window];
        Gammat(:,:,iEv) = X_subj(t_sample_g,:);
        
    end
    
    GamBlock(:,:,iSj) = nanmean(Gammat,3);
    
    % also save the replay evoked state prob, per subject, for later use:
    FO_replay(:,:,iSj) = nanmean(Gammat+repmat(Gamma_mu,[2*t_window+1,1,length(t_g)]),3);
    %FO_replay(:,:,iSj) = nanmean(Gammat);
end

FO_replay = squeeze(nanmean(reshape(FO_replay,[2*t_window+1,K,2,nSj/2]),3));
GamBlock = squeeze(nanmean(reshape(GamBlock,[2*t_window+1,K,2,nSj/2]),3));
t_stats = nanmean(GamBlock,3)./(nanstd(GamBlock,[],3)./sqrt(size(GamBlock,3)));

fig1 = figure('Position',[440 235 549 563])
for k=K:-1:1
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
    mu=nanmean(GamBlock(:,k,:),3);
    ste=nanstd(GamBlock(:,k,:),[],3)./sqrt(size(GamBlock,3));
    shadedErrorBar(t,mu,ste,{'Color',set1_cols{k},'LineWidth',lw,'LineStyle',ls},1);
    hold on;
    h(k) = plot(NaN,NaN,'Color',set1_cols{k},'LineWidth',lw,'LineStyle',ls);
 end
%leg=[];for k=1:K;leg{k}=int2str(k);end


% compute and plot significance bars:
%subplot(2,length(K_options),length(K_options)+iK);
%significant_points=t_stats;

% for k=1:K; for j=1:length(t)
%     [hyp(j,k),pvalstest(j,k)] = ttest(GamBlock(j,k,:));
% end;end

% alternatively manual 1 sided ttest (sanity check: this gives the same as
% ttest with argument 'Tail' ,'right'
n=size(GamBlock,3);
hyp = tcdf((t_stats),n-1) > (1-pthreshold);

% OR do cluster sign flipping test:
for k=1:K
    dat = permute(GamBlock(:,k,:),[3,2,1]);
    thresh = 3;
    nP=5000;
    [corrp(:,k),tstats_cluster(:,k)] = osl_clustertf(dat,thresh,nP);
end
hyp = corrp>(1-1e-3);
ylow = ylim;
if method==1
    ylow(1) = min(squash(nanmean(GamBlock,3)))*1.15;
else
    ylow(1) = -0.0689;
end
yvals = ylow(1) * [1.05:0.05:2];%,1.1,1.15,1.15,1.15,1.20,1.25,1.3,1.35];
lengththresh = 2; % remove segments less than this length
for iSt=1:K
    if 1%iSt~=7
        sig_line = find(diff([0;hyp(:,iSt);0]));
       S = reshape(sig_line,[2,length(sig_line)/2]);
       S(:,diff(S)<lengththresh) = [];
       sig_line = S(:);
       sig_line(sig_line>length(hyp))=length(hyp);
        
        
%         smooth_onpoints = ([0;conv(hyp(:,iSt),ones(lengththresh,1),'same')>0;0]);
%         smooth_offpoints= ([1;conv(~hyp(:,iSt),ones(lengththresh,1),'same')>0;1]);
%         % remove short on segments:
%         sig_line_smoothed = find(diff([0;smooth_onpoints;0]));
%         sig_line_smoothed(1:2:end) = sig_line_smoothed(1:2:end)+lengththresh/2-1;
%         sig_line_smoothed(2:2:end) = sig_line_smoothed(2:2:end)-lengththresh/2;
%         
%         sig_line = sig_line_smoothed;
%         if length(sig_line)~=length(sig_line_smoothed)
%             new_hyp = zeros(size(hyp,1),1);
%             currentstate = false;
%             for i=1:length(sig_line)-1
%                 if all(hyp(sig_line(i):sig_line(i)+lengththresh,iSt)) == ~currentstate
%                     currentstate=~currentstate;
%                 end
%                 new_hyp(sig_line(i):sig_line(i+1)-1)=currentstate; 
%             end
%         end
%         for i=2:length(sig_line)
%             if sig_line(i)-sig_line(i-1)<lengththresh
%                 if all(hyp(sig_line(i-1):sig_line(i)-1,iSt))
%                     hyp(sig_line(i-1):sig_line(i),iSt)=0;
%                 else
%                     hyp(sig_line(i-1):sig_line(i),iSt)=1;
%                 end
%             end
%         end
%        sig_line = find(diff([0;hyp(:,iSt);0]));
%        sig_line(sig_line==(2*t_window+2))=2*t_window;
        
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
ylim([yval*1.05,ylow(2)]);
for k=1:K,h(k).DisplayName=['RSN-State ',int2str(k)];end
leg=legend(h,'Location','EastOutside');
plot4paper('Time (sec) from replay event','State probability')
grid on

end