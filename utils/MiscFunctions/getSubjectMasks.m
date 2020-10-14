function [mask_AtoB,mask_BtoA,mask_boundaries,goodSamples] = getSubjectMasks(dir)
% For the first session first resting state sequence, this returns a mask
% for data that omits bad segments and out of trial sections, plus a window
% as specified by tforward and tback
rejectshortsegments=false;
if contains(dir,'Study2_FLI')
    whichstudy=2;
    nSj=66;
elseif contains(dir,'Study2')
    whichstudy=2;
    nSj=44;
elseif contains(dir,'CanonicalRS')
    whichstudy=3;
    nSj=55;
else
    whichstudy=1;
    nSj=42;
end
for iSj=1:nSj    
    if contains(dir,'600Hz')
        RST = spm_eeg_load([dir,'giles_symmetric_f_session',int2str(iSj)]);
    else
        RST = spm_eeg_load([dir,'sfold_giles_symmetric_f_session',int2str(iSj)]);
    end
    iRun = 2-mod(iSj,2);
    
    % Good timepoints/trials
    good_samplepoints = good_samples(RST);
    
    if rejectshortsegments
        ontimes = find(diff([0,good_samplepoints,0])==1);
        offtimes= find(diff([0,good_samplepoints,0])==-1);
        segmentlength= offtimes-ontimes;
        to_reject = segmentlength<50;
        for i=1:length(to_reject)
            if offtimes(i)>length(good_samplepoints)
                offtimes(i)=length(good_samplepoints);
            end
            if to_reject(i)
                good_samplepoints(ontimes(i):offtimes(i))=0;
            end
        end
    end
    
    %out of trial segments:
    evtypes = RST.events; evtypes = {evtypes(:).type}';
    evvals = RST.events; evvals = {evvals(:).value}';
    evvals(all(cellfun(@ischar,evvals),2),:) = {0}; % replace string with 0
    evvals=cell2mat(evvals); %convert cell to double
    evtimes = RST.events; evtimes = [evtimes(:).time]';
    if whichstudy==1
        RestingTriggerVals = [77, 88];  % 77 for the resting state BEFORE reward learning, 88 for the resting state AFTER reward learning
        RestStmInds = strcmp('UPPT001_down', evtypes) & ismember(evvals, RestingTriggerVals(iRun));    % onset of Resting
        RestStmTimes = evtimes(RestStmInds);

        EndTriggerVals = [99]; 
        EndStmInds = strcmp('UPPT001_down', evtypes) & ismember(evvals, EndTriggerVals);  % end of Resting
        EndStmTimes = evtimes(EndStmInds);

        wholetimeindex=zeros(size(RST,2),1);

        Fs = RST.fsample;
        if isempty(RestStmTimes) % for Subject X
            wholetimeindex(floor(EndStmTimes(end)*Fs)-(60*Fs*2)+1:floor(EndStmTimes(end)*Fs))=1;
        end


        if ~isempty(RestStmTimes)
            if isempty(EndStmTimes) || EndStmTimes(end)*Fs < (300*Fs) % for Subjects Y and Z
                wholetimeindex(floor(RestStmTimes(end)*Fs):floor(RestStmTimes(end)*Fs)+(300*Fs)-1)=1; 
            else
                wholetimeindex(floor(EndStmTimes(end)*Fs)-(300*Fs)+1:floor(EndStmTimes(end)*Fs))=1; % we can also use the start trigger
            end
        end
    else  % denotes study II:
        RestingTriggerVals = [87, 88];  % 87 for the resting state BEFORE reward learning, 88 for the resting state AFTER reward learning
        RestStmInds = strcmp('frontpanel trigger', evtypes) & ismember(evvals, RestingTriggerVals(iRun));    % onset of Resting
        RestStmTimes = evtimes(RestStmInds);

        EndTriggerVals = [99]; 
        EndStmInds = strcmp('frontpanel trigger', evtypes) & ismember(evvals, EndTriggerVals);  % end of Resting
        EndStmTimes = evtimes(EndStmInds);

        wholetimeindex=zeros(size(RST,2),1);
        Fs = RST.fsample;
        if ~isempty(RestStmTimes)
            wholetimeindex(floor(RestStmTimes(end)*Fs):floor(RestStmTimes(end)*Fs)+300*Fs-1)=1; 
        end

    end
    mask_boundaries{iSj}=logical(wholetimeindex);
    
    mask1{iSj} = logical(wholetimeindex) & good_samplepoints';
    goodSamples{iSj} = good_samplepoints;
    mask2{iSj} = logical(wholetimeindex);
    mask_GamtoScore = mask1{iSj}(goodSamples{iSj});
    goodSampleBoundaries = [0,diff(goodSamples{iSj})] + [diff(goodSamples{iSj}), 0];
    goodSampleBoundaries_masked = goodSampleBoundaries(goodSamples{iSj});
    mask_ScoreToGam = [mask1{iSj}(mask2{iSj})]; %note that mask2 just eliminates boundaries, not bad points
    NGam = sum(mask_ScoreToGam);
    mask_AtoB{iSj} = mask_ScoreToGam;
    mask_BtoA{iSj} = mask_GamtoScore;
end


end