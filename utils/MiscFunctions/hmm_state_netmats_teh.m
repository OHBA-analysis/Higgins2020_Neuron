function [ state_netmats ] = hmm_state_netmats( hmm, S )

% [ state_netmats ] = hmm_state_netmats( hmm, S )

try S.global_only=S.global_only; catch S.global_only=false; end
try parcellated_filenames=S.parcellated_filenames; catch error('Need to specify parcellated_filenames'); end

logtrans=0;

normalisation = S.normalisation;
embed=S.embed;
assignment=S.assignment;

D_fnames=S.parcellated_filenames;    

nsubs=length(D_fnames);

for subnum = 1:nsubs

    disp(['Computing for subj num ' num2str(subnum)]);
       
    try
        Dp = spm_eeg_load(prefix(D_fnames{subnum},'p'));
    catch
        Dp = spm_eeg_load((D_fnames{subnum}));
    end
        
    %state_netmats{subnum}.parcelAssignments=parcelAssignments;
    %state_netmats{subnum}.parcellation=Dp.parcellation;
    
    embed.tres=1/Dp.fsample;
        
    %     if subnum==3,
    %         Dp(2,:,:)=-Dp(2,:,:);
    %         Dp(37,:,:)=-Dp(37,:,:);
    %         Dp.save;
    %     end;

    % returns data as num_nodes x num_embeddings x ntpts
    %datap = prepare_data(Dp,normalisation,logtrans,embed,roinets_protocol);        
    datap = osl_teh_prepare_data(Dp,normalisation,logtrans,[],embed);
    
    if ~S.global_only && hmm.K>1
        
        hmm_sub = hmm; 
        hmm_sub.statepath = hmm.statepath(hmm.subj_inds==subnum); 
        hmm_sub.gamma = hmm.gamma(hmm.subj_inds==subnum,:);

%        hmm_sub = rmfield(hmm_sub,'MixingMatrix');    

        for k = 1:hmm_sub.K

             disp(['Computing for state ' num2str(k)]);

             switch assignment
                 case 'hard' 

                    inds = logical(hmm_sub.statepath == k);

                    if size(datap,2) ~= length(inds)
                        error('not implemented');
                    end

                    %normconstant = sqrt(size(inds,1) / sum(double(inds)) );

                    datap_in=datap(:,inds);% * normconstant;
                                  
                 case 'soft'
                    x = hmm_sub.gamma(:,k);

                    if size(datap,2) ~= length(x)
                        error('not implemented');
                    end

                    %normconstant = sqrt(size(x,1) / sum(x) ); 
                    %x = x * normconstant; 
                    x2 = permute(repmat(x,[1, size(datap,1)]),[2 1]);
                    datap_in=datap.*x2;
             end

            try
                [state_netmats{subnum}.state{k}]=feval(S.netmat_method,datap_in,S.netmat_method_options);
                state_netmats{subnum}.state{k}.ntpts=size(datap_in,2);
            catch
                warning(['State ' num2str(k) ' is not visited in subject ' num2str(subnum)]);
                state_netmats{subnum}.state{k}=[];
                state_netmats{subnum}.state{k}.ntpts=0;
            end
            
        end
    end
    
    [state_netmats{subnum}.global]=feval(S.netmat_method,datap,S.netmat_method_options);
    state_netmats{subnum}.global.ntpts=size(datap,3);

    %sum(state_netmats{subnum}.global.spectramt.psd(:,9,9),1)
    %sum(state_netmats{subnum}.state{k}.spectramt.psd(:,9,9),1)
    %std(squeeze(datap_in(9,1,:))),std(squeeze(datap(9,1,:)))
    %sum(state_netmats{subnum}.global.netmat(9,9,:))
    %sum(state_netmats{subnum}.state{k}.netmat(9,9,:))
    
    
    state_netmats{subnum}.netmat_method=S.netmat_method;
end                        
