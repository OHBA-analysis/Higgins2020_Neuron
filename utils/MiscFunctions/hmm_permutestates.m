function hmm_new = hmm_permutestates(hmm,new_state_ordering);
hmm_new = hmm;
for i=1:hmm.K
    hmm_new.state(i) = hmm.state(new_state_ordering(i));
    hmm_new.Dir_alpha(i) = hmm.Dir_alpha(new_state_ordering(i));
    hmm_new.Pi(i) = hmm.Pi(new_state_ordering(i));
    hmm_new.gamma(:,i) = hmm.gamma(:,new_state_ordering(i));
    if isfield(hmm,'statepath')
        hmm_new.statepath(hmm.statepath==new_state_ordering(i)) =i;
    end
    if isfield(hmm,'statemap_parcel_vectors_persubj')
        hmm_new.statemap_parcel_vectors_persubj(:,:,i) = hmm.statemap_parcel_vectors_persubj(:,:,new_state_ordering(i));
        hmm_new.statemap_parcel_vectors(:,i) = hmm.statemap_parcel_vectors(:,new_state_ordering(i));
    end
    for j=1:hmm.K
        hmm_new.P(i,j) = hmm.P(new_state_ordering(i),new_state_ordering(j));
        hmm_new.Dir2d_alpha(i,j) = hmm.Dir2d_alpha(new_state_ordering(i),new_state_ordering(j));
    end
end

end