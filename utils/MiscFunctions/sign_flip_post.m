function [energy, inds]=sign_flip_post(netmat,flips,group_netmat,sigma,num_embeddings,pinv_group_netmat,method,num_lags,inds)

energy=1e32;

flips=sparse(diag(kron(flips',ones(1,num_embeddings))));

switch method
    case 'matrix_distance_metric'
        [metric, inds]=matrix_distance_metric(flips*netmat*flips,group_netmat,num_lags,'none',inds);
        energy = -metric;
    otherwise
        if sigma==0
            %Sigmab=flips*group_netmat*flips;
            
            pinv_Sigmab=flips*pinv_group_netmat*flips;
            energy=trace(netmat*pinv_Sigmab);%+logdet(Sigmab);

        else
            Sigmab=eye(size(netmat))*sigma^2 + flips*group_netmat*flips;
            pinv_Sigmab=pinv(Sigmab);
            energy=trace(netmat*pinv_Sigmab)+logdet(Sigmab);
        end;
        
end;

end

