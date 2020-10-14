function disttoplot = plotMDS_states(hmm)
% Computes the 2D MDS distances for a given network of hmm states using
% the probability matrix distribution to define distances.

P_network = hmm.P;
% remove diagonal (ie self connections) for visualisation:
P_network(logical(eye(hmm.K)))=zeros(hmm.K,1);
P_network = P_network./repmat(sum(P_network,2),1,hmm.K);

% dist_graph = zeros(hmm.K);
% MG = mean(hmm.gamma);
% for i=1:hmm.K-1
%     for j=i+1:hmm.K
%         % marginal distribution - note this doesn't work well and not
%         % reflective of actual proximity between states (more weighted by
%         % overall state probability)
% %        dist_graph(i,j) = sum(P_network(i,setdiff([1:hmm.K],[i,j])))*MG(i) + ...
% %            sum(P_network(j,setdiff([1:hmm.K],[i,j])))*MG(j);
%         % conditional distribution:
%         dist_graph(i,j) = sum(P_network(i,setdiff([1:hmm.K],[i,j])))*1 + ...
%             sum(P_network(j,setdiff([1:hmm.K],[i,j])))*1;
%         %dist_graph(i,j) = 1 - sum(P_network(i,[i,j]));
%     end
% end
% dist_graph = dist_graph + dist_graph';


% alternative as sanity check:
P_network = 1-P_network;
P_network(logical(eye(hmm.K)))=zeros(hmm.K,1);
P_network = P_network + P_network';
dist_graph = P_network;

[disttoplot,eps] = cmdscale(dist_graph,2);



% alternative method based on windowed correlation (achieves similar results):

% % Compute distances based on windowed correlation matrix:
% windowlength = 1*sample_rate;
% T_tot = length(hmm.gamma);
% gam_windowed = reshape(hmm.gamma(1:(T_tot-mod(T_tot,windowlength)),:),[windowlength,floor(T_tot/windowlength),hmm.K]);
% gam_windowed = squeeze(sum(gam_windowed,1));
% dist_graph = dist(gam_windowed);
% 
% [disttoplot,eps] = cmdscale(dist_graph,2);

end