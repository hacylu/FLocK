%%----------------------%%
% Input
%     R: image channel (current implementation using Red Channel)
%    C_mask3: the foreground is nuclei which clump together or have
%    concave shape
% Output:
%    bs4: the binary image with nuclei centers as foreground
% Written by Hongming Xu
% ECE U of A
% feel free to use it
%%---------------------%%


function [ns]=XNucleiCenD_Clustering(R,C_mask3,largeSigma,smallSigma)


%% Method (iii) gGoG

R_hat=double(255-R);
thetaStep=pi/9;    % thetaStep
sigmaStep=-1;              % Sigma step
kerSize=largeSigma*4;    % Kernel size
alpha = 1;


%% fast old version
% [aggregated_response] = Xaggregate_gLoG_filter_response2(R_hat, largeSigma, smallSigma, sigmaStep,thetaStep, kerSize, alpha);
% bfs=zeros(size(R));
% for i=1:size(aggregated_response,3)
%     aggregated_response1=aggregated_response(:,:,i);
%     aggregated_response_norm = scaleNormalization(aggregated_response1, 0, 255);                                                                                       
%     bt=imregionalmax(aggregated_response_norm);
%     bt(~C_mask3)=0;
%     
%     if i==1
%         bfs=bt;
%     else
%         bt_dilate=imdilate(bfs,strel('disk',round(smallSigma*1.5)));
%         bt_and=bt_dilate&bt;
%         bt_temp=bt-bt_and;
%                
%         bfs=bfs|bt_temp;
%        end
%     
% end

%% modified version
[aggregated_response] = Xaggregate_gLoG_filter_response2_V2(R_hat, largeSigma, smallSigma, sigmaStep,thetaStep, kerSize, alpha); %% summation with each direction
X=[;];
Y=[];
for i=1:size(aggregated_response,3)
    aggregated_response1=aggregated_response(:,:,i);
    aggregated_response_norm = scaleNormalization(aggregated_response1, 0, 255);                                                                                       
    bt=imregionalmax(aggregated_response_norm);
    
%    bcur=bcur|bt;   %% for generating figures
    
    bt(~C_mask3)=0;
    [r,c]=find(bt);
    X=[X,[r';c']];
    ind=sub2ind(size(R),r,c);
    Y=[Y,aggregated_response_norm(ind)'];    
end

% show(aggregated_response_norm)
% [r,c]=find(bcur);
% X=[r';c'];
bandwidth=6;
[clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X,bandwidth);

ns=[];
for i=1:length(clustMembsCell)
    temp=clustMembsCell{i,1};
    [~,ind]=max(Y(temp));
    index=temp(ind);
    ns=[ns,X(:,index)];
end

end