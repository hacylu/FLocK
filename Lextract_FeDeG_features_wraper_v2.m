% by cheng lu, on 30rd Jan. 2017
% v2 are using for FeDeG v2
function [feats,feats_description] = Lextract_FeDeG_features_wraper_v2(bounds,I,properties,nuclei,para)
% compute the mean and standard deviation of CGTs across all bounds
% c = co-occurrence matrix
warning('off','all');

% para.bandWidth_space=200;% higher the bigger of FeDeG
% %             para.bandWidth_features=[20];% bandwidth of in the feature space
%
% para.bandWidth_features=[60;10];% bandwidth of in the feature space
% % para.bandWidth_features=[20;5];% bandwidth of in the feature space
% para.num_fixed_types=3;

feats=[];
feats_description=[];
set_alpha=[para.bandwidth_min:para.bandwidth_res:para.bandwidth_max];% bandWidth for space
% set_alpha_features=[para.bandwidth_min_features:para.bandwidth_res_features:para.bandwidth_max_features];% bandWidth for feature
for ai=1:size(para.bandwidth_min_features,1)
    set_alpha_features(:,ai)=[para.bandwidth_min_features(ai):para.bandwidth_res_features(ai):para.bandwidth_max_features(ai)]';% bandWidth for feature
end
% set_FeDeG_diff_attribute={'Centroid-Area','Centroid-Area-Eccentricity','Centroid-Area-MeanIntensity','Centroid-Area-Eccentricity-MeanIntensity'};
% set_FeDeG_diff_attribute={'Centroid-Area','Centroid-MeanIntensity','Centroid-Longness','Centroid-Circularity','Centroid-Eccentricity',...
%     'Centroid-Area-Eccentricity','Centroid-Area-MeanIntensity','Centroid-Area-Eccentricity-MeanIntensity'};
% set_FeDeG_diff_attribute=para.feature_space;
curpara.num_fixed_types=para.num_fixed_types;

for f=1:length(set_alpha)
    %for a=1:length(set_FeDeG_diff_attribute)
    curpara.bandWidth_space=set_alpha(f); % higher the bigger of FeDeG
%     curpara.bandWidth_features=para.bandWidth_features;
    curpara.bandWidth_features=set_alpha_features(f,:)';    
    %% computing FeDeG using meanshift clustering
    % clear dataPts;
    %  para.bwGT=bwGT;
    %  para.show=flag_showMSinterResult;
    %         para.Feature='Centroid';
    curpara.feature_space=para.feature_space;
    
    fprintf('use MS in feature space (%s) with spacial bandwidth %d to build FeDeG...\n',curpara.feature_space,curpara.bandWidth_space);
    
    %         para.Feature='Centroid-Area-Eccent-Solid';
    curpara.debug=0;
    % don't put I if don't want to debug
    if curpara.debug
        curpara.nuclei=1;
    else
        curpara.nuclei=nuclei;
    end
    
    curpara.properties=properties;
    %         para.
    [clustCent,data2cluster,cluster2dataCell,data_other_attribute,clust2types,typeCent]=Lconstruct_FeDeG_v2(I,curpara);
    %% derive features from FeDeG
    curpara.debug=0;
    % don't put I if don't want to debug
    if curpara.debug
        curpara.I=1;
    else
        curpara.I=I;
    end
    curpara.data_other_attribute=data_other_attribute;
    
    
    curpara.debug=0;
    curpara.clust2types=clust2types;
    curpara.typeCent=typeCent;
    
%     [feat,feat_description]=L_get_FeDeG_features_v2(clustCent,data2cluster,cluster2dataCell,curpara);
    [feat,feat_description]=L_get_FeDeG_features_v2(clustCent,cluster2dataCell,curpara);
    
    %%% get feature description
    temp= feat_description;
    for i=1:length(temp)
        cur=temp{i};
        str=sprintf('(%s)-bandwidth=%d',curpara.feature_space,curpara.bandWidth_space);
        cur(end+1:end+length(str))=str;
        temp{i}=cur;
    end
    feats=cat(2,feats,feat); feats_description=cat(2,feats_description,temp);
    %     fprintf('%d features at %f\n', length(feat),curpara.alpha);
end
end
% end


