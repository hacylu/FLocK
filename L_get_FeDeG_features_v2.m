%% %%%% this function is written by Cheng LU, from Case Western Reserve University. July, 2018
%%% any question please send email to hacylu@gmail.com
%% given a feature driven cell graphs (FeDeG) or FlocK, we calculate the feature here
%%% this is the second version of FeDeG feature extration
% written by Cheng Lu at July 2018. Cleveland.
% v2 added  1) more features;
%           2) provide the flag for breaking the object into k-groups by specify the k or threshold in featuer space
%           3) use different bandwidthes for space and features
% note, we use FeDeG(cell graph) to represent FeDeG/FLocK
%% features are described below
%% A. features that do not consider cluster phenotypes
%%% 1. Overlapped features: overlaped area portion between nearby FeDeG (need to turn FeDeG into binary component, then cal the overlapped area)
%%% 2. Variation within FeDeG: variation of the size of FeDeG/number of nuclei within a FeDeG (include the single cell as a FeDeG) across the whole image
%%% 3. Variation wrt to centroid: variation of the distance/size/color/shape of nucleus to the centroid in a FeDeG, so
%%%    that each FeDeG has a variation value, then use statistics across the whole image
%%% 4. Spatial arrangement features: use the centroids of FeDeG or intersected region of FeDeG as node, build global graph, extract the standard global graph feature
%%% 5. Density features: that capture the density of the FeDeG
%% B. features that consider cluster phenotypes
%%% 1. intersection properties/measurement of FeDeGs from different phenoytype
%%% 2. NNearest features: enrichement of the neibouhood FeDeGs type
%%% 3. Spatial arrangement features: use the centroids of different types of FeDeG as node, build global graph, extract the standard global graph feature

%%% input:
%          clustCent - nuclei cluster centroid that based on nuclei
%          centroid, we named a FeDeG/FlocK a cluster in the script
%          data2clusterXX - for every data point(nuclei location) which cluster it belongs to (numPts)
%          cluster2dataCellXX - for every cluster which points are in it (numClust)
%% note we assume that the first two dimensions of clustCent are location dimension
%% note that cluster2dataCell may contain empty elements
function [set_feature,set_feature_name]=L_get_FeDeG_features_v2(clustCent,cluster2dataCell,para)
set_feature_name=[]; idx_feat=1;
set_feature=[];
%% turn each FeDeG into polygon by using the centroid of nuclei as polygon vertices
% 1) collect the polygon(FeDeG) in the image. 
% 2) also record the FeDeG info. 
% 3) cal the variation of the distance/size/color/shape of nucleus to the centroid in a FeDeG

set_polygon=[]; % a cell structure, each element contain a n-by-2 vector, where n is the number of vertices, 2 is the location x and y of vertices
%%% note that for FeDeG with only one or two nuclei, we treat them seperately
labels_single_and_two_cell_cluster=zeros(1,size(clustCent,2));

% for recording the FeDeG info.
set_FeDeG_size=ones(1,size(clustCent,2));
set_nuclei_num_in_FeDeG=ones(1,size(clustCent,2));
% only consider the FeDeG has more than 3 nuclei, the others will put 0, and will be ignored when compute the actual feature
set_variation_distance_2_centroid_of_FeDeG=zeros(1,size(clustCent,2));
% this is accounting for the variation of other attributes, the first dimension is the number of FeDeG, seconde dimension is the feature(nuclear attribute)
set_variation_other_attribute_2_centroid_of_FeDeG=zeros(size(clustCent,2),size(para.data_other_attribute,1));

for i=1:size(clustCent,2)
    cur_cluster2data=cluster2dataCell{i};
    cur_centroid=[];
    cur_attribute=[];
    for j=1:length(cur_cluster2data)
        cur_centroid(j,:)=para.properties(cur_cluster2data(j)).Centroid;
        if ~isempty(para.data_other_attribute)
            cur_attribute(j,:)=para.data_other_attribute(:,cur_cluster2data(j))';
        end
    end
    % use the convex hull algorithm to determine the outter bound of a FeDeG
    if size(cur_centroid,1)>3
        K=convhull (cur_centroid(:,1),cur_centroid(:,2));
        set_polygon{i}=[cur_centroid(K,1),cur_centroid(K,2)];
        set_FeDeG_size(i)=polyarea(cur_centroid(K,1),cur_centroid(K,2)); % area
        set_nuclei_num_in_FeDeG(i)=size(cur_centroid,1);
        
        cur_mean_centroid=mean([cur_centroid(:,1),cur_centroid(:,2)]);
        set_variation_distance_2_centroid_of_FeDeG(i)= std(sqrt((cur_centroid(:,1)-cur_mean_centroid(1)).^2+(cur_centroid(:,2)-cur_mean_centroid(2)).^2));
        
        if ~isempty(para.data_other_attribute)
            cur_mean_attribute=mean(cur_attribute);
            set_variation_other_attribute_2_centroid_of_FeDeG(i,:)=std(abs(cur_attribute-repmat(cur_mean_attribute,size(cur_attribute,1),1)));
        end
    else %for FeDeGs with only one or two nuclei
        set_polygon{i}=cur_centroid;
        labels_single_and_two_cell_cluster(i)=1;
        % use a proximate nuclei size as the FeDeG size
        if size(cur_centroid,1)==1
            tmp=para.nuclei{cur_cluster2data};
            set_polygon{i}=[tmp(:,2) tmp(:,1)];% use the nuclei boundary as polygon for FeDeG with singel cell
            set_FeDeG_size(i)=para.properties(cur_cluster2data).Area;
        end
        if size(cur_centroid,1)==2
            tmp1=para.nuclei{cur_cluster2data(1)};
            tmp2=para.nuclei{cur_cluster2data(2)};
            tmp=[tmp1;tmp2];
            cur_centroid=[tmp(:,2) tmp(:,1)];
            K=convhull (cur_centroid(:,1),cur_centroid(:,2));
            set_polygon{i}=[cur_centroid(K,1),cur_centroid(K,2)];
            set_FeDeG_size(i)=polyarea(cur_centroid(K,1),cur_centroid(K,2));
        end
    end
end
%%% prepare the cell location information and store them in variable FM
AlldataPts=[];
tempC=[para.properties.Centroid];
AlldataPts(1:length(tempC)/2,1)=tempC(1:2:end);
AlldataPts(1:length(tempC)/2,2)=tempC(2:2:end);
dataPts=AlldataPts(:,1:2);
FM=dataPts(:,1:2)'; % nuclear centroid

%%  calculate the intersection and record the intersection information between each pair of FeDeG
set_intersection_FeDeG_idx_pair=[]; idx_IFeDeG=1;
% check the intersection information, the relative portion depends on which FeDeG to be in the denominator, have two possible values, for the smaller one we put in XX_small
set_intersection_FeDeG_relative_portion_small=[];
set_intersection_FeDeG_relative_portion_large=[];
set_intersection_FeDeG_absolute_value=[];
set_intersection_FeDeG_centroid=[];

d=squareform(pdist(clustCent(1:2,:)'));
d(d==0)=Inf;

for i=1:length(set_polygon)-1
    cur_polygon_i=set_polygon{i};
    [~,ind]=sort(d(i,:),'ascend');
    set_polygon_NN30_idx=ind(1:30);%check it's 30 neighbours
    set_polygon_NN30=set_polygon(set_polygon_NN30_idx);
    for j=1:length(set_polygon_NN30)
        cur_polygon_j=set_polygon_NN30{j};
        % calculate the intersection
        if ~isempty(cur_polygon_i)&&~isempty(cur_polygon_j)
            [cur_polygon_intersect_x,cur_polygon_intersect_y]= polybool('intersection',cur_polygon_i(:,1),cur_polygon_i(:,2),cur_polygon_j(:,1),cur_polygon_j(:,2));
            cur_polygon_intersect=[cur_polygon_intersect_x,cur_polygon_intersect_y];
            cur_polygon_intersect(sum(isnan(cur_polygon_intersect),2)>0,:)=[];
         end
        
        %%% if we found FeDeGs that have intersection, record the information
        if ~isempty(cur_polygon_intersect)
            set_intersection_FeDeG_idx_pair(idx_IFeDeG,:)=[i set_polygon_NN30_idx(j)];
            % record the intersection information, the relative portion depends on which FeDeG to be in the denominator, have two possible values, smaller one we put in XX_small
            area_intersect=polyarea(cur_polygon_intersect(:,1),cur_polygon_intersect(:,2));
            set_intersection_FeDeG_relative_portion_small(idx_IFeDeG)=min(area_intersect/polyarea(cur_polygon_i(:,1),cur_polygon_i(:,2)),area_intersect/polyarea(cur_polygon_j(:,1),cur_polygon_j(:,2)));
            set_intersection_FeDeG_relative_portion_large(idx_IFeDeG)=max(area_intersect/polyarea(cur_polygon_i(:,1),cur_polygon_i(:,2)),area_intersect/polyarea(cur_polygon_j(:,1),cur_polygon_j(:,2)));
            set_intersection_FeDeG_absolute_value(idx_IFeDeG)=area_intersect;
            
            set_intersection_FeDeG_centroid(idx_IFeDeG,:)=mean(cur_polygon_intersect);
            idx_IFeDeG=idx_IFeDeG+1;
        end 
    end
end
%% A. features that do not consider phenotypes
%% 1. intersection properties/measurement of FeDeGs or overlaped area portion between nearby FeDeG
%     the feature name should be informative so that we know what feature it is
set_feature(idx_feat)=(idx_IFeDeG-1)/size(clustCent,2);
set_feature_name{idx_feat}='FeDeG-portion of intersected FeDeG';% note that this feature can be greater than 1 since one FeDeG can have intersection with more than one FeDeG
idx_feat=idx_feat+1;

set_feature(idx_feat)=(idx_IFeDeG-1);
set_feature_name{idx_feat}='FeDeG-abs. number of intersected FeDeG';
idx_feat=idx_feat+1;

T=0.1:0.1:0.9;% threshold for the 'highly intersected'
for i=1:length(T)
    set_feature(idx_feat)=sum((set_intersection_FeDeG_relative_portion_large>T(i)))/size(clustCent,2);
    set_feature_name{idx_feat}=sprintf('FeDeG-portion of highly intersected FeDeG:T=%.1f',T(i));
    idx_feat=idx_feat+1;
end

for i=1:length(T)
    set_feature(idx_feat)=sum((set_intersection_FeDeG_relative_portion_large>T(i)));
    set_feature_name{idx_feat}=sprintf('FeDeG-number of highly intersected FeDeG:T=%.1f',T(i));
    idx_feat=idx_feat+1;
end

%%% statistics on the intersection area
statisti_name = [{'mean'} {'median'} {'std'} {'range'} {'min'} {'max'} {'kurtosis'} {'skewness'} ];
%the format is like this: FeDeG-statistcname(intersection area)


for mi = 1:numel(statisti_name)
    if isempty(set_intersection_FeDeG_relative_portion_small)
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s(set_intersection_FeDeG_relative_portion_small)',statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'intersection area portion:small' ')'  ];
    idx_feat=idx_feat+1;
end

for mi = 1:numel(statisti_name)
    if isempty(set_intersection_FeDeG_relative_portion_large)
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s(set_intersection_FeDeG_relative_portion_large)',statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'intersection area portion:large' ')'  ];
    idx_feat=idx_feat+1;
end

for mi = 1:numel(statisti_name)
    if isempty(set_intersection_FeDeG_absolute_value)
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s(set_intersection_FeDeG_absolute_value)',statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'intersection area abs val' ')'  ];
    idx_feat=idx_feat+1;
end

for mi = 1:numel(statisti_name)
    if isempty(set_intersection_FeDeG_relative_portion_small)
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s((set_intersection_FeDeG_relative_portion_large+set_intersection_FeDeG_relative_portion_small)/2)',...
            statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'intersection area portion:mean' ')'  ];
    idx_feat=idx_feat+1;
end
%% 2. statistics of the size of FeDeG/number of nuclei within a FeDeG (include the single cell as a FeDeG) across the whole image
for mi = 1:numel(statisti_name)
    set_feature(idx_feat)=eval(sprintf('%s(set_FeDeG_size)',statisti_name{mi}));
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'size of all FeDeG' ')'  ];
    idx_feat=idx_feat+1;
end

% statistics of FeDeG that exclude the FeDeG has one or two cells
for mi = 1:numel(statisti_name)
    if isempty(find(~labels_single_and_two_cell_cluster, 1))
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s(set_FeDeG_size(~labels_single_and_two_cell_cluster))',...
            statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'size of non-1-2cell FeDeG' ')'  ];
    idx_feat=idx_feat+1;
end

set_feature(idx_feat)=sum(~labels_single_and_two_cell_cluster)/length(labels_single_and_two_cell_cluster);
set_feature_name{idx_feat}=sprintf('FeDeG-non-1-2cell FeDeG ratio');
idx_feat=idx_feat+1;

for mi = 1:numel(statisti_name)
    set_feature(idx_feat)=eval(sprintf('%s(set_nuclei_num_in_FeDeG)',...
        statisti_name{mi}));
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'nuclei no. in FeDeG' ')'  ];
    idx_feat=idx_feat+1;
end
%% 3. variation of the distance/size/color/shape of nucleus to the centroid in a FeDeG (FeDeG has >=3 cells), so
%%% that each FeDeG has a variation value, then use statistics across the whole image, this group of features maybe trivial, since the FeDeG are form by the predefined bandwidth.
if ~isempty(para.data_other_attribute)
    for t=1:size(set_variation_other_attribute_2_centroid_of_FeDeG,2)
        cur_att=set_variation_other_attribute_2_centroid_of_FeDeG(:,t);
        for mi = 1:numel(statisti_name)
            if isempty(eval(sprintf('%s(cur_att(cur_att>0))',...
                    statisti_name{mi})))
                set_feature(idx_feat)=0;
            else
                set_feature(idx_feat)=eval(sprintf('%s(cur_att(cur_att>0))',...
                    statisti_name{mi}));
            end
            set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'var of nuclear feat' num2str(t) 'to the centroid'  ')'  ];
            idx_feat=idx_feat+1;
        end
    end
end

for mi = 1:numel(statisti_name)
    if isempty(eval(sprintf('%s(set_variation_distance_2_centroid_of_FeDeG(set_variation_distance_2_centroid_of_FeDeG>0))',...
            statisti_name{mi})))
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s(set_variation_distance_2_centroid_of_FeDeG(set_variation_distance_2_centroid_of_FeDeG>0))',...
            statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'var of cell to FeDeG centroid' ')'  ];
    idx_feat=idx_feat+1;
end
%% 4. Spatial arrangement features: use the centroids of FeDeG or intersected region of FeDeG as node, build global graph, extract the standard global graph feature
% use the centroids of FeDeG as node, build global graph, extract the standard global graph feature
if length(set_intersection_FeDeG_centroid)>8
    [vfeature,GraphFeatureDescription] = get_graph_features(clustCent(1,:)',clustCent(2,:)');
else
    % can't compute the features, and put all zeros
    vfeature=zeros(1,51);
    load('GraphFeatureDescription.mat')
end
set_feature=cat(2,set_feature,vfeature);

for i=1:length(GraphFeatureDescription)
    curname=GraphFeatureDescription{i};
    curnew=['FeDeG-FeDeG centroid-' curname];
    set_feature_name{idx_feat}=curnew;
    idx_feat=idx_feat+1;
end

%%% use the centroids of intersected region of FeDeG as node, build global graph, extract the standard global graph feature
if length(set_intersection_FeDeG_centroid)>8
    [vfeature,GraphFeatureDescription] = get_graph_features(set_intersection_FeDeG_centroid(:,1),set_intersection_FeDeG_centroid(:,2));
else
    % can't compute the features, and put all zeros
    vfeature=zeros(1,51);
    load('GraphFeatureDescription.mat')
end

set_feature=cat(2,set_feature,vfeature);

for i=1:length(GraphFeatureDescription)
    curname=GraphFeatureDescription{i};
    curnew=['FeDeG-intersected FeDeG centroid-' curname];
    set_feature_name{idx_feat}=curnew;
    idx_feat=idx_feat+1;
end
%% 5. Density features: that capture the density of the FeDeG,
%%% cal the density for FeDeGs have >3 cells and all FeDeGs
%%% the density will be the number of object(nuclei) / number of pixels the FeDeG taken

% statistics of density of FeDeG that exclude the FeDeG has one or two cells
for mi = 1:numel(statisti_name)
    if isempty(find(~labels_single_and_two_cell_cluster, 1))
        set_feature(idx_feat)=0;
    else
        set_feature(idx_feat)=eval(sprintf('%s(set_nuclei_num_in_FeDeG(~labels_single_and_two_cell_cluster)./set_FeDeG_size(~labels_single_and_two_cell_cluster))',...
            statisti_name{mi}));
    end
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'density of non-1-2cell FeDeGs' ')'  ];
    idx_feat=idx_feat+1;
end

% statistics of density of all FeDeGs
for mi = 1:numel(statisti_name)
    set_feature(idx_feat)=eval(sprintf('%s(set_nuclei_num_in_FeDeG./set_FeDeG_size)',...
        statisti_name{mi}));
    set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'density of all FeDeGs' ')'  ];
    idx_feat=idx_feat+1;
end
%% B. features that consider cluster types
%%% 1. intersection properties/measurement of FeDeGs from different FeDeG Type
%%% 2. NNearest features/Distance to other clusters features:
%%% 3. Encompass features:

if para.num_fixed_types>1
    %% 1. intersection properties/measurement of FeDeGs from different FeDeG Type
    %     set_intersection_FeDeG_idx_pair
    set_intersection_FeDeG_idx_pair_inFeDeGType=para.clust2types(set_intersection_FeDeG_idx_pair);
    for iType=1:para.num_fixed_types
        if ~isempty(find(sum(set_intersection_FeDeG_idx_pair_inFeDeGType==iType,2)==2, 1))
            % intra-type
            set_feature(idx_feat)=sum(sum(set_intersection_FeDeG_idx_pair_inFeDeGType==iType,2)==2)/sum(para.clust2types==iType);
        else
            set_feature(idx_feat)=0;
        end
        set_feature_name{idx_feat}=['FeDeG-portion of intra-type FeDeG intersction-Type' num2str(iType)];
        idx_feat=idx_feat+1;
        
        if ~isempty(find(sum(set_intersection_FeDeG_idx_pair_inFeDeGType~=iType,2)==2, 1))
            % inter-type
            set_feature(idx_feat)=sum(sum(set_intersection_FeDeG_idx_pair_inFeDeGType~=iType,2)==2)/sum(para.clust2types==iType);
        else
            set_feature(idx_feat)=0;
        end
        set_feature_name{idx_feat}=['FeDeG-portion of inter-type FeDeG intersction-Type' num2str(iType)];
        idx_feat=idx_feat+1;
    end
    
    for iType=1:para.num_fixed_types
        if ~isempty(find(sum(set_intersection_FeDeG_idx_pair_inFeDeGType==iType,2)==2, 1))
            % intra-type
            set_feature(idx_feat)=sum(sum(set_intersection_FeDeG_idx_pair_inFeDeGType==iType,2)==2);
        else
            set_feature(idx_feat)=0;
        end
        set_feature_name{idx_feat}=['FeDeG-number of intra-type FeDeG intersction-Type' num2str(iType)];
        idx_feat=idx_feat+1;
        
        if ~isempty(find(sum(set_intersection_FeDeG_idx_pair_inFeDeGType~=iType,2)==2, 1))
            % inter-type
            set_feature(idx_feat)=sum(sum(set_intersection_FeDeG_idx_pair_inFeDeGType~=iType,2)==2);
        else
            set_feature(idx_feat)=0;
        end
        set_feature_name{idx_feat}=['FeDeG-number of inter-type FeDeG intersction-Type' num2str(iType)];
        idx_feat=idx_feat+1;
    end
    
    %% 2. NNearest features: enrichement of the neibouhood FeDeGs type
    nset=[5;10;15];
    d=squareform(pdist(clustCent(1:2,:)'));
    d(d==0)=Inf;
    numClust=size(clustCent,2);
    for n=1:length(nset)
        curClust=[];
        if length(d)>nset(n)
            for i=1:numClust
                [~,ind]=sort(d(i,:),'ascend');
                curClust(i)=sum(para.clust2types(ind(1:nset(n)))~=para.clust2types(i))/nset(n);
            end
            for mi = 1:numel(statisti_name)
                set_feature(idx_feat)=eval(sprintf('%s(curClust)',...
                    statisti_name{mi}));
                set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'portion of other FeDeG types in' num2str(nset(n)) 'nearest neighbors' ')'  ];
                idx_feat=idx_feat+1;
            end
        else
            for mi = 1:numel(statisti_name)
                set_feature(idx_feat)=NaN;
                set_feature_name{idx_feat} = ['FeDeG-' statisti_name{mi} '(' 'portion of other FeDeG types in' num2str(nset(n)) 'nearest neighbors' ')'  ];
                idx_feat=idx_feat+1;
            end
        end
    end
    %% 3. Spatial arrangement features: use the centroids of different types of FeDeG as node, build global graph, extract the standard global graph feature
    for i=1:para.num_fixed_types
        if sum(para.clust2types==i)>8
            [vfeature,GraphFeatureDescription] = get_graph_features(clustCent(1,para.clust2types==i)',...
                clustCent(2,para.clust2types==i)');
        else
            % can't compute the features, and put all zeros
            vfeature=zeros(1,51);
        end
        set_feature=cat(2,set_feature,vfeature);
        
        for j=1:length(GraphFeatureDescription)
            curname=GraphFeatureDescription{j};
            curnew=['FeDeG-Type' num2str(i) 'FeDeG centroid-' curname];
            set_feature_name{idx_feat}=curnew;
            idx_feat=idx_feat+1;
        end
    end
end