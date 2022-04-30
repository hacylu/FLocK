
%% %%%% this function is written by Cheng LU, from Case Western Reserve University. July, 2018
%%% any question please send email to hacylu@gmail.com
% this function uses mean shif on pre-seg nuclei to construct FeDeGs/FlocK in the image


% -input:
%         para.feature_space: feature space you would like to explore, see
%                             below for more information, or you can add any feature combination there as well
%         para.bandWidth_space: bandwidth in the spacial space, higher of the value, bigger of the FeDeG
%         para.bandWidth_features: bandwidth in the corresponding feature space
%         para.num_fixed_types - if provide a fixed number of types, 0 means no subtype number predefined. 2
%                                means group all the FeDeG into two types. can set to any number which maybe
%                                meaningful for the type of image you have at hand. also need to consider the
%                                bandwith in feature space at the same time.
%         para.debug=1; % turn the debug mode on
%         para.nuclei=nuclei; % assign pre-segmented nuclei
%         para.properties=properties; % you can check out the nuclear properites 
%         para.shownuclei_thiness=1; % the visual effect of nuclei 
%         para.shownucleiconnection_thiness=3; % the visual effect of FeDeG

% -ouput:
%          clustCent - nuclei cluster centroid that based on nuclei centroid
%          data2clusterXX - for every data point which cluster it belongs to (numPts)
%          cluster2dataCellXX  - for every cluster which points are in it (numClust)
%          XXXXXX_filtered - for the filtered data (not complete)
%          data_other_attribute - return the attributes used for constructing FeDeG other than location, 
%                                 e.g., if you used 'Centroid+Area', then the Area will be returned. this 
%                                 will be used in cal the feature in L_get_FeDeG_features function.
%          clust2types - nuclei cluster centroids to the type corresponding while para.num_fixed_types>1
%          typeCent - the type center on feature sapce
% v2 can specify different bandwidth for coordination space Hs, and feature space Hf
% v3 added one section before the first loop to double check if the
%    cluster2dataCell contains empty elements

%%% an example on how to use this function please see show_example_v2.m


function  [clustCent,data2cluster,cluster2dataCell,data_other_attribute,clust2types,typeCent]=Lconstruct_FeDeG_v3(curIM,para)
ss=para.properties;
%% use location to form FeDeG
if strcmp(para.feature_space,'Centroid')
    AlldataPts=[];
    tempC=[ss.Centroid];
    AlldataPts(1:length(ss),1)=tempC(1:2:end);
    AlldataPts(1:length(ss),2)=tempC(2:2:end);
    dataPts=AlldataPts(:,1:2);
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v2(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=[];
end
%% use location with Area to form FeDeG 
if strcmp(para.feature_space,'Centroid-Area')
    AlldataPts=[];
    tempC=[ss.Centroid];
    AlldataPts(1:length(ss),1)=tempC(1:2:end);
    AlldataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    AlldataPts(:,idxF)=[ss.Area];
    dataPts=AlldataPts(:,1:3);
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% use location with Area and eccentricity to form FeDeG with shape
if strcmp(para.feature_space,'Centroid-Area-Eccentricity')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.Area];   
    idxF=idxF+1;
    dataPts(:,idxF)=[ss.Eccentricity];
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% use location with Area, mean intensity to form FeDeG 
if strcmp(para.feature_space,'Centroid-Area-MeanIntensity')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.Area];   
    idxF=idxF+1;
    dataPts(:,idxF)=[ss.MeanIntensity]; 
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% use location with mean intensity to form FeDeG 
if strcmp(para.feature_space,'Centroid-MeanIntensity')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.MeanIntensity]; 
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% with eccentricity, note that eccentricity may not be a good measure for shape
if strcmp(para.feature_space,'Centroid-Eccentricity')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.Eccentricity];    
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% with longness
if strcmp(para.feature_space,'Centroid-Longness')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.MajorAxisLength]./[ss.MinorAxisLength];     
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% with Circularity
if strcmp(para.feature_space,'Centroid-Circularity')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.Circularity];  
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% with shape and Intensity
if strcmp(para.feature_space,'Centroid-Area-Eccentricity-MeanIntensity')
    dataPts=[];
    tempC=[ss.Centroid];
    dataPts(1:length(ss),1)=tempC(1:2:end);
    dataPts(1:length(ss),2)=tempC(2:2:end);
    idxF=3;
    dataPts(:,idxF)=[ss.Area];   idxF=idxF+1;
    dataPts(:,idxF)=[ss.Eccentricity]; idxF=idxF+1;
    dataPts(:,idxF)=[ss.MeanIntensity]; 
    
    [clustCent,data2cluster,cluster2dataCell,clust2types,typeCent] = MeanShiftCluster_v3(dataPts',para.bandWidth_space, para.bandWidth_features,para.num_fixed_types);
    FM=dataPts(:,1:2)';
    data_other_attribute=dataPts(:,3:end)';
end
%% ... you can add more feature space combination here
% %% visulize result
if para.debug
    % show the FeDeG with the different color
    coordinate_in_space=zeros(2,length(clustCent));
    for i=1:size(coordinate_in_space,2)
        coordinate_in_space(:,i)=mean(FM(:,data2cluster==i),2);
    end
    % colorized FeDeG
    cmap=colormap('hsv');
    show(curIM,8,'FeDeG without pre-defined phenotype number');
    hold on
    numObj=size(FM,2);
    for k = 1 : numObj
        if isfield(para,'nuclei')
            plot(para.nuclei{k}(:,2), para.nuclei{k}(:,1), 'Color',cmap(mod(data2cluster(k),size(cmap,1))+1,:), 'LineWidth', para.shownuclei_thiness);            
        end
        % connecting line
        plot([coordinate_in_space(1,data2cluster(k)), FM(1,k)],[coordinate_in_space(2,data2cluster(k)), FM(2,k)],'Color',cmap(mod(data2cluster(k),size(cmap,1))+1,:), 'LineWidth',para.shownucleiconnection_thiness);        
    end
    hold off

    if para.num_fixed_types>1
        % show the same type of FeDeG with the same color
        coordinate_in_space=zeros(2,length(clustCent));
        for i=1:size(coordinate_in_space,2)
            coordinate_in_space(:,i)=mean(FM(:,data2cluster==i),2);
        end
        % colorized FeDeG
        cmap=[1 0 0; 0 1 0; 0 0 1;1 1 0; 1 0 1; 0 1 1; 1 1 1];
        cmap=cmap(1:para.num_fixed_types,:);
        
        show(curIM,9,sprintf('FeDeG with pre-defined phenotype number = %d',para.num_fixed_types));
        hold on
        numObj=size(FM,2);
        for k = 1 : numObj
            if isfield(para,'nuclei')
                plot(para.nuclei{k}(:,2), para.nuclei{k}(:,1), '-', 'Color',cmap(clust2types(data2cluster(k)),:),'LineWidth', para.shownuclei_thiness);
            end
            % connecting line
            plot([coordinate_in_space(1,data2cluster(k)), FM(1,k)],[coordinate_in_space(2,data2cluster(k)), FM(2,k)],'Color',cmap(clust2types(data2cluster(k)),:), 'LineWidth',para.shownucleiconnection_thiness);
        end
        hold off
    end
end
%% add this line in v3, to double check if the cluster2dataCell contains empty element
label_empty=zeros(length(cluster2dataCell),1);
for i=1:length(cluster2dataCell)
    cur=cluster2dataCell{i};
    if isempty(cur)
        label_empty(i)=1;
    end
end
label_empty=logical(label_empty);

clustCent(:,logical(label_empty))=[];%sum(label_empty)
cluster2dataCell(logical(label_empty))=[];
if para.num_fixed_types~=0
    clust2types(logical(label_empty))=[];
end
