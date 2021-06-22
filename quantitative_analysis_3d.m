curpath=pwd;
addpath(genpath([curpath '/nuclei_seg']));
I = imread([pwd '/img/TMA 002-G6.png']);
% crop image
I=imcrop(I,round(round([623.5 1132.5 510 414])));
figure; imshow(I);
%% **** 1 neclei segmentation using multi-resolution watershed 
[I_norm, ~, ~] = normalizeStaining(I);
I_normRed=I_norm(:,:,1);
p.scales=3:2:10;% the scale of nuclei
[nuclei, properties] = nucleiSegmentationV2(I_normRed,p);

figure; imshow(I);hold on;
for k = 1:length(nuclei)
    plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
end
hold off;
%% 
feat_values=[];
para.feature_space='Centroid-Area-MeanIntensity';

for i_space = 2:10    
    for i_area = 1:10
        for i_MeanIntensity=1:10
            para.bandWidth_space = 5*i_space ;            
            para.bandWidth_features = [5*i_area; 5*i_MeanIntensity];
            
            para.debug=0; % turn the debug mode on
            para.nuclei=nuclei; % assign pre-segmented nuclei
            para.properties=properties; % you can check out the nuclear properites
            para.shownuclei_thiness=1; % the visual effect of nuclei
            para.shownucleiconnection_thiness=3; % the visual effect of FeDeG
            
            para.num_fixed_types=0; % specify the phenotype numbers in the image
            [clustCent,data2cluster,cluster2dataCell,data_other_attribute,clust2types,typeCent]=Lconstruct_FeDeG_v3(I,para);
            %% **** extracting FeDeG features
            para.I=I;
            para.data_other_attribute=data_other_attribute;
            para.debug=0;
            para.clust2types=clust2types;
            para.typeCent=typeCent;
            
            [set_feature,set_feature_name]=L_get_FeDeG_features_v3(clustCent,cluster2dataCell,para);
            idx1 = strcmp(set_feature_name,'FeDeG-portion of intersected FeDeG');
            idx4 = strcmp(set_feature_name,'FeDeG-mean(size of non-1-2cell FeDeG)');
            idx7 = strcmp(set_feature_name,'FeDeG-mean(var of nuclear feat1to the centroid)');
            
            feat_values(i_space,i_area,i_MeanIntensity,1) = set_feature(idx1);
            feat_values(i_space,i_area,i_MeanIntensity,2) = set_feature(idx4);
            feat_values(i_space,i_area,i_MeanIntensity,3) = set_feature(idx7);
        end
    end
end

%% plot them out, the y axis is the space, x axis is the feature idx, (f1, f4, f7)

% % for i_space = 1:5    
% %     for i_area = 1:5
% %         for i_MeanIntensity=1:5
% %             para.bandWidth_space = 10*i_space ;            
% %             para.bandWidth_features = [10*i_area; 5*i_MeanIntensity];
            
figure(1);
for i_f=1:3
    curF=squeeze(feat_values(:,:,:,i_f));
    for i_space=4:10        
        subplot(7,3,(i_space-4)*3+i_f);
        curF_cur_space=squeeze(curF(i_space,:,:));
        [xx,yy]=meshgrid(5:5:50,5:5:50);
        surf(xx,yy,curF_cur_space);
%         surface(xx,yy,curF_cur_space);
        xlabel('area');
        ylabel('mean intensity');
        zlabel('feature value');
        set(gca,'FontSize',11);
%         if i_f==1
%             title('portion of intersected FeDeG');
%         end
%         if i_f==2
%             title('mean(size of non-1-2cell FeDeG)');
%         end
%          if i_f==3
%             title('mean(var of nuclear feat1 to the centroid)');
%         end
    end
end

