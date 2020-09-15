addpath(genpath('./veta_watershed'));
addpath(genpath('./GeneralLoG'));
addpath(genpath('./staining_normalization'));
dbstop if error
% IMData4SymmetricVoting_1={'8913.tif'};
% IMData4SymmetricVoting_1={'29_38.png'};
% IMData4SymmetricVoting_1={'36_44.png'};
% I=imread('C:\Users\cxl884\Google Drive\SachethImages\images\392245-1.tif');
% IMData4SymmetricVoting_1={'C:\Users\cxl884\Google Drive\SachethImages\images\392245-1.tif'};
% IMData4SymmetricVoting_1={'new_1.png'};
% strFolderName='C:\Users\cxl884\Dropbox\DL_model_nuclei_lung_examples_priya\poor_tumor_detection\';
% strFolderName='F:\Users\cxl884\Dropbox\DL_model_nuclei_lung_examples_priya\poor_tumor_detection\';
% % corr
% strFolderName='F:\Users\cxl884\Dropbox\DL_model_nuclei_lung_examples_priya\good_lymphocyte_detection_poor_tumor_detection\';

% strFolderName='F:\Nutstore\Nutstore\4PublicProgram\nuclei_seg\';
strFolderName=[pwd '\'];
dirIM=dir([strFolderName '1_Cluster#3.tif']);
for i=1:length(dirIM)
    curIMName=dirIM(i).name;
    curIM=imread([strFolderName curIMName]);
    curIMsize=size(curIM);
    [curIM_norm] = normalizeStaining(curIM);
    curIM_normRed=curIM_norm(:,:,1);%show(curIM(:,:,3))
    curIM_normRed=curIM(:,:,3);%show(curIM(:,:,3))

    %% using gLoG-Watershed 
    close all
    p.scales=[6:8];
    disp('begin nuclei segmentation using watershed');
    [nuclei, properties] = nucleiSegmentationV2_gLoG(curIM_normRed,p);
    
    figure;imshow(curIM);hold on;
    for k = 1:length(nuclei)
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
    end
    hold off;
    
    %% using multi resolution watershed, speed up veta
%     
%     p.scales=[6];
%     disp('begin nuclei segmentation using watershed');
%     [nuclei, properties] = nucleiSegmentationV2(curIM_normRed,p);
%     
%     figure;imshow(curIM);hold on;
%     for k = 1:length(nuclei)
%         plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
%     end
%     hold off;
%     %% using general LoG, in folder testIM_German/GeneralLoG
%     %%% initial segmentation
%     R=curIM_normRed; I=curIM;
%     ac=30;    % remove isolated pixels(non-nuclei pixels)
% %     show(bwSpot)
% 
%     [Nmask,cs,rs,A3]=XNucleiSeg_GL_Thresholding(R,ac);      %% thresholding based binarization % show(R)
% %     show(Nmask);
%     %%% gLoG seeds detection
%     largeSigma=10;smallSigma=8; % Sigma value
%     ns=XNucleiCenD_Clustering(R,Nmask,largeSigma,smallSigma);  %% To detect nuclei clumps
% %     53 s
%     figure(1),imshow(I(:,:,1));
%     hold on,plot(ns(2,:),ns(1,:),'y+');
%      hold on,plot(cs,rs,'y+');
% %     hold on,plot(cs,rs,'g*');
%     
%     %%% marker-controlled watershed segmentation
% %     ind=sub2ind(size(R),ns(1,:),ns(2,:));
% %     bs4=zeros(size(R));
% %     bs4(ind)=1;
% %     [bnf,blm]=XWaterShed(Nmask,bs4);
% %     
% %     %%% show segmentations
% %     bb=bwperim(Nmask);  % for debugging finial segmentations
% %     ind5=find(bb);
% %     blm(ind5)=1;
% %     bb2=bwperim(A3);
% %     ind6=find(bb2);
% %     blm(ind6)=1;
% %     overlay1=imoverlay(I,blm,[0 1 0]);
% %     figure(2),imshow(overlay1);
end