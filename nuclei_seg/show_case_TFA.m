addpath(genpath('./veta_watershed'));
addpath(genpath('./GeneralLoG'));
addpath(genpath('./staining_normalization'));
dbstop if error
% IMData4SymmetricVoting_1={'8913.tif','8914.tif'};
% IMData4SymmetricVoting_1={'29_38.png'};
% IMData4SymmetricVoting_1={'36_44.png'};
% I=imread('C:\Users\cxl884\Google Drive\SachethImages\images\392245-1.tif');
% IMData4SymmetricVoting_1={'C:\Users\cxl884\Google Drive\SachethImages\images\392245-1.tif'};
IMData4SymmetricVoting_1={'CLu.tif'};
for i=1:length(IMData4SymmetricVoting_1)
    curIMName=IMData4SymmetricVoting_1{i};
    curIM=imread(curIMName);
    curIMsize=size(curIM);
    [curIM_norm] = normalizeStaining(curIM);
    curIM_normRed=curIM_norm(:,:,3);

    %% using gLoG-Watershed 
    
%     p.scales=[4:2:6];
%     disp('begin nuclei segmentation using watershed');
%     [nuclei, properties] = nucleiSegmentationV2_gLoG(curIM_normRed,p);
%     
%     figure;imshow(curIM);hold on;
%     for k = 1:length(nuclei)
%         plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
%     end
%     hold off;
    
    %% using multi resolution watershed, speed up veta
    
%     p.scales=[4:2:6];
%     disp('begin nuclei segmentation using watershed');
%     [nuclei, properties] = nucleiSegmentationV2(curIM_normRed,p);
%     
%     figure;imshow(curIM);hold on;
%     for k = 1:length(nuclei)
%         plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
%     end
%     hold off;
    %% using general LoG, in folder testIM_German/GeneralLoG
    %%% initial segmentation
    R=curIM_normRed; I=curIM;%show(R)
    ac=10;    % remove isolated pixels(non-nuclei pixels)
%     show(bwSpot)

    [Nmask]=XNucleiSeg_GL_Thresholding_TFA(R,ac);      %% thresholding based binarization % show(R)
%     show(Nmask);
    %%% gLoG seeds detection
    largeSigma=7;smallSigma=5; % Sigma value

%     largeSigma=10;smallSigma=6; % Sigma value
    ns=XNucleiCenD_Clustering(R,Nmask,largeSigma,smallSigma);  %% To detect nuclei clumps
%     53 s

%%% remove the seed point with bright intensity

R=I(:,:,1);

idx_remove=[];
for i=1:length(ns)
    if R(ns(1,i),ns(2,i))>140;
        idx_remove=[idx_remove i];
    end
end
ns(:,idx_remove)=[];

figure(1),imshow(I(:,:,1));
hold on,plot(ns(2,:),ns(1,:),'y+');

end