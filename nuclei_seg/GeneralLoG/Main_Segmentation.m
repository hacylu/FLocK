function Main_Segmentation
clc;clear all;

load('I.mat');

R=I(:,:,1);

%% initial segmentation
ac=30;    % remove isolated pixels(non-nuclei pixels)
[Nmask,cs,rs,A3]=XNucleiSeg_GL_Thresholding(R,ac);      %% thresholding based binarization

%% gLoG seeds detection
 largeSigma=8;smallSigma=4; % Sigma value
 ns=XNucleiCenD_Clustering(R,Nmask,largeSigma,smallSigma);  %% To detect nuclei clumps
 
 figure(1),imshow(I);
 hold on,plot(ns(2,:),ns(1,:),'y+');
 hold on,plot(cs,rs,'g*');
 
%% marker-controlled watershed segmentation
 ind=sub2ind(size(R),ns(1,:),ns(2,:));
 bs4=zeros(size(R));
 bs4(ind)=1;
[bnf,blm]=XWaterShed(Nmask,bs4);

%% show segmentations
bb=bwperim(Nmask);  % for debugging finial segmentations
ind5=find(bb);
blm(ind5)=1;
bb2=bwperim(A3);
ind6=find(bb2);
blm(ind6)=1;
overlay1=imoverlay(I,blm,[0 1 0]);
figure(2),imshow(overlay1);
end