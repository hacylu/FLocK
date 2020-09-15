% conversation regarding now to do modify nuclear segmentation for H&E stain
%
% so to do HE, its actually very simple
% i changed line 7 of the example
% to be this
% [c1,c2,c3]=colour_deconvolution(io,'HE');
%    I=c1;
%    tic; [nuclei properties] = nucleiSegmentation(I); toc;
% and it works perfectly
% the only other thing i had to change was a parameter or two, depending on how big the nuclei you want to segment are


% line 693 of nucleiSegmentation specifices a range of how big the *radius of the nuclei you want to segment are
