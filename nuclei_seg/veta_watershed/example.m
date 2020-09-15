%% segment the nuclei
io = imread('0682_A_1_3_.tif');

I = imread('0682_A_1_3_.tif');

if ~exist('0682_A_1_3_.mat', 'file'); %took out ~. code wasn't running image through nucleiSegmentation since the .mat file already existed
    %but .mat file already contains information about segmentation. revise
    %code to create .mat file if it doesn't exist
    [c1] = color_deconvolution(io, 'HE');
    I = c1;
    tic; [nuclei properties] = nucleiSegmentation(I(:,:,1)); toc;
    
    
    save('0682_A_1_3_.mat', 'nuclei', 'properties');
    
else
    
    load('0682_A_1_3_.mat');
    
end

%% heuristic to remove non-epithelial regions

dilateRadius = 5;
areaOpenThreshold = 8000;

img = double(I(:,:,1))/255;

B = img < graythresh(img);

B = imdilate(B, strel('disk', dilateRadius));
B = bwareaopen(B, areaOpenThreshold );

b = false(1, length(properties));

for l = 1:length(properties)
    
    c = round( properties(l).Centroid );
    
    b(l) = B(c(2), c(1));
    
end

nuclei = nuclei(b);
properties = properties(b);

%% display nuclei and nuclei area

displayNucleiProperties(I, nuclei, properties, 'Area');
