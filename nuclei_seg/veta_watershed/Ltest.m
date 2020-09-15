%% segment the nuclei
curIM = imread('1.png');

[curIM_norm, curH, curE] = normalizeStaining(curIM);

tic; [nuclei properties] = nucleiSegmentation(curH(:,:,1)); toc;

show(curIM);

save('0682_A_1_3_.mat', 'nuclei', 'properties');


%% heuristic to remove non-epithelial regions

dilateRadius = 5;
areaOpenThreshold = 8000;

img = double(curH(:,:,1))/255;

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

displayNucleiProperties(curH, nuclei, properties, 'Area');
