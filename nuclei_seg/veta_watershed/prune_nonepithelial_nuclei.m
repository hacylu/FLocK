function [ep_nuclei, ep_properties, all_seg_info] = prune_nonepithelial_nuclei(I, nuclei, properties)
% b: index for epithelial nuclei


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

ep_nuclei = nuclei(b);
ep_properties = properties(b);

all_seg_info.nuclei = nuclei;
all_seg_info.properties = nuclei;
all_seg_info.is_epithelial = b;
