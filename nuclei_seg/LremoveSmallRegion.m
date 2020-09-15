% for remove the small regions from the bw by specified area thrshold 
% T is the specified area thrshold. 
% bw is the binary mask contains many regions.
function bw_remove=LremoveSmallRegion(bw,T)

cc= bwconncomp(bw);

stats = regionprops(cc, 'Area');
idx = find([stats.Area] > T);
bw_remove = ismember(labelmatrix(cc), idx);

