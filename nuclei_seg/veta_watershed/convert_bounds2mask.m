function seg_mask = convert_bounds2mask(img,bounds)
%seg_mask = convert_bounds2mask(img,bounds)

seg_mask = zeros(size(img,1), size(img,2));
for i = 1:numel(bounds)
    seg_mask = seg_mask + poly2mask(bounds(i).c,bounds(i).r,size(img,1),size(img,2));
end

% alternative slower implementation
% mask = false(img);
% for i = 1:numel(bounds)
%     mask = mask | poly2mask(bounds(i).r, bounds(i).c, size(img,1), size(img,2));
% end