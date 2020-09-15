function bounds = nuclei2bounds(nuclei)

bounds = [];
parfor i = 1:numel(nuclei)
    x = nuclei{i};
    bounds(i).r = x(:,1);
    bounds(i).c = x(:,2);
    [bounds(i).centroid_r,bounds(i).centroid_c] = poly_centroid(x(:,1),x(:,2));
end