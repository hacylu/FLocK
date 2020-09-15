function a = polygonClipperArea(P)
% polygonClipperArea: Calculates the area of the contour given as a
% PolygonClipper structure.

a = zeros(1, length(P));

for k = 1:length(P)
    
    a(k) = polyarea(P(k).x, P(k).y);
    
end

end