function P = polygonClipperConvert(C)
% polygonClipperConvert: Converts the contours to structures that can be
% used with PolygonClipper.
%
% The contours are of the same format as the output of bwboundaries().
%

for k = length(C):-1:1
    
    P(k).x = C{k}(:,2);
    P(k).y = C{k}(:,1);
    P(k).hole = 0;
    
end

end