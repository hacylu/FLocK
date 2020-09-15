function [normalizedImg] = scaleNormalization(img, lowBound, uppBound)

%%% img is a grey-level one

maxV = max(max(img));
minV = min(min(img));

normalizedImg = lowBound + (img - minV)*(uppBound-lowBound)/(maxV-minV);
