function [nuclei properties] = nucleiSegmentationV2(I, varargin)
% nucleiSegmentation: Performs nuclei segmentation in histopathology
% images using multi-scale and multi-marker watershed segmentation.
%
% Usage:
% ---------
% [nuclei properties] = nucleiSegmentation(I);
%
% or
%
% [nuclei properties] = nucleiSegmentation(I, 'OptionalParameter', value);
%
% Input:
% ---------
% I - Grayscale input image with intensities scaled in the [0 255] range.
%
% For optional input parameters please have a look at the getParameters()
% subroutine of this function.
%
% Output:
% ---------
% nuclei     - Nuclei contours in the same format as returned by
%              bwboundaries().
% properties - Properties of the segmented regions in the same format as
%              returned by regionprops().
%
% Copyright (c) 2012, Mitko Veta
% Image Sciences Institute
% University Medical Center
% Utrecht, The Netherlands
% All rights reserved.
%
% Cheng Lu from CCIDP speeded up the merging procedure on Oct. 3, 2015

p = getParameters(varargin{:});

% make sure the input image is in double format
I = double(I);

nuclei = cell(1, length(p.scales));
properties = cell(1, length(p.scales));
keep =[]; accepted = [];

%if isempty(gcp('nocreate')), pinfo = parcluster('local'); numpools = pinfo.NumWorkers; parpool(numpools-1); end
for n = 1:length(p.scales)
    
    disp(['Processing at scale ' num2str(p.scales(n)) ' pixels.']);
    
    Ip = morphologicalProcessing(I, p.scales(n));
    
    S = fastRadialSymmetryTransform(Ip, p.scales(n):2*p.scales(n), p.alpha, p.beta, p.kappa);
    
    LRadialSymmetry = markerControlledWatershedRadialSymmetry(Ip, S, p.h, p.scales(n), p.noBackgroundMarkers);
    
    [ERadialSymmetry LRadialSymmetry] = bwboundaries(LRadialSymmetry, 'noholes');
    
    PRadialSymmetry = regionProperties(I, LRadialSymmetry);
    
    LRegionalMinima = markerControlledWatershedRegionalMinima(Ip, p.scales(n), p.noBackgroundMarkers);
    
    [ERegionalMinima LRegionalMinima] = bwboundaries(LRegionalMinima, 'noholes');
    
    PRegionalMinima = regionProperties(I, LRegionalMinima);
    
    try
        P = cat(1, PRadialSymmetry, PRegionalMinima);
        E = cat(1, ERadialSymmetry, ERegionalMinima);
    catch ME
        fprintf('\ncatch statement added by George Lee - CCIPD, CWRU - 1/18/15: \n%s',ME.message);
        
        P1 = numel(PRadialSymmetry);
        P2 = numel(PRegionalMinima);
        
        E1 = numel(ERadialSymmetry); 
        E2 = numel(ERegionalMinima);
        
        if P1>P2,		P = PRadialSymmetry;
        else		P = PRegionalMinima;
        end
        P(1)
        
        if E1>E2,		E = ERadialSymmetry;
        else		E = ERegionalMinima;
        end
        E(1)
        
    end
    
    
    if exist('P','var')
        if ~isempty(P)
            % remove regions with area that is too small or too large for current
            % scale
            accepted = ruleBasedClassification(P, {'Area', p.scales(n)^2*pi * [1 4]});
            
            % rule-based classification (rejection of false nuclei)
            accepted = accepted & ruleBasedClassification(P, p.featureRanges);
            
            nuclei{n} = E(accepted);
            properties{n} = P(accepted);
        else
            fprintf('\nP not found for scale %i from scales [%i %i] ', p.scales(n), p.scales(1), p.scales(end))
        end
    else
        fprintf('\nP not found for scale %i from scales [%i %i] ', p.scales(n), p.scales(1), p.scales(end))
    end
    
    if ~isempty(accepted)
        keep = [keep n];
    end
    
end
%% changed by George Lee - CCIPD, CWRU - 1/18/2015
% reconciling potential issues with differing fields
% nuclei{1}
% properties{1}
% updated_nuclei = {};
% updated_properties = {};
% for i = 1:numel(nuclei)
% % 	if isstruct(nuclei{i})
% % 	if isempty(fields(nuclei{i})) || isempty(fields(properities{i}))
% % 		keep = [keep i];
% % 	end
%
% updated_nuclei(end+1) = nuclei{n};
% updated_properties(end+1) = properties{n};
% end

% nuclei{keep};
% properties{keep};

nuclei = cat(1, nuclei{keep});
properties = cat(1, properties{keep});
%% changed by George Lee - CCIPD, CWRU - 1/18/2015 ^^^

%
% nuclei = cat(1, nuclei{:});
% properties = cat(1, properties{:});


fitness = cat(1, properties.Solidity);
centroids = cat(1, properties.Centroid);

disp('Merging contours.');

% tic
idx = mergeContours_s(nuclei, fitness, p.Th, centroids);
% toc

nuclei = nuclei(idx);
properties = properties(idx);

end

%==========================================================================
% Main processing blocks

function Ip = morphologicalProcessing(I, n)

% imreconstruct is faster with single (also logical and uint8)
I = single(I);

% disk-shaped flat structuring elements
B = strel('disk', n, 0);
B2 = strel('disk', floor(n/2), 0);

% stage 1 (opening by reconstruction): remove bright unconnected objects
% smaller than the structuring element
marker1 = imerode(I, B);
stage1 = imreconstruct(marker1, I);

t = imcomplement(stage1);

% stage 2 (closing by reconstruction): remove dark unconnected objects
% smaller than the structuring element
marker2 = imerode(t, B);
stage2 = imreconstruct(marker2, t);

t = imcomplement(stage2);

% simplify the shape of the remaining dark objects, disconnect loosly
% connected objects
Ip = imclose(t, B2);

% convert (back) to double
Ip = double(Ip);

end

function S = fastRadialSymmetryTransform(I, radii, alpha, beta, kappa)

[rows cols] = size(I);

[gm dx dy] = imSobel(I);

gm = gm + eps;

dx = dx./gm;
dy = dy./gm;

S = zeros(rows,cols);

[x,y] = meshgrid(1:cols, 1:rows);

% voting coefficients: all set at -1 since the orientation only
% dark fast radial symmetry transform is needed
val = -(gm(:) > beta);

for n = 1:numel(radii);
    
    r = radii(n);
    
    negx = x - round(r*dx);
    negy = y - round(r*dy);
    
    negx( negx<1 ) = 1;
    negx( negx>cols ) = cols;
    negy( negy<1 ) = 1;
    negy( negy>rows ) = rows;
    
    O = accumarray([negy(:), negx(:)], val, size(I));
    
    if r < numel(kappa)
        kappa = kappa(r);
    else
        kappa = kappa(end);
    end
    
    O(O > kappa) = kappa;
    O(O < -kappa) = -kappa;
    
    Fo = sign(O) .* (abs(O) ./ kappa) .^ alpha;
    
    % accumulate over all radii
    S = S + r*imGauss(Fo, 0.25*r);
    
end

% normalize
S = S ./ numel(radii);

end

function L = markerControlledWatershedRadialSymmetry(I, S, h, n, noBackgroundMarkers)

% foreground markers
fg = imextendedmin(S, h);%show(fg);show(fg);
fg = imclearborder(fg);

% provisional nuclei map
t = imdilate(fg, strel('disk', 2*n));

% background markers
if noBackgroundMarkers
    bg = false(size(fg));
else
    bg = bwmorph(~t, 'thin', Inf);
end

% segmentation function
g = imSobel(I);

% modified segmentation function
D = imimposemin(g, fg | bg);

% marker-controlled watershed
L = double(watershed(D));

% clear background regions
L = L .* ismember(L, unique(L(fg)));

% clear regions on the border
L = imclearborder(L);

% reorder the labels
L = orderLabels(L);

end

function L = markerControlledWatershedRegionalMinima(I, n, noBackgroundMarkers)

% foreground markers
fg = imregionalmin(I);
fg = imclearborder(fg);
fg = imerode(fg, strel('disk', 1));

% provisional nuclei map
t = imdilate(fg, strel('disk', 2*n));

% background markers
if noBackgroundMarkers
    bg = false(size(fg));
else
    bg = bwmorph(~t, 'thin', Inf);
end

% segmentation function
g = imSobel(I);

% modified segmentation function
D = imimposemin(g, fg | bg);

% marker-controlled watershed
L = double(watershed(D));

% clear background regions
L = L .* ismember(L, unique(L(fg)));

% clear regions on the border
L = imclearborder(L);

% reorder the labels
L = orderLabels(L);

end

function accepted = ruleBasedClassification(properties, ranges)

accepted = true(length(properties), 1);

for k = 1:2:length(ranges)
    
    v = cat(1, properties.(ranges{k}));
    r = ranges{k+1};
    
    accepted = accepted & ( v >= r(1) ) & ( v <= r(2) );
    
end

end


%==========================================================================
% Additional utilities

function P = regionProperties(I, L)

PADDING = 1.5;

props = {
    'Area'
    'Perimeter'
    'Eccentricity'
    'MajorAxisLength'
    'MinorAxisLength'
    'EquivDiameter'
    'Solidity'
    'Orientation'
    'Centroid'
    'WeightedCentroid'
    };

% NOTE: the image is inverted
P = regionprops(L, 255-I, props);

% NOTE: the image is inverted
% extract separate images and label images (masks) for each region
[R M D] = labelmatrixToROI(255-I, L, PADDING, [], P);

% calculate additional properties
for k = 1:length(P)
    
    c = P(k).Centroid;
    w = P(k).WeightedCentroid;
    d = P(k).EquivDiameter;
    a = P(k).Area;
    l = P(k).Perimeter;
    
    region = R{k};
    regionMask = M{k};
    distance = D{k};
    
    pixelValues = region(regionMask);
    
    % circularity:
    circularity = 4 * pi * a / l^2;
    
    % elliptical deviation:
    % difference between the region and the elliptical approximation of the
    % region
    ellipse = distance <= 1;
    intersection = ellipse & regionMask;
    ellipticalDeviation = 1 - 2*sum(intersection(:)) / (P(k).Area + sum(ellipse(:)));
    
    % mass displacement:
    % distance betweend centroid and weighted centroid normalized with
    % equivalent diameter
    massDisplacement = sqrt(sum((c - w).^2))/d;
    
    % integrated intensity:
    integratedIntensity = sum(pixelValues);
    
    % mean intensity:
    meanIntensity = mean(pixelValues);
    
    % intensity deviation:
    intensityDeviation = std(pixelValues);
    
    % intensity range:
    intensityRange = prctile(pixelValues, 97.5)...
        - prctile(pixelValues, 2.5);
    
    % the inside boundary is defined as the residual between the image
    % region and its erosion with an isotropic strucutring element with
    % radius 1/8 of the equivalent diameter of the region
    se = strel('disk', round(d/8), 0);
    insideBoundary = xor(regionMask, imerode(regionMask, se));
    outsideBoundary = xor(regionMask, imdilate(regionMask, se));
    
    insideBoundaryValues = region(insideBoundary);
    outsideBoundaryValues = region(outsideBoundary);
    
    % inside boundary intensity statistics:
    meanInsideBoundaryIntensity = mean(insideBoundaryValues );
    insideBoundaryIntensityDeviation = std(insideBoundaryValues );
    insideBoundaryIntensityRange = prctile(insideBoundaryValues , 97.5) - prctile(insideBoundaryValues , 2.5);
    normalizedInsideBoundaryIntensity = meanInsideBoundaryIntensity / meanIntensity;
    
    % outside boundary intensuty statistics:
    meanOutsideBoundaryIntensity = mean(outsideBoundaryValues );
    outsideBoundaryIntensityDeviation = std(outsideBoundaryValues );
    outsideBoundaryIntensityRange = prctile(outsideBoundaryValues , 97.5) - prctile(outsideBoundaryValues , 2.5);
    normalizedOutsideBoundaryIntensity = meanOutsideBoundaryIntensity / meanIntensity;
    
    % boundary saliency:
    boundarySaliency = meanInsideBoundaryIntensity - meanOutsideBoundaryIntensity;
    normalizedBoundarySaliency = normalizedInsideBoundaryIntensity - normalizedOutsideBoundaryIntensity;
    
    % add to structure
    P(k).Circularity = circularity;
    P(k).EllipticalDeviation = ellipticalDeviation;
    P(k).MassDisplacement = massDisplacement;
    P(k).IntegratedIntensity =  integratedIntensity;
    P(k).MeanIntensity = meanIntensity;
    P(k).IntensityDeviation = intensityDeviation;
    P(k).IntensityRange = intensityRange;
    P(k).MeanInsideBoundaryIntensity = meanInsideBoundaryIntensity;
    P(k).InsideBoundaryIntensityDeviation = insideBoundaryIntensityDeviation;
    P(k).InsideBoundaryIntensityRange = insideBoundaryIntensityRange;
    P(k).NormalizedInsideBoundaryIntensity = normalizedInsideBoundaryIntensity;
    P(k).MeanOutsideBoundaryIntensity = meanOutsideBoundaryIntensity;
    P(k).OutsideBoundaryIntensityDeviation = outsideBoundaryIntensityDeviation;
    P(k).OutsideBoundaryIntensityRange = outsideBoundaryIntensityRange;
    P(k).NormalizedOutsideBoundaryIntensity = normalizedOutsideBoundaryIntensity;
    P(k).BoundarySaliency = boundarySaliency;
    P(k).NormalizedBoundarySaliency = normalizedBoundarySaliency;
    
end

end

function [R M D] = labelmatrixToROI(I, L, padding, interpType, P)
% labelmatrixToROI: Extracts regions from an image given with a labelmatrix
% as separate images such that the major and minor axes of the regions are
% aligned with the image axes.

if  ~exist('padding', 'var') || isempty(padding)
    padding = 1;
end

if  ~exist('interpType', 'var') || isempty(interpType)
    interpType = 'linear';
end

% the properties matrix can be given as an input if already computed
if  ~exist('P', 'var') || isempty(P)
    props = {
        'Centroid'
        'Orientation'
        'MajorAxisLength'
        'MinorAxisLength'
        };
    
    P = regionprops(L, props);
end

R = cell(length(P), 1);
M = cell(length(P), 1);
D = cell(length(P), 1);

iI = griddedInterpolant(I, interpType);
iL = griddedInterpolant(L, 'nearest');

warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

for k = 1:length(P)
    
    % location, size and orientation of the image
    c = P(k).Centroid;
    a = P(k).MajorAxisLength;
    b = P(k).MinorAxisLength;
    t = P(k).Orientation;
    
    w = round(a * padding);
    h = round(b * padding);
    phi = degToRad(t);
    
    % form a sampling grid
    x = -(h/2-0.5):(h/2-0.5);
    y = -(w/2-0.5):(w/2-0.5);
    
    [X Y] = ndgrid(x, y);
    
    % rotate and center
    A = rotationTranslationMatrix(phi, c);
    T = maketform('affine', A);
    
    [V U] = tformfwd(T, Y, X);
    
    % sample image and label matrix
    %region = interp2(I, U, V, interpType, 0);
    %regionMask = interp2(L, U, V, '*nearest', 0) == k;
    region = iI(U, V);
    regionMask = iL(U, V)==k;
    
    % NOTE: The sampling of the region mask assumes propper region
    % ordering: k-th region has label k.
    
    % normalized distance from center
    distance = sqrt((2*X/a).^2 + (2*Y/b).^2);
    
    region(isnan(region)) = 0;
    regionMask(isnan(regionMask)) = false;
    
    R{k} = region;
    M{k} = regionMask;
    D{k} = distance;
    
end

warning('on', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

end

%==========================================================================
% Helpers

function b = checkNumericBounds(x, lowerBound, higherBound)

if ~exist('lowerBound', 'var')
    lowerBound = -Inf;
end

if ~exist('higherBound', 'var')
    higherBound = Inf;
end

x = x(:);

b = isnumeric(x) && all(x >= lowerBound) && all(x <= higherBound);

end

function radians = degToRad(degrees)
% degToRad: Converts degrees to radians.

radians = degrees * pi / 180;

end

function G = imGauss(I, s)
% imGauss: Gaussian image filtering.

sze = roundToOdd(4*s);

f = fspecial('gaussian', sze, s);

G = imfilter(I, f, 'symmetric', 'conv');

end

function [gm dx dy] = imSobel(I)
% imSobel: Sobel edge detection.

s = [1 0 -1; 2 0 -2; 1 0 -1];

dx = imfilter(I, s , 'symmetric', 'conv');
dy = imfilter(I, s', 'symmetric', 'conv');

gm = sqrt(dx.^2 + dy.^2);

end

function L = orderLabels(L)
% orderLabels: Ensures the labels of the labelmatrix go from 1 to the
% number of labels.

sz = size(L);

L = L(:);

% in this case: L = b(n)
[b, ~, n] = unique(L);

% the new labels
l = 0:length(b);

L = l(n);

L = reshape(L, sz);

end


function A = rotationTranslationMatrix(theta, c)
% rotationTranslationMatrix: Returns a rotation ans translation matrix for
% a given angle and offset.

if ~exist('theta', 'var') || isempty(theta)
    theta = 0;
end

if ~exist('c', 'var') || isempty(c)
    c = zeros(1, 2);
end

A = [ cos(theta) -sin(theta) 0  ;
    sin(theta)  cos(theta) 0  ;
    c(1)        c(2)       1 ];

end

function x = roundToOdd(x)
% roundToOdd: Rounds the number to the nearest odd integer.

x = ceil(x);
x = x + mod(x-1,2);

end

%==========================================================================
% Parameters parsing

function p = getParameters(varargin)

ip = inputParser;

% scales for multi-scale processing in micrometers:
ip.addOptional('scales', 6:13, @(x)checkNumericBounds(x, 0));
%1:5
%6:13

% fast radial symmetry transform properties (see Loy & Zelinski paper):
ip.addOptional('alpha', 1, @(x)checkNumericBounds(x, 0));

ip.addOptional('beta', 60, @(x)checkNumericBounds(x, 0));

ip.addOptional('kappa', 10, @(x)checkNumericBounds(x, 0));

% h parameter for extended regional minima (see imextendedmin.m):
ip.addOptional('h', 0.6, @(x)checkNumericBounds(x, 0));

% do not calculate background markers (faster):
ip.addOptional('noBackgroundMarkers', false, @islogical);

% properties ranges for rule based rejection:
ranges = {'Solidity', [0.875 1],...
    'MassDisplacement', [0 0.08],...
    'BoundarySaliency', [20 255]
    };

ip.addOptional('featureRanges', ranges, @iscell);

% threshold value for contour merging:
ip.addOptional('Th', 0.2,...
    @(x)checkNumericBounds(x, 0));

% parse properties
ip.parse(varargin{:});
p = ip.Results;

end
