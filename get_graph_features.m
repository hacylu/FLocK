function [vfeature,GraphFeatureDescription] = get_graph_features(x,y)

% graphfeats    Calculates graph-based features for nuclear centroids
% located at (x,y) in the image. 
% 
% Necessary input:
% x,y: x and y coordinates of points that will be used for graph construction
% (i.e. nuclear centroids). 

% Output Description: vfeature contains the following:

% Voronoi Features
% 1: Area Standard Deviation
% 2: Area Average
% 3: Area Minimum / Maximum
% 4: Area Disorder
% 5: Perimeter Standard Deviation
% 6: Perimeter Average
% 7: Perimeter Minimum / Maximum
% 8: Perimeter Disorder
% 9: Chord Standard Deviation
% 10: Chord Average
% 11: Chord Minimum / Maximum
% 12: Chord Disorder

% Delaunay Triangulation
% 13: Side Length Minimum / Maximum
% 14: Side Length Standard Deviation
% 15: Side Length Average
% 16: Side Length Disorder
% 17: Triangle Area Minimum / Maximum
% 18: Triangle Area Standard Deviation
% 19: Triangle Area Average
% 20: Triangle Area Disorder

% Minimum Spanning Tree
% 21: MST Edge Length Average
% 22: MST Edge Length Standard Deviation
% 23: MST Edge Length Minimum / Maximum
% 24: MST Edge Length Disorder

% Nuclear Features
% 25: Area of polygons
% 26: Number of nuclei
% 27: Density of Nuclei
% 28: Average distance to 3 Nearest Neighbors
% 29: Average distance to 5 Nearest Neighbors
% 30: Average distance to 7 Nearest Neighbors
% 31: Standard Deviation distance to 3 Nearest Neighbors
% 32: Standard Deviation distance to 5 Nearest Neighbors
% 33: Standard Deviation distance to 7 Nearest Neighbors
% 34: Disorder of distance to 3 Nearest Neighbors
% 35: Disorder of distance to 5 Nearest Neighbors
% 36: Disorder of distance to 7 Nearest Neighbors
% 37: Avg. Nearest Neighbors in a 10 Pixel Radius
% 38: Avg. Nearest Neighbors in a 20 Pixel Radius
% 39: Avg. Nearest Neighbors in a 30 Pixel Radius
% 40: Avg. Nearest Neighbors in a 40 Pixel Radius
% 41: Avg. Nearest Neighbors in a 50 Pixel Radius
% 42: Standard Deviation Nearest Neighbors in a 10 Pixel Radius
% 43: Standard Deviation Nearest Neighbors in a 20 Pixel Radius
% 44: Standard Deviation Nearest Neighbors in a 30 Pixel Radius
% 45: Standard Deviation Nearest Neighbors in a 40 Pixel Radius
% 46: Standard Deviation Nearest Neighbors in a 50 Pixel Radius
% 47: Disorder of Nearest Neighbors in a 10 Pixel Radius
% 48: Disorder of Nearest Neighbors in a 20 Pixel Radius
% 49: Disorder of Nearest Neighbors in a 30 Pixel Radius
% 50: Disorder of Nearest Neighbors in a 40 Pixel Radius
% 51: Disorder of Nearest Neighbors in a 50 Pixel Radius

load('GraphFeatureDescription.mat')

% Calculate the Voronoi diagram.
[VX,VY] = voronoi(x,y);
[V, C] = voronoin([x(:),y(:)]);

% Okay, so:
% VX, VY - These guys contain the vertices in a way such that
% plot(VX,VY,'-',x,y,'.') creates the voronoi diagram. I don't actually
% think these are used later on, but I'm keeping them here just in case.

% C - This variable is an m by 1 cell array, where m is the number of cell
% centroids in your image. Each element in C is a vector
% with the coordinates of the vertices of that row's voronoi polygon. 
% V - This is a q by 2 matrix, where q is the number of vertices and 2 is
% the number of dimensions of your image. Each element in V contains the
% location of the vertex in 2D space.
% The idea here is that if you want to see the coordinates for the vertices
% of polygon 5, for example, you would go:
%     X = V(C{5},:)
% which would display five rows, each with the 2D coordinates of the vertex
% of polygon 5.

% Get the delaunay triangulation...
del = delaunay(x,y);    

% Returns a set of triangles such that no data points are contained in any 
% triangle's circumcircle. Each row of the numt-by-3 matrix TRI defines 
% one such triangle and contains indices into the vectors X and Y. When 
% the triangles cannot be computed (such as when the original data is 
% collinear, or X is empty), an empty matrix is returned.
    
% Get the Minimum Spanning Tree (MST) (optional: plot)
[mst.edges mst.edgelen mst.totlen] = mstree([x,y],[],0);

% Record indices of inf and extreme values to skip these cells later
Vnew        = V;
Vnew(1,:)   = [];

% Find the data points that lie far outside the range of the data
[Vsorted,I]     = sort([Vnew(:,1);Vnew(:,2)]);
N               = length(Vsorted);
Q1              = round(0.25*(N+1));
Q3              = round(0.75*(N+1));
IQR             = Q3 - Q1;
highrange       = Q3 + 1.5*IQR;
lowrange        = Q1 - 1.5*IQR;
Vextreme        = [];
Vextreme        = [Vextreme; V(find(V > highrange))];
Vextreme        = [Vextreme; V(find(V < lowrange))];

banned = [];
for i = 1:length(C)
    if(~isempty(C{i}))
        
    if(max(any(isinf(V(C{i},:)))) == 1 || max(max(ismember(V(C{i},:),Vextreme))) == 1)
        banned = [banned, i];
    end
    end
end
% If you've eliminated the whole thing (or most of it), then only ban 
% indices that are infinity (leave the outliers)
if(length(banned) > length(C)-2)
    banned = [];
    for i = 1:length(C)
        if(max(any(isinf(V(C{i},:)))) == 1)
            banned = [banned, i];
        end
    end
end

% Voronoi Diagram Features
% Area
c = 1;
d = 1;
e = d;
for i = 1:length(C)

    if(~ismember(i,banned) && ~isempty(C{i}))
        X = V(C{i},:);
        chord(1,:) = X(:,1);
        chord(2,:) = X(:,2);
        % Calculate the chord lengths (each point to each other point)
        for ii = 1:size(chord,2)
            for jj = ii+1:size(chord,2)
                chorddist(d) = sqrt((chord(1,ii) - chord(1,jj))^2 + (chord(2,ii) - chord(2,jj))^2);
                d = d + 1;
            end
        end

        % Calculate perimeter distance (each point to each nearby point)
        for ii = 1:size(X,1)-1
            perimdist(e) = sqrt((X(ii,1) - X(ii+1,1))^2 + (X(ii,2) - X(ii+1,2))^2);
            e = e + 1;
        end
        perimdist(size(X,1)) = sqrt((X(size(X,1),1) - X(1,1))^2 + (X(size(X,1),2) - X(1,2))^2);
        
        % Calculate the area of the polygon
        area(c) = polyarea(X(:,1),X(:,2));
        c = c + 1;
        clear chord X
    end
end
if(~exist('area','var'))
    vfeature = zeros(1,51);
    return; 
end
vfeature(1) = std(area); 
vfeature(2) = mean(area);
vfeature(3) = min(area) / max(area);
vfeature(4) = 1 - ( 1 / (1 + (vfeature(1) / vfeature(2))) );

vfeature(5) = std(perimdist);
vfeature(6) = mean(perimdist);
vfeature(7) = min(perimdist) / max(perimdist);
vfeature(8) = 1 - ( 1 / (1 + (vfeature(5) / vfeature(6))) );

vfeature(9) = std(chorddist);
vfeature(10) = mean(chorddist);
vfeature(11) = min(chorddist) / max(chorddist);
vfeature(12) = 1 - ( 1 / (1 + (vfeature(9) / vfeature(10))) );

% Delaunay 
% Edge length and area
c = 1;
d = 1;
for i = 1:size(del,1)
    t = [x(del(i,:)),y(del(i,:))];
    
    sidelen(c:c+2) = [sqrt( ( t(1,1) - t(2,1) )^2 + (t(1,2) - t(2,2))^2 ), ...
        sqrt( ( t(1,1) - t(3,1) )^2 + (t(1,2) - t(3,2))^2 ), ...
        sqrt( ( t(2,1) - t(3,1) )^2 + (t(2,2) - t(3,2))^2 )];
    dis(i,1:3) = sum( sidelen(c:c+2) );
    c = c + 3;
    triarea(d) = polyarea(t(:,1),t(:,2));
    d = d + 1;
end

vfeature(13) = min(sidelen) / max(sidelen);
vfeature(14) = std(sidelen); 
vfeature(15) = mean(sidelen);
vfeature(16) = 1 - (1 / (1 + (vfeature(14) / vfeature(15)) ) );

vfeature(17) = min(triarea) / max(triarea);
vfeature(18) = std(triarea);
vfeature(19) = mean(triarea);
vfeature(20) = 1 - (1 / (1 + (vfeature(18) / vfeature(19))) );


% MST: Average MST Edge Length
% The MST is a tree that spans the entire population in such a way that the
% sum of the Euclidian edge length is minimal.

vfeature(21) = mean(mst.edgelen);
vfeature(22) = std(mst.edgelen);
vfeature(23) = min(mst.edgelen) / max(mst.edgelen);
vfeature(24) = 1 - ( 1 / ( 1 + (vfeature(22)/vfeature(21)) ) ); 

% Nuclear Features
% Density
vfeature(25) = sum(area); 
vfeature(26) = size(C,1);
vfeature(27) = vfeature(26) / vfeature(25);

% Average Distance to K-NN
% Construct N x N distance matrix:
for i = 1:size(x,1)
    for j = 1:size(x,1)
        distmat(i,j) = sqrt( (x(i) - x(j))^2 + (y(i) - y(j))^2 );
    end
end
DKNN = zeros(3,size(distmat,1));
kcount = 1;
for K = [3,5,7]
    % Calculate the summed distance of each point to it's K nearest neighbours

    for i = 1:size(distmat,1)
        tmp = sort(distmat(i,:),'ascend');

        % NOTE: when finding the summed distance, throw out the first result,
        % since it's the zero value at distmat(x,x). Add 1 to K to compensate.
        DKNN(kcount,i) = sum(tmp(2:K+1));
    end
    kcount = kcount + 1;
end
vfeature(28) = mean(DKNN(1,:));
vfeature(29) = mean(DKNN(2,:));
vfeature(30) = mean(DKNN(3,:));

vfeature(31) = std(DKNN(1,:));
vfeature(32) = std(DKNN(2,:));
vfeature(33) = std(DKNN(3,:));

vfeature(34) = 1 - (1 / ( 1 + (vfeature(31) / (vfeature(28)+eps)) ));
vfeature(35) = 1 - (1 / ( 1 + (vfeature(32) / (vfeature(29)+eps)) ));
vfeature(36) = 1 - (1 / ( 1 + (vfeature(33) / (vfeature(30)+eps)) ));

% NNRR_av: Average Number of Neighbors in a Restricted Radius
% Set the number of pixels within which to search
rcount = 1;
for R = [10:10:50]

    % For each point, find the number of neighbors within R pixels
    for i = 1:size(distmat,1)

        % NOTE: as above with the K-NN calculation, we subtract 1 from the
        % number of pixels found, because this corresponds to the diagonal of
        % the distance matrix, which is always 0 (i.e. distmat(x,x) = 0 for all
        % x).
        NNRR(rcount,i) = length( find( distmat(i,:) <= R ) ) - 1;
    end
    if(sum(NNRR(rcount,:)) == 0)
        eval(['NNRR_av_' num2str(R) ' = 0;']);
        eval(['NNRR_sd_' num2str(R) ' = 0;']);
        eval(['NNRR_dis_' num2str(R) ' = 0;']);
    else
        eval(['NNRR_av_' num2str(R) ' = mean(NNRR(rcount,:));']);
        eval(['NNRR_sd_' num2str(R) ' = std(NNRR(rcount,:));']);
        eval(['NNRR_dis_' num2str(R) ' = 1 - (1 / (1 + (NNRR_sd_' num2str(R) '/NNRR_av_' num2str(R) ')));']);
    end
    
    rcount = rcount + 1;
end

vfeature(37) = NNRR_av_10;
vfeature(38) = NNRR_av_20;
vfeature(39) = NNRR_av_30;
vfeature(40) = NNRR_av_40;
vfeature(41) = NNRR_av_50;
vfeature(42) = NNRR_sd_10;
vfeature(43) = NNRR_sd_20;
vfeature(44) = NNRR_sd_30;
vfeature(45) = NNRR_sd_40;
vfeature(46) = NNRR_sd_50;
vfeature(47) = NNRR_dis_10;
vfeature(48) = NNRR_dis_20;
vfeature(49) = NNRR_dis_30;
vfeature(50) = NNRR_dis_40;
vfeature(51) = NNRR_dis_50;


%--------------------------------------------------------------------------
% MSTREE: Minimum spanning tree
%
%     Syntax: [edges,edgelen,totlen] = mstree(crds,{labels},{doplot})
%
%         crds =    [N x p] matrix of point coordinates.
%         labels =  optional [N x 1] vector of numeric or character labels.
%         doplot =  optional boolean flag indicating whether (=1) or not (=0) 
%                     to plot the points and minimum spanning tree in the space 
%                     of the first two axes [default = 0].
%         ---------------------------------------------------------------------
%         edges =   [(N-1) x 2] list of points defining edges.
%         edgelen = [(N-1) x 1] vector of edge lengths corresponding to
%                     'edges'.
%         totlen =  total length of tree (sum of edge lengths).
%

% RE Strauss, 7/9/96
%   2/9/99 -   produces plot only for 2D coordinates, but will find the tree
%                for any number of dimensions.
%   9/7/99 -   changed plot colors for Matlab v5.
%   10/8/03 -  added point numbers to plot.
%   10/14/03 - sort edges by increasing point id.
%   12/5/03 -  if p>2, produce plot for first two axes; suppress plot for p=1;
%                use rowtoval() to sort edges.
%   12/9/03 -  added optional plot labels.

function [edges,edgelen,totlen] = mstree(crds,labels,doplot)
  if (nargin < 2) labels = []; end;
  if (nargin < 3) doplot = []; end;

  [n,p] = size(crds);

  if (isempty(doplot) | p<2)
    doplot = 0;
  end;
  
  if (~isempty(labels))
    if (~ischar(labels))
      labels = tostr(labels(:));
    end;
  end;

  totlen = 0;
  edges = zeros(n-1,2);
  edgelen = zeros(n-1,1);

  dist = eucl(crds);            % Pairwise euclidean distances
  highval = max(max(dist))+1;


  e1 = 1*ones(n-1,1);           % Initialize potential-edge list
  e2 = [2:n]';
  ed =  dist(2:n,1);

  for edge = 1:(n-1)            % Find shortest edges
    [mindist,i] = min(ed);        % New edge
    t = e1(i);
    u = e2(i);
    totlen = totlen + mindist;
    if (t<u)                      % Add to actual-edge list
      edges(edge,:) = [t u];
    else
      edges(edge,:) = [u t];
    end;
    edgelen(edge) = mindist;

    if (edge < n-1)
      i = find(e2==u);              % Remove new vertex from
      e1(i) = 0;                    %   potential-edge list
      e2(i) = 0;
      ed(i) = highval;

      indx = find(e1>0);
      for i = 1:length(indx)        % Update potential-edge list
        j = indx(i);
        t = e1(j);
        v = e2(j);
        if (dist(u,v) < dist(t,v))
          e1(j) = u;
          ed(j) = dist(u,v);
        end;
      end;
    end;
  end;
  
  v = rowtoval(edges);
  [v,edges,edgelen] = sortmat(v,edges,edgelen);  % Sort by increasing pt identifier

  if (doplot)
    figure;
    plot(crds(:,1),crds(:,2),'ok');
    putbnd(crds(:,1),crds(:,2));
    
    deltax = 0.018*range(crds(:,1));
    deltay = 0.02*range(crds(:,2));
    for i = 1:n
      if (isempty(labels))
        lab = int2str(i);
      else
        lab = labels(i,:);
      end;
%       text(crds(i,1)+deltax,crds(i,2)+deltay,lab);
    end;
    
    hold on;
    for i = 1:(n-1)
      t = edges(i,1);
      u = edges(i,2);
      x = [crds(t,1); crds(u,1)];
      y = [crds(t,2); crds(u,2)];
      plot(x,y,'-k');
    end;
    hold off;
  end;

  return;
% EUCL: Calculates the euclidean distances among a set of points, or between a
%       reference point and a set of points, or among all possible pairs of two 
%       sets of points, in P dimensions.  Returns a single distance for two points.
%
%     Syntax: dists = eucl(crds1,crds2)
%
%        crds1 = [N1 x P] matrix of point coordinates.  If N=1, it is taken to
%                   be the reference point.
%        crds2 = [N2 x P] matrix of point coordinates.  If N=1, it is taken to
%                   be the reference point.
%        -----------------------------------------------------------------------
%        dists = [N1 x N1] symmetric matrix of pairwise distances (if only crds1
%                   is specified);
%                [N1 x 1]  col vector of euclidean distances (if crds1 & ref
%                   are specified);
%                [1 x N2]  row vector of euclidean distances (if ref & crds2
%                   are specified);
%                [N1 x N2] rectangular matrix of pairwise distances (if crds1
%                   & crds2 are specified);
%                [1 x 1]   scalar (if crds1 is a [2 x P] matrix or ref1 & ref2
%                   are specified).
%

% RE Strauss, 5/4/94
%   10/28/95 - output row (rather than column) vector for the (reference
%               point)-(set of points) case; still outputs column vector for the
%               (set of points)-(reference point) case.
%   10/30/95 - for double for-loops, put one matrix-col access in outer loop
%               to increase speed.
%   10/12/96 - vectorize inner loop to increase speed.
%    6/12/98 - allow for P=1.
%   11/11/03 - initialize dists to NaN for error return.

function dists = eucl(crds1,crds2)
  if (~nargin) help eucl; return; end;
  dists = NaN;

  if (nargin < 2)                     % If only crds1 provided,
    [N,P] = size(crds1);
    if (N<2)
      error('  EUCL: need at least two points');
    end;

    crds1 = crds1';                     % Transpose crds
    dists = zeros(N,N);                 % Calculate pairwise distances

    for i = 1:N-1
      c1 = crds1(:,i) * ones(1,N-i);
      if (P>1)
        d = sqrt(sum((c1-crds1(:,(i+1:N))).^2));
      else
        d = abs(c1-crds1(:,(i+1:N)));
      end;
      dists(i,(i+1:N)) = d;
      dists((i+1:N),i) = d';
    end;
    if (N==2)                            % Single distance for two points
      dists = dists(1,2);
    end;

  else                                % If crds1 & crds2 provided,
    [N1,P1] = size(crds1);
    [N2,P2] = size(crds2);
    if (P1~=P2)
      error('  EUCL: sets of coordinates must be of same dimension');
    else
      P = P1;
    end;

    crds1 = crds1';                     % Transpose crds
    crds2 = crds2';

    if (N1>1 & N2>1)                    % If two matrices provided,
      dists = zeros(N1,N2);               % Calc all pairwise distances between them
      for i = 1:N1
        c1 = crds1(:,i) * ones(1,N2);
        if (P>1)
          d = sqrt(sum((c1-crds2).^2));
        else
          d = abs(c1-crds2);
        end;
        dists(i,:) = d;
      end;
    end;

    if (N1==1 & N2==1)                  % If two vectors provided,
      dists = sqrt(sum((crds1-crds2).^2));  % Calc scalar
    end;

    if (N1>1 & N2==1)                   % If matrix & reference point provided,
      crds1 = crds1 - (ones(N1,1)*crds2')'; % Center points on reference point
      if (P>1)                              % Calc euclidean distances in P-space
         dists = sqrt(sum(crds1.^2))';
      else
         dists = abs(crds1)';
      end;
    end;                                    % Return column vector

    if (N1==1 & N2>1)                   % If reference point & matrix provided,
      crds2 = crds2 - (ones(N2,1)*crds1')'; % Center points on reference point
      if (P>1)                              % Calc euclidean distances in P-space
        dists = sqrt(sum(crds2.^2));
      else
        dists = abs(crds2);
      end;
    end;                                    % Return row vector
  end;

  return;

% PUTBND: Changes the [min,max] settings for both axes to allow a buffer 
%         beyond the range of the data.  If a single matrix of point coordinates 
%         is given, columns beyond the second are ignored.
%
%     Syntax: v = putbnd(x,y,{buffer},{nocall}) 
%                   OR 
%             v = putbnd([x y],buffer,{nocall})
%
%          x =      vector of x coordinates.
%          y =      vector of y coordinates.
%          buffer = optional buffer size, as proportion of ranges 
%                     [default = 0.05].
%          nocall = optional boolean flag indicating, if true, that the axis 
%                     settings are to be returned but that the current plot is 
%                     to be left unaltered [default = 0].
%          -------------------------------------------------------------------
%          v =      axis bounds: [xmin xmax ymin ymax].
%

% See also putbnds() if change arguments.

% RE Strauss, 9/20/97
%   12/17/99 - allow for mixed input row/col vectors.
%   3/21/00 -  allow specification of spatial buffer size.
%   3/23/00 -  avoid initial call to axis().
%   5/23/01 -  add 'nocall' option.

function v = putbnd(x,y,buffer,nocall)
  if (nargin < 2) y = []; end;
  if (nargin < 3) buffer = []; end;
  if (nargin < 4) nocall = []; end;

  if (isempty(y) | isscalar(y))
    nocall = buffer;
    buffer = y;
    if (size(x,2)>=2)
      y = x(:,2);
      x = x(:,1);
    else
      error('  PUTBND: invalid point coordinates');
    end;
  end;

  x = x(:);
  y = y(:);
  if (length(x) ~= length(y))
    error('  PUTBND: lengths of coordinate vectors are incompatible.');
  end;

  if (isempty(buffer))
    buffer = 0.05;
  end;
  if (isempty(nocall))
    nocall = 0;
  end;

  % Remove NaN's
  indx = (isfinite(x) & isfinite(y));
  x = x(indx);
  y = y(indx);

  v = zeros(1,4);
  v(1) = min(x)-buffer*range(x);
  v(2) = max(x)+buffer*range(x);
  v(3) = min(y)-buffer*range(y);
  v(4) = max(y)+buffer*range(y);

  if (v(2)-v(1) < eps)      % No variation in x
    v(1) = x(1)-1;
    v(2) = x(1)+1;
  end;
  if (v(4)-v(3) < eps)      % No variation in y
    v(3) = y(1)-1;
    v(4) = y(1)+1;
  end;

  if (~nocall)
    axis(v);
  end;

  return;
  
function y = range(x)
% Calculates the range of x
y = abs(max(x) - min(x));

% Rowtoval: Converts rows of a matrix (numeric or character) to a vector of 
%           scalars.  Useful for checking for identical rows or for sorting 
%           rows into lexicological sequence [but see sortrows() for the latter].
%
%     Usage: [vals,base] = rowtoval(x,{base})
%
%           x =     numeric or character matrix.
%           base =  optional base for transformation; useful for finding values
%                     for rows not included in a previous transformation.
%           -------------------------------------------------------------------
%           vals =  column vector of representative scalar values.
%

% RE Strauss, 8/31/99
%   3/22/02 - added option of optional base.
%   5/5/06 -  convert input matrix to numeric before processing.

function [vals,base] = rowtoval(x,base)
  if (nargin < 2) base = []; end;

  x = double(x);
  [r,c] = size(x);
  xmin = min(min(x));
  x = x - xmin;
  if (isempty(base))
    base = max(max(x))+1;
  end;

  vals = zeros(r,1);
  for i = 1:c
    vals = vals + x(:,i)*base^(c-i);
  end;

  return;
  
  
  % SORTMAT:  Sorts a primary "key" column vector, and re-sequences the rows of 
%           one or more secondary matrices into the corresponding sequence.
%
%     Usage: [outvect,outmat1,...,outmat9] = sortmat(invect,inmat1,{inmat2},...,{inmat9})
%
%           invect = primary (column) vector to be sorted.
%           inmat1 = first secondary matrix to be resequenced (can be null).
%           inmat2 = optional second matrix to be resequenced (can be null).
%             .
%             .
%             .
%           inmat9 = optional ninth matrix to be resequenced (can be null).
%           --------------------------------------------------------
%           outvect = sorted primary vector.
%           outmat1 = resequenced first matrix.
%             .
%             .
%             .
%           outmat9 = resequenced ninth matrix.
%

function [outvect,outmat1,outmat2,outmat3,outmat4,outmat5,outmat6,outmat7,outmat8,outmat9] ...
            = sortmat(invect,inmat1,inmat2,inmat3,inmat4,inmat5,inmat6,inmat7,inmat8,inmat9)

  if (nargin ~= nargout)
    error('SORTMAT: number of output arguments must match number of input arguments');
  end;

  [r,c] = size(invect);
  if (c>1)
    if (r>1)
      error('SORTMAT: first input argument must be a column vector');
    else
      invect = invect';
      r = c;
    end;
  end;

  outmat1 = [];
  outmat2 = [];
  outmat3 = [];
  outmat4 = [];
  outmat5 = [];
  outmat6 = [];
  outmat7 = [];
  outmat8 = [];
  outmat9 = [];

  [outvect,seq] = sort(invect);

  if (nargin > 1)
    if (~isempty(inmat1))
      if (size(inmat1,1) ~= r)
        error('SORTMAT: <inmat1> must have same number of rows as <invect>');
      end;
      outmat1 = inmat1(seq,:);
    end;
  end;

  if (nargin > 2)
    if (~isempty(inmat2))
      if (size(inmat2,1) ~= r)
        error('SORTMAT: <inmat2> must have same number of rows as <invect>');
      end;
      outmat2 = inmat2(seq,:);
    end;
  end;

  if (nargin > 3)
    if (~isempty(inmat3))
      if (size(inmat3,1) ~= r)
        error('SORTMAT: <inmat3> must have same number of rows as <invect>');
      end;
      outmat3 = inmat3(seq,:);
    end;
  end;

  if (nargin > 4)
    if (~isempty(inmat4))
      if (size(inmat4,1) ~= r)
        error('SORTMAT: <inmat4> must have same number of rows as <invect>');
      end;
      outmat4 = inmat4(seq,:);
    end;
  end;

  if (nargin > 5)
    if (~isempty(inmat5))
      if (size(inmat5,1) ~= r)
        error('SORTMAT: <inmat5> must have same number of rows as <invect>');
      end;
      outmat5 = inmat5(seq,:);
    end;
  end;

  if (nargin > 6)
    if (~isempty(inmat6))
      if (size(inmat6,1) ~= r)
        error('SORTMAT: <inmat6> must have same number of rows as <invect>');
      end;
      outmat6 = inmat6(seq,:);
    end;
  end;

  if (nargin > 7)
    if (~isempty(inmat7))
      if (size(inmat7,1) ~= r)
        error('SORTMAT: <inmat7> must have same number of rows as <invect>');
      end;
      outmat7 = inmat7(seq,:);
    end;
  end;

  if (nargin > 8)
    if (~isempty(inmat8))
      if (size(inmat8,1) ~= r)
        error('SORTMAT: <inmat8> must have same number of rows as <invect>');
      end;
      outmat8 = inmat8(seq,:);
    end;
  end;

  if (nargin > 9)
    if (~isempty(inmat9))
      if (size(inmat9,1) ~= r)
        error('SORTMAT: <inmat9> must have same number of rows as <invect>');
      end;
      outmat9 = inmat9(seq,:);
    end;
  end;

  return;


