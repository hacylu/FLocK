% Oct. 2, 2015. Cheng LU  speeded up the merging procedure which will takes 
% more than a whole day for a large image( 10k by 10k and up)
function idx = mergeContours_s(C, fitness, Th, centroids)

DTHRESH = 100;
DTHRESH2 = DTHRESH^2;

P = polygonClipperConvert(C);

s = length(P);

disp(' - Allocate memory for calculating distance matrix.');
disp(' - Calculating distance matrix.');

%% original method
% % distance (overlap) matrix
% % allocate memory for 1% of the size
% display('Original speed');
% tic
% D = spalloc(s, s, round(s^2/100));%D(,k,l) store the overlap percentage
% 
% for k = 1:s
%     
%     Pk = P(k);
%     
%     centroid = centroids(k,:);
%     
%     
%     
%     for l = k+1:s
%         
%         d = centroid  - centroids(l,:);
%         
%         % centroids too far, intersection improbable (speed improvement)
%         if d*d' > DTHRESH2
%             continue;
%         end
%         
%         Pl = P(l);
%         
%         intersect = polyclip(Pk, Pl, 1);
%         
%         if isempty(intersect)
%             continue;
%         end
%         
%         areaIntersect = sum(polygonClipperArea(intersect));
%         
%         overlap = max(...
%             areaIntersect / polygonClipperArea(Pk), ...
%             areaIntersect / polygonClipperArea(Pl)...
%             );
%         
%         D(k,l) = overlap;
%         
%     end
%     
% end
% 
% D = D + D';
% 
% % adjacency matrix
% A = D > Th;
% toc
%% speed up method
% distance (overlap) matrix
% display('New speed');
tic
Dn = spalloc(s, s, round(s^2/100));%D(,k,l) store the overlap percentage
% Dn=sparse(s,s);

% cal overlap matrix
for k = 1:s
    Pk = P(k);
    centroid = centroids(k,:);
    
    dC=bsxfun(@minus,centroid,centroids);
    dCM=bsxfun(@times,dC,dC);
    d2_k=sum(dCM,2);

%     d2_k=diag(dC*dC'); % take lots of time
    
    temp=(d2_k<DTHRESH2)';
    temp(1)=0;
    
    idxCur=find(temp);
    for l=1:length(idxCur) 
        
        Pl = P(idxCur(l));
        intersect = polyclip(Pk, Pl, 1);
        
        if ~isempty(intersect)
            areaIntersect = sum(polygonClipperArea(intersect));
            
            overlap = max(...
                areaIntersect / polygonClipperArea(Pk), ...
                areaIntersect / polygonClipperArea(Pl)...
                );
            
            Dn(k,idxCur(l)) = overlap;
        end
        

    end
end
% adjacency matrix
An = Dn > Th;
toc
%% resolve overlaps
list = 1:s;
c = 0;

disp(' - Resolving concurrence.');

while ~isempty(list)
    
    % find fittest contour
    [~, idx_mx] = max(fitness(list));
    
    % add it to the list of accepted contours
    c = c + 1;
    idx(c) = list(idx_mx);
    
    % remove it from the waiting list
    list(idx_mx) = [];
    
    % remove all adjacent contours to the fittest
    remove = false(size(list));
    
    for l = 1:length(list)
        remove(l) = An(idx(c), list(l));
    end
    
    list = list(~remove);
    
end

end