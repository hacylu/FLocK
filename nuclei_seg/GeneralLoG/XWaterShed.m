function [bnf,blm] = XWaterShed(bw,bm)

% Input:
%    -bw the binary image
%    -bm the binary marker image
% Output:
%    -bwSplitted the result
% Program written by Hongming Xu
% ECE, University of Alberta
% Modified Marker-Controlled Watershed Algorithm
bw=imdilate(bw,strel('disk',1));
filled=imfill(bw,'holes');
holes=filled&~bw;
bigholes=bwareaopen(holes,10);
smallholes=holes&~bigholes;
bw=bw|smallholes;

CC = bwconncomp(bw);
LL = labelmatrix(CC); %% assume that background is labelled as 0
[m,n] = size(bw);
bw_dist = zeros(m,n);
bm_temp=logical(zeros(m,n));
bw_temp=logical(zeros(m,n));
bw_dist0=-bwdist(~bw);
tic
for obj = 1:CC.NumObjects
    AA = zeros(m,n);
    AA(find(LL == obj)) = 1;
    BB = AA & bm;
    S = sum(BB(:));
    if S > 1
         bw_temp=bw_temp|AA;
         bm_temp=bm_temp|BB;
%         bw_dist1 = -bwdist(~AA);
%        [r,c]=find(bw_dist1~=0);

        [row,col]=find(BB);
        [r,c]=find(AA);
        X=[r c];Y=[row col];
        D=sqrt(sum(abs(repmat(permute(X,[1 3 2]),[1 size(Y,1) 1])...
            -repmat(permute(Y,[3 1 2]),[size(X,1) 1 1])).^2,3));
        ind1=sub2ind([m,n],Y(:,1),Y(:,2));
%        d=bw_dist1(ind1);
        d=bw_dist0(ind1);
        temp=repmat(d',size(D,1),1);
        Dmin=min(D+temp,[],2);
        ind2=sub2ind([m,n],X(:,1),X(:,2));
%        bw_dist1(ind2)=Dmin;
         bw_dist(ind2)=Dmin;
         
%% original method         
%         X = [row,col]; 
%         ind1=sub2ind(size(bm),X(:,1),X(:,2));
%         [r,c] = find(bw_dist1 ~= 0);
%         for k = 1:length(r)
%             i = r(k);
%             j = c(k);
%             y = [i j];
%             d = sqrt(sum(abs(X - repmat(y,size(X,1),1)).^2,2));
%             s_value = d + bw_dist1(ind1); 
%             bw_dist1(i,j) = min(s_value);
%         end
%% end original

%       bw_dist = bw_dist + bw_dist1;    
    end
    
end
toc
bw_dist_minImposed = imimposemin(bw_dist,bm_temp);
bw_dist_minImposed(~bw_temp) = -Inf;
L = watershed(bw_dist_minImposed,8);
stats=regionprops(L,'Area');
area=[stats(:).Area];
ind=find(max(area));
bn=bw_temp;
indm=find(L==ind | L==0);
bn(indm)=0;
bnf=bn | (bw-bw_temp);

blm=bwperim(bw-bw_temp);
ind3=find(L==0);
blm(ind3)=1;
end

