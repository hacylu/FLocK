function [stain_img] = color_deconvolution(img,stain_type,input_stain_mtx)

% if ~strcmp(class(img),'uint8')
%     error('IMG must be uint8.');
% end
img = double(img);

switch stain_type
    case {'HE','H&E'}
        stain_mtx = [0.644211 0.716556 0.266844; 
                     0.092789 0.954111 0.283111;
                     0        0        0]';
    case {'RGB'}   
%         stain_mtx = eye(3);  
        stain_mtx = [0 1 1; 
                     1 0 1;
                     1 1 0]';             
    case {'MTX','mtx'}
        stain_mtx = input_stain_mtx;
    otherwise 
        error('Unsupported STAIN_TYPE.');
end

stain_mtx = mult_col(stain_mtx,1./nonzero(sqrt(sum(stain_mtx.^2))));

if all(stain_mtx(:,2)==0)
    stain_mtx(:,2) =  stain_mtx([3 1 2],1);
end

if all(stain_mtx(:,3)==0)
    stain_mtx(:,3) = sqrt(1 - min(sum(stain_mtx.^2,2),1));
    stain_mtx(:,3) = stain_mtx(:,3)./sqrt(sum(stain_mtx(:,3).^2));
end

inv_stain_mtx = inv(stain_mtx);

A = -log(nonzero(img)/255);
A = reshape(A,size(img,1)*size(img,2),3)';
C = stain_mtx\A;          %=> inv(stain_mtx)*T;
C = reshape(C',size(img));
stain_img = uint8(min(round(255*exp(-C)),255));





