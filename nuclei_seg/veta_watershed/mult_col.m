function y = mult_col(mtx,vec)
% Y = MULT_COL(MTX,VEC)
%   MULT_COL multiplies the column vector VEC with each column of the
%   matrix MTX, producing Y.
%
% MTX: (MxN)
% VEC: (Mx1)
% Y:   (MxN)
%
% See also add_col, add_row, mult_row
%y = mult_col_c(mtx,vec);
y = bsxfun(@times, mtx, vec);