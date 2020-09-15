function y = mult_row(mtx,vec)
% Y = MULT_ROW(MTX,VEC)
%   MULT_ROW multiplies the row vector VEC with each row of the matrix MTX,
%   producing Y.
%
% MTX: (MxN)
% VEC: (1xN)
% Y:   (MxN)
%
% See also add_col, add_row, mult_col
y = mult_col(mtx',vec')';