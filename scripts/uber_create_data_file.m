function fid = uber_create_data_file(istr)
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Copyright 2021, Liheng Zheng
%
% This file is part of UBER.
%
%    UBER is free software: you can redistribute it and/or modify it under the
%    terms of the MIT License as published by Massachusetts Institute of
%    Technology. UBER is distributed in the hope that it will be useful, but
%    WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
%    FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
%
%    You should have received a copy of the MIT License along with UBER. If not,
%    see <https://opensource.org/licenses/MIT>.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%
% function fid = uber_create_data_file(istr)
% this function creates the gridded data file for UBER from a user specified
% structure.
%
% input:
%    istr, contains:
%    .fname (file name)
%    .ndim (= 1, 2, or 3)
%    $    (1)    |    (2)    |    (3)    $
%    $   .nx1    |   .nx1    |   .nx1    $
%    $           |   .nx2    |   .nx2    $
%    $           |           |   .nx3    $
%    $ .xx1(nx1) | .xx1(nx1) | .xx1(nx1) $
%    $           | .xx2(nx2) | .xx2(nx2) $
%    $           |           | .xx3(nx3) $
%    .darr
%
%    the size of data array darr is [nx1], [nx1, nx2], or [nx1, nx2, nx3] depending
%    on ndim. N.B., both MATLAB and FORTRAN arrays are column major.
%
% output:
%    fid (file id)

fid = fopen(istr.fname, 'w');
fwrite(fid, istr.ndim, 'integer*4');

switch istr.ndim
   case 1
      fwrite(fid, istr.nx1, 'integer*4');
      fwrite(fid, istr.xx1, 'double');
   case 2
      fwrite(fid, istr.nx1, 'integer*4');
      fwrite(fid, istr.nx2, 'integer*4');
      fwrite(fid, istr.xx1, 'double');
      fwrite(fid, istr.xx2, 'double');
   case 3
      fwrite(fid, istr.nx1, 'integer*4');
      fwrite(fid, istr.nx2, 'integer*4');
      fwrite(fid, istr.nx3, 'integer*4');
      fwrite(fid, istr.xx1, 'double');
      fwrite(fid, istr.xx2, 'double');
      fwrite(fid, istr.xx3, 'double');
   otherwise
      error('function uber_create_data_file: wrong istr.dim (%i)', istr.dim);
end

fwrite(fid, istr.darr, 'double');

fclose(fid);

end% function uber_create_data_file



