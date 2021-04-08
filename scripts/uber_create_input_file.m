function fid = uber_create_input_file(istr)
% function fid = uber_create_input_file(istr)
% this function creates the solution input file for UBER from a user specified
% structure.
%
% input:
%    istr, contains:
%    .fname (file name)
%    .mode  (input mode, = 1 or 2)
%    $    (1)    $                (2)                $
%    $    .np    $        .ndim (= 1, 2, or 3)       $
%    $           $    (1)    |    (2)    |    (3)    $
%    $           $   .nx1    |   .nx1    |   .nx1    $
%    $           $           |   .nx2    |   .nx2    $
%    $           $           |           |   .nx3    $
%    $           $               .nt                 $
%    $ .xx1(np)  $ .xx1(nx1) | .xx1(nx1) | .xx1(nx1) $
%    $ .xx2(np)  $           | .xx2(nx2) | .xx2(nx2) $
%    $ .xx3(np)  $           |           | .xx3(nx3) $
%    $  .tt(np)  $              .tt(nt)              $
%
% output:
%    fid (file id)

fid = fopen(istr.fname, 'w');
fwrite(fid, istr.mode, 'integer*4');

switch istr.mode
   case 1
      fwrite(fid, istr.np, 'integer*4');
      fwrite(fid, istr.xx1, 'double');
      fwrite(fid, istr.xx2, 'double');
      fwrite(fid, istr.xx3, 'double');
      fwrite(fid, istr.tt, 'double');
   case 2
      fwrite(fid, istr.ndim, 'integer*4');
      switch istr.ndim
         case 1
            fwrite(fid, istr.nx1, 'integer*4');
         case 2
            fwrite(fid, istr.nx1, 'integer*4');
            fwrite(fid, istr.nx2, 'integer*4');
         case 3
            fwrite(fid, istr.nx1, 'integer*4');
            fwrite(fid, istr.nx2, 'integer*4');
            fwrite(fid, istr.nx3, 'integer*4');
         otherwise
            error('function uber_create_input_file: wrong istr.ndim (%i)', istr.ndim);
      end
      fwrite(fid, istr.nt, 'integer*4');
      switch istr.ndim
         case 1
            fwrite(fid, istr.xx1, 'double');
         case 2
            fwrite(fid, istr.xx1, 'double');
            fwrite(fid, istr.xx2, 'double');
         case 3
            fwrite(fid, istr.xx1, 'double');
            fwrite(fid, istr.xx2, 'double');
            fwrite(fid, istr.xx3, 'double');
      end
      fwrite(fid, istr.tt, 'double');
   otherwise
      error('function uber_create_input_file: wrong istr.mode (%i)', istr.mode);
end

fclose(fid);

end% function uber_create_input_file


