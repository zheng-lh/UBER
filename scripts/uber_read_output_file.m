function ostr = uber_read_output_file(fname)
% function ostr = uber_read_output_file(fname)
% this function reads the UBER solution output file and stores its content into a
% structure.
%
% input:
%    fname (file name)
%
% output:
%    ostr, contains:
%    .fname
%    .mode (input mode, = 1 or 2)
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
%    $                     .sol                      $
%    $                     .err                      $
%
%    the sizes of the solution array sol and error array err depend on the input
%    mode and dimensions. in mode 1, they are [np, 2]. in mode 2, they are
%    [size_of_space, nt], where size_of_space is [nx1], [nx1, nx2], or [nx1, nx2,
%    nx3] depending on ndim.

fid = fopen(fname, 'r');
ostr = [];
ostr.fname = fname;
io_mode = fread(fid, 1, 'integer*4');
ostr.mode = io_mode;

switch io_mode
   case 1
      np = fread(fid, 1, 'integer*4');
      xx1 = fread(fid, [1, np], 'double');
      xx2 = fread(fid, [1, np], 'double');
      xx3 = fread(fid, [1, np], 'double');
      tt = fread(fid, [1, np], 'double');
      sol = fread(fid, [1, 2*np], 'double');
      err = fread(fid, [1, 2*np], 'double');

      sol = reshape(sol, [np, 2]);
      err = reshape(err, [np, 2]);

      ostr.np = np;
      ostr.xx1 = xx1;
      ostr.xx2 = xx2;
      ostr.xx3 = xx3;
      ostr.tt = tt;
      ostr.sol = sol;
      ostr.err = err;
   case 2
      ndim = fread(fid, 1, 'integer*4');
      ostr.ndim = ndim;
      nx1 = 1;
      nx2 = 1;
      nx3 = 1;
      switch ndim
         case 1
            nx1 = fread(fid, 1, 'integer*4');
            ostr.nx1 = nx1;
         case 2
            nx1 = fread(fid, 1, 'integer*4');
            nx2 = fread(fid, 1, 'integer*4');
            ostr.nx1 = nx1;
            ostr.nx2 = nx2;
         case 3
            nx1 = fread(fid, 1, 'integer*4');
            nx2 = fread(fid, 1, 'integer*4');
            nx3 = fread(fid, 1, 'integer*4');
            ostr.nx1 = nx1;
            ostr.nx2 = nx2;
            ostr.nx3 = nx3;
         otherwise
            error('function uber_read_output_file: wrong ndim value (%i)', ndim);
      end
      nt = fread(fid, 1, 'integer*4');
      ostr.nt = nt;
      switch ndim
         case 1
            xx1 = fread(fid, [1, nx1], 'double');
            ostr.xx1 = xx1;
         case 2
            xx1 = fread(fid, [1, nx1], 'double');
            xx2 = fread(fid, [1, nx2], 'double');
            ostr.xx1 = xx1;
            ostr.xx2 = xx2;
         case 3
            xx1 = fread(fid, [1, nx1], 'double');
            xx2 = fread(fid, [1, nx2], 'double');
            xx3 = fread(fid, [1, nx3], 'double');
            ostr.xx1 = xx1;
            ostr.xx2 = xx2;
            ostr.xx3 = xx3;
      end
      tt = fread(fid, [1, nt], 'double');
      ostr.tt = tt;
      np = nx1*nx2*nx3;
      sol = fread(fid, [1, np*nt], 'double');
      err = fread(fid, [1, np*nt], 'double');

      switch ndim
         case 1
            sol = reshape(sol, [nx1, nt]);
            err = reshape(err, [nx1, nt]);
         case 2
            sol = reshape(sol, [nx1, nx2, nt]);
            err = reshape(err, [nx1, nx2, nt]);
         case 3
            sol = reshape(sol, [nx1, nx2, nx3, nt]);
            err = reshape(err, [nx1, nx2, nx3, nt]);
      end

      ostr.sol = sol;
      ostr.err = err;
   otherwise
      error('function uber_read_output_file: wrong io_mode value (%i)', io_mode);
end

fclose(fid);

end% function uber_read_output_file



            





