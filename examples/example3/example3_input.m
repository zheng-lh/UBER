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

clear all;

Po180 = pi/180.d0;
aa0 = [30.5:1:89.5]*Po180;
na0 = numel(aa0);
pphi = [-360.d0+1.5d0:3.d0:360.d0-1.5d0]*Po180;
nph = numel(pphi);
delta_t = 7200.d0;

istr = [];
istr.fname = 'example3_in.dat';
istr.mode = 2;
istr.ndim = 2;
istr.nx1 = na0;
istr.nx2 = nph;
istr.nt = 2;
istr.xx1 = aa0;
istr.xx2 = pphi;
istr.tt = [1:istr.nt]*delta_t;

uber_create_input_file(istr);


