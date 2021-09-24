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

istr = [];
istr.fname = 'example2_in.dat';
istr.mode = 2;
istr.ndim = 1;
istr.nx1 = 101;
istr.nt = 20;
istr.xx1 = linspace(0.d0, 1.d0, istr.nx1);
tt = linspace(0.d0, 5.d-2, istr.nt+1);
istr.tt = tt(2:end);

uber_create_input_file(istr);



