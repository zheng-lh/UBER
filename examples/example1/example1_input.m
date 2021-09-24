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
istr.fname = 'example1_in.dat';
istr.mode = 2;
istr.ndim = 1;
istr.nx1 = 25;
istr.nt = 2;
istr.xx1 = linspace(2.d-2, 0.98d0, istr.nx1);
istr.tt = [5.d-2, 2.d-1];

uber_create_input_file(istr);



