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



