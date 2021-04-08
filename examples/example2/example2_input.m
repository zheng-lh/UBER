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



