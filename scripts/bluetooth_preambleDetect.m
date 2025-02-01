clear
clc

% Sample Rate from USRP should be 4 MHz. We have 13 octets of preamble and
% which comes out to 104 bits and matched filtered to 104*10 (size of match
% h)and then sized to 104*10*4 for a 4 MHz sampling rate. 
a_old = read_complex_binary('bluetooth_preamble.dat');
a_new = read_complex_binary('bluetooth_preamble_change.dat');
a_new_short = a_new(3e6+1:4e6);
load('new_preamble_bluetooth.mat');
load('preamble_bluetooth.mat');

plot(abs(xcorr(a_new_short, new_preamble_better_bluetooth)));