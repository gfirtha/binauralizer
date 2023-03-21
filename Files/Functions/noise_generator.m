clear
close all

fs = 48e3;
N = 5*fs;
out = pinknoise(N);
audiowrite('noise.wav',out,fs)