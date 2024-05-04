clear all, close all; clc;

u0     = im2double(imread('texmos3.s512.tiff'))*1;
[m,n] = size(u0);
[x,y] = meshgrid(1:n,1:m);
v0     = zeros(m,n);
a = 0.4;
v0(1:m/2,1:n/2)     = a*cos(2*pi*128/m*x(1:m/2,1:n/2)).*cos(2*pi*128/n*y(1:m/2,1:n/2));%sum(v0(:));
v0(m/2+1:end,1:n/2) = a*cos(2*pi*64/m*x(m/2+1:end,1:n/2));%sum(v0(:));
v0(1:m/2,n/2+1:end) = a*cos(2*pi*64*(x(1:m/2,n/2+1:end)/m+y(1:m/2,n/2+1:end)/n)) ;%sum(v0(:));
v0(m/2+1:end,m/2+1:end) = a*cos((2*pi*128)/m*y(m/2+1:end,1:n/2));

Im     = u0+v0;

sigma  = 3; 
lambda = 1e4;
Param.Reglambda = lambda;
Param.Sigma     = sigma;
Im              = im2double(Im);

tic; [uu,OutPut] = ImSmoothL0TVQP(Im, Param); t=toc;
figure(90); imshow(uu);
figure(91); imshow((Im-uu)+0.5);
