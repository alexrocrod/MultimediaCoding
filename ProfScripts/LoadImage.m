% Multimedia Coding
%
%
% LoadImage.m:  load image file 

clear all; close all;


clear all; close all;

% filename = 'barbara.tiff';
filename = 'DataSet\Images\barbara.tiff';

info = imfinfo(filename);
x = imread(filename);  

figure(1)
imshow(x); 

fprintf('\n'); 
fprintf('Format:          %s',info.Format); fprintf('\n'); 
fprintf('Size:            Width = %d,  Height = %d',info.Width,info.Height); fprintf('     [pixels] \n');

