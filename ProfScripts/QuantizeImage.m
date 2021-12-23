% Multimedia Coding
%
%
% QuantizeImage.m:  quantize image file using a midrise uniform quantizer

clear all; close all;

% filename = 'barbara.tiff';
filename = 'DataSet\Images\barbara.tiff';

info = imfinfo(filename);
x = imread(filename);  

figure(1)
imshow(x); colormap('gray'); axis('square')

Delta = 2^4;  % quantization with quantization step  2^4
x4 = uint8(floor(single(x)/Delta)*Delta + Delta/2);
figure(2)
imshow(x4); colormap('gray'); axis('square')

Delta = 2^6;  % quantization with quantization step  2^6
x2 = uint8(floor(single(x)/Delta)*Delta + Delta/2);
figure(3)
imshow(x2); colormap('gray'); axis('square')

Delta = 2^7;  % quantization with quantization step  2^7
x1 = uint8(floor(single(x)/Delta)*Delta + Delta/2);
figure(4)
imshow(x1); colormap('gray'); axis('square')
