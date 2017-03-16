%--------------------------------------------------------------------------
%  Examples of usage of the method proposed in:
%
%  [1] C. Chen, J. Ni, and Z. Shen, "Effective estimation of image rotation
% angle using spectral method," Signal Processing Letters, IEEE, vol. 21,
% no. 7, pp. 890--894, July 2014.
%
% The 'Lena' image used in this example can be downloaded from the USC-SIPI
% image database on:
% http://sipi.usc.edu/database/download.php?vol=misc&img=4.2.04
% For copyright information, please go to:
% http://sipi.usc.edu/database/copyright.php
%
%--------------------------------------------------------------------------
% This code is provided only for research purposes.
%--------------------------------------------------------------------------


% Clear all variables and close all figures
clear all; close all;

% add current path to the search path
addpath(genpath(cd));
% Read the image
RGB = imread('Lena.tiff');
% Take the Green channel
img = RGB(:,:,2);
% parameter settings
angle = 15; kernel = 'bicubic';
fprintf('true angle: %d\n', angle);
% rotate image
img_r = imrotate(img,angle,kernel);
% compute Txx
block_size = 128;
Txx = CalculateTxx(img_r,block_size);
% estimate rotation angle
estimatedAngle = EstimateAngle(Txx,0);
% the proposed method
fprintf('angle estimated using the proposed method: %f\n', estimatedAngle(1));
% Padin et al.'s method
fprintf('angle estimated using Padin et al.''s method: %f\n', estimatedAngle(4));
