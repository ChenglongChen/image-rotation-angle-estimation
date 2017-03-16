function subimage=CentralCrop(img,B) 
% -------------------------------------------------------------------------
%   Written Feb 10, 2012 by Chenglong Chen
%   Copyright 2012 by Chenglong Chen
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The end-user understands that the program was
% developed for research purposes and is advised not to rely exclusively on
% the program for any reason.
% -------------------------------------------------------------------------
% Contact: c.chenglong@gmail.com
% -------------------------------------------------------------------------
% Input:       img .... input image (one or more channel)
%                B .... size of the center portion to be cropped, must be
%                       even
% Output: subimage .... the center portion
% -------------------------------------------------------------------------
% A central cropping operation of size B x B applied to the image.
% -------------------------------------------------------------------------
 

% claculate the coordinates of the center portion
[row,col] = size(img(:,:,1));
i = floor(row/2)-floor(B/2);
j = floor(col/2)-floor(B/2);
X = j+1:j+B;
Y = i+1:i+B;

% crop the center portion
if(length(size(img))==3)
    subimage = img(Y,X,:);
else
    subimage = img(Y,X);
end

end