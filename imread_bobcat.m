function [A,X,Y,err] = imread_bobcat(filename,verbose)
% Reads image from Bobcat CCD camera with pixel locations
%   * Infrared imaging camera
% DKS
% 2018-05-26


% check input vars
if ~exist('verbose','var')
    verbose=0;
end


% config camera config
pixsize=20e-6;      % camera pixel pitch [m]


% load raw image
A=imread(filename);
A=double(A);

% error check
err=0;      % no err

if max(A(:))==2^16-1
    if verbose>0
        warning('Image is saturated.');
    end
    err=1;      % saturated image
end


% pixel positions
npixels=size(A);
x=pixsize*(1:npixels(1));
y=pixsize*(1:npixels(2));
[X,Y]=ndgrid(x,y);


end