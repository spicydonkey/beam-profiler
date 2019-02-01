function F = gauss2rot(p,x)
%% Gaussian function parameterisation
% 2D (bi-variate) Normal/Gaussian function with rotated principal axes with
% offset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter vector: p = [amp, mu_x, sig_x, mu_y, sig_y, theta, offset]
%
%
% DKS

%% Evaluate gaussian function
% transform to rotated coord system
%   Mrot = [cos(phi) -sin(phi); sin(phi) cos(phi)]
X(:,:,1)= (x(:,:,1)-p(2))*cos(p(6)) - (x(:,:,2)-p(4))*sin(p(6));
X(:,:,2)= (x(:,:,1)-p(2))*sin(p(6)) + (x(:,:,2)-p(4))*cos(p(6));

F=p(1)*exp(-( X(:,:,1).^2/(2*p(3)^2) + X(:,:,2).^2/(2*p(5)^2))) + p(7);

% figure(3)
% alpha(0)
% imagesc(F)
% colormap('gray')
% figure(gcf)%bring current figure to front
% drawnow
% beep
% pause %Wait for keystroke
