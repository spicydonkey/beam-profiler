% demo script for beam profile analysis from camera image
%
% DK Shin
% 20180218
%

clc;
close all;

%% config
% camera
pixsize=20e-6;      % camera pixel pitch [m]

% filter
gfilt_sig=3;        % gaussian width

% data file
fname='F0_f79+81_a70.png';
% fname='f80.png';


%% load raw image
Iraw=imread(fname);

if max(Iraw(:))==2^16-1
    warning('Image is saturated.');
end

% beam intensity - scale pixel values to [0,1]
Ibeam=double(Iraw)/double(max(Iraw(:)));

% get image axis
npixels=size(Ibeam);
x=pixsize*(1:npixels(1));
y=pixsize*(1:npixels(2));
[X,Y]=ndgrid(x,y);

% display beam profile
h_raw=figure('Name','beam profile');
s=surf(1e3*X,1e3*Y,Ibeam,'EdgeColor','none','FaceColor','interp');
cbar=colorbar;
cbar.Title.String='Intensity (a.u.)';
axis tight;
view(2);
xlabel('x [mm]');
ylabel('y [mm]');


%% filter
Ifilt=imgaussfilt(Ibeam,gfilt_sig);

% display beam profile
h_filt=figure('Name','beam profile (filtered)');
s=surf(1e3*X,1e3*Y,Ifilt,'EdgeColor','none','FaceColor','interp');
cbar=colorbar;
cbar.Title.String='Intensity (a.u.)';
axis tight;
view(2);
xlabel('x [mm]');
ylabel('y [mm]');


%% approximate beam profile
% centre
x0_approx=sum(X(:).*Ifilt(:))/sum(Ifilt(:));
y0_approx=sum(Y(:).*Ifilt(:))/sum(Ifilt(:));

% rms width
idx_x0=round(x0_approx/pixsize);
idx_y0=round(y0_approx/pixsize);
wx_approx=sqrt(sum((X(:,idx_y0).^2).*Ifilt(:,idx_y0))/sum(Ifilt(:,idx_y0))-x0_approx^2);
wy_approx=sqrt(sum((Y(idx_x0,:).^2).*Ifilt(idx_x0,:))/sum(Ifilt(idx_x0,:))-y0_approx^2);

amp_approx=Ifilt(idx_x0,idx_y0);


%% fit beam profile
% 2D Gaussian with rotated axis
Z=cat(3,X,Y);       % format indep data array
p0=[amp_approx,x0_approx,wx_approx,y0_approx,wy_approx,0,0];

[p_fit,resnorm,residual,exitflag] = lsqcurvefit(@gauss2rot,p0,Z,Ifilt);

% visualise fit
xfit=linspace(min(x),max(x));
yfit=linspace(min(y),max(y));
[Xfit,Yfit]=ndgrid(xfit,yfit);
Zfit=cat(3,Xfit,Yfit);

Ifit=gauss2rot(p_fit,Zfit);

% display beam profile
h_fit=figure('Name','fitted beam profile');
s=surf(1e3*Xfit,1e3*Yfit,Ifit,'EdgeColor','none','FaceColor','interp');
cbar=colorbar;
cbar.Title.String='Intensity (a.u.)';
axis tight;
view(2);
xlabel('x [mm]');
ylabel('y [mm]');


%% summarise beam profile
amp=p_fit(1);
x0=[p_fit(2),p_fit(4)];
w=[p_fit(3),p_fit(5)]
theta=p_fit(6);
c=p_fit(7);


% %% display fitted cross-section
% % indicate major-minor axis
% mm=tan(theta);      % gradient
% cc=x0(2)-mm*x0(1);  % y-intercept
% 
% xq=linspace(min(X(:)),max(X(:)));
% yq=x0(2);
% vq=interpn(X,Y,Ifilt,xq,yq,'linear');
% 
% figure;
% plot(xq,vq);


% display
fprintf('%0.3e\n',x0);
