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
fname='L1_F0_f80_a70_t350.png';
% fname='L3_raman3_00.png';


%%% vis: publication
config_fig.units='centimeters';
config_fig.pos_2col=[0,0,17.2,6];      % 2-column
% config_fig.pos_2col=[0,0,8.6,3.2];      % 2-column
config_fig.rend='painters';
config_fig.ax_fontsize=9;
config_fig.ax_lwid=0.5;
config_fig.mark_siz=4;
config_fig.line_wid=0.5;


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
xlabel('x (mm)');
ylabel('y (mm)');


%% filter
Ifilt=imgaussfilt(Ibeam,gfilt_sig);

% display beam profile
h_filt=figure('Name','beam profile (filtered)');
s=surf(1e3*X,1e3*Y,Ifilt,'EdgeColor','none','FaceColor','interp');
cbar=colorbar;
cbar.Title.String='Intensity (a.u.)';
axis tight;
view(2);
xlabel('x (mm)');
ylabel('y (mm)');


%% approximate beam profile
% centre
x0_approx=sum(X(:).*Ifilt(:))/sum(Ifilt(:));
y0_approx=sum(Y(:).*Ifilt(:))/sum(Ifilt(:));

% rms width
idx_x0=round(x0_approx/pixsize);
idx_y0=round(y0_approx/pixsize);
sigx_approx=sqrt(sum((X(:,idx_y0).^2).*Ifilt(:,idx_y0))/sum(Ifilt(:,idx_y0))-x0_approx^2);
sigy_approx=sqrt(sum((Y(idx_x0,:).^2).*Ifilt(idx_x0,:))/sum(Ifilt(idx_x0,:))-y0_approx^2);

amp_approx=Ifilt(idx_x0,idx_y0);


%% fit beam profile
% 2D Gaussian with rotated axis
Z=cat(3,X,Y);       % format indep data array
p0=[amp_approx,x0_approx,sigx_approx,y0_approx,sigy_approx,0,0];

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
xlabel('x (mm)');
ylabel('y (mm)');


%% summarise beam profile
beamprof.amp=p_fit(1);                   % normalised amplitude
beamprof.x0=[p_fit(2),p_fit(4)];         % beam centre
beamprof.w=2*abs([p_fit(3),p_fit(5)]);        % beam width (e^-2 radius of Intensity)
beamprof.theta=p_fit(6);                 % principal axes tilt wrt camera
beamprof.c=p_fit(7);                     % constant background intensity

printstruct(beamprof);

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
% 
% 
% % display
% fprintf('%0.3e\n',x0);


%% VIS: 
figname='beam_profile';
h=figure('Name',figname,'Units',config_fig.units,'Position',config_fig.pos_2col,...
    'Renderer','opengl');

% raw
subplot(1,2,1);

s=surf(1e3*X,1e3*Y,Ibeam,'EdgeColor','none','FaceColor','interp');
colormap('gray');
cbar=colorbar;
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Title.String='Intensity (a.u.)';
cbar.Label.FontSize=config_fig.ax_fontsize;
cbar.FontSize=config_fig.ax_fontsize;

ax=gca;
axis tight;
set(ax,'Layer','Top');
grid off;
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

view(2);
title('raw');
xlabel('x (mm)');
ylabel('y (mm)');


% fitted
subplot(1,2,2);

s=surf(1e3*Xfit,1e3*Yfit,Ifit,'EdgeColor','none','FaceColor','interp');
colormap('gray');
cbar=colorbar;
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Title.String='Intensity (a.u.)';
cbar.Label.FontSize=config_fig.ax_fontsize;
cbar.FontSize=config_fig.ax_fontsize;

ax=gca;
axis tight;
set(ax,'Layer','Top');
grid off;
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

view(2);
title('fitted');
xlabel('x (mm)');
ylabel('y (mm)');

%%% SAVE
% [~,fname_main,~]=fileparts(fname)
% vecrast(h,strjoin({get(gcf,'name'),fname_main},'_'),300,'top','pdf')