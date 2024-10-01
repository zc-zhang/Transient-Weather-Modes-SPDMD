% test Residual DMD for SCALE data simulation;
% %% Ref Residual-Dynamic-Mode-Decomposition/Examples_gallery_4/Cylinder_example
% Ref: Another look at residual DMD
%close all

%load('Vfull.mat')
% load('EnsVfull') with 
% d=3880 (dimension), M1 =11 (Ensemble initial), M2-1=120 (snapshots)


%DATA =Vfull;
DATA =Tdata;
M1 = 1;                % the number of 
M2 = size(DATA,2)-1;  % the number of snapshots 
M= M1*M2;
dim = size(DATA,1);    % the size of the state dimension
N = 120;              % the number of basis functions  N<=M1*M2
ind1=1:M2; 
ind2=1:M2;
Xa = DATA(:,ind1);
Ya = DATA(:,ind1+1);

% Kernel ResDMD
[G,K,L,PX,PY] = kernel_ResDMD(Xa,Ya,'N',N,'type',"Gaussian");
%--------------------------------
[V,LAM,W2] = eig(K,'vector');      
 
 E=diag(LAM); 
R = abs(sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(LAM).^2)));
[~,Index_res]=sort(R,'descend');



figure
scatter(real(LAM),imag(LAM),300,R,'.','LineWidth',1);
hold on
% scatter(real(LAM),imag(LAM),250,R,'.');
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
% load('cmap.mat')
% colormap(cmap2); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title(sprintf('Residuals ($N=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%exportgraphics(gcf,sprintf('SCALE_res_N%d.jpg',N),'ContentType','vector','BackgroundColor','none')

figure
loglog([0.001,1],[0.001,1],'k','linewidth',2)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'b.','markersize',20)
xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
ylabel('residual','interpreter','latex','fontsize',18)
title(sprintf('Residuals ($N=%d$)',N),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%exportgraphics(gcf,sprintf('SCALE_res2_N%d.jpg',N),'ContentType','vector','BackgroundColor','none')

%%


x_pts=-1.2:0.02:1.2;    y_pts=-0.02:0.02:1.2;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES0 = KoopPseudoSpecQR(PX,PY,1/M,z_pts);
RES0=reshape(RES0,length(y_pts),length(x_pts));

RES = KoopPseudoSpec(double(G),double(K),double(L),z_pts,'Parallel','off');	% compute pseudospectra
%RES = KoopPseudoSpec(G,K,L,z_pts,'Parallel','off');	% compute pseudospectra
RES=reshape(RES,length(y_pts),length(x_pts));

%% Plot pseudospectra

figure
hold on
v=(10.^(-10:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES0)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES0)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
axis equal;

title(sprintf('Naive Residual ($N=%d$)',N),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12);
box on
%exportgraphics(gcf,sprintf('SCALE_pseudoW_N%d.jpg',N),'ContentType','vector','BackgroundColor','none')




figure
hold on
v=(10.^(-10:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
axis equal;

title(sprintf('Pseudospectrum ($M=%d$)',M),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12);
box on
%exportgraphics(gcf,sprintf('SCALE_pseudo_N%d.jpg',M),'ContentType','vector','BackgroundColor','none')

%% Plot Koopman Modes 
% %Koopman mode decomposition
sigu=DATA(1:3880,ind1);       sigu=sigu(:,1:M2);
sigv=DATA(3881:end,ind1);     sigv=sigv(:,1:M2);
sig_vel=sqrt(sigu.^2+sigv.^2);

I_p=find(R<0.5);
K_modes_res_u=pinv(PX(1:M2,1:N)*V(:,I_p))*transpose(sigu);
K_modes_res_v=pinv(PX(1:M2,1:N)*V(:,I_p))*transpose(sigv);
K_modes_res_vel=pinv(PX(1:M2,1:N)*V(:,I_p))*transpose(sig_vel);

K_modes_u=pinv(PX(1:M2,1:N)*V)*transpose(sigu);  
K_modes_v=pinv(PX(1:M2,1:N)*V)*transpose(sigv);

% Plot specific modes
th=0.01;0.05;
close all

figure;
h=tiledlayout(1,3);

lam=0.9439+0.2458i; % spectral parameter
lam=0.888040789931911 + 0.213101183194327i;
mode=find(abs(E-lam)==min(abs(E-lam)));
lam=E(mode);
clc
R(mode)
tt=norm(PSI_x(1:M2,1:N)*V(:,I_p==mode))/sqrt(M2);

XIu=reshape((K_modes_res_u(I_p==mode,:)),[40,97])*tt;
myContours = linspace(min(real(XIu(:))),max(real(XIu(:))), 21);
nexttile
contourf((real(XIu)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))

XIv=reshape((K_modes_res_v(I_p==mode,:)),[40,97])*tt;
myContours = linspace(min(real(XIv(:))),max(real(XIv(:))), 21);
nexttile
contourf((real(XIv)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))

XIvel=reshape((K_modes_res_vel(I_p==mode,:)),[40,97])*tt;
myContours = linspace(min(real(XIvel(:))),max(real(XIvel(:))), 21);
nexttile
contourf((real(XIvel)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))


% Energy test
sige=sum(sigu.^2+sigv.^2);
K_modes_energy=pinv(PSI_x(1:M2,1:N)*V)*transpose(sige(1:M2));
% 
% 
% 
% 
% 
