clc;clear all; close all;

% Grid
Nx = 200;
Ny = Nx;
Lx = 1;
Ly = 1;
dx = Lx/Nx;
dy = Ly/Ny;
x = dx/2:dx:Lx-dx/2;
y = dy/2:dy:Ly-dy/2;
[xx,yy] = meshgrid(x,y);

% Parameters
var_lnk = 2;
corr_lenx= 5*dx;
corr_leny= 5*dy;

[perm,var_lnk_actual]= random_perm(var_lnk,corr_lenx,corr_leny,Nx,Ny,Lx,Ly);

surf(xx,yy,log10(perm),'facecolor','interp','edgecolor','none','facelighting','phong');
view([0,0,1]); axis equal tight; 
colorbar; colormap(gray); %caxis([0 1]); 
