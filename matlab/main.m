clear, clc
%% Read Image
I = double(imread('Data/I1.bmp'));

%% Automatic segmentation
[cx,cy,bounding_box,area,orientation,fftspectrum] = ...
    get_principal_component(I,80,'false');

%% Optimization with PSO

% filter parameter
wavelength = 20;

% objective function
f = @(inputs) funcObj(fftspectrum,cx,cy,wavelength,inputs(1),...
    inputs(2),inputs(3),inputs(4),'false');

% PSO initialization
n = 10; % Number of particles
xmin = [2000 0.75 0.1 -20*pi/180];
xmax = [4000 0.85 0.4  20*pi/180];
vmax = 0.1;
niter = 10; % Number of iterations
d = length(xmin); % Number of dimensions

x = (xmax-xmin).*rand(n,d)+xmin; % x initialization
v = (vmax/3+vmax/3)*rand(n,d) - vmax/3; % v inizialitation


pb = x; % Individual best inizialitation
for i = 1:n
    if i == 1
        gb = pb(i,:);
    end
    
    if f(pb(i,:)) < f(gb)
        gb = pb(i,:);
    end
end


iters = 0;
C1 = 0.72984*2.05; % Cognitive component
C2 = 0.72984*2.05; % Solcial componet
while 1
    iters = iters + 1;
    
    if iters ~= 1
        for i = 1:n
            if f(x(i,:)) < f(pb(i,:))
                pb(i,:) = x(i,:);
            end

            if f(pb(i,:)) < f(gb)
                gb = pb(i,:);
            end
        end
    end
    
    disp(['iter = ',num2str(iters),', gb(1) = ',num2str(gb(1)),...
        ', gb(2) = ',num2str(gb(2)),', gb(3) = ',num2str(gb(3)),...
        ', f(gb) = ',num2str(f(gb))])

    
    for i = 1:n
        for j = 1:d
            v(i,j) = v(i,j) + C1*rand*(pb(i,j) - x(i,j)) + ...
                C2*rand*(gb(j)-x(i,j));
            x(i,j) = x(i,j) + v(i,j);
            if x(i,j) > xmax(j)
                x(i,j) = xmax(j);
            elseif x(i,j) < xmin(j)
                x(i,j) = xmin(j);
            end
        end
    end
    
    if iters == niter
        break
    end
end

%% comparative study

% Using Hanning window filter
[phased_H,hanning_filter] = AutoFiltre(I,80,1);

% Using optimun gabor filter
s = gb(1);
sigma_r = gb(2);
theta_sigma = gb(3);
theta_0 = gb(4);

[residues_G,phased_G,gabor_filter] = funcObj(fftspectrum,cx,cy,...
    wavelength,s,sigma_r,theta_sigma,theta_0,'false');



% show results
% ------------filters
figure(1),imagesc(ifftshift(gabor_filter))
figure(2),imagesc(hanning_filter)
%-------------phase
% figure(3),imagesc(phased_H),colormap gray
% figure(4),imagesc(phased_G),colormap gray

% residue results
%-------------hanning
% figure(5),imagesc(phased_H),colormap gray
%load rec.mat
%[phased_H_roi] = imcrop(phased_H,rec);
[residue_charge_H,rnum_H]=PhaseResidues(phased_H);
[row,col] = find(residue_charge_H);
figure(6),imagesc(phased_H),colormap gray,hold on;
plot(col,row,'ro','LineWidth',2,'MarkerSize',6),hold off

%-------------gabor
% figure(7),imagesc(phased_G),colormap gray
%[phased_G_roi] = imcrop(phased_G,rec);
[residue_charge_G,rnum_G]=PhaseResidues(phased_G);
[row,col] = find(residue_charge_G);
figure(8),imagesc(phased_G),colormap gray,hold on;
plot(col,row,'ro','LineWidth',2,'MarkerSize',6),hold off

disp(['Hanning residues: ',num2str(rnum_H)])
disp(['Gabor residues: ',num2str(rnum_G)])