function [gabor_filter] = gaborfilter(I,cx,cy,s,wavelength,sigmaOnf,thetaSigma,angl,D)

% % This code is part of research work done by Jesus Pineda Castro, Jackson Meza and 
% Juan Dominguez. Universidad Tecnológica de Bolívar , Cartagena, Colombia.
% parameter description:
% --------------------------
% I ----- fourier spectrum image
% cx ---- x coordinate
% cy ---- y coordinate
% s  ---- temporal image size
% wavelength ---- wavelength
% sigmaOnf ---- sst radial component
% thetaSigma ---- angular component
% angl ------ angle
% D ---- display results? ('true'/'false')
% return : gabor filter
%===================================================================================


% s ranges
if mod(s,2)
    srange = [-(s-1)/2:(s-1)/2]/(s-1);
else
    srange = [-s/2:(s/2-1)]/s;
end

[x,y] = meshgrid(srange);

% Radial log-Gabor component of the filter.

radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters

radius(1,1) = 1;                  % Get rid of the 0 radius value at the 0

lp = lowpassfilter([s,s],.45,15);   % Radius .45, 'sharpness' 15

fo = 1.0/wavelength;                  % Centre frequency of filter.

% The following implements the log-gabor transfer function.
logGabor = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));
logGabor(round(s/2+1), round(s/2+1)) = 0;     % Set the value at the 0 frequency point
% of the filter back to zero
% (undo the radius fudge).
logGabor = logGabor.*lp;  


% Angular component of the filter.

theta = atan2(y,x);               % Matrix values contain polar angle.
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.

sintheta = sin(theta);
costheta = cos(theta);

% For each point in the filter matrix calculate the angular distance from the
% specified filter orientation.  To overcome the angular wrap-around problem
% sine difference and cosine difference values are first computed and then
% the atan2 function is used to determine angular distance.

ds = sintheta * cos(angl) - costheta * sin(angl);    % Difference in sine.
dc = costheta * cos(angl) + sintheta * sin(angl);    % Difference in cosine.
dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
spread = exp((-dtheta.^2) / (2 * thetaSigma^2));     % Calculate the angular


% filter component.
filter = spread.*logGabor;                           % Product of the two components.
filter = ifftshift(filter);
filter_I = zeros(size(I,1),size(I,2));
filter_m = ~(filter < 1e-2);

% filter info
stats = regionprops(filter_m,'BoundingBox');
b = round(stats.BoundingBox);
filter_subregion = filter(b(2):b(2)+b(4),b(1):b(1)+b(3));
[~,idx] = max(filter_subregion(:));
[rows_r,cols_r] = ind2sub([size(filter_subregion,1),size(filter_subregion,2)],idx);
centroides = [cols_r + b(1)-1,rows_r + b(2)-1];

if strcmp(D,'true')
    figure(2),imagesc(filter),hold on,
    rectangle('Position',b,'EdgeColor','r')
    plot(centroides(1),centroides(2),'+b')
    hold off;
end

% filter centering
delta = [centroides(1)-b(1), b(3) + 1, centroides(2)-b(2), b(4) + 1];
filter_I(cy-delta(3): cy-delta(3) + delta(4) - 1,cx-delta(1): cx-delta(1) + delta(2) - 1) = filter_subregion;
filter_I = filter_I(1:size(I,1),1:size(I,2));

% output
gabor_filter = ifftshift(filter_I);

end

