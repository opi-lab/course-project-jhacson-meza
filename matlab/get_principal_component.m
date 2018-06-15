function [cx,cy,bounding_box,area,orientation,F] = get_principal_component(I,T,D)

% % This code is part of research work done by Jesus Pineda Castro, Jackson Meza and 
% Juan Dominguez. Universidad Tecnológica de Bolívar , Cartagena, Colombia.
%--------------------------------------------------------------------------------------
% parameter description:
% --------------------------
% I ----------- input fringe image
% T ----------- segmentation threshold
% D ----------- display results? ('true'/'false')
% return : x - y coordinates
%=====================================================================================

% preprocessing

Gx=hanning(size(I,2));
Gy=hanning(size(I,1));
AA=Gy*Gx';
AA=(AA-min(AA(:)))/(max(AA(:))-min(AA(:)));

F = fft2(double(I));
Fim = fft2(double(I).*AA);
Ih = fftshift(abs(Fim));
Ih(:,end/2-5:end/2+5)=0;
In = round(1024*imnormalize(Ih));

if strcmp(D,'true')
    figure(700), imagesc(In), colormap gray, hold on
end
[h,w] = size(In);
offset_w = round(w*0.04);
offset_h = round(h*0.04);

%  Seeds - regiongrow segmentation

seed_tl = [offset_h offset_w];
seed_tr = [offset_h w-offset_w];
seed_bl = [h-offset_h offset_w];
seed_br = [h-offset_h w-offset_w];
if strcmp(D,'true')
    plot(offset_w,offset_h,'sr')
    plot(w-offset_w,offset_h,'sr')
    plot(offset_w,h-offset_h,'sr')
    plot(w-offset_w,h-offset_h,'sr')
end

S = zeros(h,w);
S(seed_tl(1),seed_tl(2)) = 1;
S(seed_tr(1),seed_tr(2)) = 1;
S(seed_bl(1),seed_bl(2)) = 1;
S(seed_br(1),seed_br(2)) = 1;

mask = ~regiongrow(In,S,T);
if strcmp(D,'true')
    figure(799), imagesc(mask),colormap gray
    figure(800), imagesc(In), hold on
end

I2 = In.*mask;
marco=regionprops(mask,I2,'BoundingBox','Area','Orientation');
talla=size(marco,1);
Areas = zeros([talla 1]); Rec=zeros([talla 4]); centroides=zeros([talla,2]);
orientation = zeros([talla 1]);
for i=1:talla
    Areas(i) = marco(i).Area;
    orientation(i) = marco(i).Orientation;
    Rec(i,:) = marco(i).BoundingBox; b = round(Rec(i,:));
    subregion = I2(b(2):b(2)+b(4),b(1):b(1)+b(3));
    [~,idx] = max(subregion(:));
    [rows_r,cols_r] = ind2sub([size(subregion,1),size(subregion,2)],idx);
    centroides(i,:) = [cols_r + b(1)-1,rows_r + b(2)-1];
    if strcmp(D,'true')
        rectangle('Position',Rec(i,:),'EdgeColor','r')
        plot(centroides(i,1),centroides(i,2),'+b')
    end
end

Areas_analis = Areas;
for i=1:4
    [~,p] = max(Areas_analis(:));
    pos(i) = p;
    Areas_analis(p)=0;
end

[~,p]=max(centroides(pos,1));

% outputs
cx = round(centroides(pos(p),1));
cy = round(centroides(pos(p),2));
bounding_box = Rec(pos(p),:);
area = Areas(pos(p));
orientation = orientation(pos(p));

end

