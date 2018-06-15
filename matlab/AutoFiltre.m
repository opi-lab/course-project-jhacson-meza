function [phased, FILTRE, F,mod] = AutoFiltre(I,T,pi)
% I: Imagen de entrada
% T: Umbral para el crecimiento de regiones
% pi=1 Se muestran las graficas
% pi=[] no se muestran

if nargin<3
    pi=0;
end

tic
Gx=hanning(size(I,2));
Gy=hanning(size(I,1));
AA=Gy*Gx';
AA=(AA-min(AA(:)))/(max(AA(:))-min(AA(:)));

F = fft2(double(I));
Fim = fft2(double(I).*AA);
Ih = fftshift(abs(Fim));
Ih(:,end/2-5:end/2+5)=0;
In = round(1024*imnormalize(Ih));
if pi
figure(700), imagesc(In), colormap gray, hold on
end
[h,w] = size(In);
offset_w = round(w*0.04);
offset_h = round(h*0.04);


%  Seeds
seed_tl = [offset_h offset_w];
seed_tr = [offset_h w-offset_w];
seed_bl = [h-offset_h offset_w];
seed_br = [h-offset_h w-offset_w];
if(pi)
plot(offset_w,offset_h,'sr')
plot(w-offset_w,offset_h,'sr')
plot(offset_w,h-offset_h,'sr')
plot(w-offset_w,h-offset_h,'sr')
end
% figure(701), imhist(uint8(In)) 
% C= imhist(uint8(In));


S = zeros(h,w);
S(seed_tl(1),seed_tl(2)) = 1;
S(seed_tr(1),seed_tr(2)) = 1;
S(seed_bl(1),seed_bl(2)) = 1;
S(seed_br(1),seed_br(2)) = 1;

%T =80;% Threshold
%[rel,T]=fr_ac(C,P);
%T=T+1;
% mask es la mascara a resolucion [256 nan]
%%
mask = ~regiongrow(In,S,T);
if(pi==1)
    figure(799), imagesc(mask),colormap gray
    figure(800), imagesc(In), hold on
end

marco=regionprops(mask,'BoundingBox','Area','Centroid');
talla=size(marco,1);
Areas = zeros([talla 1]); Rec=zeros([talla 4]); centroides=zeros([talla,2]);
for i=1:talla
Areas(i) = marco(i).Area;
Rec(i,:) = marco(i).BoundingBox;
centroides(i,:) = marco(i).Centroid;
if(pi==1)
    rectangle('Position',Rec(i,:),'EdgeColor','r')
    plot(centroides(i,1),centroides(i,2),'+b')
end
end

%distancia = sqrt((w/2-centroides(:,1)).^2+(h/2-centroides(:,2)).^2);
%%
Areas_analis = Areas;

for i=1:4
[m,p] = max(Areas_analis(:));
pos(i) = p;  
Areas_analis(p)=0;
end

[m,p]=max(centroides(pos,1));

cx = round(centroides(pos(p),1));
cy = round(centroides(pos(p),2));
c = round(cx-size(In,2)/2);
w = round(2/3*c);
%cy = dy-size(In,1)/2;
if(pi==1)
figure(804), imagesc(In), hold on 
rectangle('Position',[cx-w/2,cy-2*w,w,4*w],'EdgeColor','r')
plot(cx,cy,'+b')
end
%%
%toc

Gx=hanning(numel(cx-w/2:cx+w/2));
Gy=hanning(numel(cy-2*w:cy+2*w));
AA=Gy*Gx';
AA=(AA-min(AA(:)))/(max(AA(:))-min(AA(:)));

FILTRE = zeros(size(In));
FILTRE(cy-2*w:cy+2*w,cx-w/2:cx+w/2)=AA;

TF1=F.*fftshift(FILTRE);
TFI=ifft2(TF1);
IMA=imag(TFI);
REA=real(TFI);
phased=atan2(IMA,REA);
mod = sqrt(IMA.*IMA + REA.*REA);
end

