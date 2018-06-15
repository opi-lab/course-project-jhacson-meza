function [In] = imnormalize(I)

if ~isa(I,'double')
    In = double(I);
    In = In./max(In(:));
elseif any(I(:)<0)
    I = I + abs(min(I(:)));
    In = I./max(I(:));
else
    In = I./max(I(:));
end