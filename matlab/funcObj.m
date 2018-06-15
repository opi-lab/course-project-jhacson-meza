function [residue_sum,phased,gabor_filter] = funcObj(I,cx,cy,wavelength,s,sigmaOnf,thetaSigma,angl,D)
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
% return : number of residues
%===================================================================================

% gabor filter
[gabor_filter] = gaborfilter(I,cx,cy,s,wavelength,sigmaOnf,thetaSigma,angl,D);

% filtering
TF1 = I.* gabor_filter;
TFI=ifft2(TF1);
IMA=imag(TFI);
REA=real(TFI);
phased=atan2(IMA,REA);

%residue analysis
[~,residue_sum]=PhaseResidues(phased);

end

