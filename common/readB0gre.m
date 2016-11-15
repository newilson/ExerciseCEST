function [B0maps, hdr1, ampl, phases, pathname1, pathname2] = readB0gre(pathname1,pathname2)

if nargin<2 || isempty(pathname1) || isempty(pathname2)
    [hdr1,ampl,~,pathname1] = readdicomfiles2d;
    [hdr2,phases,~,pathname2] = readdicomfiles2d;
else
    [hdr1,ampl] = readdicomfiles2d(pathname1);
    [hdr2,phases] = readdicomfiles2d(pathname2);
end
[nx,ny,nacq] = size(ampl);
ampl = reshape(ampl,[nx,ny,3,nacq/3]);
phases = reshape(phases,[nx,ny,3,nacq/3]);

phases = pi * ( (phases-2048.0)/2048.0 );  % Scale phase to be between -pi to +pi

B0maps = NWcalcB0gre([],phases,hdr1(1))/hdr1(1).sf;
for ii=1:nacq/3
    B0maps(:,:,ii) = anisodiff(B0maps(:,:,ii),20,50,0.03,1);
end

return
