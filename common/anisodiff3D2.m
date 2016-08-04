% ANISODIFF - Anisotropic diffusion.
%
% Usage:
%  im = anisodiff(im, niter, kappa, lambda, option)
%
% Arguments:
%         im     - input image
%         niter  - number of iterations.
%         kappa  - conduction coefficient 20-100 ?
%         lambda - max value of .25 for stability
%         option - 1 Perona Malik diffusion equation No 1
%                  2 Perona Malik diffusion equation No 2
%
% Returns:
%         im   - diffused image.
%
% kappa controls conduction as a function of gradient.  If kappa is low
% small intensity gradients are able to block conduction and hence diffusion
% across step edges.  A large value reduces the influence of intensity
% gradients on conduction.
%
% lambda controls speed of diffusion (you usually want it at a maximum of
% 0.25)
%
% Diffusion equation 1 favours high contrast edges over low contrast ones.
% Diffusion equation 2 favours wide regions over smaller ones.

% Reference:
% P. Perona and J. Malik.
% Scale-space and edge detection using ansotropic diffusion.
% IEEE Transactions on Pattern Analysis and Machine Intelligence,
% 12(7):629-639, July 1990.
%
% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk @ csse uwa edu au
% http://www.csse.uwa.edu.au
%
% June 2000  original version.
% March 2002 corrected diffusion eqn No 2.

function diff2 = anisodiff3D2(im, niter, kappa, lambda, option)

im = double(im);
[rows,cols,slices] = size(im);


diffl = zeros(rows+2, cols+2,slices+2);
diffl(2:rows+1, 2:cols+1,2:slices+1) = im;
for i = 1:niter
    % Construct diffl which is the same as im but
    % has an extra padding of zeros around it.
    
    % North, South, East and West differences
    deltaN = diffl(1:rows,2:cols+1,2:slices+1)   - im;
    deltaS = diffl(3:rows+2,2:cols+1,2:slices+1) - im;
    deltaE = diffl(2:rows+1,3:cols+2,2:slices+1) - im;
    deltaW = diffl(2:rows+1,1:cols,2:slices+1)   - im;
    deltaF = diffl(2:rows+1,2:cols+1,1:slices) - im;
    deltaB = diffl(2:rows+1,2:cols+1,3:slices+2) - im;
    
    % Conduction
    
    if option == 1
        cN = exp(-(deltaN/kappa).^2);
        cS = exp(-(deltaS/kappa).^2);
        cE = exp(-(deltaE/kappa).^2);
        cW = exp(-(deltaW/kappa).^2);
        cF = exp(-(deltaF/kappa).^2);
        cB = exp(-(deltaB/kappa).^2);
    elseif option == 2
        cN = 1./(1 + (deltaN/kappa).^2);
        cS = 1./(1 + (deltaS/kappa).^2);
        cE = 1./(1 + (deltaE/kappa).^2);
        cW = 1./(1 + (deltaW/kappa).^2);
        cF = 1./(1 + (deltaF/kappa).^2);
        cB = 1./(1 + (deltaB/kappa).^2);
    end
    
    im = im + lambda*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW+ cF.*deltaF + cB.*deltaB);
end
diff2 = im;
end

