function [Ar, Ah, Av] = cart2lud3(Ax, Ay, Az, theta, phi)
% CART2LUD3(Ax, Ay, Az, theta, phi)  
% Convert vector from cartesian to Ludwig 3 coordinates 
% compute elements in spherical coordinates using standard transformations
Ar = Ax.*sind(theta).*cosd(phi)+Ay.*sind(theta).*sind(phi)+Az.*cosd(theta);
Ah = Ax.*( 1 + (cosd(theta)-1).*(cosd(phi)).^2 ) + Ay.*( (cosd(theta)-1).*cosd(phi).*sind(phi) ) - Az.*( -sind(theta).*cosd(phi) );
Av = Ax.*( (cosd(theta)-1).*cosd(phi).*sind(phi) ) + Ay.*( 1 + (cosd(theta)-1).*(sind(phi)).^2 ) - Az.*( -sind(theta).*sind(phi) );
end
