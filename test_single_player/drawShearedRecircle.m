function [lat_coords,long_coords] = drawShearedRecircle(long, lat, shear_angle, half_width, half_height)
    m = shear_angle;
    w = half_width;
    h = half_height;
    X = linspace(-w, w,50);
    Y1 = (-(m*w^2).*X-( (h*w)^4 *(w^4-X.^4)).^(1/4))/w^2;
    Y2 = (-(m*w^2).*X+( (h*w)^4 *(w^4-X.^4)).^(1/4))/w^2;
    lat_coords = [X, X(end:-1:1)]+lat;
    long_coords = [Y1, Y2(end:-1:1)]+long;
end

