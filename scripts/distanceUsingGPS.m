function d = distanceUsingGPS(GPS1, GPS2)
R = 6371e3; % metres
diffLat_rad = (GPS2(1)-GPS1(1))*pi/180;
diffLong_rad = (GPS2(2)-GPS1(2))*pi/180;
GPS1_lat_rad = GPS1(1)*pi/180;
GPS2_lat_rad = GPS2(1)*pi/180;

a = sin(diffLat_rad/2) * sin(diffLat_rad/2) + cos(GPS2_lat_rad) * cos(GPS1_lat_rad) * sin(diffLong_rad/2) * sin(diffLong_rad/2);
c = 2 * atan2(sqrt(a), sqrt(1-a));
d = R * c;