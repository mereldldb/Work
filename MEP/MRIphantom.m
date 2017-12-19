function p=MRIphantom(n)
%               A    a     b    x0    y0    phi
%        ---------------------------------
ellipse = [     0.8000    0.7200    0.9500         0         0         0
                0.1200    0.6900    0.9200         0         0         0
                0.9800    0.6624    0.8740         0   -0.0184         0
                0.7450    0.6524    0.8640         0   -0.0184         0
                0.9800    0.4100    0.1600   -0.2200         0  -72.0000
                0.9800    0.3100    0.1100    0.2200         0   72.0000
                0.6170    0.2100    0.2500         0    0.3500         0
                0.9500    0.0460    0.0460         0    0.1000         0
                0.9500    0.0460    0.0230   -0.0800   -0.6050         0
                0.9500    0.0460    0.0230    0.0600   -0.6050  -90.0000
                0.9500    0.0460    0.0460         0   -0.1000         0
                0.9500    0.0230    0.0230         0   -0.6050         0
];

overlap = [0 0 0 0 7 0 8 0 0 0 5 0];

p = zeros(n);

xax =  ( (0:n-1)-(n-1)/2 ) / ((n-1)/2); 
xg = repmat(xax, n, 1);   % x coordinates, the y coordinates are rot90(xg)

for k = 1:size(ellipse,1)    
   asq = ellipse(k,2)^2;       % a^2
   bsq = ellipse(k,3)^2;       % b^2
   phi = ellipse(k,6)*pi/180;  % rotation angle in radians
   x0 = ellipse(k,4);          % x offset
   y0 = ellipse(k,5);          % y offset
   A = ellipse(k,1);           % Amplitude change for this ellipse
   x=xg-x0;                    % Center the ellipse
   y=rot90(xg)-y0;  
   cosp = cos(phi); 
   sinp = sin(phi);
   idx=find(((x.*cosp + y.*sinp).^2)./asq + ((y.*cosp - x.*sinp).^2)./bsq <= 1); 
   p(idx) = A;

end

for j = 1:size(ellipse,1)
  if overlap(j) ~= 0
        e1 = j;
        asq = ellipse(e1,2)^2;       
        bsq = ellipse(e1,3)^2;       
        phi = ellipse(e1,6)*pi/180;  
        x0 = ellipse(e1,4);          
        y0 = ellipse(e1,5);          
        A = ellipse(e1,1);           
        x=xg-x0;                    
        y=rot90(xg)-y0;  
        cosp = cos(phi); 
        sinp = sin(phi);
        e2 = overlap(j); 
        cosp2 = cos(ellipse(e2,6)*pi/180);
        sinp2 = sin(ellipse(e2,6)*pi/180);
        asq2 = ellipse(e2,2)^2;       
        bsq2 = ellipse(e2,3)^2; 
        x0_2 = ellipse(e2,4);         
        y0_2 = ellipse(e2,5);          
        A2 = ellipse(e2,1);           
        x2=xg-x0_2;                    
        y2=rot90(xg)-y0_2;  

        idx2 = find((((x.*cosp + y.*sinp).^2)./asq + ((y.*cosp - x.*sinp).^2)./bsq <= 1) & ...
            (((x2.*cosp2 + y2.*sinp2).^2)./asq2 + ((y2.*cosp2 - x2.*sinp2).^2)./bsq2 <= 1)); 
        p(idx2) = .5*(A+A2);    
  end 
end   
