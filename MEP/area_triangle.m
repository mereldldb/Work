function ar = area_triangle(a,b,c)
ar = .5*abs(a(1)*(b(2)-c(2))+b(1)*(c(2)-a(2))+c(1)*(a(2)-b(2)));