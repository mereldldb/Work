function coord = calc_trans(a, b, oth_coord, typ)
if typ == 'x' %we want to calculate x
    coord = abs((b(2)-oth_coord)/(a(2)-b(2))) * a(1) + abs((a(2)-oth_coord)/(a(2)-b(2))) * b(1);
elseif typ == 'y'
    coord = abs((b(1)-oth_coord)/(a(1)-b(1))) * a(2) + abs((a(1)-oth_coord)/(a(1)-b(1))) * b(2);
end