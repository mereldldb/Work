function G = generate_rotation_matrix(npx_original, npy_original, ratxy, rot_angle)
 
    %npx_original = 32; %no. of pixels x-direction
    %npy_original = 32; %no. of pixels y-direction
    %ratxy = 1; %pixel height (y)/pixel length (x)
 
    ratxy2 = ratxy; 
    no_of_divisions = 0;
    while ratxy2 > 1.5
        no_of_divisions = no_of_divisions + 1;
        ratxy2 = ratxy2/2;
    end
 
    no_of_multiplications = 0;
    while ratxy2 < .75
        no_of_multiplications = no_of_multiplications + 1;
        ratxy2 = ratxy2*2;
    end
 
    npx = npx_original*2^no_of_multiplications;
    npy = npy_original*2^no_of_divisions;
 
    length_x = 1;
    length_y = ratxy2;
 
    theta = degtorad(rot_angle);
 
    rot_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
 
    pixel_numbers = reshape(1:npx*npy,npy,npx);
 
    center_pixels_x = ((-npx/2+.5):(npx/2-.5))*length_x;
    vertices_pixels_x = (-npx/2:npx/2)*length_x;
 
    center_pixels_y = ((npy/2-.5):-1:(-npy/2+.5))*length_y;    
    vertices_pixels_y = (npy/2:-1:-npy/2)*length_y;
 
    [center_pixels_x, center_pixels_y] = meshgrid(center_pixels_x, center_pixels_y);
    [vertices_pixels_x, vertices_pixels_y] = meshgrid(vertices_pixels_x, vertices_pixels_y);
 
    %% how much padding do we need?
    padding_radius = sqrt((npx/2*length_x)^2+(npy/2*length_y)^2);
    npx_padding_one_side = ceil((padding_radius-npx/2*length_x)/length_x)+1; %we have to add the same amount of pixels on the left and on the right
    npy_padding_one_side = ceil((padding_radius-npy/2*length_y)/length_y)+1; %we have to add the same amount of pixels at the top and at the bottom
 
    npx_padded = npx+npx_padding_one_side*2;
    npy_padded = npy+npy_padding_one_side*2;
 
    %% with padding
 
    padded_pixel_numbers = reshape(1:npx_padded*npy_padded,npy_padded,npx_padded);
    %which pixels in the padded array are the real pixels?
    real_pixels_in_padded = padded_pixel_numbers(npy_padding_one_side+1:npy+npy_padding_one_side,npx_padding_one_side+1:npx+npx_padding_one_side);
 
    center_padded_x = ((-npx_padded/2+.5):(npx_padded/2-.5))*length_x;
    vertices_padded_x = (-npx_padded/2:npx_padded/2)*length_x;
 
    center_padded_y = ((npy_padded/2-.5):-1:(-npy_padded/2+.5))*length_y;    
    vertices_padded_y = (npy_padded/2:-1:-npy_padded/2)*length_y;
 
    [center_padded_x, center_padded_y] = meshgrid(center_padded_x, center_padded_y);
    [vertices_padded_x, vertices_padded_y] = meshgrid(vertices_padded_x, vertices_padded_y);
 
    %% special cases: rotation by 0 or 180 degrees
    if rot_angle == 0
        G = speye(npx*npy);
        return
    elseif rot_angle == 180
        G = rot90(speye(npx*npy));
        return
    end
 
    %%
    cent_rot_x = zeros(npy,npx);
    cent_rot_y = zeros(npy,npx);
    vert_rot_x = zeros(npy+1,npx+1);
    vert_rot_y = zeros(npy+1,npx+1);
    center_which_pixel = zeros(npy,npx);
    vert_which_pixel = zeros(npy+1,npx+1);
 
    for i_pixels = 1:npx*npy
        cent_rot = rot_matrix*[center_pixels_x(i_pixels); center_pixels_y(i_pixels)];
        cent_rot_x(i_pixels) = cent_rot(1);
        cent_rot_y(i_pixels) = cent_rot(2);
        % in which "old" pixel is the new center located?
        [~,center_which_pixel(i_pixels)] = min(abs((cent_rot_x(i_pixels)-reshape(center_padded_x,[],1)).^2+(cent_rot_y(i_pixels)-reshape(center_padded_y,[],1)).^2));
    end
 
    for i_vertices = 1:(npx+1)*(npy+1)
        vert_rot = rot_matrix*[vertices_pixels_x(i_vertices); vertices_pixels_y(i_vertices)];
        vert_rot_x(i_vertices) = vert_rot(1);
        vert_rot_y(i_vertices) = vert_rot(2);
        % in which "old" pixel is the new vertex located?
        %[~,vert_which_pixel(i_vertices)] = min(abs((vert_rot_x(i_vertices)-reshape(center_padded_x,[],1)).^2+(vert_rot_y(i_vertices)-reshape(center_padded_y,[],1)).^2));
    end
 
%     %%
%     figure(2)
%     hold on
%     plot(reshape(vertices_padded_x,[],1),reshape(vertices_padded_y,[],1),'b.')
%     plot(reshape(vert_rot_x,[],1),reshape(vert_rot_y,[],1),'r.')
%     plot(reshape(cent_rot_x,[],1),reshape(cent_rot_y,[],1),'k.')
%     axis square
 
    %%
    G = sparse(npx*npy,npx_padded*npy_padded);
    i_pixels = 0;
    count_0_vertex = 0;
 
    correction_x = 0;
    correction_y = 0;
 
    if rem(npx,2) == 1
        correction_x = .5*length_x;
    end
    if rem(npy,2) == 1
        correction_y = .5*length_y;
    end
 
    for i_x = 1:npx
        for i_y = 1:npy
            i_pixels = i_pixels + 1;
    %         if i_y == 3 && i_x == 1
    %             break
    %         end
 
            cp = center_which_pixel(i_pixels);
            pixels_of_interest = [cp-npy_padded-1, cp-npy_padded, cp-npy_padded+1; cp-1 cp cp+1; cp+npy_padded-1 cp+npy_padded cp+npy_padded+1]';
 
            vertices_x = vert_rot_x(i_y:i_y+1,i_x:i_x+1);
            vertices_y = vert_rot_y(i_y:i_y+1,i_x:i_x+1);
            for l_count = 1:4
               dis_to_centers = round(10^13*abs((vertices_x(l_count)-reshape(center_padded_x,[],1)).^2+(vertices_y(l_count)-reshape(center_padded_y,[],1)).^2))/10^13;
               %dis_to_centers = abs((vertices_x(l_count)-reshape(center_padded_x,[],1)).^2+(vertices_y(l_count)-reshape(center_padded_y,[],1)).^2);
               centers_closest = find(dis_to_centers == min(dis_to_centers));
 
               if length(centers_closest) == 1
                   vert_pixels(l_count) = centers_closest;
               else
                   if ismember(cp,centers_closest)
                       vert_pixels(l_count) = cp;
                   elseif ismember(pixels_of_interest(2),centers_closest)
                        vert_pixels(l_count) = pixels_of_interest(2);
                        vertices_y(l_count) = round((vertices_y(l_count)+correction_y)/length_y)*length_y;
                   elseif ismember(pixels_of_interest(4),centers_closest)
                        vert_pixels(l_count) = pixels_of_interest(4);
                        vertices_x(l_count) = round((vertices_x(l_count)+correction_x)/length_x)*length_x;
                   elseif ismember(pixels_of_interest(6),centers_closest)
                        vert_pixels(l_count) = pixels_of_interest(6);
                        vertices_x(l_count) = round((vertices_x(l_count)+correction_x)/length_x)*length_x;
                   elseif ismember(pixels_of_interest(8),centers_closest)
                        vert_pixels(l_count) = pixels_of_interest(8);
                        vertices_y(l_count) = round((vertices_y(l_count)+correction_y)/length_y)*length_y;
                   end
               end
            end
            vert_pixels = reshape(vert_pixels,2,2);
    %         vert_pixels = vert_which_pixel(i_y:i_y+1,i_x:i_x+1);
 
 
            vertices_location_grid = cell(3);
 
            for o_count = 1:4
                vertices_location_grid{find(pixels_of_interest == vert_pixels(o_count))} = [vertices_location_grid{find(pixels_of_interest == vert_pixels(o_count))}, [vertices_x(o_count); vertices_y(o_count)]];
            end
            
            comb = [1 2; 1 3; 2 4; 3 4]; %combinations of vertices: left top bottom right
 
            x_trans = cell(1,4);
            y_trans = cell(1,4);
            cent_pixel1 = cell(1,4);
            cent_pixel2 = cell(1,4);
 
            contributing_pixels = [];
            transitions_to_other = cell(3); %transition from pixel to what other pixels?
            transitions_grid = cell(3);
            matr = zeros(3);        
 
            [~,~,vert_indices] = intersect(vert_pixels,pixels_of_interest);
            vert_logical_grid = zeros(3);
            vert_logical_grid(vert_indices) = pixels_of_interest(vert_indices) > 0;
 
            if rot_angle == 90 || rot_angle == 270
                if length_x == length_y
                    switch rot_angle
                        case 90
                            target_pixel =  padded_pixel_numbers(npx_padding_one_side+npx+1-i_x,npy_padding_one_side+i_y);
                        case 270
                            target_pixel = padded_pixel_numbers(npx_padding_one_side+i_x,npy_padding_one_side+npy+1-i_y);
                    end
                    G(i_pixels,target_pixel) = 1;  
                else
                cont_pixels_grid = (vert_logical_grid > 0);
                switch rot_angle
                    case 90
                        vertices_x = rot90(vertices_x);
                        vertices_y = rot90(vertices_y);
                    case 270
                        vertices_x = rot90(rot90(rot90(vertices_x)));
                        vertices_y = rot90(rot90(rot90(vertices_y)));
                end
 
                if cont_pixels_grid(2,1) == 1 && cont_pixels_grid(2,3) == 1
                    cont_pixels_grid(2,2) = 1;
                end
                if cont_pixels_grid(1,2) == 1 && cont_pixels_grid(3,2) == 1
                    cont_pixels_grid(2,2) = 1;
                end
                if cont_pixels_grid(1,1) == 1 && cont_pixels_grid(1,3) == 1
                    cont_pixels_grid(1,2) = 1;
                end
                if cont_pixels_grid(3,1) == 1 && cont_pixels_grid(3,3) == 1
                    cont_pixels_grid(3,2) = 1;
                end
                if cont_pixels_grid(1,3) == 1 && cont_pixels_grid(3,3) == 1
                    cont_pixels_grid(2,3) = 1;
                end
                if cont_pixels_grid(1,1) == 1 && cont_pixels_grid(3,1) == 1
                    cont_pixels_grid(2,1) = 1;
                end
 
                pixel_counter = 0;
 
                if (mod(vertices_x(1)+correction_x,length_x) < 1e-12) || (mod(vertices_x(1)+correction_x,length_x) > length_x-1e-12)
                    cont_pixels_grid(:,1) = 0;
                end
                if (mod(vertices_x(3)+correction_x,length_x) < 1e-12) || (mod(vertices_x(3)+correction_x,length_x) > length_x-1e-12)
                    cont_pixels_grid(:,3) = 0;
                end
                if (mod(vertices_y(1)+correction_y,length_y) < 1e-12) || (mod(vertices_y(1)+correction_y,length_y) > length_y-1e-12)
                    cont_pixels_grid(1,:) = 0;
                end
                if (mod(vertices_y(2)+correction_y,length_y) < 1e-12) || (mod(vertices_y(2)+correction_y,length_y) > length_y-1e-12)
                    cont_pixels_grid(3,:) = 0;
                end
 
                for m_count = 1:9
                    eps = 1e-6;
                    %sum(sum(cont_pixels_grid))
                    if cont_pixels_grid(m_count) ~= 0
                        if sum(sum(cont_pixels_grid)) == 2
                            pixel_counter = pixel_counter + 1;
                            if length_x < length_y 
                                switch pixel_counter
                                    case 1
                                        point1 = [vertices_x(1);  vertices_y(1)];
                                        point2 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];
                                        point3 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                        point4 = [vertices_x(2);  vertices_y(2)];
                                    case 2
                                        point1 = [vertices_x(3);  vertices_y(3)];
                                        point2 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                        point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];
                                        point4 = [vertices_x(4);  vertices_y(4)];
                                end
                            else
                                switch pixel_counter
                                    case 1
                                        point1 = [vertices_x(1);  vertices_y(1)];
                                        point2 = [vertices_x(1); floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point4 = [vertices_x(3);  vertices_y(3)];
                                    case 2
                                        point1 = [vertices_x(2);  vertices_y(2)];
                                        point2 = [vertices_x(2); ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [vertices_x(4); ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                        point4 = [vertices_x(4);  vertices_y(4)];
                                end
                            end
 
                        elseif sum(sum(cont_pixels_grid)) == 3
                            pixel_counter = pixel_counter + 1;
                            if length_x < length_y 
                                switch pixel_counter
                                    case 1
                                        point1 = [vertices_x(1);  vertices_y(1)];
                                        point2 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];
                                        point3 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                        point4 = [vertices_x(2);  vertices_y(2)];
                                    case 2 
                                        point1 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];
                                        point2 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                        point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];
                                        point4 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                    case 3
                                        point1 = [vertices_x(3);  vertices_y(3)];
                                        point2 = [ceil((vertices_x(3)+correction_x+eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                        point3 = [ceil((vertices_x(4)+correction_x+eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];
                                        point4 = [vertices_x(4);  vertices_y(4)];                               
                                end
                            else
                                switch pixel_counter
                                    case 1
                                        point1 = [vertices_x(1);  vertices_y(1)];
                                        point2 = [vertices_x(1); floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point4 = [vertices_x(3);  vertices_y(3)];
                                    case 2 
                                        point1 = [vertices_x(1); floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point2 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [vertices_x(4); floor((vertices_y(4)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point4 = [vertices_x(2); floor((vertices_y(2)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                    case 3
                                        point1 = [vertices_x(2);  vertices_y(2)];
                                        point2 = [vertices_x(2);  floor((vertices_y(2)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [vertices_x(4);  floor((vertices_y(4)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point4 = [vertices_x(4);  vertices_y(4)];
                                end
                            end
                        elseif sum(sum(cont_pixels_grid)) == 4
                            pixel_counter = pixel_counter + 1;
                            switch pixel_counter
                                case 1
                                    point1 = [vertices_x(1);  vertices_y(1)];
                                    point2 = [vertices_x(1); floor((vertices_y(1)-eps)/length_y)*length_y] - [0; correction_y];
                                    point3 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                    point4 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];
                                case 2
                                    point1 = [vertices_x(2);  vertices_y(2)];
                                    point2 = [vertices_x(2); ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                    point3 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                    point4 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                case 3
                                    point1 = [vertices_x(3);  vertices_y(3)];
                                    point2 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                    point3 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                    point4 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                case 4
                                    point1 = [vertices_x(4);  vertices_y(4)];
                                    point2 = [vertices_x(4); ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                    point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                    point4 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];
                            end
                        elseif sum(sum(cont_pixels_grid)) == 6
                            pixel_counter = pixel_counter + 1;
                            if length_x < length_y                  
                                switch pixel_counter
                                    case 1
                                        point1 = [vertices_x(1);  vertices_y(1)];
                                        point2 = [vertices_x(1); floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];                            
                                    case 2
                                        point1 = [vertices_x(2);  vertices_y(2)];
                                        point2 = [vertices_x(2); ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                    case 3
                                        point1 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];
                                        point2 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point3 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                    case 4
                                        point1 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                        point2 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; ceil((vertices_y(4)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];
                                    case 5
                                        point1 = [vertices_x(3);  vertices_y(3)];
                                        point2 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                    case 6
                                        point1 = [vertices_x(4);  vertices_y(4)];
                                        point2 = [vertices_x(4); ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];  
                                end
                            else
                                switch pixel_counter
                                    case 1
                                        point1 = [vertices_x(1);  vertices_y(1)];
                                        point2 = [vertices_x(1); floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; vertices_y(1)] - [correction_x; 0];                           
                                    case 2
                                        point1 = [vertices_x(1); floor((vertices_y(1)-eps)/length_y)*length_y] - [0; correction_y];
                                        point2 = [ceil((vertices_x(1)+correction_x+eps)/length_x)*length_x; floor((vertices_y(1)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point3 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [vertices_x(2); ceil((vertices_y(2)+eps)/length_y)*length_y] - [0; correction_y];
                                    case 3
                                        point1 = [vertices_x(2);  vertices_y(2)];
                                        point2 = [vertices_x(2); ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; ceil((vertices_y(2)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [ceil((vertices_x(2)+correction_x+eps)/length_x)*length_x; vertices_y(2)] - [correction_x; 0];
                                    case 4
                                        point1 = [vertices_x(3);  vertices_y(3)];
                                        point2 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; vertices_y(3)] - [correction_x; 0];
                                    case 5
                                        point1 = [vertices_x(3); floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [0; correction_y];
                                        point2 = [floor((vertices_x(3)+correction_x-eps)/length_x)*length_x; floor((vertices_y(3)+correction_y-eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [vertices_x(4); ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                    case 6
                                        point1 = [vertices_x(4);  vertices_y(4)];
                                        point2 = [vertices_x(4); ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [0; correction_y];
                                        point3 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; ceil((vertices_y(4)+correction_y+eps)/length_y)*length_y] - [correction_x; correction_y];
                                        point4 = [floor((vertices_x(4)+correction_x-eps)/length_x)*length_x; vertices_y(4)] - [correction_x; 0];                      
                                end
                            end
                        end
                        XY = [point1, point2, point3, point4];
                        matr(m_count) = polyarea(XY(1,:),XY(2,:));
                    end
                end
                G(i_pixels,pixels_of_interest) = reshape(matr/(length_x*length_y),1,[]);
                end
            else
            if sum(sum(vert_logical_grid)) ~= 4
                [vert_pixels_vec_count,vert_pixels_vec] = hist(reshape(vert_pixels,[],1),unique(vert_pixels));
                vert_logical_grid(vert_indices) = vert_pixels_vec_count;
            end
 
            if sum((reshape(vertices_x,[],1) == 0) & (reshape(vertices_y,[],1) == 0)) ~= 0
                %do something that allows us to calculate areas for pixels with vertices (0,0)
                where_is_0 = find((reshape(vertices_x,[],1) == 0) & (reshape(vertices_y,[],1) == 0));
                vert_pixels(where_is_0) = cp;
                %we have to add this in case one of the vertices is in the center pixel --> then (0,0) is not counted as a transition point
                switch where_is_0
                    case 1
                        x_trans{1} = 0;
                        y_trans{1} = 0;
                        x_trans{2} = 0;
                        y_trans{2} = 0;
                    case 2
                        x_trans{1} = 0;
                        y_trans{1} = 0;
                        x_trans{3} = 0;
                        y_trans{3} = 0;
                    case 3
                        x_trans{4} = 0;
                        y_trans{4} = 0;
                        x_trans{2} = 0;
                        y_trans{2} = 0;
                    case 4
                        x_trans{4} = 0;
                        y_trans{4} = 0;
                        x_trans{3} = 0;
                        y_trans{3} = 0;
                end
            end
            %else
            x_trans_tot = [];
            y_trans_tot = [];
                for j_count = 1:4
                    if vert_pixels(comb(j_count,1)) ~= vert_pixels(comb(j_count,2))
                        if ismember(abs(vert_pixels(comb(j_count,1))-vert_pixels(comb(j_count,2))),[2*npy_padded+1,npy_padded+1,1,npy_padded-1, 2*npy_padded-1]) %one top/bottom boundary is crossed --> know y
                            y_trans{j_count} = [y_trans{j_count}, mean([center_padded_y(vert_pixels(comb(j_count,1))),center_padded_y(vert_pixels(comb(j_count,2)))])];
                            x_trans{j_count} = [x_trans{j_count}, calc_trans([vertices_x(comb(j_count,1)); vertices_y(comb(j_count,1))],[vertices_x(comb(j_count,2)); vertices_y(comb(j_count,2))],y_trans{j_count}(end),'x')];
                        elseif ismember(abs(vert_pixels(comb(j_count,1))-vert_pixels(comb(j_count,2))),[npy_padded-2,2,npy_padded+2]) %two top/bottom boundaries are crossed
                            y_trans{j_count} = [y_trans{j_count}, mean([center_padded_y(vert_pixels(comb(j_count,1))),center_padded_y(vert_pixels(comb(j_count,2)))])-.5*length_y, mean([center_padded_y(vert_pixels(comb(j_count,1))),center_padded_y(vert_pixels(comb(j_count,2)))])+0.5*length_y];
                            x_trans{j_count} = [x_trans{j_count}, calc_trans([vertices_x(comb(j_count,1)); vertices_y(comb(j_count,1))],[vertices_x(comb(j_count,2)); vertices_y(comb(j_count,2))],y_trans{j_count}(end-1),'x'), calc_trans([vertices_x(comb(j_count,1)); vertices_y(comb(j_count,1))],[vertices_x(comb(j_count,2)); vertices_y(comb(j_count,2))],y_trans{j_count}(end),'x')];                        
                        end
                        if  ismember(abs(vert_pixels(comb(j_count,1))-vert_pixels(comb(j_count,2))),[npy_padded-2,npy_padded-1,npy_padded,npy_padded+1,npy_padded+2]) %one left/right boundary is crossed --> know x
                            x_trans{j_count} =  [x_trans{j_count}, mean([center_padded_x(vert_pixels(comb(j_count,1))),center_padded_x(vert_pixels(comb(j_count,2)))])];
                            y_trans{j_count} = [y_trans{j_count}, calc_trans([vertices_x(comb(j_count,1)); vertices_y(comb(j_count,1))],[vertices_x(comb(j_count,2)); vertices_y(comb(j_count,2))],x_trans{j_count}(end),'y')];
                        elseif ismember(abs(vert_pixels(comb(j_count,1))-vert_pixels(comb(j_count,2))),[2*npy_padded-1,2*npy_padded,2*npy_padded+1]) %two left/right boundaries
                            x_trans{j_count} =  [x_trans{j_count}, mean([center_padded_x(vert_pixels(comb(j_count,1))),center_padded_x(vert_pixels(comb(j_count,2)))])-0.5*length_x, mean([center_padded_x(vert_pixels(comb(j_count,1))),center_padded_x(vert_pixels(comb(j_count,2)))])+0.5*length_x];
                            y_trans{j_count} = [y_trans{j_count}, calc_trans([vertices_x(comb(j_count,1)); vertices_y(comb(j_count,1))],[vertices_x(comb(j_count,2)); vertices_y(comb(j_count,2))],x_trans{j_count}(end-1),'y'), calc_trans([vertices_x(comb(j_count,1)); vertices_y(comb(j_count,1))],[vertices_x(comb(j_count,2)); vertices_y(comb(j_count,2))],x_trans{j_count}(end),'y')];
 
                        end
                    end
                    x_trans_tot = [x_trans_tot, x_trans{j_count}];
                    y_trans_tot = [y_trans_tot, y_trans{j_count}];
                
                end
                if sum((abs(x_trans_tot) < 1e-13) & (abs(y_trans_tot) < 1e-13)) > 1
                    %(0,0) is counted twice when one of the vertices is (0,0) and the other is not in the center pixel
                    indices_trans_0 = find(abs(x_trans_tot) < 1e-13 & abs(y_trans_tot) < 1e-13);
                    x_trans_tot(indices_trans_0(2:end)) = [];
                    x_trans_tot(indices_trans_0(1)) = 0;                
                    y_trans_tot(indices_trans_0(2:end)) = [];
                    y_trans_tot(indices_trans_0(1)) = 0;                
                end
                x_trans_old = x_trans_tot;
                y_trans_old = y_trans_tot;
 
                trans_tot = [x_trans_tot', y_trans_tot'];
                trans_tot = unique(trans_tot,'rows');
                x_trans_tot = trans_tot(:,1)';
                y_trans_tot = trans_tot(:,2)';
 
                k_count_on_vertex = [];
                for k_count = 1:length(x_trans_tot)
                    % by determining where the transition between pixels occurs, we know which two pixels contribute to the final value
 
                    distance_to_center_pixels = (x_trans_tot(k_count)-reshape(center_padded_x(pixels_of_interest),[],1)).^2+(y_trans_tot(k_count)-reshape(center_padded_y(pixels_of_interest),[],1)).^2;
                    min_distance = min(distance_to_center_pixels);
 
                    contributing_pixels = [contributing_pixels, distance_to_center_pixels==min_distance];
                    indices = find(contributing_pixels(:,end) == 1);
                    if length(indices) == 4 
                        if x_trans_tot(k_count) == 0 && y_trans_tot(k_count) == 0
                            indices = sort(indices);
                            count_0_vertex = count_0_vertex + 1;
                            switch count_0_vertex
                                case 1
                                    if (0 < theta) && (theta < pi/2)
                                       indices = indices([1,2]);
                                    elseif (pi/2 < theta) && (theta < pi)
                                       indices = indices([2,4]);
                                    elseif (pi < theta) && (theta < 3*pi/2)
                                        indices = indices([3,4]);
                                    elseif (3*pi/2 < theta) && (theta < 2*pi)
                                        indices = indices([1,3]);
                                    end
                                case 2
                                    if (0 < theta) && (theta < pi/2)
                                       indices = indices([2,4]);
                                    elseif (pi/2 < theta) && (theta < pi)
                                       indices = indices([3,4]);
                                    elseif (pi < theta) && (theta < 3*pi/2)
                                        indices = indices([1,3]);
                                    elseif (3*pi/2 < theta) && (theta < 2*pi)
                                        indices = indices([1,2]);
                                    end
                                case 3
                                    if (0 < theta) && (theta < pi/2)
                                       indices = indices([1,3]);
                                    elseif (pi/2 < theta) && (theta < pi)
                                       indices = indices([1,2]);
                                    elseif (pi < theta) && (theta < 3*pi/2)
                                        indices = indices([2,4]);
                                    elseif (3*pi/2 < theta) && (theta < 2*pi)
                                        indices = indices([3,4]);
                                    end
                                case 4
                                    if (0 < theta) && (theta < pi/2)
                                       indices = indices([3,4]);
                                    elseif (pi/2 < theta) && (theta < pi)
                                       indices = indices([1,3]);
                                    elseif (pi < theta) && (theta < 3*pi/2)
                                        indices = indices([1,2]);
                                    elseif (3*pi/2 < theta) && (theta < 2*pi)
                                        indices = indices([2,4]);
                                    end
                            end
                            contributing_pixels(:,end) = zeros(9,1);
                            contributing_pixels(indices,end) = 1;
                        else
                             if ismember(1,indices)
                                if vert_logical_grid(1) == 1
                                     indices = [1 2 2 5];
                                     contributing_pixels(4,end) = 0;
                                else
                                    indices = [2 5 4 5]';
                                    contributing_pixels(1,end) = 0;
                                 end
                             elseif ismember(3,indices)
                                 if vert_logical_grid(3) == 1
                                     indices = [2 3 2 5];
                                     contributing_pixels(6,end) = 0;                                 
                                 else
                                     indices = [2 5 6 5]';
                                     contributing_pixels(3,end) = 0;
                                 end
                             elseif ismember(7,indices)
                                 if vert_logical_grid(7) == 1
                                     indices = [7 8 5 8];
                                     contributing_pixels(4,end) = 0;                                 
                                 else
                                     indices = [4 5 8 5]';
                                     contributing_pixels(7,end) = 0;    
                                 end
                             elseif ismember(9,indices)
                                 if vert_logical_grid(9) == 1
                                     indices = [8 9 5 8];
                                     contributing_pixels(6,end) = 0;                                 
                                 else   
                                     indices = [6 5 8 5]';
                                     contributing_pixels(9,end) = 0;
                                 end
                             end
                        end
                    end
 
                    transitions_to_other{indices(1)} = [transitions_to_other{indices(1)}, indices(2)];
                    transitions_to_other{indices(2)} = [transitions_to_other{indices(2)}, indices(1)];
 
                    transitions_grid{indices(1)} = [transitions_grid{indices(1)}, [x_trans_tot(k_count); y_trans_tot(k_count)]];
                    transitions_grid{indices(2)} = [transitions_grid{indices(2)}, [x_trans_tot(k_count); y_trans_tot(k_count)]];
 
                    if length(indices) == 4                   
                        transitions_to_other{indices(3)} = [transitions_to_other{indices(3)}, indices(4)];
                        transitions_to_other{indices(4)} = [transitions_to_other{indices(4)}, indices(3)];
 
                        transitions_grid{indices(3)} = [transitions_grid{indices(3)}, [x_trans_tot(k_count); y_trans_tot(k_count)]];
                        transitions_grid{indices(4)} = [transitions_grid{indices(4)}, [x_trans_tot(k_count); y_trans_tot(k_count)]];
 
                        k_count_on_vertex = [k_count_on_vertex, k_count];
                    end
                end
                          
                cont_pixels_grid = reshape(sum(contributing_pixels,2) > 0,3,3);
                for kk_count = k_count_on_vertex
                    if sum((x_trans_old == x_trans_tot(kk_count)) & (y_trans_old == y_trans_tot(kk_count))) == 1 
                        x_trans_sub = x_trans_tot;
                        x_trans_sub(kk_count) = Inf;
                        y_trans_sub = y_trans_tot;
                        y_trans_sub(kk_count) = Inf;
                        [~,ind_wrong_trans] = min((x_trans_sub-x_trans_tot(kk_count)).^2+(y_trans_sub-y_trans_tot(kk_count)).^2);
                        x_wrong_trans = x_trans_tot(ind_wrong_trans);
                        y_wrong_trans = y_trans_tot(ind_wrong_trans);
                        for s_count = 1:9
                            if cont_pixels_grid(s_count) ~= 0 
                                ind2_wrong_trans = find(transitions_grid{s_count}(1,:) == x_wrong_trans & transitions_grid{s_count}(2,:) == y_wrong_trans);
                                %ind2_wrong_trans = find(sum(transitions_grid{s_count} == [x_wrong_trans; y_wrong_trans]) == 2);
                                transitions_grid{s_count}(:,ind2_wrong_trans) = [];  
                                transitions_to_other{s_count}(ind2_wrong_trans) = [];
                            end
                        end
                    end
 
                end
 
                for m_count = [1:4, 6:9]
                    if isempty(transitions_to_other{m_count}) == 0              
                        trans_points_this_pixel = transitions_grid{m_count};
                        vert_no_this_pixel = vert_logical_grid(m_count);
                        trans_to_other_this_pixel = transitions_to_other{m_count};
                        vert_location_this_pixel = vertices_location_grid{m_count};
 
                        if vert_no_this_pixel == 0 %no vertex in pixel
                            %v
                            if ismember(m_count,[2,4,6,8])
                                switch m_count
                                    case 2
                                        XY = [trans_points_this_pixel(:,1),[ceil((trans_points_this_pixel(1,1)+correction_x)/length_x)*length_x - correction_x; trans_points_this_pixel(2,1)], [ceil((trans_points_this_pixel(1,2)+correction_x)/length_x)*length_x-correction_x; trans_points_this_pixel(2,2)], trans_points_this_pixel(:,2)];
                                        matr(m_count) = polyarea(XY(1,:),XY(2,:));
                                    case 4
                                        XY = [trans_points_this_pixel(:,1),[trans_points_this_pixel(1,1); floor((trans_points_this_pixel(2,1)+correction_y)/length_y)*length_y-correction_y], [trans_points_this_pixel(1,2); floor((trans_points_this_pixel(2,2)+correction_y)/length_y)*length_y-correction_y], trans_points_this_pixel(:,2)];
                                        matr(m_count) = polyarea(XY(1,:),XY(2,:));
                                    case 6
                                        XY = [trans_points_this_pixel(:,1),[trans_points_this_pixel(1,1); ceil((trans_points_this_pixel(2,1)+correction_x)/length_y)*length_y-correction_x], [trans_points_this_pixel(1,2); ceil((trans_points_this_pixel(2,2)+correction_y)/length_y)*length_y-correction_y], trans_points_this_pixel(:,2)];
                                        matr(m_count) = polyarea(XY(1,:),XY(2,:));
                                    case 8
                                        XY = [trans_points_this_pixel(:,1),[floor((trans_points_this_pixel(1,1)+correction_x)/length_x)*length_x-correction_x; trans_points_this_pixel(2,1)], [floor((trans_points_this_pixel(1,2)+correction_x)/length_x)*length_x-correction_x; trans_points_this_pixel(2,2)], trans_points_this_pixel(:,2)];
                                        matr(m_count) = polyarea(XY(1,:),XY(2,:));
                                end
                            else
                                matr(m_count) = .5*abs(trans_points_this_pixel(1,1)-trans_points_this_pixel(1,2))*abs(trans_points_this_pixel(2,1)-trans_points_this_pixel(2,2));
                            end
                        elseif vert_no_this_pixel == 1 %one vertex in pixel
%                             i_pixels
%                             i_x
%                             i_y
                            if trans_to_other_this_pixel(1) == trans_to_other_this_pixel(2) %the transitions occur on the same boundary
                                matr(m_count) = area_triangle(vert_location_this_pixel, trans_points_this_pixel(:,1),trans_points_this_pixel(:,2));
                            else %two different boundaries
                                triangle1 = area_triangle(vert_location_this_pixel,trans_points_this_pixel(:,1),trans_points_this_pixel(:,2));
                                triangle2 = .5*abs(trans_points_this_pixel(1,1)-trans_points_this_pixel(1,2))*abs(trans_points_this_pixel(2,1)-trans_points_this_pixel(2,2));
                                rectangle = 0;
                                if (rem(trans_points_this_pixel(1,1),length_x) == 0) && (rem(trans_points_this_pixel(1,2),length_x) == 0)
                                    switch m_count
                                        case 4
                                            rectangle = length_x*abs(min(trans_points_this_pixel(2,:))- (floor((min(trans_points_this_pixel(2,:))+correction_y)/length_y)*length_y-correction_y));
                                        case 6
                                            rectangle = length_x*abs(max(trans_points_this_pixel(2,:))- (ceil((max(trans_points_this_pixel(2,:))+correction_y)/length_y)*length_y-correction_y));
                                    end
                                elseif (rem(trans_points_this_pixel(2,1),length_y) == 0) && (rem(trans_points_this_pixel(2,2),length_y) == 0)
                                    switch m_count
                                        case 2
                                            rectangle = length_y*abs(max(trans_points_this_pixel(1,:))- (ceil((max(trans_points_this_pixel(1,:))+correction_x)/length_x)*length_x-correction_x));
                                        case 8
                                            rectangle = length_y*abs(min(trans_points_this_pixel(1,:))- (floor((min(trans_points_this_pixel(1,:))+correction_x)/length_x)*length_x-correction_x));
                                    end
                                end
                                matr(m_count) = triangle1+triangle2+rectangle;
                            end
                        elseif vert_no_this_pixel == 2 %two vertices in one pixel
                            if trans_to_other_this_pixel(1) == trans_to_other_this_pixel(2) %same boundary
                                [~,closest_to_vert_1] = min([norm(vert_location_this_pixel(:,1)-trans_points_this_pixel(:,1)),norm(vert_location_this_pixel(:,1)-trans_points_this_pixel(:,2))]);
            %                     triangle1 = area_triangle(vert_this_pixel(:,1),transitions_grid(:,furthest_from_vert_1),vert_this_pixel(:,2));
            %                     triangle2 = area_triangle(vert_this_pixel(:,1),transitions_grid(:,1),transitions_grid(:,2));
                                XY = [vert_location_this_pixel(:,2)'; vert_location_this_pixel(:,1)'; trans_points_this_pixel(:,closest_to_vert_1)'; trans_points_this_pixel(:,3-closest_to_vert_1)'];    
                                matr(m_count) = polyarea(XY(:,1), XY(:,2));
                            else
                                other_point = zeros(2,1);
                                padded_pixel_number = pixels_of_interest(m_count);
                                center_of_rotated_pixel_x = cent_rot_x(i_y,i_x);
                                center_of_rotated_pixel_y = cent_rot_y(i_y,i_x);
                                vertices_this_pixel_x = repmat([floor(center_of_rotated_pixel_x), ceil(center_of_rotated_pixel_x)],2,1);
                                vertices_this_pixel_y = repmat([floor(center_of_rotated_pixel_y); ceil(center_of_rotated_pixel_y)],1,2);
                                [~,index_vertex_closest_to_center] = min((reshape(vertices_this_pixel_x,[],1)-center_of_rotated_pixel_x).^2+(reshape(vertices_this_pixel_y,[],1)-center_of_rotated_pixel_y).^2);
                                vertex_closest_to_center = [vertices_this_pixel_x(index_vertex_closest_to_center); vertices_this_pixel_y(index_vertex_closest_to_center)];
                                [~,closest_to_vert_1] = min([norm(vert_location_this_pixel(:,1)-trans_points_this_pixel(:,1)),norm(vert_location_this_pixel(:,1)-trans_points_this_pixel(:,2))]);
                                XY = [vert_location_this_pixel(:,2)'; vert_location_this_pixel(:,1)'; trans_points_this_pixel(:,closest_to_vert_1)'; vertex_closest_to_center'; trans_points_this_pixel(:,3-closest_to_vert_1)']; 
                                matr(m_count) = polyarea(XY(:,1), XY(:,2));
                            end
                        else
                            disp('?');
                        end
                    end
                    if matr(m_count) < 1e-13
                        matr(m_count) = 0;
                    end
                end
                matr = matr/(length_x*length_y);
                matr(5) = 1- sum(sum(matr));
                G(i_pixels,pixels_of_interest) = reshape(matr,1,[]);
 
            end
        end
    end
 
    G = G(:,reshape(real_pixels_in_padded,[],1));
    %%
%     col_counter = 0;
%     if no_of_divisions > 0
%         G_temp = G;
%         for c_count = 1:npy_original*npx_original
%             G(c_count,:) = 1/2^no_of_divisions*sum(G_temp(c_count*2^no_of_divisions-2^no_of_divisions+1:c_count*2^no_of_divisions,:));
%         end
%     elseif no_of_multiplications > 0
%         G_temp = G;
%         for d_count = 1:npx_original*npy_original
%             G(d_count,:) = 1/2^no_of_multiplications*sum(G_temp(d_count+npy*col_counter):npy:(d_count+npy*col_counter+(2^no_of_multiplications-1)*npy),:);
%             if mod(d_count,npy_original) == 0 
%                 col_counter = col_counter+1;
%             end
%         end
%     end
     
end
