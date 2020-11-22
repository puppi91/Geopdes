% TEST_STOKES_SQUARE_NEUMANN_BC: Neumann bc

function u = test_stokes_square_neumann_bc (x, y, iside)


% syms x y
% u1 =  (sin(x) .* cos(y));
% u2 =  (-sin(y) .* cos(x));
% p = x - 5/6*2/3;
% u1x = diff(u1,x);
% u1y = diff(u1,y);
% u2x = diff(u2,x);
% u2y = diff(u2,y);
% J = [u1x, u1y; u2x, u2y];
% % Neumann b.c. Dun-pn
% % on 1 
% n = [-1,0]';
% J*n-p*n
% % on 2 
% n = [1,0]';
% J*n-p*n
% % on 3 
% n = [0,-1]';
% J*n-p*n
% % on 4 
% n = [0,1]';
% J*n-p*n

uxx = cos(x).*cos(y);
uxy = -sin(x).*sin(y);
uyx = sin(x).*sin(y);
uyy = -cos(x).*cos(y);

p = x - 5/6*2/3;

switch(iside)
    case 1
        u1 = -uxx+p;
        u2 = -uyx;
        u = cat(1, ...
            reshape (u1, [1, size(x)]), ...
            reshape (u2, [1, size(x)]));
    case 2
        u1 = uxx-p;
        u2 = uyx;
        u = cat(1, ...
            reshape (u1, [1, size(x)]), ...
            reshape (u2, [1, size(x)]));
    case 3
        u1 = -uxy;
        u2 = -uyy+p;
        u = cat(1, ...
            reshape (u1, [1, size(x)]), ...
            reshape (u2, [1, size(x)])); 
    case 4
        u1 = uxy;
        u2 = uyy-p;
        u = cat(1, ...
            reshape (u1, [1, size(x)]), ...
            reshape (u2, [1, size(x)]));
end

end


