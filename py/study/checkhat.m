% CHECKHAT  Check my understanding of hat functions on two levels.

x = 0:.001:1;

% build lambda on the coarse mesh
lc = zeros(size(x));
m1 = (x>0.25) & (x<0.5);
m2 = (x>=0.5) & (x<0.75);
lc(m1) = 4.0 * (x(m1)-0.25);
lc(m2) = 4.0 * (0.75-x(m2));

% build three lambdas on fine mesh
lf = zeros(3,length(x));
x8 = 0:.125:1;
for j = 0:2
    m1 = (x>x8(j+3)) & (x<x8(j+4));
    m2 = (x>=x8(j+4)) & (x<x8(j+5));
    lf(j+1,m1) = 8.0 * (x(m1)-x8(j+3));
    lf(j+1,m2) = 8.0 * (x8(j+5)-x(m2));
end
lcomb = 0.5*lf(1,:) + lf(2,:) + 0.5*lf(3,:);
norm(lc-lcomb)

% plot
subplot(2,1,1)
plot(x,lc)
title('coarse-mesh hat')
subplot(2,1,2)
plot(x,lcomb)
title('combination of fine-mesh hats')

