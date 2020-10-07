function z=distance(x1,x2,theda,T)
y=2.^(abs(log2((x1+1)./(x2+1))));
z=1./(1+exp(-theda*(y-T)));


