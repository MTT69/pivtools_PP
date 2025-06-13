function G1=gamma1(x,y,u,v,d)

    s=size(u);
    SE = strel('disk', d, 0);
    theSE=getnhood(SE);
    [I,J]=find(theSE);
    I=I-d-1;
    J=J-d-1;
    i0=round(s(1)/2); j0=round(s(2)/2); ind0=sub2ind(s,i0,j0);  
    IN=sub2ind(s,I+i0,J+j0);
    
    A=cart2pol(x(IN)-x(ind0), y(IN)-y(ind0));
    
    U=nansurround(u,d);
    V=nansurround(v,d);
    S=size(U);
    
    localS=[];
    for k=1:length(I)
        ANGLE=ones(S(1),S(2)).*A(k);
        [T,~]=cart2pol(circshift(U,[I(k) J(k)]), circshift(V,[I(k) J(k)]));
        localS=cat(3,localS,sin(ANGLE-T));
    end
    G1=-sum(localS,3, 'omitnan')/length(I);
    G1=unsurround(G1,d);