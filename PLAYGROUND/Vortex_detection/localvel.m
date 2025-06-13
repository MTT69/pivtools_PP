function Vp=localvel(u,v,d)
    s=size(u);

    SE = strel('disk', d, 0);
    theSE=getnhood(SE);
    [I,J]=find(theSE);
    I=I-d-1;
    J=J-d-1;
    
    localU=[];
    localV=[];
    for k=1:length(I)
        localU=cat(3,localU,circshift(u,[I(k) J(k)]));
        localV=cat(3,localV,circshift(v,[I(k) J(k)]));
    end
    Vp=cat(3,nanmean(localU,3),nanmean(localV,3));
    
    ROI=ones(s(1),s(2));
    IND=(isnan(u) & isnan(v));
    ROI(IND)=nan;
    ROI([1:d,end-d:end],:)=nan; ROI(:,[1:d,end-d:end])=nan;
    ROI=cat(3,ROI,ROI);
    Vp=Vp.*ROI;