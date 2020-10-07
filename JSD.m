function JSD_S=JSD(S)
S=S+1;
X1=(1:2)';
JSD_S=zeros(size(S));
for i=1:size(S,1)
    ideal=zeros(2,1)+1e-10;
    ideal(1)=1;
    for j=1:size(S,2)
        c=S(:,j);
        if max(c)<2 || S(i,j)<1
            JSD_S(i,j)=1;
            continue
        end
        a=prctile(S([1:i-1, i+1:size(S,1)],j),75);
        b=[ S(i,j);a];
        observed=b./sum(b)+eps;
        M=(observed+ideal)./2+eps;
        JSD_S(i,j)=1/2*kldiv(X1,ideal,M)+1/2*kldiv(X1,observed,M);
     end
end    


