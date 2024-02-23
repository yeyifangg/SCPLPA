function F = LPA(SIMdd, SIMnnc, Y,gama,beta)
% SIMdd£∫the disease-disease similarity matrix
% SIMnnc: the miRNA-miRNA similarity matrix
% Y : the ground truth (the known disease-miRNA associations)
% F : the result predicted by our method
%gama,beta 0.1-0.9
[nnc,nd]=size(Y);
%  πÈ“ªªØ
M=sum(SIMnnc);
for i=1:nnc
    for j=1:nnc
        SIMnnc(i,j)=SIMnnc(i,j)/(((M(i)*M(j))^0.5));
    end
end

D=sum(SIMdd);
for i=1:nd
    for j=1:nd
        SIMdd(i,j)=SIMdd(i,j)/(((D(i)*D(j))^0.5));
    end
end

%prediction from miRNA space
if nargin < 3
    gama=0.95;
    beta = 0.6; 
end


% PT=Y';
% PT0=Y;
% P0=Y';

PT=Y;
PT0=Y';
P0=Y;


k= 0;
delta = 1;

while  (delta > 1e-6)
    PT1 = (1-gama)*SIMnnc*PT+gama*P0;
    delta =abs(sum(sum((abs(PT1)-abs(PT)))));
    PT = PT1;
    k= k + 1;
end

%prediction from Disease space
delta = 1;
while  (delta > 1e-6)
    DD = (1-gama)*SIMdd*PT0+gama*P0';
    delta =abs(sum(sum((abs(DD)-abs(PT0)))));
    PT0 =DD;
    k= k + 1;
end
F=(beta*PT1+(1-beta)*DD');
end


