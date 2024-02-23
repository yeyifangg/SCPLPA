function sim=getSimilarityRNA(interaction,MM)
%用在global local验证时高斯相似性的重计算，testi可取1:5340,表示knowndiseasemirnainteraction.txt中的第testi行association视作0
[nd, nl] = size(interaction);
sim=zeros(nl, nl);
      %calculate gamam for Gaussian kernel calculation
for i=1:nl
sm(i)=norm(interaction(:,i))^2;
   end
    gamam=nl/sum(sm');

    %calculate Gaussian kernel for the similarity between miRNA: km
    for i=1:nl
        for j=1:nl
    km(i,j)=exp(-gamam*(norm(interaction(:,i)-interaction(:,j)))^2);
       end
    end

 for i=1:nl
     for j=1:nl
        
         
         if (MM(i,j)==0)
             sim(i,j)=km(i,j);
        else
             sim(i,j)=MM(i,j);
         end
       
        
         %sim(i,j)=(MM(i,j)+km(i,j))/2;
                     
     end
 end
    
    
    
end
