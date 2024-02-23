function sim=getSimilarityDisease(interaction,DD)
%用在global local验证时高斯相似性的重计算，testi可取1:5340,表示knowndiseasemirnainteraction.txt中的第testi行association视作0
[nd, nl] = size(interaction);
sim=zeros(nd, nd);
for i=1:nd
sd(i)=norm(interaction(i,:))^2;
   end
    gamad=nd/sum(sd');
    
    %calculate Gaussian kernel for the similarity between disease: kd
    for i=1:nd
        for j=1:nd
    kd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
       end
    end

  
     for i=1:nd
     for j=1:nd
        
         
%          if (DD(i,j)==0)
%              sim(i,j)=kd(i,j);
%         else
%              sim(i,j)=(kd(i,j)+DD(i,j))/2;
%          end
%        
        
        sim(i,j)=(DD(i,j)+kd(i,j))/2;
                     
     end
 end
    
    

end
