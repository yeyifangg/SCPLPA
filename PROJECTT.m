function [score]=PROJECTT(SIM,LD)


 projectT=SIM*LD; 
 [nnc,nd]=size(LD);
      for i=1:nd
  p=LD(:,i);
   if(norm(p)~=0)

    projectT(:,i)=projectT(:,i)/norm(p);  

   end
      end 
score=projectT;

end


