clear all
S = diag([6 0.5 0.4 0.2 0.1 0.1 0.1]);

for i = 1:100
   x = mvnrnd(zeros(size(S,1),1),S,100);
   y = mvnrnd(zeros(size(S,1),1),S,100);
   
   %[s(i),a(:,i)]=dim.krzsim(x,y,3);
   p(i) = dim.krztest(x,y,5);
end
