S=[0 0; 0 1; 0 2; 0 3; 1 0; 1 1; 1 2; 1 3; 2 0; 2 1; 2 2; 2 3; 3 0; 3 1; 3 2; 3 3];
N=size(S,1);
for j=1:N
    Y(j)=S(j,1)^2+exp(S(j,2));
end
Y=Y';
theta=[10 10];
lob=[0.1 0.1];
upb=[20 20];
[dmodel,perf]=dacefit(S,Y,@regpoly0,@corrgauss,theta,lob,upb);
X=gridsamp([0 0; 3 3],10);
[YX MSE]=predictor(X,dmodel);
X1=reshape(X(:,1),10,10);
X2=reshape(X(:,2),10,10);
YX=reshape(YX,size(X1));
figure(1);
mesh(X1,X2,YX);
hold on;
plot3(S(:,1),S(:,2),Y,'.k','MarkerSize',10);
hold off;