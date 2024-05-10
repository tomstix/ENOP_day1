
% Assignment 7 Excercise 1
clear all;

it=20;      % Number of iterations for each run 
ftest=[];  % Best function output values
Select=3;   % Selection of the function cases
L_it=20;    % Number of algorithm run cycles
switch Select 

    case 1      % Rosenbrock for n=2
    fr=@rosenbrock;
    fr_1=@(x,y) rosenbrock([x y]);

    for i=1:L_it
       %ftest(i)=MRGS(fun,n,x_min,x_max,y_min,y_max,z_min,z_max,iterations)
        ftest(i)=MRGS(fr,2,-2,2,-2,2,0,0,it);
    end
    figure;
    fsurf(fr_1,[-2,2,-2,2]); 
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    title('Rosenbrock Surface Plot')
    figure;
    fcontour(fr_1,[-2,2,-2,2],'LevelStep',7);
    xlabel('x-axis');
    ylabel('y-axis');
    title('Rosenbrock Contour with start,optim., and best points')
    hold on;
       ftest(1)=MRGS(fr,2,-2,2,-2,2,0,0,it);
    hold off;

    case 2      % f_p for n=1
    fp=@f_p;
    for i=1:L_it
        ftest(i)=MRGS(fp,1,-2,2,0,0,0,0,it);
    end   

    case 3      % f_p for n=2
    fp=@f_p;
    fp_1=@(x,y) f_p([x y]);
    for i=1:L_it
        ftest(i)=MRGS(fp,2,-2,2,-2,2,0,0,it);
    end

    figure;
    fsurf(fp_1,[-2,2,-2,2]);
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    title('fp Surface Plot');

    figure;
    fcontour(fp_1,[-2,2,-2,2]);
    xlabel('x-axis');
    ylabel('y-axis');
    title('fp countour with start, optim., and best points');
    hold on;
        ftest(1)=MRGS(fp,2,-2,2,-2,2,0,0,it);
    hold off;

    case 4      % fa for n=1
    fa=@auckley;
    fa_1=@(x,y) auckley([x,y]);
    for i=1:L_it
        ftest(i)=MRGS(fa,1,-10,10,0,0,0,0,it);
    end
 
    case 5      % fa for n=2
    fa=@auckley;
    fa_1=@(x,y) auckley([x,y]);
    for i=1:L_it
        ftest(i)=MRGS(fa,2,-10,10,-10,10,0,0,it);
    end
   
    figure;
    fsurf(fa_1,[-10,10,-10,10])
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    title('Auckley Surface Plot')   
       
    figure;
    fcontour(fa_1,[-10,10,-10,10]);
    xlabel('x-axis');
    ylabel('y-axis');
    title('Auckley countour with start, optim., and best points');
    hold on;
    ftest(1)=MRGS(fa,2,-10,10,-10,10,0,0,it);
    hold off;
    
    case 6      % fa for n=3
    fa=@auckley;
    for i=1:L_it
        ftest(i)=MRGS(fa,3,-10,10,-10,10,-10,10,it);
    end
    
end
   disp("Max:    "+max(ftest))  % Display results
   disp("Min:    "+min(ftest))
   disp("Average:"+mean(ftest))
   disp("STEDV:  "+std(ftest))

%function MRGS
function f_best = MRGS(f,n,x_min,x_max,y_min,y_max,z_min,z_max,k_max)
Maxfunev=50;
options = optimset('MaxFunEvals',Maxfunev,'Display','off'); %
k=0; 
f_best=inf;
 
while k<k_max
    x_new=GSP(n,x_min,x_max,y_min,y_max,z_min,z_max);   % Generates a random a starting point

    if n==1
        [x_new,f_new]=fmincon(f,x_new,[],[],[],[],[x_min]',[x_max]','',options);
    end

    if n==2
       plot(x_new(1),x_new(2),'o','MarkerSize',5,'MarkerFaceColor', 'k');%,LineWidth=3) %plots the starting point. 
       hold on;
       [x_new,f_new]=fmincon(f,x_new,[],[],[],[],[x_min y_min]',[x_max y_max]','',options);
       plot(x_new(1),x_new(2),'o','MarkerSize',5,'MarkerFaceColor', 'm');%,LineWidth=3) %plots the new resultant point)
    end

    if n==3
        [x_new,f_new]=fmincon(f,x_new,[],[],[],[],[x_min y_min z_min]',[x_max y_max z_max]','',options);
    end
    if f_new<f_best
        x_best=x_new;
        f_best=f_new;
        
        if n==2
        plot(x_new(1),x_new(2),'x','MarkerSize',10,'LineWidth',2,'Color','r');  % Adds the new best point       
        end
        
    end
      k=k+1;    % Updates the k
      disp("Run Results_________________________")
      disp("Iteration No:   "+(k));
      disp("Best Fun value: "+(f_best));
      disp("Fun evaluations: "+(k*Maxfunev));
end

end

%Function to generate a random starting point
function x = GSP(n, x_min, x_max, y_min, y_max, z_min, z_max)
    x = zeros(1, n);

    % Define the ranges for each of the dimensions
    ranges = {x_min, x_max; y_min, y_max; z_min, z_max};

    for i = 1:n
        minVal = ranges{i, 1};  % Min
        maxVal = ranges{i, 2};  % Max
        x(i) = rand(1) * (maxVal - minVal) + minVal;  % Compute the coordinate
    end
end


%Rosenbrock function
function fr=rosenbrock(x)
%n=length(x);
fr = (1-x(1)).^2+100*(x(2)-x(1).^2).^2;
end

%Function fp
function fp=f_p(x)
n=length(x);
if n==1
    fp=-exp(-((norm(x)).^2)/2)*cos(10.*x);
elseif n==2
    fp=-exp(-((norm(x)).^2)/2)*prod(cos(10.*x));
end
end

%Auckley function
function fa = auckley(x)
    n = length(x);
    fa = -20*exp( -0.2*sqrt((1/n) * sum(x.^2)) ) - exp((1/n) * sum(cos(2*pi.*x))) + 21;
end

