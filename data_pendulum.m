function [dy]= rigid (t,y)     
L1=input('Input the length to the first point mass: ');   %this is the length from the pivot to the first point mass
L2=input('Input the length to the second point mass: ');  %This is the length from the first point mass to the second one
Theta1=input('Input the angle from the vector <0,-1> for mass 1 to start: ');  %This is the angle between L1 and the vector <0,-1>
Theta2=input('Input the angle from the vector <0,-1> for mass 2 to start: ');     %This is the angle between L2 and the vector <0,-1>
w1=input('Input the angular velocity for mass 1 to start: ');
w2=input('Input the angular velocity for mass 1 to start: ');
m1=input('Input the mass of the point mass 1: ');        %Mass of 1st point mass
m2=input('Input the mass of the point mass 2: ');         %Mass of 2nd point mass
g=9.8;  %This is the value of gravity 
Initial_values=[Theta1; w1; Theta2; w2]; %set intial values for ode45 and sets x(1)=Theta1 x(2)=w1 x(3)=Theta2 x(4)=w2 to be passed to ode45
    function x_solution = subfunction(t,x)
    x_solution=zeros(4,1);
    x_solution(1)=x(2); %derivative of Theta1 which is w1
    x_solution(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
    x_solution(3)=x(4); %derivative of Theta2 which is w2
    x_solution(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
    end
[t,x]=ode45(@subfunction,[0 20],Initial_values) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
figure
plot(t,x(:,1))
title('Theta1 vs time')
xlabel('Time')
ylabel('Angle')
figure
plot(t,x(:,3))
title('Theta2 vs time')
xlabel('Time')
ylabel('Angle')
figure
plot(t,x(:,2))
title('w1 vs time')
xlabel('Time')
ylabel('Angular Velocity')
figure
plot(t,x(:,4))
title('w2 vs time')
xlabel('Time')
ylabel('Angular Velocity')
figure 
plot(x(:,1),x(:,3))
title('Theta1 vs Theta2')
xlabel('Theta1 Angle')
ylabel('Theta2 Angle')
x1=L1*sin(x(:,1));           %this and the next couple lines takes thetas and angular velocity and converts them to cartesian coordinates
y1=-L1*cos(x(:,1));
x2=L2*sin(x(:,3))+L1*sin(x(:,1));
y2=-L2*cos(x(:,3))-L1*cos(x(:,1));
x0=zeros(569);
y0=zeros(569);
figure
hold on
for i=1:569                                 %after running this I found 569 to be the length of the solution vector thus 569 is used
    plot([x0(i) x1(i)],[y0(i) y1(i)],'r')    %still cant get this to not show the previous data
    plot([x1(i) x2(i)],[y1(i) y2(i)],'b')     %This is looping through the rows of the solution vector and plotting them somewhat real time
    plot(x1(i),y1(i),'.k','MarkerSize',30)    %the 2 above plots, plot the vectors connecting pivots to points and this line and the one below plot the point masses with large ball markers
    plot(x2(i),y2(i),'.k','MarkerSize',30)
    pause(.00351)                             %pause gives animated effect
end
hold off
%run the program again except for add small random number to the initial
%angular velocity of point mass 2
Input_values2=[Theta1; w1; Theta2; w2+rand];
[t,z]=ode45(@subfunction,[0 20],Input_values2)
x1_2=L1*sin(z(:,1));                    %this and the next couple lines takes thetas and angular velocity and converts them to cartesian coordinates
y1_2=-L1*cos(z(:,1));
x2_2=L2*sin(z(:,3))+L1*sin(z(:,1));
y2_2=-L2*cos(z(:,3))-L1*cos(z(:,1));
x0=zeros(569);  %setting up origin for the plotting
y0=zeros(569);
figure
hold on
for i=1:569                                      %after running this I found 569 to be the length of the solution vector thus 569 is used
    plot([x0(i) x1_2(i)],[y0(i) y1_2(i)],'r')    %still cant get this to not show the previous data
    plot([x1_2(i) x2_2(i)],[y1_2(i) y2_2(i)],'b')   %This is looping through the rows of the solution vector and plotting them somewhat real time
    plot(x1_2(i),y1_2(i),'.k','MarkerSize',30)      %the 2 above plots, plot the vectors connecting pivots to points and this line and the one below plot the point masses with large ball markers
    plot(x2_2(i),y2_2(i),'.k','MarkerSize',30)
    pause(.00351)                                 %pause just allows this to seem animated.
end
hold off
figure
comet(x2,y2)    %comet functions plot "animated" coordinates of 2nd point mass
figure
comet(x2_2,y2_2)   %same as previous but is for the slightly adjusted one (added rand to w2)


%Bifurcation
% takes a ton of lines but is really the same thing repeated over and over,
% since I couldnt put the function inside a loop
%This is the bifurcation for a system with m1 m2 L1 L2 all equalling 1 and
%the initial conditions theta1 w1 and w2 all set to 0 with varying theta2
%to start
 %set intial values for ode45 and sets x(1)=Theta1 x(2)=w1 x(3)=Theta2 x(4)=w2 to be passed to ode45

 g=9.8
 m1=1
 m2=1
 L1=1
 L2=1
 %note the initial values of [0; 0; 0; 0] are ommited as they yield zero,
 %they are reintroduced at the end during plotting
        Initial_values2=[0; 0; .2; 0];
     function x_solution2 = subfunction2(t,x);
        x_solution2=zeros(4,1);
        x_solution2(1)=x(2); %derivative of Theta1 which is w1
        x_solution2(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution2(3)=x(4); %derivative of Theta2 which is w2
        x_solution2(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x2]=ode45(@subfunction,[0 20],Initial_values2) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
    Initial_values3=[0; 0; .4; 0];
     function x_solution3 = subfunction3(t,x);

        x_solution3=zeros(4,1);
        x_solution3(1)=x(2); %derivative of Theta1 which is w1
        x_solution3(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution3(3)=x(4); %derivative of Theta2 which is w2
        x_solution3(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x3]=ode45(@subfunction,[0 20],Initial_values3) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2

     Initial_values4=[0; 0; .6; 0];
     function x_solution4 = subfunction4(t,x);
    
        x_solution4=zeros(4,1);
        x_solution4(1)=x(2); %derivative of Theta1 which is w1
        x_solution4(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution4(3)=x(4); %derivative of Theta2 which is w2
        x_solution4(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x4]=ode45(@subfunction,[0 20],Initial_values4) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
Initial_values5=[0; 0; .8; 0];
     function x_solution5 = subfunction5(t,x);
    
        x_solution5=zeros(4,1);
        x_solution5(1)=x(2); %derivative of Theta1 which is w1
        x_solution5(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution5(3)=x(4); %derivative of Theta2 which is w2
        x_solution5(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x5]=ode45(@subfunction,[0 20],Initial_values5) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
Initial_values6=[0; 0; 1; 0];
     function x_solution6 = subfunction6(t,x);
    
        x_solution6=zeros(4,1);
        x_solution6(1)=x(2); %derivative of Theta1 which is w1
        x_solution6(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution6(3)=x(4); %derivative of Theta2 which is w2
        x_solution6(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x6]=ode45(@subfunction,[0 20],Initial_values6) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
Initial_values7=[0; 0; 1.2; 0];
     function x_solution7 = subfunction7(t,x);
    
        x_solution7=zeros(4,1);
        x_solution7(1)=x(2); %derivative of Theta1 which is w1
        x_solution7(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution7(3)=x(4); %derivative of Theta2 which is w2
        x_solution7(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x7]=ode45(@subfunction,[0 20],Initial_values7) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
 Initial_values8=[0; 0; 1.4; 0];
     function x_solution8 = subfunction8(t,x);
   
        x_solution8=zeros(4,1);
        x_solution8(1)=x(2); %derivative of Theta1 which is w1
        x_solution8(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution8(3)=x(4); %derivative of Theta2 which is w2
        x_solution8(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x8]=ode45(@subfunction,[0 20],Initial_values8) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2
 Initial_values9=[0; 0; 1.6; 0];
     function x_solution9 = subfunction9(t,x);
   
        x_solution9=zeros(4,1);
        x_solution9(1)=x(2); %derivative of Theta1 which is w1
        x_solution9(2)=(-g*(2*m1+m2)*sin(x(1))-m2*g*sin(x(1)-2*x(3))-2*sin(x(1)-x(3))*m2*(x(4)^2*L2+x(2)^2*L1*cos(x(1)-x(3))))/(L1*(2*m1+m2-m2*cos(2*x(1)-2*x(3)))); %derivative of w1 or 2nd derivative of Theta1 expressed in terms of the 1st derivative of theta1 
        x_solution9(3)=x(4); %derivative of Theta2 which is w2
        x_solution9(4)=(2*sin(x(1)-x(3))*(x(2)^2*L1*(m1+m2)+g*(m1+m2)*cos(x(1))+x(4)^2*L2*m2*cos(x(1)-x(3))))/(L2*(2*m1+m2-m2*cos(2*x(1)-2*x(3))))
end
    [t,x9]=ode45(@subfunction,[0 20],Initial_values9) %solves the system of differential equations x(:,1)=theta1 ,2=w1  ,3=theta2  ,4=w2

%b1 will be zeros anyways %the rest of this puts all the solutions to
%theta2 into one big vector to be plotted. This is how the previous
%bifurcaiton diagram was constructed (the one we went over in class)
b2=x2(:,3)'
b3=x3(:,3)'
b4=x4(:,3)'
b5=x5(:,3)'
b6=x6(:,3)'
b7=x7(:,3)'
b8=x8(:,3)'
b9=x9(:,3)'
Theta2=(0:.2:1.6)
b_total=[zeros(1,300);b2(1,1:300);b3(1,1:300);b4(1,1:300);b5(1,1:300);b6(1,1:300);b7(1,1:300);b8(1,1:300);b9(1,1:300)]
figure
plot(Theta2, b_total, '-k', 'markersize', 100);
grid on;
end
