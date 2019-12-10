#plotting the movement of leapfrogging vortex rings
#author: Melanie Markham
#Date: 09/23/2019
#g=gamma, variable for strength of vortex 
#r is distance from one vortex to the other
#rv is the vector form of r
using Plots
using LinearAlgebra
function velocity(g,rv,r)
    cross(g,rv)/(2*pi*r^2)
end
function position(v,t)
    v*t
end
function r(x2,x1,y2,y1)
    sqrt((x2-x1)^2+(y2-y1)^2)
end
function rv(x2,x1,y2,y1)
    [x2-x1,y2-y1,0]
end
g=[0.0;0.0;1.0];
t=collect(0.0:0.01:40.0);
p1=zeros(3,length(t));
p2=zeros(3,length(t));
p3=zeros(3,length(t));
p4=zeros(3,length(t));
v1=zeros(3,length(t));
v2=zeros(3,length(t));
v3=zeros(3,length(t));
v4=zeros(3,length(t));
r12=zeros(1,length(t));
r13=zeros(1,length(t));
r14=zeros(1,length(t));
r23=zeros(1,length(t));
r24=zeros(1,length(t));
r34=zeros(1,length(t));
r12v=zeros(3,length(t));
r13v=zeros(3,length(t));
r14v=zeros(3,length(t));
r23v=zeros(3,length(t));
r24v=zeros(3,length(t));
r34v=zeros(3,length(t));

v1[:,1]=[0.0;0.0;0.0];
p1[:,1]=[0.0;-0.5;0.0];
v2[:,1]=[0.0;0.0;0.0];
p2[:,1]=[0.0;0.5;0.0];
v3[:,1]=[0.0;0.0;0.0];
p3[:,1]=[1.0;0.5;0.0];
v4[:,1]=[0.0;0.0;0.0];
p4[:,1]=[1.0;-0.5;0.0];
r12[1]=r(p2[1,1],p1[1,1],p2[2,1],p1[2,1]);
r13[1]=r(p3[1,1],p1[1,1],p3[2,1],p1[2,1]);
r14[1]=r(p4[1,1],p1[1,1],p4[2,1],p1[2,1]);
r23[1]=r(p3[1,1],p2[1,1],p3[2,1],p2[2,1]);
r24[1]=r(p4[1,1],p2[1,1],p4[2,1],p2[2,1]);
r34[1]=r(p4[1,1],p3[1,1],p4[2,1],p3[2,1]);
r12v[:,1]=rv(p2[1,1],p1[1,1],p2[2,1],p1[2,1]);
r13v[:,1]=rv(p3[1,1],p1[1,1],p3[2,1],p1[2,1]);
r14v[:,1]=rv(p4[1,1],p1[1,1],p4[2,1],p1[2,1]);
r23v[:,1]=rv(p3[1,1],p2[1,1],p3[2,1],p2[2,1]);
r24v[:,1]=rv(p4[1,1],p2[1,1],p4[2,1],p2[2,1]);
r34v[:,1]=rv(p4[1,1],p3[1,1],p4[2,1],p3[2,1]);

for i=2:length(t)
    v1[:,i]=velocity(-g,r12v[:,i-1],r12[i-1])+velocity(-g,r13v[:,i-1],r13[i-1])+velocity(g,r14v[:,i-1],r14[i-1])
    p1[:,i]=p1[:,i-1]+position(v1[:,i],(t[i]-t[i-1]))
    v2[:,i]=velocity(g,-r12v[:,i-1],r12[i-1])+velocity(-g,r23v[:,i-1],r23[i-1])+velocity(g,r24v[:,i-1],r24[i-1])
    p2[:,i]=p2[:,i-1]+position(v2[:,i],(t[i]-t[i-1]))
    v3[:,i]=velocity(g,-r13v[:,i-1],r13[i-1])+velocity(-g,-r23v[:,i-1],r23[i-1])+velocity(g,r34v[:,i-1],r34[i-1])
    p3[:,i]=p3[:,i-1]+position(v3[:,i],(t[i]-t[i-1]))
    v4[:,i]=velocity(g,-r14v[:,i-1],r14[i-1])+velocity(-g,-r24v[:,i-1],r24[i-1])+velocity(-g,-r34v[:,i-1],r34[i-1])
    p4[:,i]=p4[:,i-1]+position(v4[:,i],(t[i]-t[i-1]))
    r12[i]=r(p2[1,i],p1[1,i],p2[2,i],p1[2,i])
    r13[i]=r(p3[1,i],p1[1,i],p3[2,i],p1[2,i])
    r14[i]=r(p4[1,i],p1[1,i],p4[2,i],p1[2,i])
    r23[i]=r(p3[1,i],p2[1,i],p3[2,i],p2[2,i])
    r24[i]=r(p4[1,i],p2[1,i],p4[2,i],p2[2,i])
    r34[i]=r(p4[1,i],p3[1,i],p4[2,i],p3[2,i])
    r12v[:,i]=rv(p2[1,i],p1[1,i],p2[2,i],p1[2,i])
    r13v[:,i]=rv(p3[1,i],p1[1,i],p3[2,i],p1[2,i])
    r14v[:,i]=rv(p4[1,i],p1[1,i],p4[2,i],p1[2,i])
    r23v[:,i]=rv(p3[1,i],p2[1,i],p3[2,i],p2[2,i])
    r24v[:,i]=rv(p4[1,i],p2[1,i],p4[2,i],p2[2,i])
    r34v[:,i]=rv(p4[1,i],p3[1,i],p4[2,i],p3[2,i])
end

plot(p1[1,:],p1[2,:],label="p1",linecolor=:blue,ylims=(-2,2),yticks=0,xticks=0:2.5:13)
plot!(p2[1,:],p2[2,:],label="p2",linecolor=:blue)
plot!(p3[1,:],p3[2,:],label="p3",linecolor=:red)
plot!(p4[1,:],p4[2,:],label="p4",linecolor=:red)
