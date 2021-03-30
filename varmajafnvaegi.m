function V=varmajafnvaegi(a,b,beta1,beta2,h);

%hja okkur er q=0

xb=0;
yb=0;

N=a/h;
M=b/h;

u=a/N;
k=b/M;%ma orugglega fara

%r=@(x,y) 0;

x=xb+(0:N)*u;
y=yb+(0:M)*k;

A=zeros((N+1)*(M+1),(N+1)*(M+1));
b_vigur=zeros((N+1)*(M+1),1);

%innri punktar
for i=2:N
    for j=2:M
        %r_sum=r(x(i)-2*h/3,y(j)-h/3)+r(x(i)-h/3,y(j)-2*h/3)+r(x(i)+h/3,y(j)-h/3);
        %r_sum=r_sum+r(x(i)+2*h/3,y(j)+h/3)+r(x(i)+h/3,y(j)+2*h/3)+r(x(i)-h/3,y(j)+h/3);
        
        pkt=i+(j-1)*(N+1);
        
        %A(pkt,pkt)=4-h^2*r_sum/18;
        %A(pkt,pkt-1)=-1-h^2*(r(x(i)-h/3,y(j)+h/3)+r(x(i)-2*h/3,y(j)-h/3))/18;
        %A(pkt,i-1+(j-2)*(N+1))=-h^2*(r(x(i)-2*h/3,y(j)-h/3)+r(x(i)-h/3,y(j)-2*h/3))/18;
        %A(pkt,i+(j-2)*(N+1))=-1-h^2*(r(x(i)-h/3,y(j)-2*h/3)+r(x(i)+h/3,y(j)-h/3))/18;
        %A(pkt,pkt+1)=-1-h^2*(r(x(i)+h/3,y(j)-h/3)+r(x(i)+2*h/3,y(j)+h/3))/18;
        %A(pkt,i+1+j*(N+1))=-h^2*(r(x(i)+2*h/3,y(j)+h/3)+r(x(i)+h/3,y(j)+2*h/3))/18;
        %A(pkt,i+j*(N+1))=-1-h^2*(r(x(i)+h/3,y(j)+2*h/3)+r(x(i)-h/3,y(j)+h/3))/18;
        
        A(pkt,pkt)=2*(u^2+k^2)/(u*k);
        A(pkt,pkt-1)=-k/u;
        A(pkt,i-1+(j-2)*(N+1))=0;
        A(pkt,i+(j-2)*(N+1))=-u/k;
        A(pkt,pkt+1)=-k/u;
        A(pkt,i+1+j*(N+1))=0;
        A(pkt,i+j*(N+1))=-u/k;
        
        %f_sum=f(x(i)-2*h/3,y(j)-h/3)+f(x(i)-h/3,y(j)-2*h/3)+f(x(i)+h/3,y(j)-h/3);
        %f_sum=f_sum+f(x(i)+2*h/3,y(j)+h/3)+f(x(i)+h/3,y(j)+2*h/3)+f(x(i)-h/3,y(j)+h/3);
        
        b_vigur(pkt)=0;%-h^2*f_sum/6;
    end
end

%nedri jadar
for i=1:N+1
    j=1;
    pkt=i+(j-1)*(N+1);
    A(pkt,pkt)=1;
    b_vigur(pkt)=psi_1(x(i),a,beta1);
end

%efri jadar
for i=1:N+1
    j=M+1;
    pkt=i+(j-1)*(N+1);
    A(pkt,pkt)=1;
    b_vigur(pkt)=psi_2(x(i),b,beta2);
end
%vinstri jadar
for j=2:M
    i=1;
    pkt=i+(j-1)*(N+1);
    A(pkt,pkt)=(u^2+k^2)/(u*k);%pkt
    A(pkt,i+j*(N+1))=-k/(2*u);%upp
    A(pkt,i+(j-2)*(N+1))=-k/(2*u);%nidur
    A(pkt,pkt+1)=-h/k;%til haegri
    A(pkt,i+1+j*(N+1))=0;%ska
    b_vigur(pkt)=0;
end

%haegri jadar
for j=2:M
    i=N+1;
    pkt=i+(j-1)*(N+1);
    A(pkt,pkt)=(u^2+k^2)/(u*k);%pkt
    A(pkt,i+j*(N+1))=-k/(2*u);%upp
    A(pkt,i+(j-2)*(N+1))=-k/(2*u);%nidur
    A(pkt,pkt-1)=-h/k;%til vinstri
    A(pkt,i-1+(j-2)*(N+1))=0;%ska
    b_vigur(pkt)=0;
end    

A = sparse(A);
c=A\b_vigur;
size(c)
mn=(M+1)*(N+1);
HZ=(reshape(c(1:mn),N+1,M+1))'
