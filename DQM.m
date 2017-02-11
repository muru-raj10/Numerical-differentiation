function Coef=DQM(X,m) %input a column vectors X and required derivative
N= length(X);
C=zeros(N,N,N-1); %Lagrange polynomial of deg N-1 can be differentiated at most N-1 times. 

Val=zeros(N,N);
for j = 1 : N
    Val(:,j)=X(j)-X;
    Val(j,j)=1;
end
l1 = prod(Val);  %l(1)(xj) for the denominator

for j = 1 : N
    for i = 1 : N 
        if j ~= i
            C(i,j,1)=l1(i)/((X(i)-X(j))*l1(j));
        end
    end
end

Temp=-sum(C,2);% a column vector 
for i = 1 : N
    C(i,i,1)=Temp(i);
end

for k = 2 : N-1
    for j = 1 : N
        for i = 1 : N
            if i ~=j
                C(i,j,k)= k*(C(i,j,1)*C(i,i,k-1)-(C(i,j,k-1)/(X(i)-X(j))));
            end
        end
    end
    Temp2=-sum(C,2);
    for i = 1 : N
        C(i,i,k)=Temp2(i,k);
    end
end

Coef=zeros(N,N,m);
for k = 1 : m
Coef(:,:,k) = C(:,:,k);
end