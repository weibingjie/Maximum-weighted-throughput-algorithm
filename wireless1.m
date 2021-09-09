clc;clear;
n=30;  %方便改变随机点的个数
T=500;  %迭代次数
L=rand(n);  %产生n*n的邻接矩阵
for i=1:n
    L(i,i)=0;  %将对角线上的数置为0
end
L=10*L;
L=floor(L);  %向下去整 
L=mod(L,2);
for i=1:n
    for j=1:i
        L(j,i)=L(i,j);
    end  
end 
k=100;  %数据对象个数
V=100;
a=1;
C=1000*ones(n,n);  %链路容量先设置为定值
for i=1:n
    C(i,i)=0;  %将对角线上的数置为0
end
B=30;   %每个节点的缓存大小，假设每个节点的缓存大小相等
C_max=1000;
lamda=30;
d_max=2*lamda+C_max; %根据2A_max+C_max
h=zeros(n,k); %缓存初始化
h(n,:)=1;  %假设源节点在最后一个节点
h_MWT=h;
h_DFC=h;
h_VIP=h;
h_LRU=h;
h_LCE=h;
Q_MWT_1=zeros(n,k);
W_MWT_1=V*a*ones(n,k);
D_MWT_1=zeros(n,k);
Z_MWT_1=zeros(n,k);
v_MWT=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
r_MWT=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
u_MWT=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
u_average=zeros(1,T);%记录平均值
Q_MWT_average=zeros(T,1);
Q_DFC_average=zeros(T,1);
q_MWT=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
d_MWT=zeros(n,k);
dn_MWT=zeros(n,k);
ds_MWT=zeros(n,k);
Q_DFC_1=zeros(n,k);
D_DFC_1=zeros(n,k);
Z_DFC_1=zeros(n,k);
v_DFC=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
r_DFC=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
u_DFC=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
q_DFC=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
Q_VIP_1=zeros(n,k);
Q_LRU_1=zeros(n,k);
Q_LCE_1=zeros(n,k);
t=1; %迭代变量
while t<=T
    A=poissrnd(lamda,n,k);  %请求外在到达率
    %分配请求速率
    r_VIP=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
    r_LRU=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵
    r_LCE=zeros(n,n,k); %三维矩阵，生成k个n*n的矩阵               
    for i=1:n
        for j=1:n
            k_VIP=[];k_LRU=[];k_LCE=[];
            if(L(i,j)==1)
                for m=1:k
                    v_MWT(i,j,m)=max((1-h_MWT(i,m))*Q_MWT_1(i,m)-(1-h_MWT(j,m))*Q_MWT_1(j,m),0);
                    v_DFC(i,j,m)=max((1-h_DFC(i,m))*Q_DFC_1(i,m)-(1-h_DFC(j,m))*Q_DFC_1(j,m),0);
                    k_VIP=[k_VIP,max((1-h_VIP(i,m))*Q_VIP_1(i,m)-(1-h_VIP(j,m))*Q_VIP_1(j,m),0)];
                    k_LRU=[k_LRU,max((1-h_LRU(i,m))*Q_LRU_1(i,m)-(1-h_LRU(j,m))*Q_LRU_1(j,m),0)];
                    k_LCE=[k_LCE,max((1-h_LCE(i,m))*Q_LCE_1(i,m)-(1-h_LCE(j,m))*Q_LCE_1(j,m),0)];
                end
            end
            [m_VIP,index_VIP]=max(k_VIP);
            [m_LRU,index_LRU]=max(k_LRU);
            [m_LCE,index_LCE]=max(k_LCE);
            r_VIP(i,j,index_VIP)=C_max;
            r_LRU(i,j,index_LRU)=C_max;
            r_LCE(i,j,index_LCE)=C_max;
        end
    end
    for i=1:n
        for j=1:n
            if(sum(v_MWT(i,j,:))<=C(i,j))
                for m=1:k
                    r_MWT(i,j,m)=v_MWT(i,j,m);
                end
            else
                e=0;%辅助变量
                for m=1:k
                    if(v_MWT(i,j,m)>0)
                        e=e+1;
                    end
                end
                x=(sum(v_MWT(i,j,:))-C(i,j))/e;  %简化最优值
                for m=1:k
                    r_MWT(i,j,m)=max(v_MWT(i,j,m)-x,0);
                end
            end
            if(sum(v_DFC(i,j,:))<=C(i,j))
                for m=1:k
                    r_DFC(i,j,m)=v_DFC(i,j,m);
                end
            else
                e=0;%辅助变量
                for m=1:k
                    if(v_DFC(i,j,m)>0)
                        e=e+1;
                    end
                end
                x=(sum(v_DFC(i,j,:))-C(i,j))/e;  %简化最优值
                for m=1:k
                    r_DFC(i,j,m)=max(v_DFC(i,j,m)-x,0);
                end
            end
        end
    end
    %请求丢弃速率，最终丢弃速率，数据和虚拟数据传输速率
    for j=1:n
        for m=1:k
            if(Q_MWT_1(j,m)>W_MWT_1(j,m))
                d_MWT(j,m)=d_max;
            else
                d_MWT(j,m)=0;
            end
            if(W_MWT_1(j,m)>V*a)
                ds_MWT(j,m)=d_max;
            else
                ds_MWT(j,m)=0;
            end
            if(sum(r_MWT(:,j,m))<=D_MWT_1(j,m))
                for i=1:n
                    u_MWT(j,i,m)=r_MWT(i,j,m);
                    q_MWT(j,i,m)=0;
                end
            else
                e=0;%辅助变量
                for i=1:n
                    if(r_MWT(i,j,m)>0)
                        e=e+1;
                    end
                end
                b=(sum(r_MWT(:,j,m))-D_MWT_1(j,m))/e;%辅助变量
                for i=1:n
                    u_MWT(j,i,m)=max(r_MWT(i,j,m)-b,0);
                    q_MWT(j,i,m)=r_MWT(i,j,m)-u_MWT(j,i,m);
                end
            end
            if(sum(r_DFC(:,j,m))<=D_DFC_1(j,m))
                for i=1:n
                    u_DFC(j,i,m)=r_DFC(i,j,m);
                    q_DFC(j,i,m)=0;
                end
            else
                e=0;%辅助变量
                for i=1:n
                    if(r_DFC(i,j,m)>0)
                        e=e+1;
                    end
                end
                b=(sum(r_DFC(:,j,m))-D_DFC_1(j,m))/e;%辅助变量
                for i=1:n
                    u_DFC(j,i,m)=max(r_DFC(i,j,m)-b,0);
                    q_DFC(j,i,m)=r_DFC(i,j,m)-u_DFC(j,i,m);
                end
            end
        end
    end
    u_average(t)=sum(u_MWT(:))/n/k;
    Q_MWT_average(t)=mean(Q_MWT_1(:));
    Q_DFC_average(t)=mean(Q_DFC_1(:));
    %更新队列
    for i=1:n
        for m=1:k
            Q_MWT_2(i,m)=(1-h_MWT(i,m))*max(max(Q_MWT_1(i,m)-sum(r_MWT(:,i,m)),0)-d_MWT(i,m),0)+...
                A(i,m)+min(Z_MWT_1(i,m),A(i,m))+sum(r_MWT(i,:,m));
            dn_MWT(i,m)=min(max(Q_MWT_1(i,m)-sum(r_MWT(i,:,m)),0),(1-h_MWT(i,m))*d_MWT(i,m));
            W_MWT_2(i,m)=max(W_MWT_1(i,m)-(1-h_MWT(i,m))*ds_MWT(i,m),0)+dn_MWT(i,m);
            D_MWT_2(i,m)=max(D_MWT_1(i,m)-sum(u_MWT(i,:,m))-max(A(i,m)-Z_MWT_1(i,m),0),0)+...
                h_MWT(i,m)*Q_MWT_1(i,m)+sum(u_MWT(:,i,m));
            Z_MWT_2(i,m)=max(Z_MWT_1(i,m)-sum(q_MWT(i,:,m))-A(i,m),0)+sum(q_MWT(:,i,m))+...
                (1-h_MWT(i,m))*max(sum(r_MWT(:,i,m))-Z_MWT_1(i,m)-D_MWT_1(i,m),0);
            Q_DFC_2(i,m)=(1-h_DFC(i,m))*max(Q_DFC_1(i,m)-sum(r_DFC(:,i,m)),0)+A(i,m)+...
                min(Z_DFC_1(i,m),A(i,m))+sum(r_DFC(i,:,m));
            D_DFC_2(i,m)=max(D_DFC_1(i,m)-sum(u_DFC(i,:,m))-max(A(i,m)-Z_DFC_1(i,m),0),0)+...
                h_DFC(i,m)*Q_DFC_1(i,m)+sum(u_DFC(:,i,m));
            Z_DFC_2(i,m)=max(Z_DFC_1(i,m)-sum(q_DFC(i,:,m))-A(i,m),0)+sum(q_DFC(:,i,m))+...
                (1-h_DFC(i,m))*max(sum(r_DFC(:,i,m))-Z_DFC_1(i,m)-D_DFC_1(i,m),0);
            Q_VIP_2(i,m)=max(max(Q_VIP_1(i,m)-sum(r_VIP(:,i,m)),0)+A(i,m)+sum(r_VIP(i,:,m))-h_VIP(i,m)*80000,0);
            Q_LRU_2(i,m)=(1-h_LRU(i,m))*max(Q_LRU_1(i,m)-sum(r_LRU(:,i,m)),0)+A(i,m)+sum(r_LRU(i,:,m));
            Q_LCE_2(i,m)=(1-h_LCE(i,m))*max(Q_LCE_1(i,m)-sum(r_LCE(:,i,m)),0)+A(i,m)+sum(r_LCE(i,:,m));
        end
    end
    %缓存算法
    h_MWT=h;
    h_DFC=h;
    h_VIP=h;
    h_LRU=h;
    h_LCE=h;
    for i=1:n-1
        h_LCE(i,:)=randperm(k)<=B;
        k_MWT=[];F_MWT=[];
        k_DFC=[];F_DFC=[];
        k_VIP=[];F_VIP=[];
        k_LRU=[];F_LRU=[];
        for m=1:k
            for j=1:n
                if(u_MWT(j,i,m)~=0&&ismember(m,k_MWT)==0)
                    k_MWT=[k_MWT;m];
                end
                if(u_DFC(j,i,m)~=0&&ismember(m,k_DFC)==0)
                    k_DFC=[k_DFC;m];
                end
                if(ismember(m,k_VIP)==0)
                    k_VIP=[k_VIP;m];
                end
                if(r_LRU(i,j,m)~=0&&ismember(m,k_LRU)==0)
                    k_LRU=[k_LRU;m];
                end
            end
        end
        l_MWT=length(k_MWT);
        l_DFC=length(k_DFC);
        l_VIP=length(k_VIP);
        l_LRU=length(k_LRU);
        if(l_LRU<=B)
            for j=1:l_LRU
                h_LRU(i,k_LRU(j))=1;
            end
        else
           F_LRU=k_LRU(randperm(l_LRU,B));
           for j=1:length(F_LRU)
               h_LRU(i,j)=1;
           end
        end  
        for j=1:l_MWT
            F_MWT=[F_MWT;Q_MWT_1(i,k_MWT(j))];
        end
        for j=1:l_DFC
            F_DFC=[F_DFC;max(Q_DFC_1(i,k_DFC(j))-sum(r_DFC(i,:,m)))];
        end
        for j=1:l_VIP
            F_VIP=[F_VIP;Q_VIP_1(i,k_VIP(j))];
        end
        if(l_MWT<=B)
            for j=1:l_MWT
                h_MWT(i,k_MWT(j))=1;
            end
        else
           [F_sort,I]=sort(F_MWT);
           for j=l_MWT-B+1:l_MWT
               h_MWT(i,k_MWT(I(j)))=1;
           end
        end 
        if(l_DFC<=B)
            for j=1:l_DFC
                h_DFC(i,k_DFC(j))=1;
            end
        else
           [F_sort,I]=sort(F_DFC);
           for j=l_DFC-B+1:l_DFC
               h_DFC(i,k_DFC(I(j)))=1;
           end
        end  
        if(l_VIP<=B)
            for j=1:l_VIP
                h_VIP(i,k_VIP(j))=1;
            end
        else
           [F_sort,I]=sort(F_VIP);
           for j=l_VIP-B+1:l_VIP
               h_VIP(i,k_VIP(I(j)))=1;
           end
        end 
    end
    Q_MWT_1=Q_MWT_2;
    W_MWT_1=W_MWT_2;
    D_MWT_1=D_MWT_2;
    Z_MWT_1=Z_MWT_2;
    Q_DFC_1=Q_DFC_2;
    D_DFC_1=D_DFC_2;
    Z_DFC_1=Z_DFC_2;
    Q_VIP_1=Q_VIP_2;
    Q_LRU_1=Q_LRU_2;
    Q_LCE_1=Q_LCE_2;
    t=t+1;
end
%求请求队列平均值
Q_MWT=mean(Q_MWT_1(:))
Q_DFC=mean(Q_DFC_1(:))
Q_VIP=mean(Q_VIP_1(:))
Q_LRU=mean(Q_LRU_1(:))
Q_LCE=mean(Q_LCE_1(:))
r_k=mean(u_average)
x=0:1:T-1;
plot(x,Q_MWT_average,'-b',x,Q_DFC_average,'-r');
set(gca,'XTick',[0:100:500]) %x轴范围0-500，间隔100
set(gca,'YTick',[0:100:500]) %y轴范围0-500，间隔100
legend('with request dropping queue','without request dropping queue');   %右上角标注
xlabel('Time')  %x轴坐标描述
ylabel('Average Request Queue Backlog') %y轴坐标描述