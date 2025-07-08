
% II,II II,IS II,SI II,SS
% IS,II IS,IS IS,SI IS,SS
% SI,II SI,IS SI,SI SI,SS
% SS,II SS,IS SS,SI SS,SS

Mu=[1 4/5 4/5 0;
    0 0   0   0;
    0 0   0   0;
    0 1/5 1/5 1];

% I,I I,S
% S,I S,S

Eq=@(q)[1 q;
        0 1-q]^2;

Eqinv=@(qa)[1 -qa/(1-qa);
           0 1/(1-qa)]^2;

Reset=[2/12 2/6;1/12 1/6];
% Reset=[1/6 1/6;1/6 1/6];
preset=0;

filename='EntBob_final_model_v2';

%

u=sparse([1 0;0 0]);d=sparse([0 0;0 1]); 
l=sparse([0 0;1 0]);r=sparse([0 1;0 0]);
id=speye(2);

Nanc=0;
Mvec=[4 6 8 10 12 16 20]+Nanc;
Nsamp_M=[5000 4000 4000 4000 4000 3000 2000];
NM=length(Mvec);
Nt_M=4*Mvec;
qvec_M=cell(NM,1);
sigmavec_M=cell(NM,1);
qbarvec_M=cell(NM,1);
qmean=.1;
for nm=1:NM
         M=Mvec(nm);
if M<=12
  sigmavec_M{nm}=(0:.05:1.4)*qmean;
elseif M==16
  sigmavec_M{nm}=[(0:.05:.6),(.6:.05:1.1)]*qmean;
elseif M==20
  sigmavec_M{nm}=(.4:.05:1)*qmean;
end

p=.9;
q1=qmean-sigmavec_M{nm}.*sqrt((1-p)/p);
q2=(qmean-p*q1)/(1-p);
qvec_M{nm}=[q1;q2];
%mmax=size(qvec_M{nm},2);
qbarvec_M{nm}=1-abs(1-q1).^p.*abs(1-q2).^(1-p);

end

%%
OP_M=cell(NM,1);
MI_M=cell(NM,1);


for nm=1:NM
    qvec=qvec_M{nm};
    qbarvec=qbarvec_M{nm};
    mmax=size(qvec_M{nm},2);
    Nsamp=Nsamp_M(nm);
    Nt=Nt_M(nm);
    t_M{nm}=[1 2:Nt];
    t=t_M{nm};
    Ntn=length(t);
    M=Mvec(nm);
    OP=zeros(Nt,Nsamp,mmax,M);
    MI=zeros(Nt,Nsamp,mmax,M*(M-1)/2);
    tic
    for m=1:mmax
        Oi=cell(M,1);
        Oi(1:M)={id};
        uM=cell(M,1);
        dM=cell(M,1);
        lM=cell(M,1);
        rM=cell(M,1);
        weight=zeros(2^M,1);
        weightA=zeros(2^M,1);
        weightAB=zeros(2^M,1);
        weightB=zeros(2^M,1);
        for h=0:(2^M-1)
            hbin=dec2bin(h,M);
            weight(h+1)=sum(hbin=='1');
            weightA(h+1)=sum(hbin([1]) == '1');
            weightAB(h+1)=sum(hbin([1 2]) == '1');
            weightB(h+1)=sum(hbin([2]) == '1');
%             weightABconj(h+1)=sum(hbin(3:M) == '1');
        end
        h=0;
        psi0=zeros(2^M,1);
        psi0(1)=2^M/(2^M+1);
        psi0(end)=1/(2^M+1);        
%         psi0=(1/3-h).^(weight).*(2/3+h).^(M-weight);



        for i=1:M
            Oitmp=Oi;
            Oitmp(i)={u};
            uM{i}=ProductStateOp(Oitmp,M);
            Oitmp(i)={d};
            dM{i}=ProductStateOp(Oitmp,M);
            Oitmp=Oi;
            Oitmp(i)={l};
            lM{i}=ProductStateOp(Oitmp,M);
            Oitmp=Oi;
            Oitmp(i)={r};
            rM{i}=ProductStateOp(Oitmp,M);
        end

        C1M=cell(M,1);
        C2M=cell(M,1);

        for i=1:M
            C1i=1;
            if i ~= 1
                C1i=uM{1}*uM{i}+dM{1}*dM{i}+lM{1}*rM{i}+rM{1}*lM{i};
            end
            C1M{i}=C1i;
        end
        for j=1:M
            C2j=1;
            if j~=2
                C2j=uM{2}*uM{j}+dM{2}*dM{j}+lM{2}*rM{j}+rM{2}*lM{j};
            end
            C2M{j}=C2j;
        end
        

        toc

        q=qvec(:,m);
        qbar=qbarvec(m);
        tic
        [M m]
        parfor nr=1:Nsamp
       	   idM2=speye(2^(M-2));
           Psi=psi0;
           layer=0;
           nnt=0;

           OP_tmp=zeros(Nt,M);
           MI_tmp=zeros(Nt,M*(M-1)/2);
           %sitedisorder=binornd(1,1-p,M,1)+1;
           %qi=q(sitedisorder);
           Eaqijt=kron(Eqinv(qbar),Eqinv(qbar));

           for nt=1:Nt
               for k=1:((M-Nanc)/2)
                   ij=randperm(M-Nanc,2);
                   ij=sort(ij);
                   i=ij(1);j=ij(2);

                   qit=q(binornd(1,1-p)+1);
                   qjt=q(binornd(1,1-p)+1);                   

                   Eqijt=kron(Eq(qit),Eq(qjt));
                   
                   Muij=kron(sparse((Eaqijt*Eqijt*Mu)),idM2);                 
                   if j ~= 2
                       Psi=C2M{j}*C1M{i}*Muij*C2M{j}*C1M{i}*Psi;
                   else
                       Psi=Muij*Psi;
                   end
               end
               Psi=Psi/sum(Psi);
               weighttmpA=weightA;
               weighttmpB=weightB;
               weighttmpAB=weightAB;
               if mod(nt,M)==0 || nt==Nt
                   OP_tmp_tmp=-log2((2.^(2*weighttmpA-1))'*Psi); %'
                   MI_tmp_tmp=OP_tmp_tmp ... %
                       -log2((2.^(2*weighttmpB-1))'*Psi) ... %'
                       +log2((2.^(2*weighttmpAB-2))'*Psi); %'


                   nnnind=1;
                   for nnn=3:M
                       weighttmpB=C2M{nnn}*weightB;
                       weighttmpAB=C2M{nnn}*weightAB;
                       nnnind=nnnind+1;
                       MI_tmp_tmp(nnnind)=OP_tmp_tmp ...
                           -log2((2.^(2*weighttmpB-1))'*Psi) ... %'
                           +log2((2.^(2*weighttmpAB-2))'*Psi); %'
                       
                   end
                   for mmm=2:(M)
                       weighttmpA=C1M{mmm}*weightA;
                       OP_tmp_tmp(mmm)=-log2((2.^(2*weighttmpA-1))'*Psi); %'
                       if mmm<M
                           weighttmpB=C2M{mmm+1}*weightB;
                           weighttmpAB=C2M{mmm+1}*C1M{mmm}*weightAB;
                           nnnind=nnnind+1;
                           MI_tmp_tmp(nnnind)=OP_tmp_tmp(mmm) ... %
                               -log2((2.^(2*weighttmpB-1))'*Psi) ... %'
                               +log2((2.^(2*weighttmpAB-2))'*Psi); %'
                           
                       end
                       for nnn=(mmm+2):M
                           weighttmpB=C2M{nnn}*weightB;
                           weighttmpAB=C2M{nnn}*C1M{mmm}*weightAB;
                           nnnind=nnnind+1;
                           MI_tmp_tmp(nnnind)=OP_tmp_tmp(mmm) ... %
                               -log2((2.^(2*weighttmpB-1))'*Psi) ... %'
                               +log2((2.^(2*weighttmpAB-2))'*Psi); %'
                       end
                   end
                   OP_tmp(nt,:)=OP_tmp_tmp;
                   MI_tmp(nt,:)=MI_tmp_tmp;
               end

           end
 	      OP(:,nr,m,:)=OP_tmp;
          MI(:,nr,m,:)=MI_tmp;

        end
        OP_M{nm}=OP;
        MI_M{nm}=MI;
        clear psi0 weightB weightA weightAB weightABconj weight C1M C2M C12M uM dM lM rM
        save(filename,'-v7.3')
        toc
    end
end





