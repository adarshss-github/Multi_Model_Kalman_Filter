[M,K] = shrblgMat(ones(20,1)*1e3,ones(20,1)*15e5) ;
[f,w,Phi] = undmodpar(M,K) ;
save('ShrBldgMat') ;
C = rayleighdampmat(M,K,[1,2],[f(1),f(2)],[1/100,1/100]) ; %Damping

load('ShrBldgMat')
L = C ;
A = [zeros(20,20) eye(20); -inv(M)*K -inv(M)*L] ;
B = [zeros(20,1); -inv(M)*M*ones(20,1)] ;
Ca = zeros(2,20) ;
Ca(1,1) = 1 ;
Ca(1,2) = 1 ;
C = Ca*[-inv(M)*K -inv(M)*L] ;
D = Ca*-inv(M)*M*ones(20,1) ;
sys1 = ss(A,B,C,D) ;
 bode(sys1) ;
 
 hold on 
 
 A = [zeros(20,20) eye(20); -inv(M)*K -inv(M)*L] ;
B = [zeros(20,1); -inv(M)*M*ones(20,1)] ;
Ca = zeros(1,20) ;
Ca(1,2) = 1 ;
C = Ca*[-inv(M)*K -inv(M)*L] ;
D = Ca*-inv(M)*M*ones(20,1) ;
sys2 = ss(A,B,C,D) ;
 bode(sys1,sys2) ;
 
  A = [zeros(20,20) eye(20); -inv(M)*K -inv(M)*L] ;
B = [zeros(20,1); -inv(M)*M*ones(20,1)] ;
Ca = zeros(1,20) ;
Ca(1,3) = 1 ;
C = Ca*[-inv(M)*K -inv(M)*L] ;
D = Ca*-inv(M)*M*ones(20,1) ;
sys3 = ss(A,B,C,D) ;

 bode(sys1,sys2,sys3) ;
 
 [Ug,tg] = gausswn(0.3962,2,100,0,200,5) ;
save('ShrBldgMat.mat')
rms(Ug(1,:))
Ug = 9.81*Ug ;
save('ShrBldgMat.mat')

infmat = ones(20,1) ;
 F = zeros(20,20001) ;
Xddum = zeros(20,20001) ;
Xvdum = zeros(20,20001) ;
Xadum = zeros(20,20001) ;
Xd = cell(5,1) ;
Xv = cell(5,1) ;
Xa = cell(5,1) ;

for i = 1:1:5
    
    F = -M*infmat*Ug(i,:) ;
  [Xddum,Xvdum,Xadum,t] = NBMdof(1/2,1/6,M,C,K,F,100,zeros(20,1),zeros(20,1),200) ;  
  Xd{i,1} = Xddum ;
   Xv{i,1} = Xvdum ; 
   Xa{i,1} = Xadum ; 
   
end

save('ShrBldgMat.mat')

qddum = zeros(20,20001) ; %Note it seems the qa, qv and qd stored have some issues.
qvdum = zeros(20,20001) ;
qadum = zeros(20,20001) ;
qd = cell(5,1) ;
qv = cell(5,1) ;
qa = cell(5,1) ;

 for i = 1:1:5
    
     Xddum = Xd{i,1} ;
     Xvdum = Xv{i,1} ;
     Xadum = Xa{i,1} ;
   qddum = inv(Phi)*Xddum ;
   qvdum = inv(Phi)*Xvdum ;
   qadum = inv(Phi)*Xadum ;
   qd{i,1} = qddum ;
   qv{i,1} = qvdum ;
   qa{i,1} = qadum ;
   
 end

 save('ShrBldgMat.mat')
 
 %%2 Modes and 2 Sensors
 qd = qd{1,1} ;
 qd = qd([1,2],:) ;
 qv = qv{1,1} ;
 qv = qv([1,2],:) ;
 qdv = [qd;qv] ;
 
 