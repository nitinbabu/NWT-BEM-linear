clc
clear 
close all

% previos work in separate files
% 27.03.2021 Kernel Integration Figured (using quadl)
% 11.04.2021 BCs and Matrix Construction Complete. Pending Validation
%  fixed Diagonal Aii = -SUM(Aij) [Becker eq 3.37]
% [PENDING] Singularity condition Placeholders fix 
% [PENDING] Validation
% [PENDING] Test lin-BEM for AddedMass problem

%Global Coordinates Input mesh

%     ^
%  y  |
%     |
%      - - - ->
%         x
% x = [0:0.5:10, 10 10 10 10 10 10 10,9.5:-0.5:0, 0 0 0 0 0 0 0 ]';
% y = [zeros(1,20),0:0.4286:3.0002,(3+zeros(1,19)),3.0002:-0.4286:0]';
% Nj = 55;
% Ej = Nj;
%******************
% Geometry
%******************
x = [0 1 2 3 3 3 3 2 1 0 0 0 0];
y = [0 0 0 0 1 2 3 3 3 3 2 1 0];
%BC =[ N N N D D D N N N D D D]; D-1 N-0 Specifies type of boundary like a
%flag (Note in Lin-BEM BC are on the nodes themselves)
bcn = [0 0 0 1 1 1 0 0 0 1 1 1];
bc = [ 0 0 0 1 2 3 0 0 0 0 0 0];
p = zeros(size(bc))';
q = zeros(size(bc))';

for c = 1:length(bc)
   
    if bcn(c) == 1        
        p(c) = bc(c);
        q(c) = NaN; 
    else
        q(c) = bc(c);
        p(c) = NaN;
    end
end


Ej = 13;
Nj = 13;
plot(x,y,'o-')
hold on
% placeholder values of phi and phin
% p  =  zeros(Nj,1);
% 
% q  = zeros(Nj,1);
% q = [0 0 0 1/2 1/sqrt(2) sqrt(3)/2   0 0 0 0 0 0];
% in natural frame over element Ej


%Influence coefficients
%  !! define nfx and nfy

for ms = 1:Nj-1
     
    A = [x(ms),y(ms)];        %intl point
    B = [x(ms+1),y(ms+1)];    %fnl point
    V = B - A;
    midV = A + 0.5 * V;
    normal2 = [V(2),  -V(1)]; % which normal change minus
    % plot([A(1), B(1)], [A(2), B(2)], 'r');
%     hold on
%     plot([midV(1), midV(1) + normal2(1)], [midV(2), midV(2) + normal2(2)], '.-b');

    xn_tail(ms) =  midV(1); xn_head(ms) =   midV(1) + normal2(1);
    yn_tail(ms) = midV(2);  yn_head(ms) = midV(2) + normal2(2);

        distn(ms) = sqrt((xn_head(ms)-xn_tail(ms)).^2+(yn_head(ms)-yn_tail(ms)).^2);
        nfx(ms) = (xn_head(ms)-xn_tail(ms))./distn(ms);
        nfy(ms) = (yn_head(ms)-yn_tail(ms))./distn(ms);

end

for s = 1:Nj-1
    for f = 1:Nj-1
% Here I am forming G  and Gn = div(G)*n 2D green's velocity potential 
if s==f % Doublecheck, currrently placeholders 11.04.21
    G(s,f) = 0.5;
    Gn(s,f) = 0;
else
    G(s,f) = -(2*pi)^-1*log(sqrt((x(s)-x(f)).^2+(y(s)-y(f)).^2));
    Gn(s,f) = -(2*pi)^-1*((x(s)-x(f))/((x(s)-x(f)).^2+(y(s)-y(f)).^2)).*nfx(f)-(2*pi)^-1*((y(s)-y(f))/((x(s)-x(f)).^2+(y(s)-y(f)).^2)).*nfy(f);
end
end
end


M = Ej;
% consider an element
%  (1)o-------------------o(2)

for P = 1:Nj-1 %Source Point
    for Q = 1:Nj-1 %Field Point
        K1 = Gn(P,Q);
        K2 = G(P,Q);
        
        %This bit chances of errors creeping
        xj = x(P);
        xj1 = x(Q);
        yj = y(P);
        yj1 = y(Q);
        
                N1 = @(eta) 0.5*(1-eta);%N1;
                N2 = @(eta) 0.5*(1+eta);%N2;
                
%                 k
                % % Jacobian of transformation
                dN1de = -0.5; dN2de = 0.5;
                dxde = dN1de.*xj + dN2de.*xj1;
                dyde = dN1de.*yj + dN2de.*yj1;

                J = @(eta) sqrt((dxde).^2 + (dyde).^2);
%                 K1 = 0.65; K2 = 0.95;
                id1 = @(eta) K1.*N1(eta).*J(eta); % [Disregard]q1 N1 refer to k-sum counter, not K1 Kernel
                id2 = @(eta) K1.*N2(eta).*J(eta);

                is1 = @(eta) K2.*N1(eta).*J(eta);
                is2 = @(eta) K2.*N2(eta).*J(eta);
                
                IS = quadl(is1,-1,1,1e-5) + quadl(is2,-1,1,1e-5);
                ID = quadl(id1,-1,1,1e-5) + quadl(id2,-1,1,1e-5);
            
        Am = ID; % X-plys with potential
        Bm = IS; % X-plys with gradient
      
      
       MatA(P,Q) = Am;
       MatB(P,Q) = Bm;
    end
end

%Better MatA
diagA = 0;
for P = 1:Nj-1 %Source Point
    for Q = 1:Nj-1 %Field Point
        if P ~= Q
        diagA = diagA-MatA(P,Q);
        else
        end
    end
    MatA(P,P) = diagA;
    diagA = 0;
end
% The previous steps kindof create the matrices A and B these need to be
% rearranged in some clever manner to get the product of known BCs and the
% corresponding As and Bs to the RHS and the Unkowns and the corresponding
% As and Bs i.e. influence coefficients whose product will get a LHS

% The  Caveat is to rearrange it with approprite sign changes.

% Constructing Ax = B then x =Ainverse B [.... A. BECKER]

sumB = 0;
A = zeros(size(MatA));
for P = 1:Nj-1 %Source Point
    for Q = 1:Nj-1 %Field Point
        
        if bcn(Q) == 1
            sumB = sumB - MatA(P,Q)*p(Q);
        else
            sumB = sumB + MatB(P,Q)*q(Q);
        end
         
        %
        if bcn(Q) == 1
            A(P,Q) = MatA(P,Q);
        else
            A(P,Q) = -MatB(P,Q);
        end
    
    end
    B(P) = sumB;
end
B = B';

z = inv(A).*B;




