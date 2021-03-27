clc
clear 
close all

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

x = [0 1 2 3 3 3 3 2 1 0 0 0 0];
y = [0 0 0 0 1 2 3 3 3 3 2 1 0];
Ej = 13;
Nj = 13;
plot(x,y,'o-')
hold on
% placeholder values of phi and phin
p  = zeros(Nj,1);
q  = zeros(Nj,1);

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
if s==f
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

q1 = 10;
p1 = 2;
q2 = 10;
p2 = 2;
for P = 1:Nj-1 %Source Point
    for Q = 1:Nj-1 %Field Point
        K1 = G(P,Q);
        K2 = Gn(P,Q);
        xj = x(P);
        xj1 = x(Q);
        yj = y(P);
        yj1 = y(Q);
        
        for m = 1:M
            
               
                    N1 = @(eta) 0.5*(1-eta);%N1;
               
                    N2 = @(eta) 0.5*(1+eta);%N2;
                
%                 k
                % % Jacobian of transformation
                dN1de = -0.5; dN2de = 0.5;
                dxde = dN1de.*xj + dN2de.*xj1;
                dyde = dN1de.*yj + dN2de.*yj1;

                J = @(eta) sqrt((dxde).^2 + (dyde).^2);
%                 K1 = 0.65; K2 = 0.95;
                is1 = @(eta) K1.*N1(eta).*J(eta); % q1 N1 refer to k-sum counter, not K1 Kernel
                is2 = @(eta) K1.*N2(eta).*J(eta);

                id1 = @(eta) K2.*N1(eta).*J(eta);
                id2 = @(eta) K2.*N2(eta).*J(eta);
                
                IS = q1.*quadl(is1,-1,1) + q2.*quadl(is2,-1,1);
                ID = p1.*quadl(id1,-1,1) + p2.*quadl(id2,-1,1);
            
            Am(m) = IS;
            Bm(m) = ID;
        end
        Amk = sum(Am);
        Bmk = sum(Bm);
       MatA(P,Q) = Amk;
       MatB(P,Q) = Bmk;
    end
end


