% Spectral Graph Learning
% Assume it is GWSS

N = 15; % Number of nodes
threshold = 0.01; % Thresholding in A

Layout = [14,40;28,-36;46,-56;60,-28;34,14;0,2;32,58;48,-44;42,14;54,20;18,68;-2,64;-48,4;60,4;-50,38];
X = Y';

Cx = (1/15)*(X*(X'));
[eVec,eVal] = eig(Cx);

% We have to minimize -> min ||Z||_{0} subject to Z = U \Lambda U^{T}

Lambda = sdpvar(N,1);
L = diag(Lambda);
Z = eVec*L*(eVec');

constraints = [(Z*ones(N,1))==zeros(N,1),sum(diag(Z))==N,Z==eVec*L*(eVec')];

for i=1:N
    for j=1:N
        if(i==j)
            constraints = [constraints,Z(i,j)>=0];
            continue;
        else
            constraints = [constraints,Z(i,j)<=0];
        end
    end
end

objective = norm(Z,1);
options = sdpsettings('solver','mosek');

optimize(constraints,objective,options);
Z_est = abs(value(Z));

D_hat = diag(Z_est);
A_hat = Z_est - D_hat;
W = A_hat;
A_hat(W>threshold) = 1;
A_hat(W<=threshold) = 0;

G_hat = grasp_struct; % Defining the struct variable
G_hat.A = A_hat; % Assigning Adjacency
% Assigning coordinates to nodes
G_hat.layout = Layout; % Assigning coordinates for each node

figure()
grasp_show_graph(gca,G_hat);
title('Learnt Graph')
xlim([-52,62]);
ylim([-58,70]);