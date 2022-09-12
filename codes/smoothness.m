N = 15; % Number of Nodes
M = 257; % Number of realizations
alpha = 4; % Optimization Parameter
beta = 1; % Sparsity Parameter

% ------------------------ part (a) -------------------------------------

data = Y';
X = data;
Layout = [14,40;28,-36;46,-56;60,-28;34,14;0,2;32,58;48,-44;42,14;54,20;18,68;-2,64;-48,4;60,4;-50,38];

% -------------------------- part (b) -----------------------------------

yalmip('clear')

% Problem is - argmin ||X-Z||_F^2 + alpha*Tr(Z'Lz)
Z_hat = sdpvar(N,M);
L_hat = sdpvar(N,N);

constraints = [(L_hat*ones(N,1))==zeros(N,1),sum(diag(L_hat))==N];
for i=1:N
    for j=1:N
        if(i==j)
            %constraints = [constraints,L_hat(i,j)>=0];
            continue;
        else
            constraints = [constraints,L_hat(i,j)<=0];
        end
    end
end

objective = norm(X-Z_hat,'fro')^2 + alpha*sum(diag(Z_hat'*L_hat*Z_hat)) + beta*norm(abs(L_hat),1);
options = sdpsettings('solver','mosek','mosek.MSK_IPAR_BI_MAX_ITERATIONS',50,'verbose',1);

optimize(constraints,objective,options);
%optimize(constraints,objective);
Z_est = value(Z_hat);
L_est = value(L_hat);

[U_hat,~] = eig(L_est);
y_hat = inv(U_hat)*Z_est; % Estimated signal

D_hat = diag(L_est);
A_hat = L_est - D_hat;

G_hat = grasp_struct; % Defining the struct variable
G_hat.A = A_hat; % Assigning Adjacency
% Assigning coordinates to nodes
G_hat.layout = Layout; % Assigning coordinates for each node

figure()
grasp_show_graph(gca,G_hat);
title('Learnt Graph')
xlim([-52,62]);
ylim([-58,70]);

RMSErr = norm(y_hat-yData,'fro')^2;
message = ['The obtained error is ', num2str(RMSErr)];
disp(message);

RMS = [];
alphaIter = 1:10;

for alpha = alphaIter
    objective = norm(X-Z_hat,'fro')^2 + alpha*sum(diag(X'*L_hat*X - Z_hat'*L_hat*epsilon - epsilon'*L_hat*Z_hat - epsilon'*L_hat*epsilon));
    optimize(constraints,objective);
    Z_est = value(Z_hat);
    L_est = value(L_hat);

    [U_hat,Lambda_hat] = eig(L_est);
    y_hat = inv(U_hat)*Z_est; % Estimated signal
    RMSErr = norm(y_hat-yData,'fro')^2;
    RMS = [RMS,RMSErr];
end

figure()
plot(alphaIter,RMS,'black',LineWidth=1.5);
xlabel('Optimization Parameter (\alpha)');
ylabel('RMS Error');
sgtitle('Effect of Changing Optimization Parameter (\alpha)')