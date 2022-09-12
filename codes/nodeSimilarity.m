clc
clear

N = 15;

x1 = load("new_data\303_data.mat").Y;
%x = x(2,:);
Layout = [14,40;28,-36;46,-56;60,-28;34,14;0,2;32,58;48,-44;42,14;54,20;18,68;-2,64;-48,4;60,4;-50,38];

% ------------------------------ Time = 2 ----------------------------

x = x1(2,:);
x = dataNorm(x);

sigma = 1;
gamma = 0.941;

W = zeros(N,N);

for i = 1:N
    for j = 1:N
        if(i==j)
            continue;
        end
        W(i,j) = exp((-(x(i)-x(j))^2)/2*(sigma^2));
    end
end

A = W;
W(A>=gamma) = 1;
W(A<gamma) = 0;

G = grasp_struct; % Defining the struct variable
G.A = W; % Assigning Adjacency
G.layout = Layout; % Assigning coordinates for each node

figure()
grasp_show_graph(gca,G,'node_values',x);
xlim([-52,62]);
ylim([-58,70]);
sgtitle('Generated Graph');
colorbar

% ------------------------------ Time = 100 ----------------------------

x = x1(100,:);
x = dataNorm(x);

sigma = 1;
gamma = 0.998;

W = zeros(N,N);

for i = 1:N
    for j = 1:N
        if(i==j)
            continue;
        end
        W(i,j) = exp((-(x(i)-x(j))^2)/2*(sigma^2));
    end
end

A = W;
W(A>=gamma) = 1;
W(A<gamma) = 0;

G = grasp_struct; % Defining the struct variable
G.A = W; % Assigning Adjacency
G.layout = Layout; % Assigning coordinates for each node

figure()
grasp_show_graph(gca,G,'node_values',x);
xlim([-52,62]);
ylim([-58,70]);
sgtitle('Generated Graph');
colorbar