S = diag([1 15])
x = mvnrnd([0 0],S,1000);
S = diag([80 3])
y = mvnrnd([0 0],S,1000);

% Rx = rotation_angle_axis(deg2rad(45),[sqrt(2)/2, 0.0, sqrt(2)/2]);
% Ry = rotation_angle_axis(deg2rad(-70),[sqrt(2)/2, 0.0, sqrt(2)/2]);
theta = deg2rad(45);
Rx = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
x = x*Rx;
y = y*Rx;


[Q,D] = dim.cpca({cov(x) cov(y)},{1000 1000},'tol',1e-10);

Q*diag(D(:,1))*Q'
cov(x)
Q*cov(x)*Q'
Q*diag(D(:,2))*Q'
cov(y)
Q*cov(y)*Q'

figure;
subplot(211); hold on
plot(x(:,1),x(:,2),'o');
plot(y(:,1),y(:,2),'ro');

xx = x*Q;
yy = y*Q;
subplot(212); hold on
plot(xx(:,1),xx(:,2),'o');
plot(yy(:,1),yy(:,2),'ro');
