S = diag([50 1])
x = mvnrnd([0 0],S,1000);
S = diag([5 1])
y = mvnrnd([0 0],S,1000);

theta = deg2rad(45);
R = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
x = x*R;
y = y*R;


[aInd, p, varToThisDimInTSpace, varToThisDimInSSpace, varToThisDimWorstCase, randCI, sigs] = ...
  varAlignment(x,y, 2)


S = diag([50 1 10])
x = mvnrnd([0 0 0],S,1000);
S = diag([5 10 1])
y = mvnrnd([0 0 0],S,1000);
R = rotation_angle_axis(deg2rad(45),[sqrt(2)/2, 0.0, sqrt(2)/2]);

x = x*R;
y = y*R;

[aInd, p, varToThisDimInTSpace, varToThisDimInSSpace, varToThisDimWorstCase, randCI, sigs] = ...
  varAlignment(x,y, 3)

p = 10;
S = diag([rand(1,10)*.2]);
x = mvnrnd(zeros(size(S,1),1),S,1000);
S = diag([rand(1,10)*.2]);
y = mvnrnd(zeros(size(S,1),1),S,1000);

[aInd, p, varToThisDimInTSpace, varToThisDimInSSpace, varToThisDimWorstCase, randCI, sigs] = ...
  varAlignment(x,y, 6)
