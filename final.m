function final(scene, scale)
if nargin < 1
    scene = 0;
end

% Set up the scene here.
% Note that links don't have to be placed exactly. The first call to
% solveLinkage() will automatically snap the links so that all the
% constraints are satisfied.    
linkages = [];
pins = [];
particles = [];
switch scene

    case 10
        fileID = fopen('dimensions.txt','wt');
        fprintf(fileID,'%6s %12s Link Name\n','Original Measurement','Scaled Measurement');

        % Extra credit!
        % 
        % Backpack link (ground)
        links(1).angle = pi/2; % rotation from the positive x-axis
        links(1).pos = [-1;0]; % position of the center of rotation
        links(1).verts = [0,-0.1;9*scale,-0.1;9*scale,0.1;0,0.1]'; % display vertices
        fprintf(fileID,'%15f %17f Backpack Link\n',[9, 9*scale]);
      

        % driver link link
        links(2).angle = 0;
        links(2).pos = [-1;0];
        links(2).verts = [0,-0.1;6.5*scale,-0.1;6.5*scale,0.1;0,0.1]';
        fprintf(fileID,'%15f %17f Driver Link\n',[6.5, 6.5*scale]);

        % top scissor
        links(3).angle = 0;
        links(3).pos = [-1;9*scale];
        links(3).verts = [0,-0.1;12.5*scale,-0.1;12.5*scale,0.1;0,0.1]';
        fprintf(fileID,'%15f %17f Top Scissor Joint\n',[12.5, 12.5*scale]);


        % Bottom scissor
        links(4).angle = pi/2;
        links(4).pos = [6.5*scale;0];
        links(4).verts = [0,-0.1;15*scale,-0.1;15*scale,0.1;0,0.1]';
        fprintf(fileID,'%15f %17f Bottom Scissor Joint\n',[15, 15*scale]);


        %long wing top
        links(5).angle = 0; 
        links(5).pos = [6.5*scale;15*scale];
        links(5).verts = [0,-0.1;18.5*scale,-0.1;18.5*scale,0.1;0,0.1]';
        fprintf(fileID,'%15f %17f Long Wing\n',[18.5, 18.5*scale]);


        %wing connect 
        links(6).angle = pi/2; 
        links(6).pos = [12.5*scale; 9*scale]; 
        links(6).verts = [0,-0.1;9.5*scale,-0.1;9.5*scale,0.1;0,0.1]';
        fprintf(fileID,'%15f %17f Wing Support\n',[9.5, 9.5*scale]);

        
        % Which link is grounded?
        grounded = 1;
        % Which link is the driver?
        % Note: the driver must be attached (with a pin) to the ground.
        driver = 2;
        
        % Bottom-left
        pins(1).links = [1,2];
        pins(1).pts = [0,0;0,0]';
        % top-left
        pins(2).links = [1,3];
        pins(2).pts = [9*scale,0;0,0]';
        % Left-top
        pins(3).links = [2,4];
        pins(3).pts = [6.5*scale,0;0,0]';
        % Right-top
        pins(4).links = [3,4];
        pins(4).pts = [8.5*scale,0;5.5*scale,0]'; % pin location on link3 is randomized

        %wing boi 
        pins(5).links = [4,5];
        pins(5).pts = [15*scale,0; 0,0]';

        %wing support 
        pins(6).links = [3, 6]; 
        pins(6).pts = [12.5*scale,0;0,0]';

        %wing support top 
        pins(7).links = [5, 6];
        pins(7).pts = [4*scale,0; 9.5*scale,0]';
       
        
        % List of tracer particles for display
        particles(1).link = 4; % which link?
        particles(1).pt = [0.5;0.1]; % tracer particle point in local coords
        particles(1).ptsWorld = zeros(2,0); % transformed points, initially empty
        particles(2).link = 4;
        particles(2).pt = [2.5;-0.1];
        particles(2).ptsWorld = zeros(2,0);

        fclose(fileID);
end

% Initialize
for i = 1 : length(links)
    links(i).grounded = (i == grounded);
    links(i).driver = (i == driver);
    % These target quantities are only used for grounded and driver links
    links(i).angleTarget = links(i).angle;
    links(i).posTarget = links(i).pos;
end
drawScene(links,pins,particles,true);

% lsqnonlin options
if verLessThan('matlab','8.1')
    opt = optimset(...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
else
    opt = optimoptions('lsqnonlin',...
        'Jacobian','on',...
        'DerivativeCheck','off',...
        'Display','off'); % final-detailed iter-detailed off
end

% Simulation loop
%t = 0; % current time
%T = 1; % final time
dt = 0.01; % time step
%angVel = 2*pi; % driver angular velocity
bounds = [links(driver).angleTarget, links(driver).angleTarget + pi/4 ];
%bounds = [-1*pi/2, pi/4];
frequency = 10; 
links(driver).angleTarget = bounds(1); 
add = true;
%while t < T
while true
    % Procedurally set the driver angle.
    % Right now, the target angle is being linearly increased, but you may
    % want to do something else.
    
    if links(driver).angleTarget > bounds(2)
        add = false;
        %links(driver).angleTarget = links(driver).angleTarget - dt*frequency; 
    end 

    if links(driver).angleTarget < bounds(1)
        add = true; 
        %links(driver).angleTarget = links(driver).angleTarget + dt*frequency; 
    end
    if add
        links(driver).angleTarget = links(driver).angleTarget + dt*frequency; 
    end

    if ~add
        links(driver).angleTarget = links(driver).angleTarget - dt*frequency; 
    end
   

    %links(driver).angleTarget = links(driver).angleTarget + dt*angVel;
    % Solve for linkage orientations and positions
    [links,feasible] = solveLinkage(links,pins,opt);
    % Update particle positions
    particles = updateParticles(links,particles);
    % Draw scene
    drawScene(links,pins,particles);
    % Quit if over-constrained
    if ~feasible
        break;
    end
    %t = t + dt;
end

end

%%
function [R,dR] = rotationMatrix(angle)
c = cos(angle);
s = sin(angle);
% Rotation matrix
R = zeros(2);
R(1,1) = c;
R(1,2) = -s;
R(2,1) = s;
R(2,2) = c;
if nargout >= 2
    % Rotation matrix derivative
    dR = zeros(2);
    dR(1,1) = -s;
    dR(1,2) = -c;
    dR(2,1) = c;
    dR(2,2) = -s;
end
end

%%
function [links,feasible] = solveLinkage(links,pins,opt)
nlinks = length(links);
% Extract the current angles and positions into a vector
angPos0 = zeros(3*nlinks,1);
for i = 1 : nlinks
    link = links(i);
    ii = (i-1)*3+(1:3);
    angPos0(ii(1)) = link.angle;
    angPos0(ii(2:3)) = link.pos;
end
% Limits
lb = -inf(size(angPos0));
ub =  inf(size(angPos0));
% Solve for angles and positions
[angPos,r2] = lsqnonlin(@(angPos)objFun(angPos,links,pins),angPos0,lb,ub,opt);
% If the mechanism is feasible, then the residual should be zero
feasible = true;
if r2 > 1e-6
    fprintf('Mechanism is over constrained!\n');
    feasible = false;
end
% Extract the angles and positions from the values in the vector
for i = 1 : length(links)
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end
end

%%
function [f,J] = objFun(angPos,links,pins)
nlinks = length(links);
npins = length(pins);
% Temporarily change angles and positions of the links. These changes will
% be undone when exiting this function.
for i = 1 : nlinks
    ii = (i-1)*3+(1:3);
    links(i).angle = angPos(ii(1));
    links(i).pos = angPos(ii(2:3));
end

% Evaluate constraints
ndof = 3*nlinks;
ncon = 3 + 3 + 2*npins; % 3 for ground, 3 for driver, 2*npins for pins
f = zeros(ncon,1);
J = zeros(ncon,ndof);
k = 0;
% Some angles and positions are fixed
for i = 1 : nlinks
    link = links(i);
    if link.grounded || link.driver
        % Grounded and driver links have their angles and positions
        % prescribed.
        f(k+1,    1) = link.angle - link.angleTarget;
        f(k+(2:3),1) = link.pos - link.posTarget;
        % The Jacobian of this constraint is the identity matrix
        J(k+(1:3),k+(1:3)) = eye(3);
        k = k + 3;
    end
end
% Pin constraints
for i = 1 : npins
    pin = pins(i);
    rows = k+(1:2); % row index of this pin constraint
    indLinkA = pin.links(1); % array index of link A
    indLinkB = pin.links(2); % array index of link B
    linkA = links(indLinkA);
    linkB = links(indLinkB);
    [Ra,dRa] = rotationMatrix(linkA.angle);
    [Rb,dRb] = rotationMatrix(linkB.angle);
    % Local positions
    ra = pin.pts(:,1);
    rb = pin.pts(:,2);
    % World positions
    xa = Ra * ra + linkA.pos;
    xb = Rb * rb + linkB.pos;
    p = xa(1:2) - xb(1:2);
    f(rows,1) = p;
    % Column indices for the angles and positions of links A and B
    colAngA = (indLinkA-1)*3 + 1;
    colPosA = (indLinkA-1)*3 + (2:3);
    colAngB = (indLinkB-1)*3 + 1;
    colPosB = (indLinkB-1)*3 + (2:3);
    % The Jacobian of this constraint is the partial derivative of f wrt
    % the angles and positions of the two links.
    J(rows,colAngA) = dRa * ra;
    J(rows,colPosA) = eye(2);
    J(rows,colAngB) = -dRb * rb;
    J(rows,colPosB) = -eye(2);
    k = k + 2;
end
end

%%
function particles = updateParticles(links,particles)
% Transform particle position from local to world
for i = 1 : length(particles)
    particle = particles(i);
    link = links(particle.link);
    R = rotationMatrix(link.angle);
    x = R * particle.pt + link.pos;
    % Append world position to the array (grows indefinitely)
    particles(i).ptsWorld(:,end+1) = x;
end
end

%%
function drawScene(links,pins,particles,initialize)
if nargin < 4
    initialize = false;
end
if initialize
    clf;
else
    cla;
end
hold on;
grid on;
% Draw links
for i = 1 : length(links)
    link = links(i);
    R = rotationMatrix(link.angle);
    % Draw frame
    p = link.pos; % frame origin
    s = 0.1; % frame display size
    px = p + s*R(:,1); % frame x-axis
    py = p + s*R(:,2); % frame y-axis
    plot([p(1),px(1)],[p(2),px(2)],'r','LineWidth',3);
    plot([p(1),py(1)],[p(2),py(2)],'g','LineWidth',3);
    % Draw link geometry
    if link.grounded
        color = [1 0 0];
    elseif link.driver
        color = [0 1 0];
    else
        color = [0 0 1];
    end
    E = [R,link.pos;0,0,1]; % transformation matrix
    vertsLocal = [link.verts;ones(1,size(link.verts,2))];
    vertsWorld = E * vertsLocal;
    plot(vertsWorld(1,[1:end,1]),vertsWorld(2,[1:end,1]),'Color',color);
end
% Draw pins
for i = 1 : length(pins)
    pin = pins(i);
    linkA = links(pin.links(1));
    linkB = links(pin.links(2));
    Ra = rotationMatrix(linkA.angle);
    Rb = rotationMatrix(linkB.angle);
    xa = Ra * pin.pts(:,1) + linkA.pos;
    xb = Rb * pin.pts(:,2) + linkB.pos;
    plot(xa(1),xa(2),'co','MarkerSize',10,'MarkerFaceColor','c');
    plot(xb(1),xb(2),'mx','MarkerSize',10,'LineWidth',2);
end
% Draw particles
for i = 1 : length(particles)
    particle = particles(i);
    if ~isempty(particle.ptsWorld)
        plot(particle.ptsWorld(1,:),particle.ptsWorld(2,:),'k');
        plot(particle.ptsWorld(1,end),particle.ptsWorld(2,end),'ko');
    end
end
axis equal;
drawnow;
end
