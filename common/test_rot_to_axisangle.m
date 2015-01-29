% Test the function rot_to_axisangle
%
% Yoel Shkolnisky, January 2015.

%% Rotation with no degeneracy
rotaxis=[1;1;1];
rotaxis=rotaxis./norm(rotaxis);
gamma=pi/3;
R=axisangle_to_rot(rotaxis,gamma);
[rotaxis2,gamma2]=rot_to_axisangle(R);
fprintf('Test1:\n');
fprintf('Error in rotation axis   %5.2e degrees\n', acosd(dot(rotaxis,rotaxis2)));
fprintf('Error in in-plane angle  %5.2e degrees\n', abs(gamma-gamma2));
fprintf('Error in frobenious norm %5.2e \n', norm(R-axisangle_to_rot(rotaxis2,gamma2),'fro'));
fprintf('*********************\n\n');

%% No rotation 1
R=eye(3);
[rotaxis2,gamma2]=rot_to_axisangle(eye(3));
fprintf('Test2:\n');
fprintf('Error in frobenious norm %5.2e \n', norm(R-axisangle_to_rot(rotaxis2,gamma2),'fro'));
fprintf('*********************\n\n');

%% No rotation 2
rotaxis=[1;1;1];
rotaxis=rotaxis./norm(rotaxis);
gamma=0;
R=axisangle_to_rot(rotaxis,gamma);
[rotaxis2,gamma2]=rot_to_axisangle(R);
fprintf('Test3:\n');
fprintf('Error in frobenious norm %5.2e \n', norm(R-axisangle_to_rot(rotaxis2,gamma2),'fro'));
fprintf('*********************\n\n');

%% In-plane rotation by pi
rotaxis=[1;1;1];
rotaxis=rotaxis./norm(rotaxis);
gamma=pi;
R=axisangle_to_rot(rotaxis,gamma);
[rotaxis2,gamma2]=rot_to_axisangle(R);
fprintf('Test4:\n');
fprintf('Error in frobenious norm %5.2e \n', norm(R-axisangle_to_rot(rotaxis2,gamma2),'fro'));
fprintf('*********************\n\n');

%% Random rotation
[R,~,~]=svd(randn(3));
if det(R)<0
    R(:,1)=-R(:,1);
end
[rotaxis2,gamma2]=rot_to_axisangle(R);
fprintf('Test5:\n');
fprintf('Error in frobenious norm %5.2e \n', norm(R-axisangle_to_rot(rotaxis2,gamma2),'fro'));
fprintf('*********************\n\n');
