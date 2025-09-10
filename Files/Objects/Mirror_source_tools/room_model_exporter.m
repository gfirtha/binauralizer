clear
% Room setup:
vertices = [-2, 2
    1, 2
    3, -2
    -1, -2,
    -2, -1]'*2;
walls = [ 1, 2
    2, 3
    3, 4
    4, 5
    5, 1]';
% vertices = [-2, 2
%     2, 2
%     2, -2
%     -2, -2]'*2.5;
% walls = [ 1, 2
%     2, 3
%     3, 4
%     4, 1]';
% 
% R = 3;
% N = 15;
% fi = linspace(0,30-360/N,N)';
% vertices = [R*cosd(fi),R*sind(fi)]';
% walls = [(1:N)' circshift((1:N)',-1)]';
save room_model.mat
%
