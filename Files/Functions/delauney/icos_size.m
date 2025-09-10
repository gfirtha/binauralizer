function [ point_num, edge_num, face_num, face_order_max ] = icos_size ( )

%*****************************************************************************80
%
%% icos_size() gives "sizes" for an icosahedron.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 July 2007
%
%  Author:
%
%    John Burkardt
%
%  Output:
%
%    integer POINT_NUM, the number of points.
%
%    integer EDGE_NUM, the number of edges.
%
%    integer FACE_NUM, the number of faces.
%
%    integer FACE_ORDER_MAX, the maximum order of any face.
%
  point_num = 12;
  edge_num = 30;
  face_num = 20;
  face_order_max = 3;

  return
end
