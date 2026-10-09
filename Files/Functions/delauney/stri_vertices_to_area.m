function area = stri_vertices_to_area ( r, v1, v2, v3 )

%*****************************************************************************80
%
%% stri_vertices_to_area() computes the area of a spherical triangle.
%
%  Discussion:
%
%    A sphere centered at 0 in 3D satisfies the equation:
%
%      X * X + Y * Y + Z * Z = R * R
%
%    A spherical triangle is specified by three points on the surface
%    of the sphere.
%
%    The area formula is known as Girard's formula.
%
%    The area of a spherical triangle is:
%
%      AREA = ( A + B + C - PI ) * R * R
%
%    where A, B and C are the (surface) angles of the triangle.
%
%    Thanks to Mark Bydder for pointing out a problem that can occur if
%    two of the vertices are essentially equal.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 October 2015
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real R, the radius of the sphere.
%
%    real V1(3), V2(3), V3(3), the vertices of the triangle.
%
%  Output:
%
%    real AREA, the area of the spherical triangle.
%

%
%  Compute the lengths of the sides of the spherical triangle.
%
  [ as, bs, cs ] = stri_vertices_to_sides ( r, v1, v2, v3 );
%
%  Set AREA to zero if one side is essentially zero.
%
  if ( min ( as, min ( bs, cs ) ) < sqrt ( eps ) * ( max ( as, max ( bs, cs ) ) + 1.0 ) )
    area = 0.0;
    return
  end
%
%  Get the spherical angles.
%
  [ a, b, c ] = stri_sides_to_angles ( r, as, bs, cs );
%
%  Get the area.
%
  area = stri_angles_to_area ( r, a, b, c );

  return
end
