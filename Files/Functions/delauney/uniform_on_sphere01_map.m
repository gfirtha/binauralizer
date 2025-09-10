function x = uniform_on_sphere01_map ( dim_num, n )

%*****************************************************************************80
%
%% uniform_on_sphere01_map() maps uniform points onto the unit sphere.
%
%  Discussion:
%
%    The sphere has center 0 and radius 1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    03 August 2005
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Russell Cheng,
%    Random Variate Generation,
%    in Handbook of Simulation,
%    edited by Jerry Banks,
%    Wiley, 1998, pages 168.
%
%    Reuven Rubinstein,
%    Monte Carlo Optimization, Simulation, and Sensitivity
%    of Queueing Networks,
%    Wiley, 1986, page 234.
%
%  Input:
%
%    integer DIM_NUM, the dimension of the space.
%
%    integer N, the number of points.
%
%  Output:
%
%    real X(DIM_NUM,N), the points.
%
  x = zeros ( dim_num, n );

  for j = 1 : n
%
%  Fill a vector with normally distributed values.
%
    v = randn ( dim_num, 1 );
%
%  Compute the length of the vector.
%
    v_norm = norm ( v );
%
%  Normalize the vector.
%
    x(1:dim_num,j) = v / v_norm;

  end

  return
end
