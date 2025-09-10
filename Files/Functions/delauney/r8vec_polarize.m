function [ a_normal, a_parallel ] = r8vec_polarize ( n, a, p )

%*****************************************************************************80
%
%% r8vec_polarize() decomposes an R8VEC into normal and parallel components.
%
%  Discussion:
%
%    The (nonzero) vector P defines a direction.
%
%    The vector A can be written as the sum
%
%      A = A_normal + A_parallel
%
%    where A_parallel is a linear multiple of P, and A_normal
%    is perpendicular to P.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 August 2010
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer N, the number of entries in the array.
%
%    real A(N), the vector to be polarized.
%
%    real P(N), the polarizing direction.
%
%  Output:
%
%    real A_NORMAL(N,1), A_PARALLEL(N,1), the normal
%    and parallel components of A.
%
  a = a(:);
  p = p(:);

  p_norm = norm ( p );

  if ( p_norm == 0.0 )
    a_normal(1:n,1) = a(1:n,1);
    a_parallel(1:n,1) = 0.0;
    return
  end

  a_dot_p =  ( a(1:n,1)' * p(1:n,1) ) / p_norm;

  a_parallel(1:n,1) = a_dot_p * p(1:n,1) / p_norm;

  a_normal(1:n,1) = a(1:n,1) - a_parallel(1:n,1);

  return
end
