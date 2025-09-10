function r8vec_print ( n, a, label )

%*****************************************************************************80
%
%% r8vec_print() prints an R8VEC.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    31 July 2018
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer N, the dimension of the vector.
%
%    real A(N), the vector to be printed.
%
%    string LABEL, a title.
%
  if ( 0 < length ( label ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, '%s\n', label );
  end
  fprintf ( 1, '\n' );

  for i = 1 : n
    fprintf ( 1, '%6d: %12g\n', i, a(i) );
  end

  return
end
