function i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, label )

%*****************************************************************************80
%
%% i4mat_transpose_print_some() prints some of an I4MAT, transposed.
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
%    integer M, N, the number of rows and columns.
%
%    integer A(M,N), an M by N matrix to be printed.
%
%    integer ILO, JLO, the first row and column to print.
%
%    integer IHI, JHI, the last row and column to print.
%
%    string LABEL, a title.
%
  incx = 10;

  if ( 0 < length ( label ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, '%s\n', label );
  end

  if ( m <= 0 || n <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  (None)\n' );
    return
  end

  for i2lo = max ( ilo, 1 ) : incx : min ( ihi, m )

    i2hi = i2lo + incx - 1;
    i2hi = min ( i2hi, m );
    i2hi = min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    fprintf ( 1, '\n' );
    fprintf ( 1, '  Row: ' );
    for i = i2lo : i2hi
      fprintf ( 1, '%7d  ', i );
    end
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Col\n' );
    fprintf ( 1, '\n' );

    j2lo = max ( jlo, 1 );
    j2hi = min ( jhi, n );

    for j = j2lo : j2hi

      fprintf ( 1, '%5d: ', j );
      for i2 = 1 : inc
        i = i2lo - 1 + i2;
        fprintf ( 1, '%7d  ', a(i,j) );
      end
      fprintf ( 1, '\n' );

    end

  end

  return
end
