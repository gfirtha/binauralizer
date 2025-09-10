function value = i4_wrap ( value, lo, hi )

%*****************************************************************************80
%
%% i4_wrap() forces an integer to lie between given limits by wrapping.
%
%  Example:
%
%    LO = 4, HI = 8
%
%    In   Out
%
%    -2     8
%    -1     4
%     0     5
%     1     6
%     2     7
%     3     8
%     4     4
%     5     5
%     6     6
%     7     7
%     8     8
%     9     4
%    10     5
%    11     6
%    12     7
%    13     8
%    14     4
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 June 2020
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer VALUE, an integer value.
%
%    integer LO, HI, the desired bounds for the integer value.
%
%  Output:
%
%    integer VALUE, a "wrapped" version of IVAL.
%
  value = lo + mod ( value - lo, hi + 1 - lo );

  return
end
