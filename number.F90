module number

  implicit none

!from BGW
  integer, parameter :: dp = kind(1.0d0)

!from icm.cpp
  integer,parameter :: NCELL = 2
  integer,parameter :: NRAND = 2500000
  integer,parameter :: RAND_MAX = 2147483647
  real(dp),parameter :: HARTREE = 27.21138505 
  real(dp),parameter :: EPS5 = 1.0E-5
  real(dp),parameter :: EPS6 = 1.0E-6
  real(dp),parameter :: EPS7 = 1.0E-7
  real(dp),parameter :: EPS8 = 1.0E-8
  real(dp),parameter :: EPS9 = 1.0E-9
  real(dp),parameter :: INF9 = 1.0E+9
  real(dp),parameter :: BOHR = 0.52917721092
  real(dp), parameter :: PI = 3.141592653589793238462643383279
  real(dp), parameter :: FPI = 4 * PI
  complex(kind(0d0)),parameter :: IMG = (0.0,1.0)

!from BGW
  real(DP), parameter :: PI_D = 3.1415926535897932384626433832795_dp
  real(DP), parameter :: TOL_Small = 1.0d-6
  real(DP), parameter :: TOL_Zero = 1.0d-12
  real(DP), parameter :: TOL_Degeneracy = 1.0d-6


end module number
