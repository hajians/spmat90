!> @file test.f90
program test

  use sparsemat

  type(spmat)             :: A, K, Id       ! define two sparse matrix
  integer, parameter      :: len = 20       ! define the length of the matrices
  real (kind=spkind_real) :: f(len)         ! define the force vector
  real (kind=spkind_real) :: u(len)         ! define the solution vector
  real (kind=spkind_real) :: h              ! mesh size
  
  Id = diag_sp(len)             ! building Identity matrix

  K  = diag_sp(len,1) + diag_sp(len,-1) ! building stiffness matrix
  K  = K + (-2.0) * diag_sp(len)        ! building stiffness matrix

  h = 1.0 / (len-1)

  A = h**(-2)*K + 0.25 * Id

  call setvalspmat(A,1,1,1.0)     ! set boundary conditions
  call setvalspmat(A,1,2,0.0)     ! set boundary conditions
  call setvalspmat(A,len,len,1.0) ! set boundary conditions
  call setvalspmat(A,len,len-1,0.0) ! set boundary conditions
  
  call showspmat(A)             ! inspecting the matrix by printing it
                                ! in terminal
  call writespmat("bin/matrix", A)

  f(1:(len/2)) = -1.0; f(1) = 0.0; f(len) = 0.0
  
  call solve_gmres(A, u, f, len, len, 1.0E-8, 1.0E-8)
  call writevector("bin/solution", u)
  
end program test
