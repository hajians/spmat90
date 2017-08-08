!> @file spmat.f90 contains subroutines that manipulate a sparse matrix
!> @author Soheil Hajian
!> @version 1.1
!! last update: August 8, 2017
!! April 5, 2016
!!
!! --------------------------------------------------
!! Available sparse matrix representations:
!! --------------------------------------------------
!!
!! 1. COO (coordinate):                  type(spmat)
!! 2. CSR (compressed storage row):      type(sp_csr)
!!
!! --------------------------------------------------
!! MAIN subroutines to manipulate sparse matrix:
!! -------------------------------------------------- 
! | buildspmat        | builds a sparse matrix in COO format                 |
! | clearspmat        | destroy a sparse matrix in COO format                |
! | writespmat        | writes a sparse matrix of COO format in a file       |
! | showspmat         | shows a sparse matrix of COO format in terminal      |
! | readspmat         | read an entry of a COO matrix and return the value   |
! | insert2spmat      | inserts an entry in a COO matrix at a given position |
! | insertblock2spmat | inserts a block COO matrix into another COO matrix   |
! | SPtranspose       | transpose a COO matrix                               |
! | alpha * Asp       | 'alph' is a real*8 and Asp is a COO matrix           |
! | Asp * vec         | multiplies a COO sparse matrix to a vector           |
! | Asp + Bsp         | sum of two sparse matrices. It might be expensive    |
! | Bsp = Asp         | Assign Asp to Bsp                                    |
! | copysp2sp         | copy a COO matrix into another COO matrix            |
! | buildsp_csr       | builds a CSR matrix                                  |
! | sp_coo2csr        | converts a COO to CSR format                         |
! 
!! @todo: 1. cleanup sp_coo2csr using sort subroutine
!!        2. readspmat using sort
!!        3. add a sort flag in the type(spmat)
!!        4. maybe a finalize subroutine to do sorting and erase blank is good

module sparsemat

  ! gmres solver
  use mgmres_module, only : mgmres_st
  
  implicit none

  !> a parameter to determine how many digits we need
  integer, parameter :: spkind = selected_int_kind (8)

  !> a parameter to set double precision on reals
  integer, parameter :: spkind_real = 8

  !> spmat is a type for sparse matrices in COO format.
  type spmat

     logical :: init = .false.

     !> size of the sparse matrix
     integer (kind=spkind) :: N_i, N_j

     !> size of the i, j, k arrays
     integer (kind=spkind) :: N_array

     !> number of non-zero elements of the i, j, k arrays
     integer (kind=spkind) :: N_nonzero

     integer (kind=spkind), allocatable :: i(:)
     integer (kind=spkind), allocatable :: j(:)

     real (kind=spkind_real), allocatable :: k(:)

  end type spmat

  !> sp_csr is a type for sparse matrices in CSR format.
  type sp_csr
     logical :: init = .false.

     !> size of the sparse matrix
     integer (kind=spkind) :: N_i, N_j

     !> size of the i, j, k arrays
     integer (kind=spkind) :: N_array

     !> size of the row index
     integer (kind=spkind) :: N_row

     !> number of non-zero elements of the j, k arrays
     integer (kind=spkind) :: N_nonzero

     !> number of non-zero in i-array
     integer (kind=spkind) :: N_nonzero_row

     !> compressed row index
     integer (kind=spkind), allocatable :: i(:)
     !> column index
     integer (kind=spkind), allocatable :: j(:)

     !> value array
     real (kind=spkind_real), allocatable :: k(:)

  end type sp_csr

  !> interface operator for multiplication
  !! at the moment there is only scalar-to-matrix exists.
  interface operator(*)
     module procedure scalarMultSp, SpMultVec
  end interface operator(*)

  !> interface operator for summation
  interface operator(+)
     module procedure sumSpSp
  end interface operator(+)

  !> interface operator for summation
  interface operator(-)
     module procedure subtrSpfromSp
  end interface operator(-)

  !> assignement interface for sparse matrices
  interface assignment(=)
     module procedure assignsp2sp, assignsp2full
  end interface assignment(=)

contains

  !> buildspmat initiates a sparse matrix of size N_i, N_j
  !! with possible size of the i, j, k arrays.
  subroutine buildspmat(Asp, N_i, N_j, N_array)
    type(spmat),                     intent(inout) :: Asp
    integer (kind=spkind),           intent(in)    :: N_i
    integer (kind=spkind),           intent(in)    :: N_j
    integer (kind=spkind), optional, intent(in)    :: N_array

    integer (kind=spkind) :: dummy

    if ( Asp%init.eqv..true. ) then
       write(6,*) "error in buildspmat: Asp is init."
       stop
    endif

    Asp%N_i = N_i
    Asp%N_j = N_j

    Asp%N_nonzero = 0

    if ( present(N_array) ) then
       dummy = N_array
    else
       dummy = 1
    endif

    allocate( Asp%i( dummy ) )
    allocate( Asp%j( dummy ) )
    allocate( Asp%k( dummy ) )

    Asp%i = 0
    Asp%j = 0
    Asp%k = 0.0d0

    Asp%N_array = dummy

    Asp%init = .true.

  end subroutine buildspmat

  !> clearspmat clears a sparse matrix type
  subroutine clearspmat(Asp)
    type(spmat), intent(inout) :: Asp

    if (Asp%init.eqv..false.) then
       write(6,*) "error in clearspmat"
       stop
    endif

    deallocate( Asp%i, Asp%j, Asp%k )

    Asp%N_i       = 0
    Asp%N_j       = 0
    Asp%N_array   = 0
    Asp%N_nonzero = 0

    Asp%init = .false.
  end subroutine clearspmat

  !> write vector into file
  subroutine writevector(filename,u)
    character(len=*),        intent(in)    :: filename
    real (kind=spkind_real), intent(in)    :: u(:)

    integer (kind=spkind) :: n
    integer, parameter :: MyUnit = 11

    open(MyUnit, file=filename)

    do n=1, size(u)
       write(MyUnit,*) u(n)
    enddo

    close(MyUnit)

  end subroutine writevector
  
  !> write sp mat in COO format
  subroutine writespmat(filename,Asp)
    character(len=*), intent(in)    :: filename
    type(spmat),      intent(in)    :: Asp

    integer (kind=spkind) :: n
    integer, parameter :: MyUnit = 11

    if ( Asp%init.eqv..false. ) then
       write(6,*) "Sparse matrix is not init in writespmat"
       stop
    endif

    open(MyUnit, file=filename)
    !! write dimension of the matrix and # of non-zero elements
    write(MyUnit,*) Asp%N_i, Asp%N_j, Asp%N_nonzero

    do n=1, Asp%N_nonzero
       write(MyUnit,*) Asp%i(n), Asp%j(n), Asp%k(n)
    enddo

    close(MyUnit)

  end subroutine writespmat

  !> show sparse matrix
  subroutine showspmat(Asp)
    type(spmat), intent(in) :: Asp

    integer (kind=spkind) :: n

    if ( Asp%init.eqv..false. ) then
       write(6,*) "Sparse matrix is not init in showspmat"
       stop
    endif

    write(6,*) "size of sparse matrix:     ", Asp%N_i, "*", Asp%N_j
    write(6,*) "number of non-zero entries:", Asp%N_nonzero
    write(6,*) "size of i, j, k arrays:     ", Asp%N_array

    do n=1, Asp%N_nonzero
       write(6,*) Asp%i(n), ",", Asp%j(n), ",", Asp%k(n)
    enddo

  end subroutine showspmat

  !> delete blank space of the arrays
  subroutine deleteblankspmat(Asp)
    type(spmat), intent(inout) :: Asp

    integer (kind=spkind), dimension( Asp%N_nonzero ) :: dummy_int
    real (kind=spkind_real),      dimension( Asp%N_nonzero ) :: dummy_real

    dummy_int = Asp%i(1:Asp%N_nonzero)
    deallocate( Asp%i )
    allocate( Asp%i( Asp%N_nonzero ) )
    Asp%i = dummy_int

    dummy_int = Asp%j(1:Asp%N_nonzero)
    deallocate( Asp%j )
    allocate( Asp%j( Asp%N_nonzero ) )
    Asp%j = dummy_int

    dummy_real = Asp%k(1:Asp%N_nonzero)
    deallocate( Asp%k )
    allocate( Asp%k( Asp%N_nonzero ) )
    Asp%k = dummy_real  

    Asp%N_array = Asp%N_nonzero
  end subroutine deleteblankspmat

  !> insert blank space into the arrays
  subroutine insertblank2spmat(Asp,N_blank)
    type(spmat), intent(inout) :: Asp
    integer (kind=spkind) , intent(in) :: N_blank

    integer (kind=spkind), dimension( Asp%N_nonzero ) :: dummy_int
    real (kind=spkind_real),      dimension( Asp%N_nonzero ) :: dummy_real

    Asp%N_array = Asp%N_array + N_blank

    dummy_int = Asp%i(1:Asp%N_nonzero)
    !
    deallocate( Asp%i )
    allocate( Asp%i(Asp%N_array) )
    Asp%i(1:Asp%N_nonzero) = dummy_int

    dummy_int = Asp%j(1:Asp%N_nonzero)
    !
    deallocate( Asp%j )
    allocate( Asp%j(Asp%N_array) )
    Asp%j(1:Asp%N_nonzero) = dummy_int

    dummy_real = Asp%k(1:Asp%N_nonzero)
    !
    deallocate( Asp%k )
    allocate( Asp%k(Asp%N_array) )
    Asp%k(1:Asp%N_nonzero) = dummy_real

  end subroutine insertblank2spmat

  !> insert an element in the sparse matrix
  subroutine insert2spmat(Asp,i,j,val)
    type(spmat),           intent(inout) :: Asp
    integer (kind=spkind), intent(in)    :: i
    integer (kind=spkind), intent(in)    :: j
    real (kind=spkind_real),      intent(in)    :: val

    integer (kind=spkind) :: N_blank

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in insert2spmat: init"
       stop
    endif

    if ( (Asp%N_i < i) .or. (Asp%N_j < j) ) then
       write(6,*) "error in insert2spmat: N_i < i or N_j < j"
       stop
    endif

    if ( Asp%N_array.eq.Asp%N_nonzero ) then
       N_blank = max( Asp%N_nonzero/2, 1 )
       call insertblank2spmat(Asp,N_blank)
    endif

    Asp%N_nonzero = Asp%N_nonzero + 1

    Asp%i( Asp%N_nonzero ) = i
    Asp%j( Asp%N_nonzero ) = j

    Asp%k( Asp%N_nonzero ) = val
  end subroutine insert2spmat

  !> insertblock2spmat adds a block matrix in the position (i,j)
  subroutine insertblock2spmat(Asp,i,j,Ablc)
    type(spmat),           intent(inout) :: Asp
    integer (kind=spkind), intent(in)    :: i
    integer (kind=spkind), intent(in)    :: j
    type(spmat),           intent(in)    :: Ablc

    integer (kind=spkind) :: id

    if ( (Asp%init.eqv..false.) .or. (Ablc%init.eqv..false.) ) then
       write(6,*) "error in insertblock2spmat: Asp or Ablc is not init"
       stop
    endif

    if ( (Ablc%N_i > (Asp%N_i-i+1)) .or. (Ablc%N_j > (Asp%N_j-j+1)) ) then
       write(6,*) "error in insertblock2spmat: Ablc cannot be inserted", &
            Ablc%N_i , (Asp%N_i-i+1), Ablc%N_j, (Asp%N_j-j+1)
       stop
    endif

    do id=1, Ablc%N_nonzero
       call insert2spmat(Asp, i+Ablc%i(id)-1, j+Ablc%j(id)-1, Ablc%k(id) )
    enddo

  end subroutine insertblock2spmat

  !> transpose a sprase matrix in COO format
  function SPtranspose(Asp) result(Atsp)
    type(spmat), intent(in) :: Asp
    type(spmat) :: Atsp

    type(spmat) :: Adummy

    call copysp2sp(Asp,Atsp)
    call copysp2sp(Asp,Adummy)

    Atsp%N_i = Adummy%N_j
    Atsp%N_j = Adummy%N_i

    Atsp%i = Adummy%j
    Atsp%j = Adummy%i

    call clearspmat(Adummy)

  end function SPtranspose

  !> scalarMultSp multiplies a real number alpha to the sparse matrix Asp
  function scalarMultSp(alpha, Asp) result(alphAsp)
    real (kind=spkind_real),      intent(in)    :: alpha
    type(spmat),           intent(in)    :: Asp

    type(spmat)                          :: alphAsp

    call copysp2sp(Asp, alphAsp)

    alphAsp%k = alpha * alphAsp%k

  end function scalarMultSp

  !> multiplies a spmat to a vector
  function SpMultVec(Asp, vecX) result(vecY)
    real (kind=spkind_real), intent(in) :: vecX(:)
    type(spmat),      intent(in) :: Asp

    real (kind=spkind_real), allocatable :: vecY(:)
    real (kind=spkind_real), dimension( size(vecX) ) :: vecDummy

    integer (kind=spkind) :: n
    integer (kind=spkind) :: i, j

    if (Asp%init.eqv..false.) then
       write(6,*) "error in SpMultVec: Asp is not init", Asp%init
       stop
    endif

    if ( Asp%N_j .ne. size(VecX) ) then
       write(6,*) "error in SpMultVec: rank mismatch", Asp%N_j, size(VecX)
       stop
    endif

    vecDummy = vecX

    allocate( vecY(Asp%N_i) );    vecY     = 0.0d0

    do n=1, Asp%N_nonzero
       i = Asp%i(n)
       j = Asp%j(n)

       vecY(i) = vecY(i) + Asp%k(n) * vecDummy(j)
    enddo

  end function SpMultVec

  !> sumSpSp sums two sparse matrices
  function sumSpSp(Asp, Bsp) result(Csp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(in)    :: Bsp
    type(spmat)                          :: Csp

    integer (kind=spkind) :: i

    if ( (Asp%init.eqv..false.) .or. (Bsp%init.eqv..false.) ) then
       write(6,*) "error in sumSpSp: Asp or Bsp is not init"
       stop
    endif

    if ( (Asp%N_i.ne.Bsp%N_i) .or. (Asp%N_j.ne.Bsp%N_j) ) then
       write(6,*) "error in sumSpSp: dim not compatible"
       stop
    endif

    call copysp2sp(Asp,Csp)

    do i=1, Bsp%N_nonzero
       call insert2spmat(Csp, Bsp%i(i), Bsp%j(i), Bsp%k(i))
    enddo

    !! this part is very slow! needs to be fixed with sorting
    call sum_redundant_fast(Csp, Csp)

    call deleteblankspmat(Csp)

  end function sumSpSp

  !> substract sparse matrix from sparse matrix: Asp - Bsp
  function subtrSpfromSp(Asp,Bsp) result(Csp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(in)    :: Bsp
    type(spmat)                          :: Csp

    Csp = sumSpSp(Asp, (-1.0_spkind_real) * Bsp)

  end function subtrSpfromSp

  !> read an entry of the sparse matrix.
  !! we just do a linear search in the arrays.
  subroutine readspmat(Asp,i,j,val,flagOut,nOut)
    type(spmat),                     intent(in)    :: Asp
    integer (kind=spkind),           intent(in)    :: i
    integer (kind=spkind),           intent(in)    :: j
    real (kind=spkind_real),         intent(inout) :: val
    integer,               optional, intent(inout) :: flagOut
    integer (kind=spkind), optional, intent(inout) :: nOut
    integer (kind=spkind) :: n
    integer               :: flag

    if ( (i>Asp%N_i).or.(j>Asp%N_j).or.(Asp%init.eqv..false.) ) then
       write(6,*) "error in readspmat: "
    endif

    val  = 0.0d0
    flag = 0
    
    do n=1, Asp%N_nonzero
       !
       if ( Asp%i(n).eq.i ) then
          if ( Asp%j(n).eq.j ) then
             !
             if (present(nOut)) then
                nOut = n
             endif
             val = val + Asp%k(n)
             flag = flag + 1
             !
          endif
       endif
       !
    enddo

    if (flag.ge.2) then
       write(6,*) "warning in readspmat: multiple entries found."
    endif

    if (present(flagOut)) then
       flagOut = flag
    endif

   
  end subroutine readspmat

  !> setvalspmat set the value of a particular entry into a given value
  subroutine setvalspmat(Asp, i, j, val)
    type(spmat),                     intent(inout) :: Asp
    integer (kind=spkind),           intent(in)    :: i
    integer (kind=spkind),           intent(in)    :: j
    real (kind=spkind_real),         intent(in)    :: val

    integer :: flag, n
    real (kind=spkind_real) :: valDummy
    
    call readspmat(Asp, i, j, valDummy, flag, n)

    if (flag > 1) then
       write(6,*) "error in setvalspmat: multiple entries found"
       stop
    elseif (flag.eq.1) then
       Asp%k(n) = val
    else
       call insert2spmat(Asp, i, j, val)
    endif
    
  end subroutine setvalspmat
  
  !> copysp2sp copies a sparse matrix into another one.
  subroutine copysp2sp(Asp,Bsp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(inout) :: Bsp

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in copysp2sp: ", Asp%init
       stop
    endif
    !
    if (Bsp%init.eqv..true. ) then
       write(6,*) "warning in copysp2sp: Bsp is erased and re-init."
       write(6,*) "if you are copying Asp to Asp: ", &
            "This is not funny! You will get error!"
       call clearspmat( Bsp )
    endif

    call buildspmat(Bsp, Asp%N_i, Asp%N_j, Asp%N_array)

    Bsp%N_nonzero = Asp%N_nonzero

    Bsp%i(1:Asp%N_nonzero) = Asp%i(1:Asp%N_nonzero)
    Bsp%j(1:Asp%N_nonzero) = Asp%j(1:Asp%N_nonzero)
    Bsp%k(1:Asp%N_nonzero) = Asp%k(1:Asp%N_nonzero)

  end subroutine copysp2sp

  !> assignsp2sp copies Asp to Bsp and is used for interface only
  subroutine assignsp2sp(Bsp,Asp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(inout) :: Bsp

    type(spmat) :: Acopy

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in copysp2sp: ", Asp%init
       stop
    endif

    call copysp2sp(Asp, Acopy)

    if ( Bsp%init.eqv..true. ) then
       call clearspmat(Bsp)
    endif

    call copysp2sp(Acopy, Bsp)

  end subroutine assignsp2sp

  !> assign sparse matrix to an allocatable full matrix
  subroutine assignsp2full(Bfull,Asp)
    type(spmat),                   intent(in)    :: Asp
    real (kind=spkind_real), allocatable, intent(inout) :: Bfull(:,:)

    integer :: n

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in assignsp2full: ", Asp%init
       stop
    endif

    if ( allocated(Bfull) ) then
       deallocate( Bfull )
    endif

    allocate( Bfull(Asp%N_i, Asp%N_j) )

    Bfull = 0.0d0

    do n=1, Asp%N_nonzero
       Bfull( Asp%i(n), Asp%j(n) ) = Asp%k(n)
    enddo

  end subroutine assignsp2full

  !> sum_redundant finds multiple entries at the same position
  !! in the sparse matrix form and sum them. This feature is used
  !! in finite element method, where local assembly is used.
  !! you can use same matrix for both arguments:
  !! call sum_redundant(A1,A1)
  subroutine sum_redundant(Asp,Bsp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(inout) :: Bsp

    logical,     allocatable, dimension(:)   :: flag
    type(spmat)                              :: Acopy

    integer (kind=spkind) :: m, n
    integer (kind=spkind) :: i, j
    real (kind=spkind_real) :: val

    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in sum_redundant: ", Asp%init
       stop
    endif

    call copysp2sp(Asp,Acopy)   ! copy Asp to Acopy and use Acopy from now on

    if (Bsp%init.eqv..true. ) then
       write(6,*) "warning in sum_redundant: Bsp is erased and re-init."
       call clearspmat( Bsp )
    endif

    call buildspmat(Bsp, Acopy%N_i, Acopy%N_j, Acopy%N_nonzero)


    allocate( flag(Acopy%N_nonzero) ) ! allocate flag
    flag = .false.                         ! -1 means not visited

    do n=1, Acopy%N_nonzero
       !
       if ( flag(n).eqv..false. ) then ! if not visited
          i   = Acopy%i(n)
          j   = Acopy%j(n)
          val = Acopy%k(n)
          !
          do m=(n+1), Acopy%N_nonzero
             if ( (Acopy%i(m).eq.i) ) then
                if ( (Acopy%j(m).eq.j) ) then
                   flag(m) = .true.
                   val = val + Acopy%k(m)
                endif
             endif
          enddo
          !
          call insert2spmat(Bsp,i,j,val) ! finally insert into Bsp
       endif
       !
    enddo

  end subroutine sum_redundant

  !> this is supposed to be a fast version of the sum_redundant
  subroutine sum_redundant_fast(Asp,Csp)
    type(spmat),           intent(in)    :: Asp
    type(spmat),           intent(inout) :: Csp

    logical,     allocatable, dimension(:)   :: flag
    type(spmat)                              :: Acopy

    integer (kind=spkind) :: m, n
    integer (kind=spkind) :: i, j
    real (kind=spkind_real)      :: val

    type(spmat) :: Bsp
    
    if ( Asp%init.eqv..false. ) then
       write(6,*) "error in sum_redundant_fast: ", Asp%init
       stop
    endif

    call copysp2sp(Asp,Acopy)   ! copy Asp to Acopy and use Acopy from now on

    ! if (Csp%init.eqv..true. ) then
    !    write(6,*) "warning in sum_redundant_fast: Bsp is erased and re-init."
    !    call clearspmat( Bsp )
    ! endif

    call buildspmat(Bsp, Acopy%N_i, Acopy%N_j, Acopy%N_nonzero)

    allocate( flag(Acopy%N_nonzero) ) ! allocate flag
    flag = .false.                    ! .false. means not-visited

    call sortspmat(Acopy)       ! sort according to the row and column index   

    do n=1, Acopy%N_nonzero
       if ( flag(n).eqv..false.) then
          !
          i   = Acopy%i(n)
          j   = Acopy%j(n)
          val = Acopy%k(n)
          do m=(n+1), Acopy%N_nonzero
             if ( (i.eq.Acopy%i(m)) .and. (j.eq.Acopy%j(m)) ) then
                val = val + Acopy%k(m)
                flag(m) = .true.
             else
                goto 12345
             endif
          enddo
          !
12345     continue
          call insert2spmat(Bsp, i, j, val) ! finally insert
          !
       endif
    enddo

    Csp = Bsp
    call clearspmat(Acopy)
    call clearspmat(Bsp)
    
    deallocate( flag )
  end subroutine sum_redundant_fast

  !> sorts the entries of Asp according to the row index first and then
  !! column index.
  subroutine sortspmat(Asp)
    type(spmat), intent(inout) :: Asp

    integer (kind=spkind), allocatable :: row_index(:), col_index(:)
    integer (kind=spkind), allocatable :: interval(:,:)

    integer (kind=spkind) :: i, j
    integer (kind=spkind) :: n, m
    integer (kind=spkind) :: counter

    if ( (Asp%init.eqv..false.) .or. Asp%N_nonzero.eq.0 ) then
       write(6,*) "error in sortspmat: ", Asp%init, Asp%N_nonzero
       stop
    endif

    allocate( row_index(Asp%N_nonzero) )

    call quick_sort( Asp%i(1:Asp%N_nonzero) , row_index )

    Asp%j(1:Asp%N_nonzero) = Asp%j( row_index )
    Asp%k(1:Asp%N_nonzero) = Asp%k( row_index )

    allocate( interval(Asp%N_i,2) )
    interval = 0

    i = Asp%i(1);    counter = 1;     interval(i,1) = counter
    !
    do n=2, Asp%N_nonzero
       j = Asp%i(n)
       if (i.ne.j) then
          interval(i,2) = counter
          i = j
          interval(i,1) = counter + 1
       endif
       counter = counter + 1
    enddo
    interval(i,2) = counter

    deallocate( row_index )

    ! sorting column
    do i=1, Asp%N_i
       n = interval(i,1)
       m = interval(i,2)
       if ( n.ne.0 ) then
          allocate( col_index(m-n+1) )
          allocate( row_index(m-n+1) ); row_index = (/ (j, j=n, m) /)
          !
          call quick_sort( Asp%j(n:m), col_index )
          Asp%k(n:m) = Asp%k( row_index(col_index) )
          !
          deallocate(col_index)
          deallocate(row_index)
       endif
    enddo

    deallocate( interval )

  end subroutine sortspmat

  !> build a sparse matrix in CSR form.
  !! N_row is an estimate for the number 
  subroutine buildsp_csr(Asp,N_i,N_j,N_array,N_row)
    type(sp_csr),                    intent(inout) :: Asp
    integer (kind=spkind),           intent(in)    :: N_i
    integer (kind=spkind),           intent(in)    :: N_j
    integer (kind=spkind), optional, intent(in)    :: N_array
    integer (kind=spkind), optional, intent(in)    :: N_row

    integer (kind=spkind) :: dummy
    integer (kind=spkind) :: dummy_row

    if ( Asp%init.eqv..true. ) then
       write(6,*) "error in buildspmat: Asp is init."
       stop
    endif

    Asp%N_i = N_i
    Asp%N_j = N_j

    Asp%N_nonzero = 0
    Asp%N_nonzero_row = 0

    if ( present(N_array) ) then
       dummy = N_array
    else
       dummy = 1
    endif

    if ( present(N_row) ) then
       dummy_row = N_row
    else
       dummy_row = 1
    endif

    allocate( Asp%i( dummy_row ) )
    allocate( Asp%j( dummy ) )
    allocate( Asp%k( dummy ) )

    Asp%i = 0
    Asp%j = 0
    Asp%k = 0.0d0

    Asp%N_array = dummy

    Asp%N_row   = dummy_row

    Asp%init = .true.    

  end subroutine buildsp_csr

  !> convert COO format to CSR format
  subroutine sp_coo2csr(Aorg,Bsp)
    type(spmat),  intent(in)    :: Aorg
    type(sp_csr), intent(inout) :: Bsp

    type(spmat) :: Asp
    integer, allocatable :: row_index(:), col_index(:), sort_index(:)

    integer :: row, dummy
    integer :: i

    if ( (Aorg%init.eqv..false.) .or. (Bsp%init.eqv..true.) ) then
       write(6,*) "error in sp_coo2csr: ", Aorg%init, Bsp%init
       stop
    endif

    if ( (Aorg%N_nonzero.eq.0) ) then
       write(6,*) "error in sp_coo2csr: Aorg is empty"
       stop
    endif

    !! copy Aorg to Asp. Asp is guaranteed to have no blank space
    call copysp2sp(Aorg,Asp)

    !! first sort the row components of the Asp

    allocate( sort_index(Asp%N_nonzero) )
    allocate( row_index(Asp%N_nonzero) ) ! this will be deallocate later on


    ! init Bsp in CSR format. Size of i-array is equal N_i+1
    call buildsp_csr(Bsp, Asp%N_i, Asp%N_j, Asp%N_nonzero, Asp%N_i+1)

    row_index = Asp%i

    call quick_sort(row_index, sort_index)

    Bsp%j = Asp%j(sort_index) ! copy the column index
    Bsp%k = Asp%k(sort_index) ! copy the values

    Bsp%N_nonzero = Asp%N_nonzero ! copy the non-zero

    dummy             = row_index(1)            ! first row entry in the Asp
    Bsp%N_nonzero_row = 1
    Bsp%i(1)          = 1

    do i=2, Asp%N_nonzero
       row = row_index(i)
       if (row.ne.dummy) then
          ! empty row check
          if ( row > (dummy+1) ) then
             write(6,*) "error in sp_coo2csr: empty row in ", dummy
             stop
          endif

          dummy = row
          Bsp%N_nonzero_row = Bsp%N_nonzero_row + 1
          Bsp%i( Bsp%N_nonzero_row ) = i
       endif
    enddo


    Bsp%N_nonzero_row = Bsp%N_nonzero_row + 1 ! last entry
    Bsp%i(Bsp%N_nonzero_row) = Bsp%N_nonzero + 1

    !! delete empty spaces in Bsp%i
    deallocate( row_index )
    allocate( row_index(Bsp%N_nonzero_row) )
    row_index = Bsp%i(1:Bsp%N_nonzero_row)


    deallocate( Bsp%i );
    allocate( Bsp%i(Bsp%N_nonzero_row) )

    Bsp%i     = row_index
    Bsp%N_row = Bsp%N_nonzero_row

    !! sorting the column index
    do i=1, (Bsp%N_nonzero_row-1)
       dummy = Bsp%i(i+1) - Bsp%i(i)
       allocate( col_index(dummy) )
       call quick_sort( Bsp%j( Bsp%i(i):(Bsp%i(i+1)-1) ), col_index )
       Bsp%k( Bsp%i(i):(Bsp%i(i+1)-1) ) = Bsp%k( col_index + Bsp%i(i) - 1 )
       deallocate( col_index )
    enddo


  end subroutine sp_coo2csr

  !> diag_sp produce a diagonal sparse matrix 
  function diag_sp(n,k) result(out)
    integer (kind=spkind)           :: n
    integer (kind=spkind), optional :: k
    type(spmat)                     :: out

    integer :: i, j

    if (present(k)) then
       j = abs(k)
    else
       j = 0
    endif

    if ( abs(j) > n ) then
       write(6,*) "error in diag_sp"
       stop
    endif

    call buildspmat(out, n, n, n)
    do i=1, n-j
       call insert2spmat(out,i,i+j,1.0_spkind_real)
    enddo

    if (present(k)) then
       if (k<0) then
          out = SPtranspose(out)
       endif
    end if
    
  end function diag_sp
  
  !> diagsp takes diagonal components of Asp
  subroutine diagsp(Asp,Dsp,flag)
    type(spmat), intent(in)              :: Asp
    type(spmat), intent(inout)           :: Dsp
    logical,     intent(in), optional    :: flag

    real (kind=spkind_real) :: val
    integer :: i

    if (Asp%init.eqv..false.) then
       write(6,*) "error in diagsp: Asp is not init", Asp%init
    endif

    ! if flag is not present check the Dsp
    if (present(flag).eqv..false.) then
       if ( Dsp%init.eqv..true. ) then
          call clearspmat( Dsp )
       endif
       !
       call buildspmat(Dsp, Asp%N_i, Asp%N_i, Asp%N_i)
       !
    endif

    do i=1, Asp%N_i
       call readspmat(Asp,i,i,val)
       call insert2spmat(Dsp,i,i,val)
    enddo

  end subroutine diagsp

  !> solve_gmres uses mgmres fortran library to solve A*x=f
  subroutine solve_gmres(Asp, x, rhs, ite_max, ite_in, tol_abs, tol_rel)
    type(spmat),      intent(in)    :: Asp
    real (kind=spkind_real), intent(inout) :: x(:)
    real (kind=spkind_real), intent(in)    :: rhs(:)
    integer,          intent(in)    :: ite_max
    integer,          intent(in)    :: ite_in
    real (kind=spkind_real), intent(in)    :: tol_abs
    real (kind=spkind_real), intent(in)    :: tol_rel

    if (Asp%init.eqv..false.) then
       write(6,*) "error in solve_gmres: Asp is not init"
       stop
    elseif (Asp%N_i .ne. Asp%N_j ) then
       write(6,*) "error in solve_gmres: Asp is not square"
       stop
    elseif ( (Asp%N_i.ne.size(x)) .or. (Asp%N_i.ne.size(rhs)) ) then
       write(6,*) "error in solve_gmres: vectors are not of system size"
       stop
    endif

    call mgmres_st(Asp%N_i, Asp%N_nonzero, Asp%i, Asp%j, Asp%k, &
         x, rhs, ite_max, ite_in, tol_abs, tol_rel)

  end subroutine solve_gmres

  !! ----------------------------------------- !!
  !!                                           !!
  !!                                           !!
  !! Here comes external subroutines needed    !!
  !! for this module:                          !!
  !! 1. quick_sort                             !!
  !!                                           !!
  !!                                           !!
  !! ----------------------------------------- !!

  RECURSIVE SUBROUTINE quick_sort(list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.

    IMPLICIT NONE
    INTEGER, DIMENSION (:), INTENT(IN OUT)  :: list
    INTEGER, DIMENSION (:), INTENT(OUT)  :: order

    ! Local variable
    INTEGER :: i

    DO i = 1, SIZE(list)
       order(i) = i
    END DO

    CALL quick_sort_1(1, SIZE(list))

  CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

      INTEGER, INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER             :: i, j, itemp
      INTEGER             :: reference, temp
      INTEGER, PARAMETER  :: max_simple_sort_size = 6

      IF (right_end < left_end + max_simple_sort_size) THEN
         ! Use interchange sort for small lists
         CALL interchange_sort(left_end, right_end)

      ELSE
         ! Use partition ("quick") sort
         reference = list((left_end + right_end)/2)
         i = left_end - 1; j = right_end + 1

         DO
            ! Scan list from left end until element >= reference is found
            DO
               i = i + 1
               IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
               j = j - 1
               IF (list(j) <= reference) EXIT
            END DO


            IF (i < j) THEN
               ! Swap two out-of-order elements
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            ELSE IF (i == j) THEN
               i = i + 1
               EXIT
            ELSE
               EXIT
            END IF
         END DO

         IF (left_end < j) CALL quick_sort_1(left_end, j)
         IF (i < right_end) CALL quick_sort_1(i, right_end)
      END IF

    END SUBROUTINE quick_sort_1


    SUBROUTINE interchange_sort(left_end, right_end)

      INTEGER, INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER             :: i, j, itemp
      INTEGER             :: temp

      DO i = left_end, right_end - 1
         DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            END IF
         END DO
      END DO

    END SUBROUTINE interchange_sort

  END SUBROUTINE quick_sort

end module sparsemat
