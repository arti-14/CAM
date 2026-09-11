!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_emulator.f90
!!
!! \brief
!! Updraft emulator module
!! 
!!
!! \author Noora Hyttinen (FMI) noora.hyttinen@fmi.fi
!!
!!
!! \details
!! Includes modules m_util, m_cov, m_noise, m_gp, m_cov_sqexp, m_cov_sqexp_param4, 
!! m_cov_lin, m_cov_linsqexp, m_cov_all, m_noise_value_only, m_noise_param2, m_noise_all, m_gp_dense
!! from https://github.com/ots22/gpf (last access 1st April 2025)
!! and mo_emulator for updraft standard deviation predictions
!!
!! \belongs_to
!!
!! \copyright
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module m_util
  implicit none
  integer, parameter ::  dp = selected_real_kind(13,300)

  integer, parameter :: max_name_len=100

contains
  ! Solve a system of linear equations A*x = b for a matrix A and
  ! vectors x and b.  Implemented as a wrapper around LAPACK dgesv.
  function solve(M,b) result(x)
    real(dp), intent(in) :: M(:,:)
    real(dp) :: b(:)
    real(dp), dimension(size(M,1),size(M,2)) :: A
    real(dp), dimension(size(b,1)) :: x
    integer ipiv(size(b,1)), N, info

    ! check M is square and b conforms
    if (size(M,1) /= size(M,2)) then
       write (0,*) "solve: Array M passed to solve should be square.  Got ", &
            size(M,1), "x", size(M,2)
       error stop
    end if

    N = size(M, 1)

    ! A and x are overwritten on output of dgesv, so copy
    A = M
    x = b
    call dgesv(N, 1, A, N, ipiv, x, N, info)

    ! check for success
    if (info /= 0) then
       write (0,*) "solve: dgesv returned an error code (", info, ")"
       error stop
    end if

  end function solve

  ! in place inverse
  subroutine ninv(A)
    real(dp), dimension(:,:), intent(inout) :: A
    integer, dimension(:), allocatable :: piv
    real(dp), dimension(1) :: lwork_real
    real(dp), dimension(:), allocatable :: work
    integer :: lwork, N, info

    N = size(A,1)

    allocate(integer :: piv(N))

    call dgetrf(N, N, A, N, piv, info)
    
    call dgetri(N, A, N, piv, lwork_real, -1, info)
    lwork = nint(lwork_real(1))
    allocate(real(dp) :: work(lwork))
    call dgetri(N, A, N, piv, work, lwork, info)
    
  end subroutine ninv

  function logdet(A)
    real(dp) :: logdet
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:,:), allocatable :: tmp
    integer, dimension(:), allocatable :: ipiv
    integer :: N, info, i

    N = size(A,1)

    allocate(real(dp) :: tmp(N, N))
    allocate(integer :: ipiv(N))

    tmp = A
    call dgetrf(N, N, tmp, N, ipiv, info)

    logdet = 0.0_dp
    do i=1,N
       logdet = logdet + log(abs(tmp(i,i)))
    end do

  end function logdet
  
end module m_util


module m_cov
use m_util
implicit none

private
public cov_fn

  type, abstract :: cov_fn
   contains
     procedure(ntheta_required), deferred, nopass :: ntheta_required
     procedure(cov_val), deferred, nopass :: cov_val
     procedure(dcov_x1), deferred, nopass :: dcov_x1
     procedure(dcov_x2), deferred, nopass :: dcov_x2
     procedure(d2cov_xx), deferred, nopass :: d2cov_xx
     procedure :: cov
  end type cov_fn

  abstract interface
     pure function ntheta_required(dims)
       import cov_fn
       integer :: ntheta_required
       integer, intent(in) :: dims
     end function ntheta_required

     pure function cov_val(x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: cov_val
       real(dp), intent(in), dimension(:) :: x, y, hypers
     end function cov_val
     
     pure function dcov_x1(n,x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: dcov_x1
       real(dp), intent(in), dimension(:) :: x, y, hypers
       integer, intent(in) :: n
     end function dcov_x1
     
     pure function dcov_x2(n,x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: dcov_x2
       real(dp), intent(in), dimension(:) :: x, y, hypers
       integer, intent(in) :: n
     end function dcov_x2
     
     pure function d2cov_xx(n,m,x,y,hypers)
       use m_util, only: dp
       import cov_fn
       real(dp) :: d2cov_xx
       real(dp), intent(in), dimension(:) :: x, y, hypers
       integer, intent(in) :: n, m
     end function d2cov_xx
  end interface

contains

  pure function cov(cf,n,m,x,y,hypers)
    class(cov_fn), intent(in) :: cf
    real(dp) :: cov
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    if (n.eq.0 .and. m.eq.0) then
       cov = cf%cov_val(x,y,hypers)
    else if (n.eq.0) then
       cov = cf%dcov_x2(m,x,y,hypers)
    else if (m.eq.0) then
       cov = cf%dcov_x1(n,x,y,hypers)
    else
       cov = cf%d2cov_xx(n,m,x,y,hypers)
    end if
  end function cov

end module m_cov


module m_noise
  use m_util
  implicit none
  
  private
  public noise_model
  
  type, abstract :: noise_model
   contains
     procedure(nparams_required), deferred, nopass :: nparams_required
     ! returns a vector of values for the value and each dimension of the gradient
     procedure(noise), deferred, nopass :: noise
  end type noise_model

  abstract interface
     pure function nparams_required(dims)
       import noise_model
       integer, intent(in) :: dims
       integer nparams_required       
     end function nparams_required

     pure function noise(obs_type, params)
       use m_util, only: dp
       import noise_model
       real(dp) noise
       integer, intent(in) :: obs_type
       real(dp), intent(in) :: params(:)
     end function noise
  end interface

end module m_noise


module m_gp
  use m_util
  use m_cov
  use m_noise

  private
  public BaseGP, nlog_lik, set_hyperparams

  type, abstract :: BaseGP
     ! The noise hyperparameter(s). This is passed to the noise model
     ! (`noise_model'), which determines its precise meaning.
     real(dp), dimension(:), allocatable :: nu 
     ! covariance hyperparameters
     ! meaning depends on covariance function `covariance'
     real(dp), dimension(:), allocatable :: theta
     ! inputs
     real(dp), dimension(:,:), allocatable :: x
     ! types of the observations
     integer, dimension(:), allocatable :: obs_type
     ! observations
     real(dp), dimension(:), allocatable :: t
     ! the covariance function
     class(cov_fn), allocatable :: covariance
     ! the noise model
     class(noise_model), allocatable :: noise_model
   contains
     procedure(log_lik), deferred :: log_lik
     procedure(update_matrices), deferred :: update_matrices
     procedure(predict), deferred :: predict
     procedure(write_out), deferred :: write_out
  end type BaseGP

  abstract interface
     ! The log likelihood of the hyperparameters
     function log_lik(this)
       use m_util, only: dp
       import BaseGP
       class(BaseGP), intent(in) :: this
       real(dp) log_lik
     end function log_lik

     ! Helper routine to update the internal state.  Called when an
     ! observation or hyperparameter changes and the covariance matrix
     ! must be recomputed.
     subroutine update_matrices(this)
       import BaseGP
       class(BaseGP), intent(inout) :: this
     end subroutine update_matrices

     ! Make a prediction of the underlying function value at
     ! coordinate `xnew'.
     function predict(this, xnew, obs_type_new)
       use m_util, only: dp
       import BaseGP
       class(BaseGP), intent(in) :: this
       real(dp) predict
       real(dp), dimension(:), intent(in) :: xnew
       integer, optional, intent(in) :: obs_type_new
     end function predict

     ! Serialize to a file
     subroutine write_out(this, filename)
       import BaseGP
       class(BaseGP), intent(in) :: this
       character(len=*), intent(in) :: filename
     end subroutine write_out
  end interface

contains

  ! Set the hyperparameters of `gp' to `hypers', assuming they are
  ! ordered [nu(1:nnu), theta(1:ntheta)], where nu are the noise
  ! hyperparameters and theta are the covariance hyperparameters.
  ! Calls `update_matrices'.
  subroutine set_hyperparams(gp, hypers)
    class(BaseGP), intent(inout) :: gp
    real(dp), dimension(:) :: hypers
    integer nnu
    integer ntheta
    nnu = size(gp%nu)
    ntheta = size(gp%theta)

    if ((any(hypers(1:nnu).ne.gp%nu)) &
         & .or.any(hypers(nnu+1:nnu+ntheta).ne.gp%theta)) then
       gp%nu = hypers(1:nnu)
       gp%theta = hypers(nnu+1:nnu+ntheta)
       call gp%update_matrices
    end if
  end subroutine set_hyperparams

  ! A routine with the required interface for the NLopt library, for
  ! maximizing the log-likelihood.
  subroutine nlog_lik(val, n, hypers, grad, need_gradient, gp)
    class(BaseGP) :: gp
    integer :: n, need_gradient
    real(dp) :: val, hypers(n)
    real(dp), intent(inout) :: grad(n)

    call set_hyperparams(gp, hypers)
    val = gp%log_lik()

    call output_params(gp,val)
  end subroutine nlog_lik

  ! Helper routine for nlog_lik: print out the current hyperparameters
  ! (associated with `gp') and their corresponding log-likelihood
  ! (`val').
  subroutine output_params(gp,val)
    class(BaseGP) :: gp
    real(dp) :: val
    integer, save :: u
    logical, save :: assigned = .false.
    if (.not.assigned) open(newunit=u, file="LOG_LIK_OPTIM"); assigned=.true.
    write (u,*) "log likelihood: noise = ", gp%nu, " theta = ", & 
         gp%theta, "log(likelihood) = ", val
  end subroutine output_params

end module m_gp

module m_cov_sqexp
  use m_util
  use m_cov
  
  implicit none
  
  private
  public cov_sqexp
  
  type, extends(cov_fn) :: cov_sqexp
   contains
     procedure, nopass :: ntheta_required
     procedure, nopass :: cov_val
     procedure, nopass :: dcov_x1
     procedure, nopass :: dcov_x2
     procedure, nopass :: d2cov_xx
  end type cov_sqexp

contains
  pure function ntheta_required(dims)
    integer ntheta_required
    integer, intent(in) :: dims
 ! the scale parameter plus one `r' parameter for each dimension
    ntheta_required = dims+1
  end function ntheta_required

  pure function cov_val(x,y,hypers)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, hypers
    real(dp) :: scale, r(size(x,1))
    scale = hypers(1)
    r(:) = hypers(2:)
    cov_val = scale * exp(-0.5*sum((x-y)**2/r**2))
  end function cov_val

  pure function dcov_x1(n,x,y,hypers)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(size(x,1))
    r(:) = hypers(2:)
    dcov_x1 = -(x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers)
  end function dcov_x1

  pure function dcov_x2(n,x,y,hypers)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(size(x,1))
    r(:) = hypers(2:)
    dcov_x2 = (x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,hypers)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    real(dp) :: r(size(x,1))
    r(:) = hypers(2:)
    d2cov_xx = merge(cov_val(x,y,hypers)/r(n)**2, 0.0_dp, n.eq.m) &
         - (x(n)-y(n))/r(n)**2 * dcov_x2(m,x,y,hypers)
  end function d2cov_xx
end module m_cov_sqexp


module m_cov_sqexp_param4
  use m_util
  use m_cov
  
  implicit none
  
  private
  public cov_sqexp_param4
  
  type, extends(cov_fn) :: cov_sqexp_param4
   contains
     procedure, nopass :: ntheta_required
     procedure, nopass :: cov_val
     procedure, nopass :: dcov_x1
     procedure, nopass :: dcov_x2
     procedure, nopass :: d2cov_xx
  end type cov_sqexp_param4

contains
  pure function ntheta_required(dims)
    integer ntheta_required
    integer, intent(in) :: dims
    ntheta_required = 4 ! scale parameter, and r parameters for diag, off-diag and energy
  end function ntheta_required

  pure function cov_val(x,y,hypers)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, hypers
    real(dp) :: r(7), scale
    scale = hypers(1)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    cov_val = scale * exp(-0.5*sum((x-y)**2/r**2))
  end function cov_val

  pure function dcov_x1(n,x,y,hypers)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(7)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    dcov_x1 = -(x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers)
  end function dcov_x1

  pure function dcov_x2(n,x,y,hypers)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(7)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    dcov_x2 = (x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,hypers)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    real(dp) :: r(7)
    r(1:3) = hypers(2)
    r(4:6) = hypers(3)
    r(7) = hypers(4)
    d2cov_xx = merge(cov_val(x,y,hypers)/r(n)**2, 0.0_dp, n.eq.m) &
         - (x(n)-y(n))/r(n)**2 * dcov_x2(m,x,y,hypers)
  end function d2cov_xx
end module m_cov_sqexp_param4


module m_cov_lin
  use m_util
  use m_cov
  
  implicit none
  
  private
  public cov_lin
  
  type, extends(cov_fn) :: cov_lin
   contains
     procedure, nopass :: ntheta_required
     procedure, nopass :: cov_val
     procedure, nopass :: dcov_x1
     procedure, nopass :: dcov_x2
     procedure, nopass :: d2cov_xx
  end type cov_lin

contains
  pure function ntheta_required(dims)
    integer ntheta_required
    integer, intent(in) :: dims
 ! no parameters to be estimated
    ntheta_required = 0
  end function ntheta_required

  pure function cov_val(x,y,hypers)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, hypers
    cov_val = 10 + 10 * (sum(x * y)) ! why don't we use the crossing parameter here? It would need to be MAP estimated
  end function cov_val

  pure function dcov_x1(n,x,y,hypers)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    dcov_x1 = 10 * x(n) * y(n)
  end function dcov_x1

  pure function dcov_x2(n,x,y,hypers)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    dcov_x2 = 10 * x(n) * y(n)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,hypers)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    if (n .eq. m) then
    d2cov_xx = 10
    else
    d2cov_xx = 0
    end if
  end function d2cov_xx
end module m_cov_lin


module m_cov_linsqexp
  use m_util
  use m_cov
  
  implicit none
  
  private
  public cov_linsqexp
  
  type, extends(cov_fn) :: cov_linsqexp
   contains
     procedure, nopass :: ntheta_required
     procedure, nopass :: cov_val
     procedure, nopass :: dcov_x1
     procedure, nopass :: dcov_x2
     procedure, nopass :: d2cov_xx
  end type cov_linsqexp

contains
  pure function ntheta_required(dims)
    integer ntheta_required
    integer, intent(in) :: dims
 ! the scale parameter plus one `r' parameter for each dimension
    ntheta_required = dims+1
  end function ntheta_required

  pure function cov_val(x,y,hypers)
    real(dp) :: cov_val
    real(dp), intent(in), dimension(:) :: x, y, hypers
    real(dp) :: scale, r(size(x,1))
    scale = hypers(1)
    r(:) = hypers(2:)
    cov_val = scale * exp(-0.5*sum((x-y)**2/r**2)) + 10 + 10 * (sum(x * y))
  end function cov_val

  pure function dcov_x1(n,x,y,hypers)
    real(dp) :: dcov_x1
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(size(x,1))
    r(:) = hypers(2:)
    dcov_x1 = -(x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers) + 10 * x(n) * y(n)
  end function dcov_x1

  pure function dcov_x2(n,x,y,hypers)
    real(dp) :: dcov_x2
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n
    real(dp) :: r(size(x,1))
    r(:) = hypers(2:)
    dcov_x2 = (x(n)-y(n))/r(n)**2 * cov_val(x,y,hypers) + 10 * x(n) * y(n)
  end function dcov_x2

  pure function d2cov_xx(n,m,x,y,hypers)
    real(dp) :: d2cov_xx
    real(dp), intent(in), dimension(:) :: x, y, hypers
    integer, intent(in) :: n, m
    real(dp) :: r(size(x,1))
    r(:) = hypers(2:)
    d2cov_xx = merge(cov_val(x,y,hypers)/r(n)**2, 0.0_dp, n.eq.m) &
         - (x(n)-y(n))/r(n)**2 * dcov_x2(m,x,y,hypers)
    if (n .eq. m) then
    d2cov_xx = d2cov_xx + 10
    end if      
  end function d2cov_xx
end module m_cov_linsqexp


module m_cov_all
  use m_cov
  use m_cov_sqexp
  use m_cov_sqexp_param4
  use m_cov_lin
  use m_cov_linsqexp

  use m_util, only: max_name_len

contains

  function cov_fn_to_string(CovFunction) result (cov_fn_name)
    class(cov_fn), intent(in) :: CovFunction
    character(len=max_name_len) cov_fn_name

    select type(cf => CovFunction)
    type is (cov_sqexp) 
       cov_fn_name = 'SQEXP'
    type is (cov_sqexp_param4)
       cov_fn_name = 'SQEXP4PARAM'
    type is (cov_lin)
       cov_fn_name = 'LIN'
    type is (cov_linsqexp)
       cov_fn_name = 'LINSQEXP'
    class default
       cov_fn_name = 'UNKNOWN'
    end select
  end function cov_fn_to_string

  subroutine string_to_cov_fn(cov_fn_name, CovFunction)
    character(len=max_name_len), intent(in) :: cov_fn_name
    class(cov_fn), intent(out), allocatable :: CovFunction

    select case (cov_fn_name)
    case ('SQEXP')
       allocate(cov_sqexp :: CovFunction)
    case ('SQEXP4PARAM')
       allocate(cov_sqexp_param4 :: CovFunction)
    case ('LIN')
       allocate(cov_lin :: CovFunction)
    case ('LINSQEXP')
       allocate(cov_linsqexp :: CovFunction)
    case default
       print *, "unknown covariance function type, ", cov_fn_name
       stop 1
    end select
  end subroutine string_to_cov_fn

end module m_cov_all


module m_noise_value_only
  use m_noise
  use m_util
  implicit none
  
  private
  public noise_value_only

  type, extends(noise_model) :: noise_value_only
   contains
     procedure, nopass :: nparams_required
     procedure, nopass :: noise
  end type noise_value_only

contains
  pure function nparams_required(dims)
    integer, intent(in) :: dims
    integer nparams_required
    ! single noise level to be applied to the target value (and zero
    ! to the derivatives)
    nparams_required = 1
  end function nparams_required

  pure function noise(obs_type, params)
    integer, intent(in) :: obs_type
    real(dp), intent(in) :: params(:)
    real(dp) noise
    noise = merge(params(1), 0.0_dp, obs_type.eq.0)
  end function noise
  
end module m_noise_value_only


module m_noise_param2
  use m_noise
  use m_util
  implicit none

  private
  public noise_param2

  type, extends(noise_model) :: noise_param2
   contains
     procedure, nopass :: nparams_required
     procedure, nopass :: noise
  end type noise_param2

contains

  pure function nparams_required(dims)
    integer, intent(in) :: dims
    integer nparams_required
    ! value and derivative noise parameters
    nparams_required = 2
  end function nparams_required

  pure function noise(obs_type, params)
    integer, intent(in) :: obs_type
    real(dp), intent(in) :: params(:)
    real(dp) noise

    noise = 0

    if (obs_type.eq.0) then
       noise = params(1)
    else if (obs_type.ge.1.and.obs_type.le.6) then
       noise = params(2)
    end if
  end function noise

end module m_noise_param2


module m_noise_all
  use m_noise
  use m_noise_value_only
  use m_noise_param2

  use m_util, only: max_name_len

contains
  
  function noise_model_to_string(NoiseModel) result(noise_model_name)
    class(noise_model), intent(in) :: NoiseModel
    character(len=max_name_len) noise_model_name

    select type(nm => NoiseModel)
    type is (noise_value_only)
       noise_model_name = 'VAL'
    type is (noise_param2)
       noise_model_name = 'PARAM2'
    class default
       noise_model_name = 'UNKNOWN'
    end select
  end function noise_model_to_string

  subroutine string_to_noise_model(noise_model_name, NoiseModel)
    character(len=max_name_len), intent(in) :: noise_model_name
    class(noise_model), intent(out), allocatable :: NoiseModel

    select case (noise_model_name)
    case ('VAL')
       allocate(noise_value_only :: NoiseModel)
    case ('PARAM2')
       allocate(noise_param2 :: NoiseModel)
    case default
       print *, "unknown noise model, ", noise_model_name
       stop 1
    end select
  end subroutine string_to_noise_model

end module m_noise_all


module m_gp_dense
! A Gaussian process of full rank.  See [1, section 2], which gives an
! introduction to Gaussian processes.
! 
! [1] J. Quinonero and C. Rasmussen. Analysis of Some Methods for
! Reduced Rank Gaussian Process Regression, in Switching and Learning
! in Feedback Systems: European Summer School on Multi-Agent Control,
! Maynooth, Ireland, September 8-10, 2003, Revised Lectures and
! Selected Papers, Springer, 2005
  use m_gp
  use m_util
  use m_cov_all
  use m_noise_all
  implicit none

  private
  public DenseGP, read_DenseGP

  type, extends(BaseGP) :: DenseGP
     ! covariance matrix and inverse (C is represented by Q in the
     ! notation of ref [1])
     real(dp), dimension(:,:), allocatable :: C, invC
     ! precomuted product, used in the prediction
     real(dp), dimension(:), allocatable :: invCt
   contains
     procedure log_lik
     procedure update_matrices
     procedure predict
     procedure write_out
  end type DenseGP

  interface DenseGP
     module procedure make_DenseGP
     module procedure read_DenseGP
  end interface DenseGP

contains

  subroutine alloc_DenseGP(gp, n, ntheta, dims, CovFunction, NoiseModel)
    type(DenseGP), intent(inout) :: gp
    ! n: number of observations
    ! ntheta: number of covariance hyperparameters
    ! dims: dimension of the inputs
    integer, intent(in) :: n, ntheta, dims
    class(cov_fn) :: CovFunction
    class(noise_model) :: NoiseModel
    
    allocate(real(dp) :: gp%nu(NoiseModel%nparams_required(dims)))
    allocate(real(dp) :: gp%theta(ntheta))
    allocate(real(dp) :: gp%x(n, dims))
    allocate(integer :: gp%obs_type(n))
    allocate(real(dp) :: gp%t(n))
    allocate(real(dp) :: gp%C(n,n))
    allocate(real(dp) :: gp%invC(n,n))
    allocate(real(dp) :: gp%invCt(n))
    allocate(gp%covariance, mold=CovFunction)
    allocate(gp%noise_model, mold=NoiseModel)
  end subroutine alloc_DenseGP

  function make_DenseGP(nu, theta, x, obs_type, t, CovFunction, NoiseModel) result(gp)
    type(DenseGP) :: gp
    ! noise hyperparameters
    real(dp), dimension(:), intent(in) :: nu
    ! covariance hyperparameters
    real(dp), dimension(:), intent(in) :: theta
    ! training input coordinates
    real(dp), dimension(:,:), intent(in) :: x
    ! training input observation types
    integer,  dimension(:), intent(in) :: obs_type
    ! training outputs
    real(dp), dimension(:), intent(in) :: t

    class(cov_fn) :: CovFunction
    class(noise_model) :: NoiseModel
    
    integer :: n, d
    n = size(t)
    d = size(x,2)

    if (size(theta) /= CovFunction%ntheta_required(d)) then
       print *, "size of theta does not match number of hyperparameters required by the covariance function"
       stop 1
    end if
    
    if (size(nu) /= NoiseModel%nparams_required(d)) then
       print *, "size of nu (noise params) does not match number required by the noise model"
       stop 1
    end if

    call alloc_DenseGP(gp, n, size(theta), d, CovFunction, NoiseModel)

    gp%nu = nu
    gp%theta = theta
    gp%obs_type = obs_type
    gp%t = t
    gp%x = x
    call update_matrices(gp)    
  end function make_DenseGP

  subroutine write_out(this, filename)
    class(DenseGP), intent(in) :: this
    character(len=*), intent(in) :: filename
    character(len=max_name_len) :: cov_fn_name
    character(len=max_name_len) :: noise_model_name
    integer :: u ! unit number for output

    cov_fn_name = cov_fn_to_string(this%covariance)
    noise_model_name = noise_model_to_string(this%noise_model)

    open(newunit=u, file=filename)

    write (u,'(A)') "DenseGP"
    write (u,'(I10)') size(this%t), size(this%theta), size(this%nu), size(this%x,2)
    write (u,'(A)') trim(cov_fn_name)
    write (u,'(A)') trim(noise_model_name)
    write (u,'(es24.15)') this%nu, this%theta, this%x
    write (u,'(I4)') this%obs_type
    write (u,'(es24.15)') this%C, this%invC, this%invCt

    close(u)
  end subroutine write_out

  function read_DenseGP(filename) result(gp)
    character(len=*), intent(in) :: filename
    type(DenseGP) :: gp
    integer n, ntheta, nnu, d, u
    character(len=max_name_len) :: label
    character(len=max_name_len) :: cov_fn_name
    character(len=max_name_len) :: noise_model_name
    class(cov_fn), allocatable :: CovFunction
    class(noise_model), allocatable :: NoiseModel
    open(newunit=u, file=filename)
    read (u,'(A)') label

    if (trim(label) /= "DenseGP") then
       print *, "Incompatible data file"
       stop 1
    end if

    read (u,'(I10)') n, ntheta, nnu, d

    read (u,'(A)') cov_fn_name 
    call string_to_cov_fn(cov_fn_name, CovFunction)
    if (ntheta /= CovFunction%ntheta_required(d)) then
       print *, "ntheta does not match number required by the covariance function"
       stop 1
    end if

    read (u,'(A)') noise_model_name
    call string_to_noise_model(noise_model_name, NoiseModel)
    if (nnu /= NoiseModel%nparams_required(d)) then
       print *, "size of nu (noise params) does not match number required by the noise model"
       stop 1
    end if

    call alloc_DenseGP(gp, n, ntheta, d, CovFunction, NoiseModel)
    
    read (u,'(es24.15)') gp%nu, gp%theta, & 
         gp%x
    read (u,'(I4)') gp%obs_type
    read (u,'(es24.15)') gp%C, gp%invC, gp%invCt
    close(u)
  end function read_DenseGP

  subroutine update_matrices(this)
    class(DenseGP), intent(inout) :: this
    integer :: i,j,n
    real(dp) noise
    ! always a small amount of noise to stabilize the inversion
    real(dp), parameter :: noise_stab = 1e-9_dp

    n = size(this%t)
    
    do i=1,n
       do j=1,n
          if (i.ne.j) then
             noise = 0.0_dp
          else
             noise = this%noise_model%noise(this%obs_type(i), this%nu) + noise_stab
          end if
          ! Q in ref [1] (just below eq. (7) therein)
          this%C(i,j) = this%covariance%cov(this%obs_type(i), this%obs_type(j), &
               this%x(i,:), this%x(j,:), this%theta) + noise
       end do
    end do
    
    this%invC = this%C
    call ninv(this%invC)
    this%invCt = solve(this%C, this%t)
  end subroutine update_matrices

  function predict(this, xnew, obs_type_new)
    real(dp) predict
    class(DenseGP), intent(in) :: this
    real(dp), dimension(:), intent(in) :: xnew
    integer, optional, intent(in) :: obs_type_new
    
    integer :: obs_type_new1, i
    
    real(dp), dimension(size(this%t)) :: k
    if (present(obs_type_new)) then
       obs_type_new1 = obs_type_new
    else 
       obs_type_new1 = 0
    endif

    do i=1,size(this%t)
       k(i) = this%covariance%cov(obs_type_new1,this%obs_type(i),xnew,this%x(i,:),this%theta)
    end do
    ! Predictive mean from eq. (6) in ref [1]. Here, k corresponds to
    ! k^{*} in the reference, and t to y
    predict = dot_product(k, this%invCt)
  end function predict

  function log_lik(this)
    class(DenseGP), intent(in) :: this
    real(dp) log_lik
    ! First line of eq. (8) in ref. [1]
    log_lik = -0.5_dp * (logdet(this%C) + dot_product(this%t, this%invCt))
  end function log_lik

end module m_gp_dense


MODULE mo_emulator
  use physconst,                only: gravit
  
  ! Modules from this file, needed for GPE
  USE m_gp
  USE m_util
  USE m_gp_dense
  
  ! Modules included in OpenIFS
  !USE mo_kind,                 ONLY : dp
  !USE mo_physical_constants,   ONLY : grav
  !USE mo_netcdf,               ONLY: nf_max_name, nf_open, NF_NOERR, nf_nowrite, nf_inq_dimid, nf_inq_dimlen, &
  !                            nf_get_var_double, nf_get_vara_double, nf_close, nf_inq_varid
  !USE mo_exception,            ONLY: finish
  !USE YOMMP0,   ONLY : MYPROC  ! Number of processor for printing output 
  !USE YOMCT3, ONLY : NSTEP     ! To check step number for reading emulator only at the first step
  IMPLICIT NONE

  PRIVATE
  !Define emulator, global, one emulator for vertical wind sigma
  CLASS(BaseGP), allocatable :: gp_emu_w
  ! Define arrays needed by the standardisation of emulator input and output
  REAL(dp), DIMENSION(:,:), allocatable :: ex1_w  !Training input, used for standardation of emulator input
  REAL(dp), DIMENSION(:), allocatable   :: t1_w   !Training output

  PUBLIC :: emulator,initEmulator
CONTAINS



  SUBROUTINE  readTraining(train_filename,train_n,train_dim,train_in,train_out)

    CHARACTER(len=*), INTENT(IN)                        :: train_filename
    INTEGER                                             :: uu, train_n, i, train_dim
    INTEGER, DIMENSION(train_n)                         :: a
    REAL(dp), DIMENSION(train_n,train_dim), INTENT(OUT) :: train_in
    REAL(dp), DIMENSION(train_n), INTENT(OUT)           :: train_out
    
    ! Read in emulator training data. These are used only in the
    ! standardisation and only their means are used.
    open(newunit=uu, file=train_filename)
    read (uu,*) (train_in(i,1:train_dim), a(i), train_out(i), i=1,train_n)
    close(uu)

    
  END SUBROUTINE readTraining

  SUBROUTINE initEmulator()
    !This routine initialize emulators, reads emulator files
    REAL(dp), DIMENSION(1418,8)             :: train_in
    REAL(dp), DIMENSION(1418)               :: train_out
    logical, save :: weights_loaded = .false.
    if (weights_loaded) return
    !Emulator
    IF (allocated(gp_emu_w)) DEALLOCATE(gp_emu_w)
    allocate(gp_emu_w, source = DenseGP('gp_emu_w'))  
    
    !Training data
    CALL readTraining('DATA',1418,8,train_in,train_out)
    IF (allocated(ex1_w)) DEALLOCATE(ex1_w)
    allocate(ex1_w(1418,8), SOURCE=train_in)
    IF (allocated(t1_w)) DEALLOCATE(t1_w)
    allocate(t1_w(1418), SOURCE=train_out)
    weights_loaded = .true.
 
  END SUBROUTINE initEmulator


!>
!! This subroutine calculates vertical interacls from ground to certain level
!! used to calculate lwp inside boundary layer
  SUBROUTINE calculate_lwp(x,klev,cloud_top,zdpg,cfrac,integral) 
       INTEGER, INTENT(IN)    :: cloud_top, klev
       REAL(dp),INTENT(IN)    :: x(klev),zdpg(klev),cfrac(klev)
       REAL(dp),INTENT(OUT)   :: integral
       REAL(dp)               :: zintegral
       INTEGER                :: jk
       zintegral = 0
       DO 100 jk = cloud_top,klev
          IF (cfrac(jk) > 0.1) THEN                                    ! Divide by cloud fraction if there is > 10% cloud in that level
              zintegral  = zintegral  + x(jk) *zdpg(jk)/cfrac(jk)
          ELSE
              zintegral  = zintegral  + x(jk) *zdpg(jk)
          END IF
       100 END DO
       integral=zintegral
  END SUBROUTINE calculate_lwp



  SUBROUTINE find_cloud2(lev,indhi,x,cloud_top,cloud_base)
! Input: lev   = number of layers
!        indhi = index of the highest layer considered 
!                (e.g., the 700 hPa level)
!        x(lev) = cloud liquid water content [g/kg]
! Output: cloud_top  = full-level index for the lowermost layer of the lowest cloud 
!         cloud_base = full-level index for the uppermost layer of the lowest cloud 
   
!    INTEGER, PARAMETER :: dp = selected_real_kind(13, 300)
! This would be better:
    !USE mo_kind, ONLY : dp
    USE m_util
    INTEGER, INTENT(IN)   :: lev, indhi
    REAL(dp), INTENT(in)  :: x(lev)
    INTEGER, INTENT(out)  :: cloud_top, cloud_base

    INTEGER :: ilev

! Initialize to undefined
    cloud_top=-999
    cloud_base=-999
    DO ilev=lev,indhi,-1
      IF (x(ilev) >= 0.01) THEN
        cloud_base=ilev
        EXIT
      END IF
    ENDDO 
    IF (cloud_base >= indhi) THEN
      DO ilev=cloud_base,indhi,-1
        IF (x(ilev) < 0.01) EXIT
      ENDDO 
      cloud_top=MAX(indhi,ilev+1)
    END IF

    RETURN
  END SUBROUTINE find_cloud2

  SUBROUTINE find_cloud(x,lev,cloud_top,cloud_base)
        integer, parameter :: dp = selected_real_kind(13, 300)
    INTEGER, INTENT(IN)   :: lev
        REAL(dp),INTENT(IN)   :: x(lev)
        REAL(dp)     :: x1(lev)
    INTEGER, INTENT(OUT)  :: cloud_base,cloud_top
        INTEGER               :: i,ind
        cloud_top=0
        cloud_base=0
    !write(*,*) x,'find cloud profile'
    x1=x
    x1 = x1(size(x1):1:-1) !reverse order
    WHERE(x1 >= 0.01) x1 = 999
    WHERE(x1 < 0.01) x1=-999
    cloud_base =  SIZE(x1,1)-MAXLOC(x1,1)
        cloud_top =  SIZE(x1,1)-(MINLOC(x1(MAXLOC(x1,1):SIZE(x1,1)),1)-1+MAXLOC(x1,1)-1)-1 !change to calculate right cloud top
  END SUBROUTINE find_cloud

!>
!! This subroutine computes mask where emulator is applied for each timestep.
!! If region is inside training space emulator is also applied.
!! Returns mask and emulated value of vertical wind standard deviation.
  SUBROUTINE emulator(kidia,kfdia, kbdim, ktdia, klev, klevp1, &
                      pfull, phalf, &
              pxlm1,pxim1,pqm1,cfrac,shf, cosmu, tpot, pwsigma, pemu_mask, pemu_cb)
  ! aps = surface pressure
  ! tpot = potential temperature
  ! shf = sensible heat flux
  ! cfrac = cloud fraction

  ! Input for the emulator 

    INTEGER, INTENT(IN)    :: kidia,kfdia, kbdim, ktdia, klev, klevp1
    ! kidia = starting horizontal point
    ! kfdia = ending horizontal point
    ! ktdia = starting vertical point (highest level)
    ! klevp1 = number of half levels (klevp1=klev+1)
    ! klev = number of vertical points, highest number is the surface
    ! kbdim = number of horizontal points
    
 
    REAL(dp), INTENT(in) :: &
      pfull(kbdim,klev),    & ! Full-level pressure [Pa]  
      phalf(kbdim,klevp1),  & ! Half-level pressure [Pa]  
      pxlm1(kbdim,klev),    & ! cloud liquid water content [kg/kg]
      pxim1(kbdim,klev),    & ! cloud ice content [kg/kg]
      pqm1(kbdim,klev),     & ! specific humidity [kg/kg]
      cfrac(kbdim,klev),    & ! cloud fraction in levels [0,1] 
      shf(kbdim),           & ! surface sensible heat flux [W/m2]
      cosmu(kbdim),         & ! solar zenith angle [angle or cos of angle?]
      tpot(kbdim,klev)        ! potential temperature [K]

    REAL(dp), INTENT(INOUT) :: pwsigma(kbdim,klev) ! Standard deviation of vertical wind [m/s]
                                                   ! Contains constant values as input, stratocumulus points will be updated by emulator
    REAL(dp), INTENT(INOUT)   :: pemu_mask(kbdim)    ! emulator mask 
    REAL(dp), INTENT(INOUT)   :: pemu_cb(kbdim,klev)	    ! Mask for single layer stratocumulus cloud base (3D variable)
    !Local variables
    INTEGER                :: zfull700(kbdim),zhalf700(kbdim) ! index for nearest level at 700hpa level. 
    REAL(dp)               :: hm(kbdim,klevp1), & ! Half level altitudes [m]
                              dz(kbdim,klev)  ! Thicknesses of levels [m]

    INTEGER                :: jk,        &
                              cloud_top, &                               ! index for top of the lowest cloud layer
                              cloud_base,inv_ind_hi,inv_ind_lo           ! index for base of lowest cloud layer,index's for inversion calculations
    REAL(dp)               :: zlwp700(kbdim),ziwp700(kbdim), &           ! LWP and IWP inside 700hpa layer
                              zlwp(kbdim),ziwp(kbdim), &                 ! Total lwp and iwp
                              pdp(kbdim,klev), pdpg(kbdim,klev), &       ! pressure and altitude difference
                              g_rcp,zx(klev),zx1(klevp1),ztotalq(klev)   ! zx is just temporary vector,ztotalq is total water

    REAL(dp) :: meanlt, stdlt
    REAL(dp), DIMENSION(:) :: emuInputVec(8),emuInputVecStandard(8)
    
    ! Variables used as the emulator input, the units are ones that are needed for the emulator
    REAL(dp)            :: pcloud_top(kbdim)        ! cloud top [m]
    REAL(dp)            :: pcloud_base(kbdim)       ! cloud base [m]
    REAL(dp)            :: phpalim700(kbdim)        ! 700hpa limit debug[m]
    REAL(dp)            :: ptpot_inv(kbdim)         ! potential temperature inversion strengt [K]
    REAL(dp)            :: ptpot_pbl(kbdim)         ! potential temperature at pbl [K]
    REAL(dp)            :: ph2o_inv(kbdim)          ! h2o inversion [kg/kg]
    REAL(dp)            :: ppbl_num(kbdim)          ! particle number at pbl [m3]
    REAL(dp)            :: ppbl_h(kbdim)            ! pbl height [m] 
    REAL(dp)            :: ppres0(kbdim)            ! surface pressure [Pa]
    REAL(dp)            :: pshf(kbdim)              ! sensible heat flux [W/m2]
    REAL(dp)            :: pemu_lwp(kbdim)          ! lwp input for emulator [kg/m2]
    ! Emulator output
    REAL(dp)            :: pemu_w(kbdim)            ! emulated updraft velocity [m/s]
    REAL(dp) 		:: aps(kbdim) 		    ! surface pressure [Pa]
    
    ! CHANGE THIS TO READ THE EMULATOR ONLY ON THE FIRST STEP TO SAVE TIME. SUBROUTINE readTraining needs to be changed for this
    !IF (NSTEP <= 1) Then
        CALL initEmulator()
    !END IF
 
    
    
    g_rcp = 1._dp / gravit

    
    ! mean and std for standardizing emulator input and output
    meanlt = mean(t1_w,1418)
    stdlt  = std(t1_w,meanlt,1418)

    DO jk = kidia,kfdia
        ! For LWP calculations
        aps(jk) = phalf(jk,klevp1)
        pdp(jk,ktdia:klev) = phalf(jk,ktdia+1:klevp1) - phalf(jk,ktdia:klev)
        pdpg(jk,ktdia:klev) = g_rcp * pdp(jk,ktdia:klev)
        hm(jk,ktdia:klevp1) = 44307.69396_dp*(1.0_dp-(phalf(jk,ktdia:klevp1)/101325.0_dp)**0.190284_dp)! half levels in m
        dz(jk,ktdia:klev) = hm(jk,ktdia:klevp1-1) - hm(jk,ktdia+1:klevp1)! Level thicknesses in m
    END DO

    ! Emulator is only trained for low level clouds (clouds below 700hpa level). 
    ! 1. step find index for level 700hpa both for half levels and full levels
    zfull700 = -1
    zhalf700 = -1
    DO 101 jk = kidia,kfdia
       zx =pfull(jk,:)
            WHERE(zx .GE. 70000) zx = 999
       zx1 =phalf(jk,:)
            WHERE(zx1 .GE. 70000) zx1 = 999
           zfull700(jk) = MINLOC(MERGE(0,1,zx == 999),DIM=1)
           zhalf700(jk) = MINLOC(MERGE(0,1,zx1 == 999),DIM=1)
    101 END DO

    ! Loop over grid points
    DO 102 jk =kidia,kfdia
        pemu_mask(jk) = 0.0_dp
        pemu_lwp(jk) = -999
        pemu_w(jk) = -999
        pemu_cb(jk,:) = 0.0_dp

        ! Eliminate points where there is no low level clouds. There is cloud if there is more than 0.01 g/kg of cloud water below 700hPa
        IF(COUNT(pxlm1(jk,zfull700(jk):klev)*1000 .GE. 0.01) < 1.0_dp) THEN
            ptpot_inv(jk) = -999
            ph2o_inv(jk) = -999
            !ph2o_pbl(jk) = -999
            ptpot_pbl(jk) = -999
            ppbl_h(jk) = -999
            ppres0(jk) = -999
            pshf(jk) = -999
            pcloud_top(jk) = -999
            pcloud_base(jk) = -999
            phpalim700(jk) = -999
            pemu_w(jk) = -999
            CYCLE
        END IF
        ! 3. Eliminate all points where there is fog, check if there is cloud water on lowest layer and excluce these points
        IF(pxlm1(jk,klev)*1000 .GE. 0.01) THEN        ! First level water content > 0.01 g/kg
            ptpot_inv(jk) = -999
            ph2o_inv(jk) = -999
            !ph2o_pbl(jk) = -999
            ptpot_pbl(jk) =- 999
            ppbl_h(jk) = -999
            ppres0(jk) = -999
            pshf(jk) = -999
            pcloud_top(jk) = -999
            pcloud_base(jk) = -999
            phpalim700(jk) = -999
            pemu_w(jk) = -999            
            CYCLE
        END IF
        ! 4. step, locate lowest cloud
        CALL find_cloud2(klev,zfull700(jk),pxlm1(jk,:)*1000,cloud_top,cloud_base)
        
        ! 5. step calculate LWP inside lowest layer layer, this is used to identify that most of the cloud water in the column is in low level cloud
        ! In-cloud LWP, using cloud fraction of each level (OpenIFS has cloud fraction for each level)
        ! Output zlwp and ziwp are in g/m2
        CALL calculate_lwp(pxlm1(jk,:)*1000,klev,cloud_top,pdpg(jk,:),cfrac(jk,:),zlwp700(jk))
        CALL calculate_lwp(pxlm1(jk,:)*1000,klev,1,pdpg(jk,:),cfrac(jk,:),zlwp(jk)) 
        CALL calculate_lwp(pxim1(jk,:)*1000,klev,cloud_top,pdpg(jk,:),cfrac(jk,:),ziwp700(jk))
        CALL calculate_lwp(pxim1(jk,:)*1000,klev,1,pdpg(jk,:),cfrac(jk,:),ziwp(jk))
        
        ! Acceptable cloud conditions
        IF((zlwp700(jk) > ((zlwp(jk)+ziwp(jk))*0.5)) .AND. ( 0.1*zlwp700(jk)  > ziwp700(jk))) THEN
            ! Add cloud base to mask
            pemu_cb(jk,cloud_base)=1.0_dp
            ! Calculate emulator inputs
            inv_ind_hi = cloud_top-2
            ztotalq = pxlm1(jk,:)+pqm1(jk,:) ! [kg/kg]
            inv_ind_lo = MIN(cloud_base+2,klev)
            ptpot_inv(jk) = MAXVAL(tpot(jk,inv_ind_hi:inv_ind_lo))-MINVAL(tpot(jk,inv_ind_hi:inv_ind_lo))
            !ptpot_inv(jk) = MAXVAL(tpot(jk,inv_ind_hi:klev))-MINVAL(tpot(jk,inv_ind_hi:klev))
            ph2o_inv(jk) = MAXVAL(ztotalq(inv_ind_hi:inv_ind_lo))-MINVAL(ztotalq(inv_ind_hi:inv_ind_lo))
            !ph2o_pbl(jk) = MAXVAL(ztotalq(inv_ind_hi:inv_ind_lo))
            ptpot_pbl(jk) = MINVAL(tpot(jk,inv_ind_hi:inv_ind_lo))
            ppbl_h(jk) = 44307.69396*((aps(jk)/101325)**0.190284-(phalf(jk,cloud_top)/101325)**0.190284) ! in m, positive value
            ppres0(jk) = aps(jk)                             ! Surface level pressure [Pa]
            pshf(jk) = shf(jk)                               ! Sensible heat flux [W/m2]
            pcloud_top(jk) = cloud_top
            pcloud_base(jk) = cloud_base
            phpalim700(jk) = REAL(zfull700(jk))
            pemu_lwp(jk) = zlwp700(jk)*1e-3                  !Convert from g/m2 to kg/m2
  
            ! Set emulator input
            emuInputVec(1)=-pshf(jk)                                    ! Sensible heat flux W/m2, positive downwards in OIFS, positive upwards in UCLALES and emulator
            emuInputVec(2)=cosmu(jk)                                    ! cos of solar zenith angle
            emuInputVec(3)=ppres0(jk)                                   ! Surface pressure Pa
            emuInputVec(4)=ph2o_inv(jk)                  		! Delta qt kg/kg
            emuInputVec(5)=ptpot_pbl(jk)                                ! Theta PBL (cb) K
            emuInputVec(6)=MAX(ptpot_inv(jk),2.43457)                   ! Delta theta K, minimum value 2.43K is the lowest value in the training data
            emuInputVec(7)=pemu_lwp(jk)                                 ! In-cloud LWP kg/m2
            emuInputVec(8)=ppbl_h(jk)                                   ! Height of PBL (cloud top) m
            ! Emulate updraft for condition variables that are within the training data of the emulator

            ! Standardization returns -999 if all variables in emuInputVec are not within the training data ranges
            CALL standardize_emulator_inputs(emuInputVec(:),emuInputVecStandard(:),1418,8,ex1_w)
            ! Edit mask for the cases that don't fit the emulator cloud conditions
            IF(ANY(emuInputVecStandard(:)==-999)) THEN
                pemu_mask(jk) = -1.0_dp
                
            ELSE
                pemu_w(jk) = unstandardize_s(gp_emu_w%predict(emuInputVecStandard(:)),meanlt,stdlt)
                ! Emulator was used, update the corresponding value in pwsigma
                IF (pemu_w(jk) < 0.1) THEN
                    pwsigma(jk,cloud_top:klev) = 0.1_dp ! set a lower limit to pwsigma
                ELSE
                    pwsigma(jk,cloud_top:klev) = pemu_w(jk)
                END IF
                pemu_mask(jk) = 1.0_dp
                !WRITE(8000+MYPROC,*) 'input = ', emuInputVec, pemu_w(jk), pemu_mask(jk)
            END IF
        END IF
    102 END DO 

   END SUBROUTINE emulator

   SUBROUTINE standardize_emulator_inputs(emu_input_vec,emu_input_vec_out,train_n,train_dim,train_out)
     ! Standardise emulator inputs. If any input value is outside the training data, return -999
     INTEGER, INTENT(IN)                                  :: train_n, train_dim
     REAL(dp), DIMENSION(train_dim), INTENT(IN)           :: emu_input_vec
     REAL(dp), DIMENSION(train_dim), INTENT(out)          :: emu_input_vec_out
     REAL(dp)                                             :: stdx, meanx,min_,max_
     REAL(dp), DIMENSION(train_n,train_dim), INTENT(in)   :: train_out
     INTEGER                                              :: jj

     DO jj = 1,train_dim
        meanx = mean(train_out(:,jj),train_n)
        stdx  = std(train_out(:,jj),meanx,train_n)
    max_ = MAXVAL(train_out(:,jj))
    min_ = MINVAL(train_out(:,jj))
    IF ((emu_input_vec(jj) >= min_) .AND. (emu_input_vec(jj) <= max_)) THEN
            emu_input_vec_out(jj) = standardize(emu_input_vec(jj),meanx,stdx)
    ELSE
            emu_input_vec_out(jj) = -999
        
    END IF
     END DO
     

  END SUBROUTINE  standardize_emulator_inputs
  
  SUBROUTINE standardize_emulator_inputs2(emu_input_vec,emu_input_vec_out,train_n,train_dim,train_out)
     ! Standardise emulator inputs. Don't change the valu to -999 if some are outside the training data
     INTEGER, INTENT(IN)                                  :: train_n, train_dim
     REAL(dp), DIMENSION(train_dim), INTENT(IN)           :: emu_input_vec
     REAL(dp), DIMENSION(train_dim), INTENT(out)          :: emu_input_vec_out
     REAL(dp)                                             :: stdx, meanx,min_,max_
     REAL(dp), DIMENSION(train_n,train_dim), INTENT(in)   :: train_out
     INTEGER                                              :: jj

     DO jj = 1,train_dim
        meanx = mean(train_out(:,jj),train_n)
        stdx  = std(train_out(:,jj),meanx,train_n)
        emu_input_vec_out(jj) = standardize(emu_input_vec(jj),meanx,stdx)
     END DO
     

  END SUBROUTINE  standardize_emulator_inputs2

   ! Some helper functions 

  FUNCTION standardize(x,meanx,stdx) RESULT(res)
    REAL(dp) :: x
    REAL(dp) :: res
    REAL(dp) :: meanx 
    REAL(dp) :: stdx

    res = (x - meanx)/(stdx)
  END FUNCTION standardize

   FUNCTION unstandardize_s(x,meanx,stdx) RESULT(res)
     REAL(dp) x
     REAL(dp) :: res
     REAL(dp) :: meanx 
     REAL(dp) :: stdx

     res = (x * stdx) + meanx
   END FUNCTION unstandardize_s
   
   FUNCTION mean(x,dmn) RESULT(res)
     INTEGER dmn
     REAL(dp) x(dmn)
     REAL(dp) :: res

     res = SUM(x)/dmn
   END FUNCTION mean

   FUNCTION std(x,meanx,dmn) RESULT(res)
     INTEGER :: dmn
     REAL(dp) :: x(dmn)
     REAL(dp) :: meanx
     REAL(dp) :: res

     res = SQRT(SUM((x - meanx)**2)/dmn)
   END FUNCTION std

END MODULE mo_emulator

   
    
    

