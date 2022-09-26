MODULE NUM_REC_UTILS
  !****************************************************************************
  !** module partially taken from Numerical Recipes in Fortran 90 (nrutil)   **
  !** contains some subroutines not related to numerical calculation but     **
  !** used by subroutines of MODULE NUM_REC.                                 **
  !****************************************************************************
  IMPLICIT NONE

  INCLUDE "precision.f90"


  ! Parameters or crossover rom serial to parallel algorithms (these are
  ! used only within this NUM_REC_UTILS module):

  INTEGER(LONG),PARAMETER ::NPAR_ARTH=16,NPAR2_ARTH=8
  !Each NPAR2 must be <= the corresponding NPAR.
  INTEGER(LONG),PARAMETER ::NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER(LONG),PARAMETER ::NPAR_CUMSUM=16
  INTEGER(LONG),PARAMETER ::NPAR_CUMPROD=8
  INTEGER(LONG),PARAMETER ::NPAR_POLY=8
  INTEGER(LONG),PARAMETER ::NPAR_POLYTERM=8

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  FUNCTION assert_eq2(n1,n2,string)
    ! Report and die if integers not all equal (used or size checking).
    CHARACTER(LEN=*),INTENT(IN)::string
    INTEGER,INTENT(IN)::n1,n2
    INTEGER ::assert_eq2

    IF (n1 == n2) THEN
       assert_eq2=n1
    ELSE
       WRITE(*,*)'nrerror:an assert_eq failed with this tag:',&
            string
       STOP 'program terminated by assert_eq2'
    END IF
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*),INTENT(IN)::string
    INTEGER,INTENT(IN)::n1,n2,n3
    INTEGER ::assert_eq3

    IF (n1 == n2 .AND. n2 == n3) THEN
       assert_eq3=n1
    ELSE
       WRITE(*,*)'nrerror:an assert_eq failed with this tag:',&
            string
       STOP 'program terminated by assert_eq3'
    END IF
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*),INTENT(IN)::string
    INTEGER,INTENT(IN)::n1,n2,n3,n4
    INTEGER ::assert_eq4

    IF (n1 == n2 .AND. n2 == n3 .AND. n3 == n4) THEN
       assert_eq4=n1
    ELSE
       WRITE(*,*) 'nrerror:an assert_eq failed with this tag:',&
            string
       STOP 'program terminated by assert_eq4'
    END IF
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*),INTENT(IN)::string
    INTEGER,DIMENSION(:),INTENT(IN)::nn
    INTEGER ::assert_eqn

    IF (ALL(nn(2:) == nn(1))) THEN
       assert_eqn=nn(1)
    ELSE
       WRITE (*,*) 'nrerror:an assert_eq failed with this tag:',&
            string
       STOP 'program terminated by assert_eqn'
    END IF
  END FUNCTION assert_eqn

  SUBROUTINE nrerror(string)
    ! Report a message,then die.
    CHARACTER(LEN=*),INTENT(IN)::string

    WRITE (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

END MODULE NUM_REC_UTILS




MODULE NUM_REC
  !****************************************************************************
  !** module partially taken from Numerical Recipes in Fortran 90 (nr)       **
  !** contains some subroutines related to numerical interpolation.          **
  !****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  PRIVATE
  PUBLIC :: splin2, splie2
  PUBLIC :: splint, spline


CONTAINS

  SUBROUTINE tridag_ser(a,b,c,r,u)
    USE NUM_REC_UTILS,ONLY :assert_eq,nrerror
    IMPLICIT NONE

    REAL(DP),DIMENSION(:),INTENT(IN)::a,b,c,r
    REAL(DP),DIMENSION(:),INTENT(OUT)::u
    ! Solves for a vector u of size N the tridiagonal linear set given by
    ! equation (2.4.1) using a serial algorithm. Input vectors b (diagonal
    ! elements)and r (right-hand sides) have size N , while a and c (off-diagonal
    ! elements) are size N-1.
    REAL(DP),DIMENSION(SIZE(b))::gam ! One vector of workspace,gam is needed.
    INTEGER(LONG)::n,j
    REAL(DP)::bet

    n=assert_eq((/SIZE(a)+1,SIZE(b),SIZE(c)+1,SIZE(r),SIZE(u)/),'tridag_ser')
    bet=b(1)
    IF (bet == 0.0) CALL nrerror('tridag_ser: Error at code stage 1')
    ! If this happens then you should rewrite your equations as a set of order
    ! N-1 ,with u2 trivially eliminated.
    u(1)=r(1)/bet
    DO j=2,n ! Decomposition and forward substitution.
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j-1)*gam(j)
       IF (bet == 0.0) & !Algorithm fails ; see elow routine in Vol.1.
            CALL nrerror('tridag_ser: Error at code stage 2')
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    END DO
    DO j=n-1,1,-1 !Backsubstitution.
       u(j)=u(j)-gam(j+1)*u(j+1)
    END DO
  END SUBROUTINE tridag_ser


  RECURSIVE SUBROUTINE tridag(a,b,c,r,u)
    USE NUM_REC_UTILS,ONLY :assert_eq,nrerror

    IMPLICIT NONE
    REAL(DP),DIMENSION(:),INTENT(IN)::a,b,c,r
    REAL(DP),DIMENSION(:),INTENT(OUT)::u
    ! Solves for a vector u of size N the tridiagonal linear set given by
    ! equation (2.4.1) using a parallel algorithm. Input vectors b (diagonal
    ! elements)and r (right-hand sides) have size N , while a and c (off-diagonal
    ! elements) are size N-1.
    INTEGER(LONG),PARAMETER ::NPAR_TRIDAG=4 !Determines when serial algorithm is invoked.
    INTEGER(LONG)::n,n2,nm,nx
    REAL(DP),DIMENSION(SIZE(b)/2)::y,q,piva
    REAL(DP),DIMENSION(SIZE(b)/2-1)::x,z
    REAL(DP),DIMENSION(SIZE(a)/2)::pivc

    n=assert_eq((/SIZE(a)+1,SIZE(b),SIZE(c)+1,SIZE(r),SIZE(u)/),'tridag_par')
    IF (n < NPAR_TRIDAG) THEN
       CALL tridag_ser(a,b,c,r,u)
    ELSE
       IF (MAXVAL(ABS(b(1:n)))==0.0)& ! Algorithm fails;see elow routine in Vol.1.
            CALL nrerror('tridag_par: possible singular matrix')
       n2=SIZE(y)
       nm=SIZE(pivc)
       nx=SIZE(x)
       piva =a(1:n-1:2)/b(1:n-1:2) ! Zero the odd a's and even c's,giving x,y,z,q
       pivc =c(2:n-1:2)/b(3:n:2)
       y(1:nm)=b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
       q(1:nm)=r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
       IF (nm <n2) THEN
          y(n2)=b(n)-piva(n2)*c(n-1)
          q(n2)=r(n)-piva(n2)*r(n-1)
       END IF
       x =-piva(2:n2)*a(2:n-2:2)
       z =-pivc(1:nx)*c(3:n-1:2)
       CALL tridag(x,y,z,q,u(2:n:2)) ! Recurse and get even u's.
       u(1)=(r(1)-c(1)*u(2))/b(1)        ! Substitute and get odd u's.
       u(3:n-1:2)=(r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2)&
            -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
       IF (nm ==n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    END IF
  END SUBROUTINE tridag


  FUNCTION locate(xx,x)
    IMPLICIT NONE
    REAL(DP),DIMENSION(:),INTENT(IN)::xx
    REAL(DP),INTENT(IN)::x
    INTEGER(LONG)::locate
    ! Given an array xx(1:N),and given a value x ,returns a value j such that x
    ! is between xx(j)and xx(j+1).xx must be monotonic, either increasing or
    ! decreasing. j =0 or j =N is returned to indicate that x is out of range.
    INTEGER(LONG)::n,jl,jm,ju
    LOGICAL ::ascnd

    n=SIZE(xx)
    ascnd =(xx(n)>=xx(1)) !True if ascending order of table,false otherwise.
    jl=0 ! Initialize ower
    ju=n+1 ! and upper limits.
    DO
       IF (ju-jl <=1) EXIT ! Repeat unti this condition is satisfied.
       jm=(ju+jl)/2 ! Compute a midpoint,
       IF (ascnd .EQV.(x >=xx(jm))) THEN
          jl=jm  !and replace either the lower limit
       ELSE
          ju=jm ! or the upper limit,as appropriate.
       END IF
    END DO
    IF (x ==xx (1)) THEN ! Then set the output,being careful with the endpoints.
       locate=1
    ELSE IF (x == xx(n))THEN
       locate=n-1
    ELSE
       locate=jl
    END IF
  END FUNCTION locate

  SUBROUTINE spline(x,y,yp1,ypn,y2)
    USE NUM_REC_UTILS,ONLY :assert_eq
    IMPLICIT NONE

    REAL(DP),DIMENSION(:),INTENT(IN)::x,y
    REAL(DP),INTENT(IN)::yp1,ypn
    REAL(DP),DIMENSION(:),INTENT(OUT)::y2
    ! Given arrays x and y of ength N containing a tabulated function,i.e.,
    ! Yi =f(Xi),with X1 < X2 <...<XN ,and given values yp1 and ypn for the first
    ! derivative of the interpolating function at points 1 and N ,respectively,
    ! this routine returns an array y2 of length N that contains the second
    ! derivatives of the interpolating function at the tabulated points Xi .
    ! If yp1 and/or ypn are equal to 1 ×10^30 or larger, the routine is signaled
    ! to set the corresponding boundary condition for a natural spline,with zero
    ! second derivative on that boundary.
    INTEGER(LONG)::n
    REAL(DP),DIMENSION(SIZE(x))::a,b,c,r

    n=assert_eq(SIZE(x),SIZE(y),SIZE(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1) !Set up the tridiagonal equations.
    r(1:n-1)=6.0_DP*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_DP*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    IF (yp1 >0.99e30_DP)THEN !The lower boundary condition is set either to be natural
       r(1)=0.0
       c(1)=0.0
    ELSE ! or else to have a specified first derivative.
       r(1)=(3.0_DP/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    END IF

    IF (ypn >0.99e30_DP)THEN !The upper boundary condition is set either to be natural
       r(n)=0.0
       a(n)=0.0
    ELSE ! or else to have a specified first derivative.
       r(n)=(-3.0_DP/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    END IF
    CALL tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline

  FUNCTION splint(xa,ya,y2a,x)
    USE NUM_REC_UTILS,ONLY :assert_eq,nrerror
    IMPLICIT NONE

    REAL(DP),DIMENSION(:),INTENT(IN)::xa,ya,y2a
    REAL(DP),INTENT(IN)::x
    REAL(DP)::splint
    ! Given the arrays xa and ya ,which tabulate a function (with the XAi' s in
    ! increasing or decreasing order),and given the array y2a ,which is the
    ! output from spline above,and given a value of x ,this routine returns a
    ! cubic-spline interpolated value.The arrays xa ,ya and y2a area of the
    ! same size.
    INTEGER(LONG)::khi,klo,n
    REAL(DP)::a,b,h

    n=assert_eq(SIZE(xa),SIZE(ya),SIZE(y2a),'splint')
    klo=MAX(MIN(locate(xa,x),n-1),1)
    ! We will find the right place in the table by means of locate's bisection
    ! algorithm.This is optimal if sequential calls to this routine are at
    ! random values of x . If sequential calls are in order,and closely spaced,
    ! one would do better to store previous values of klo and khi and test if
    ! they remain appropriate on the next call.
    khi=klo+1 ! klo and khi now bracket the input value of x .
    h=xa(khi)-xa(klo)
    IF (h==0.0) CALL nrerror('bad xa input in splint') !The xa's must be distinct.
    a=(xa(khi)-x)/h ! Cubic spline polynomia is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_DP
  END FUNCTION splint


  SUBROUTINE splie2(x1a,x2a,ya,y2a)
    USE NUM_REC_UTILS,ONLY :assert_eq
    IMPLICIT NONE

    REAL(DP),DIMENSION(:),INTENT(IN)::x1a,x2a
    REAL(DP),DIMENSION(:,:),INTENT(IN)::ya
    REAL(DP),DIMENSION(:,:),INTENT(OUT)::y2a
    ! Given an M ×N tabulated function ya ,and N tabulated independent variables
    ! x2a ,this routine constructs one-dimensional natural cubic splines of the
    ! rows of ya and returns the second derivatives in the M ×N array y2a .
    ! (The array x1a is included in the argument list merely for consistency
    ! with routine splin2.)
    INTEGER(LONG)::j,m,ndum

    m=assert_eq(SIZE(x1a),SIZE(ya,1),SIZE(y2a,1),'splie2:m')
    ndum=assert_eq(SIZE(x2a),SIZE(ya,2),SIZE(y2a,2),'splie2:ndum')
    DO j=1,m
       CALL spline(x2a,ya(j,:),1.0e30_DP,1.0e30_DP,y2a(j,:)) ! Values 1.0E30 signal a natural spline.
    END DO
  END SUBROUTINE splie2

  FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
    USE NUM_REC_UTILS,ONLY :assert_eq
    IMPLICIT NONE

    REAL(DP),DIMENSION(:),INTENT(IN)::x1a,x2a
    REAL(DP),DIMENSION(:,:),INTENT(IN)::ya,y2a
    REAL(DP),INTENT(IN)::x1,x2
    REAL(DP)::splin2
    ! Given x1a ,x2a ,ya as described in splie2 and y2a as produced by that
    ! routine;and given a desired interpolating point x1 ,x2 ;this routine
    ! returns an interpolated function value by bicubic spline interpolation.
    INTEGER(LONG)::j,m,ndum
    REAL(DP),DIMENSION(SIZE(x1a))::yytmp,y2tmp2

    m=assert_eq(SIZE(x1a),SIZE(ya,1),SIZE(y2a,1),'splin2:m')
    ndum=assert_eq(SIZE(x2a),SIZE(ya,2),SIZE(y2a,2),'splin2:ndum')
    DO j=1,m
       yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
       ! Perform m evaluations of the row splines constructed by splie2 ,using
       ! the one-dimensiona spline evaluator splint.
    END DO
    CALL spline(x1a,yytmp,1.0e30_DP,1.0e30_DP,y2tmp2)
    ! Construct the one-dimensional column spline and evaluate it.
    splin2=splint(x1a,yytmp,y2tmp2,x1)
  END FUNCTION splin2


END MODULE NUM_REC
