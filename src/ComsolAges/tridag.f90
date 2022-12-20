module m_tridag
contains
  subroutine tridag(a,b,c,r,u,n)
    implicit none

! inpout + output parameter
    integer(4), intent(in) :: n

    real(8), intent(out), dimension(n) :: u

    real(8), intent(in), dimension(n) :: a, b, c, r


! local vars
    real(8) :: bet
    real(8), dimension(:), allocatable :: gam

    integer(4) :: j

    allocate(gam(n))

    if (b(1) == 0) then
       print *, "tridag.f90: b(1) == 0"
       stop
    end if

    bet = b(1)
    u(1)=r(1)/bet

    do j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       if (bet == 0) then
          print *, "tridag.f90: bet == 0"
          stop
       end if
       u(j)=(r(j)-a(j)*u(j-1))/bet
    end do

    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do

    deallocate(gam)

  end subroutine tridag
end module m_tridag
