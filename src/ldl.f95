subroutine ldl(a,n,tol,info)
  
    implicit none
  
    integer, intent(in) :: n
    integer, intent(inout) :: info
    integer :: p,k,i
    double precision, dimension(n) :: d,r
    double precision, intent(inout), dimension(n,n) :: a
    double precision, intent(in) :: tol
  
    d=0.0d0
    r=0.0d0
    do k = 1, n
        do p = 1,k-1
            r(p)=d(p)*a(k,p)
        end do
        d(k)=a(k,k) - sum(a(k,1:(k-1))*r(1:(k-1)))
        if(d(k) .GT. tol) then
            do i = k+1, n
                a(i,k)=(a(i,k)-sum(a(i,1:(k-1))*r(1:(k-1))))/d(k)
            end do
        else
            if(d(k) .LT. -0.5d0) then
                info = -1
            end if
            d(k) = 0.0d0
            a((k+1):n,k)=0.0d0
        end if
    end do
    do k = 1,n
        a(k,k) = d(k)
        a(k,(k+1):n) = 0.0d0
    end do


end subroutine ldl
