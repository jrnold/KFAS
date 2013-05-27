subroutine approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,  m, &
r, n, theta, u, ytilde, dist,maxiter,tolf,rankp,convtol)

    implicit none

    integer, intent(in) ::  p, m, r, n,dist
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,rankp
    integer ::  i, k,rankp2,d, j

    double precision, intent(in) :: tolf,convtol
    double precision, intent(in), dimension(n) :: u
    double precision, intent(in), dimension(n,1) :: yt
    double precision, intent(in), dimension(1,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(n) :: theta
    double precision, intent(inout), dimension(n,1) :: ytilde
    double precision, intent(inout), dimension(1,1,n) :: ht
    double precision, dimension(1,n) :: vt,ft,finf
    double precision, dimension(m,1,n) :: kt,kinf
    double precision  :: lik,err,mn
    double precision, dimension(m,n) :: ahat,ahat0
    double precision, dimension(m,n+1) :: at
    double precision, dimension(m,m,n+1) :: pt,pinf

    double precision, external :: ddot

    mn = dble(m*n)
    err = 1000.0d0
    ahat0=0.0d0

    k=0
    select case(dist)
        case(1)
            do i=1,n
                if(ymiss(i,1).EQ.0) then
                    ht(1,1,i) = exp(-theta(i))/u(i)
                    ytilde(i,1) =  yt(i,1)*ht(1,1,i) + theta(i) - 1.0d0
                end if
            end do
        case(2)
            do i=1,n
                if(ymiss(i,1).EQ.0) then
                    ht(1,1,i) = (1.0d0 +exp(theta(i)))**2/(u(i)*exp(theta(i)))
                    ytilde(i,1) = theta(i) + ht(1,1,i)*yt(i,1) - 1.0d0 - exp(theta(i))
                end if
            end do
        case(3)
            do i=1,n
                if(ymiss(i,1).EQ.0) then
                    ht(1,1,1:i) = (1.0d0-exp(theta(i)))**2/(u(i)*exp(theta(i)))
                    ytilde(1,1:i) = theta(i) + ht(1,1,i)*yt(i,1) + 1.0d0 - exp(-theta(i))
                end if
            end do
    end select

    do while(err .GT. convtol .AND. k .LT. maxiter)

        k=k+1

        d=0
        j=0
        rankp2 = rankp

        call kfilter(ytilde, ymiss, timevar, zt, ht,tt, rtv, qt, a1, p1, p1inf, p,n,m,r,d,j,&
        at, pt, vt, ft,kt, pinf, finf, kinf, lik, tolf,rankp2)

        call kssim(ymiss, timevar, zt, tt, at, pt, vt, ft, kt, ahat, &
        pinf,  kinf,  finf,  d, j, p, m, n, tolf)


        do i=1,n
            theta(i) = ddot(m,zt(1,1:m,(i-1)*timevar(1)+1),1,ahat(1:m,i),1)
        end do

        err= sqrt(sum((ahat-ahat0)**2))/mn

        ahat0=ahat

        select case(dist)
            case(1)
                do i=1,n
                    if(ymiss(i,1).EQ.0) then
                        ht(1,1,i) = exp(-theta(i))/u(i)
                        ytilde(i,1) =  yt(i,1)*ht(1,1,i) + theta(i) - 1.0d0
                    end if
                end do
            case(2)
                do i=1,n
                    if(ymiss(i,1).EQ.0) then
                        ht(1,1,i) = (1.0d0 +exp(theta(i)))**2/(u(i)*exp(theta(i)))
                        ytilde(i,1) = theta(i) + ht(1,1,i)*yt(i,1) - 1.0d0 - exp(theta(i))
                    end if
                end do
            case(3)
                do i=1,n
                    if(ymiss(i,1).EQ.0) then
                        ht(1,1,1:i) = (1.0d0-exp(theta(i)))**2/(u(i)*exp(theta(i)))
                        ytilde(1,1:i) = theta(i) + ht(1,1,i)*yt(i,1) + 1.0d0 - exp(-theta(i))
                    end if
                end do
        end select
    end do
    maxiter=k


end subroutine approx
