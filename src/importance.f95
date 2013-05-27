subroutine importance(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
p, n, m, r, theta, maxiter,tolf,rankp,convtol, nnd,nsim,epsplus,etaplus,&
aplus1,c,tol,info,antit,yo,n2,nsim2,w,asim,nd,ndl)


    implicit none

    integer, intent(in) ::  p, m, r, n,nnd,info,antit,nsim,n2,nsim2,dist,ndl
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(n2) :: yo
    integer, intent(in), dimension(ndl) :: nd
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,rankp
    integer ::  t, i
    double precision, intent(in) :: tolf,convtol,tol
    double precision, intent(in), dimension(n) :: u
    double precision, intent(in), dimension(n,1) :: yt
    double precision, intent(in), dimension(1,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in),dimension(nsim) :: c
    double precision, intent(inout), dimension(1,n,nsim) :: epsplus
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout), dimension(m,nsim) :: aplus1
    double precision, intent(inout), dimension(n) :: theta
    double precision, dimension(1,1,n) :: ht
    double precision, dimension(m,n,nsim2) :: asim
    double precision, dimension(n,nsim2) :: thetasim
    double precision, dimension(n,1) :: ytilde
    double precision, dimension(n2) :: dn,ueth,uethk
    double precision, dimension(nsim2) :: w

    double precision, external :: ddot

    ht=0.0d0

    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,  m, &
    r, n, theta, u, ytilde, dist,maxiter,tolf,rankp,convtol)


    if(antit .EQ. 1) then
        call alphasim(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
        p1inf, nnd, nsim, epsplus, etaplus, aplus1, p, &
        n, m, r, info, tolf,rankp,c,asim,tol,nd,ndl)
    else
        call alphasimnat(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
        p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, &
        n, m, r, info, tolf,rankp,asim,tol,nd,ndl)
    end if

    do i=1,nsim2
        do t=1,n
            thetasim(t,i) = ddot(m,zt(1,1:m,(t-1)*timevar(1)+1),1,asim(1:m,t,i),1)
        end do
    end do

    dn= (ytilde(yo,1)-theta(yo))**2

    select case(dist)
        case(1)
            ueth = exp(theta(yo))
            !dp = ueth**yt(yo,1)*exp(-ueth)
            do i=1,nsim2 !!!
                !uethk = exp(thetasim(yo,i))
                w(i) = product((exp(thetasim(yo, i)-theta(yo))**yt(yo,1)*exp(-u(yo)*(exp(thetasim(yo,i))-ueth)))/ &
                (exp(-0.5d0/ht(1,1,yo)*( (ytilde(yo,1)-thetasim(yo, i))**2 - dn ))))
            end do
        case(2)
            ueth = exp(theta(yo))/(1.0d0+exp(theta(yo)))
            do i=1,nsim2
                uethk = exp(thetasim(yo,i))/(1.0d0+exp(thetasim(yo,i)))
                w(i) = product( ((uethk/ueth)**yt(yo,1)*((1.0d0-uethk)/(1.0d0-ueth))**(u(yo)-yt(yo,1))) &
                /(exp(-0.5d0/ht(1,1,yo)*((ytilde(yo,1)-thetasim(yo, i))**2 -dn ))))
            end do
        case(3)
            ueth = 1-exp(theta(yo))
    end select


end subroutine importance
