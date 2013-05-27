subroutine ngsmooth(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, theta,&
dist, p, n, m, r, tolf,rankp, nnd,nsim,epsplus,etaplus,&
aplus1,c,tol,info,yo,n2,nsim2,alphahat,alphavar,meany,vary,maxiter,convtol,nd,ndl)


    implicit none

    integer, intent(in) ::  p, m, r, n,nnd,info,nsim,n2,nsim2,dist,maxiter,ndl
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(n2) :: yo
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(ndl) :: nd

    integer, intent(inout) :: rankp
    integer ::  t, i
    double precision, intent(in) :: tolf,tol,convtol
    double precision, intent(in), dimension(n) :: u,theta
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
    double precision, intent(inout), dimension(m,n) :: alphahat
    double precision, intent(inout), dimension(m,m,n) :: alphavar
    double precision, intent(inout), dimension(n) :: meany,vary
    double precision, dimension(n,1) :: ytilde
    double precision, dimension(1,1,n) :: htilde
    double precision, dimension(m,n,nsim2) :: asim
    double precision, dimension(n,nsim2) :: thetasim

    double precision, dimension(n2) :: dn,ueth,uethk
    double precision, dimension(nsim2) :: w
    double precision :: help

    double precision, external :: ddot

    htilde=0.0d0
    ytilde=0.0d0
    call approx(yt, ymiss, timevar, zt, tt, rtv, htilde, qt, a1, p1,p1inf, p,  m, &
    r, n, theta, u, ytilde, dist,maxiter,tolf,rankp,convtol)

    call alphasim(ymiss,timevar, ytilde, zt, htilde, tt, rtv, qt, a1, p1, &
    p1inf, nnd, nsim, epsplus, etaplus, aplus1, p, &
    n, m, r, info, tolf,rankp,c,asim,tol,nd,ndl)

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
                (exp(-0.5d0/htilde(1,1,yo)*( (ytilde(yo,1)-thetasim(yo, i))**2 - dn ))))
            end do
        case(2)
            ueth = exp(theta(yo))/(1.0d0+exp(theta(yo)))
            do i=1,nsim2
                uethk = exp(thetasim(yo,i))/(1.0d0+exp(thetasim(yo,i)))
                w(i) = product( ((uethk/ueth)**yt(yo,1)*((1.0d0-uethk)/(1.0d0-ueth))**(u(yo)-yt(yo,1))) &
                /(exp(-0.5d0/htilde(1,1,yo)*((ytilde(yo,1)-thetasim(yo, i))**2 -dn ))))
            end do
        case(3)
            ueth = 1-exp(theta(yo))
    end select


    w = w/sum(w)

    alphahat=0.0d0
    alphavar=0.0d0

    do i = 1, nsim2
        alphahat(1:m,1:n) = alphahat(1:m,1:n) + asim(1:m,1:n,i)*w(i)
    end do
    do t = 1, n
        do i = 1, nsim2
            call dsyr('u',m,w(i),asim(1:m,t,i),1,alphavar(1:m,1:m,t),m)
        end do
        call dsyr('u',m,-1.0d0,alphahat(1:m,t),1,alphavar(1:m,1:m,t),m)
        do i = 2, m
            alphavar(i,1:(i-1),t)=alphavar(1:(i-1),i,t)
        end do
    end do
    meany=0.0d0
    vary=0.0d0
    select case(dist)
        case(1)
            do t = 1, n
                do i = 1, nsim2
                    meany(t) = meany(t) + exp(ddot(m,zt(1,1:m,(t-1)*timevar(1)+1),1,asim(1:m,t,i),1))*w(i)
                end do
            end do
            meany = meany*u
            vary = meany
        case(2)
            do t = 1, n
                do i = 1, nsim2
                    help = exp(ddot(m,zt(1,1:m,(t-1)*timevar(1)+1),1,asim(1:m,t,i),1))
                    vary(t) =  vary(t) + help/((1+help)**2)*w(i)
                    meany(t) = meany(t) + help/(1+help)*w(i)
                end do
            end do
            meany = meany*u
            vary = vary*u
    end select




end subroutine ngsmooth
