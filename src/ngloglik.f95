subroutine ngloglik(yt, ymiss, timevar, zt, ht, tt, rtv, qt, a1, p1,p1inf, p,m,&
r, n, lik, theta, u, ytilde, dist,maxiter,tolf,rankp,convtol, control, &
nnd,nsim,epsplus,etaplus,aplus1,c,tol,info,antit,sim,yo,n2,nsim2,eg,nd,ndl)

    implicit none

    integer, intent(in) ::  p, m, r, n,nnd,info,antit,nsim,sim,n2,nsim2,dist,ndl
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(n2) :: yo
    integer, intent(in), dimension(ndl) :: nd
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,rankp,control
    integer ::  t, i, rankp2, d, j
    double precision, intent(in) :: tolf,convtol,tol
    double precision, intent(in), dimension(n) :: u
    double precision, intent(in), dimension(n,1) :: yt
    double precision, intent(in), dimension(1,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in), dimension(m,nsim) :: aplus1
    double precision, intent(in),dimension(nsim) :: c
    double precision, intent(inout), dimension(1,1,n) :: ht
    double precision, intent(inout), dimension(n,1) :: ytilde
    double precision, intent(inout), dimension(n) :: theta
    double precision, intent(inout), dimension(1,n,nsim) :: epsplus
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout) :: lik
    double precision, dimension(m,n,nsim2) :: asim
    double precision, dimension(n,nsim2) :: thetasim
    double precision, dimension(n2) :: dn,ueth,uethk
    double precision, dimension(nsim2) :: w
    double precision, dimension((nsim2-1)*control+1) :: wc
    double precision, dimension(1,1,n) :: epshatvar
    double precision, dimension(1,n) :: vt,ft,finf
    double precision, dimension(m,1,n) :: kt,kinf
    double precision :: eg,meanf


    double precision, external :: ddot

    rankp2 = rankp

    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,  m, &
            r, n, theta, u, ytilde, dist,maxiter,tolf,rankp2,convtol)

    d=0
    j=0
    ft=0.0d0
    kt=0.0d0
    kinf=0.0d0
    finf=0.0d0
    rankp2 = rankp

    call kfd2(ytilde, ymiss, timevar, zt, ht,tt, rtv, qt, a1, p1, p1inf, vt, ft,kt, &
    finf, kinf, d, j, p, m, r, n, lik, tolf,rankp2)

    if(control .EQ. 1) then
        call epshvar(ymiss, timevar, zt, tt, ft, kt, kinf, finf,&
        d, m, n, tolf, ht,epshatvar)
    end if

    if(sim .EQ. 1) then
        rankp2 = rankp
        if(antit .EQ. 1) then
            call alphasim(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
            p1inf, nnd, nsim, epsplus, etaplus, aplus1, p, &
            n, m, r, info, tolf,rankp2,c,asim,tol,nd,ndl)
        else
            call alphasimnat(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
            p1inf,nnd,nsim,epsplus, etaplus, aplus1, p, &
            n, m, r, info,tolf,rankp2,asim,tol,nd,ndl)
        end if
        do i=1,nsim2
            do t=1,n
                thetasim(t,i) = ddot(m,zt(1,1:m,(t-1)*timevar(1)+1),1,asim(1:m,t,i),1)
            end do
        end do
    else
        thetasim(1:n,1) = theta(1:n)
    end if


    dn= (ytilde(yo,1)-theta(yo))**2

    select case(dist)
        case(1)
            ueth = exp(theta(yo))
            do i=1,nsim2
                w(i) = product( exp( yt(yo,1)*(thetasim(yo, i)-theta(yo))-u(yo)*(exp(thetasim(yo,i))-ueth)) &
                /exp(-0.5d0/ht(1,1,yo)*( (ytilde(yo,1)-thetasim(yo, i))**2 - dn )))
            end do
            if(control .EQ. 1) then
                meanf = -0.125d0*sum(u(yo)*ueth*epshatvar(1,1,yo)**2) !l'''' = -u*exp(theta)
                do i=1,nsim2
                    wc(i) = w(i) + (sum(u(yo)*ueth*(thetasim(yo,i)-theta(yo))**3)/6.0d0 &
                    +sum(u(yo)*ueth*(thetasim(yo,i)-theta(yo))**4)/24.0d0 + meanf) !miinukset supistuneet pois
                end do
                if(sum(wc).LE.0.0d0) then
                    control = -1
                    eg = log(sum(w))
                else
                    eg = log(sum(wc))
                end if
            else
                eg = log(sum(w))
            end if


        case(2)
            ueth = log(1.0d0+exp(theta(yo)))
            do i=1,nsim2               
                w(i) = product(  exp( yt(yo,1)*(thetasim(yo, i)-theta(yo))-&
                u(yo)*(log(1.0d0+exp(thetasim(yo,i)))-ueth))/ &
                exp(-0.5d0/ht(1,1,yo)*((ytilde(yo,1)-thetasim(yo, i))**2 -dn )))
            end do
            if(control .EQ. 1) then
                ueth = u(yo)*exp(theta(yo))*(exp(2.0d0*theta(yo))-&
                4.0d0*exp(theta(yo))+1.0d0)/(1.0d0+exp(theta(yo)))**4
                uethk = u(yo)*exp(theta(yo))*(exp(theta(yo))-1.0d0)/(1.0d0+exp(theta(yo)))**3
                meanf = 0.125d0*sum(ueth*epshatvar(1,1,yo)**2)
                do i=1,nsim2
                    wc(i) = w(i) - (sum(uethk*(thetasim(yo,i)-theta(yo))**3)/6.0d0 + &
                    sum(ueth*(thetasim(yo,i)-theta(yo))**4)/24.0d0 - meanf)
                end do
                if(sum(wc).LE.0.0d0) then
                    control = -1
                    eg = log(sum(w))
                else
                    eg = log(sum(wc))
                end if
            else
                eg = log(sum(w))
            end if
        case(3)
            ueth = 1-exp(theta(yo))
      ! do i=1,nsim2
      !    uethk = 1-exp(thetasim(yo,i))
      !    w(i) = product( ((uethk/ueth)**u(yo)*((1.0d0-uethk)/(1.0d0-ueth))**(yt(yo,1))) &
      !         /(exp(-0.5d0/ht(1,1,yo)*(ytilde(yo,1)-thetasim(yo, i))**2)/dn))
      ! end do
      ! if(control .EQ. 1) then
      !    ueth = u(yo)*exp(theta(yo))*(exp(2.0d0*theta(yo))-4.0d0*exp(theta(yo))+1.0d0)/(1.0d0+exp(theta(yo)))**4
      !    meanf = 0.125d0*sum(ueth*epshatvar(1,1,yo)**2)
      !    do i=1,nsim2
      !       w(i) = w(i) + sum(ueth*(thetasim(yo,i)-theta(yo))**4)/24.0d0 - meanf
      !    end do
      ! end if
      ! eg = log(sum(w))
    end select



end subroutine ngloglik
