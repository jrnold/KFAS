subroutine distsimnata(ymiss,timevar, yt, zt, tt, rtv, qt, a1, p1, &
p1inf, nnd, nsim, etaplus, aplus1, p, n, m, r, info,tolf,rankp,etasim,tol,nd,ndl)
    !Subroutine for simulating disturbance vectors of linear gaussian state space model
    !Called by R function simDisturbance


    implicit none
  
    integer, intent(inout) :: info,rankp
    integer, intent(in) :: p, m, r, n, nnd,nsim,ndl
    integer ::  t, i, d, j,k
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
     integer, intent(in), dimension(ndl) :: nd
    double precision, intent(in) :: tolf,tol
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
        double precision, dimension(p,p,1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(inout), dimension(r,n,nsim) :: etasim
    double precision, intent(inout), dimension(r,n,nsim) :: etaplus
    double precision, intent(inout), dimension(m,nsim) :: aplus1
    double precision, dimension(n,p) :: yplus
    double precision, dimension(m,n+1) :: at,aplus
    double precision, dimension(p,n) :: vt,ft,finf
    double precision, dimension(m,m,n+1) :: pt,pinf
    double precision, dimension(m,p,n) :: kt,kinf
    double precision, dimension(r,n) :: etahat,etaplushat
    double precision, dimension(r,r,(n-1)*timevar(5)+1) :: cholqt
    double precision, dimension(m,m) :: cholp1
    double precision, dimension(r,r) :: rcholhelp
    double precision :: lik
    double precision, external :: ddot

    external kssim

    aplus=0.0d0

      ht=0.0d0
    call kfilter(yt, ymiss, timevar, zt, ht,tt, rtv, qt, a1, p1, p1inf, p,n,m,r,d,j,&
    at, pt, vt, ft,kt, pinf, finf, kinf, lik, tolf,rankp)



    call dssima(ymiss, timevar, zt, tt, rtv, qt, vt, ft, kt, kinf,  finf,&
    d, j, p, m, n, r, tolf,etahat)


    do t = 1, (n-1)*timevar(5)+1
        if(r.EQ.1) then
            cholqt(1,1,t)=sqrt(qt(1,1,t))
        else
            rcholhelp(1:r,1:r) = qt(1:r,1:r,t)
            call ldl(rcholhelp,r,tol,info)
            if(info .NE. 0) then
                info=2
                return
            end if
            do i=1,r
                cholqt(i,i,t)=sqrt(rcholhelp(i,i))
            end do
            do i=1,r-1
                cholqt((i+1):r,i,t) = rcholhelp((i+1):r,i)*cholqt(i,i,t)
            end do
        end if
    end do
  
  
    if(nnd.GT.0) then
        if(m.EQ.1) then
            cholp1(1,1)=sqrt(p1(1,1))
        else
            cholp1 = p1
            call ldl(cholp1,m,tol,info)
            if(info .NE. 0) then
                info=3
                return
            end if
            do i=1,m
                cholp1(i,i)=sqrt(cholp1(i,i))
            end do
            do i=1,m-1
                cholp1((i+1):m,i) = cholp1((i+1):m,i)*cholp1(i,i)
            end do
        end if
    end if


    do i = 1, nsim

        if(ndl.GT.0) then
            aplus(nd,1) = a1(nd)
        end if
        if(nnd.GT.0) then
            call dtrmv('l','n','n',m,cholp1,m,aplus1(1:m,i),1)
            aplus(1:m,1) = aplus(1:m,1)+aplus1(1:m,i)
        end if

        do t = 1, n
            do k = 1, p
                if(ymiss(t,k).EQ.0) then
                    yplus(t,k) = ddot(m,zt(k,1:m,(t-1)*timevar(1)+1),1,aplus(1:m,t),1)
                end if
            end do
            call dtrmv('l','n','n',r,cholqt(1:r,1:r,(t-1)*timevar(5)+1),r,etaplus(1:r,t,i),1)
            call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,aplus(1:m,t),1,0.0d0,aplus(1:m,t+1),1)
            call dgemv('n',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,etaplus(1:r,t,i),1,1.0d0,aplus(1:m,t+1),1)
        end do
     
        call kfsim(yplus, ymiss, timevar, zt, tt, a1, at, vt, ft,kt, &
        finf,  kinf,  d, j, p, m, n, tolf)


        call dssima(ymiss, timevar, zt, tt, rtv, qt, vt, ft, kt, kinf,  finf,&
        d, j, p, m, n, r, tolf,etaplushat)


        do t = 1, n
            etasim(1:r,t,i) = etahat(1:r,t) - etaplushat(1:r,t) + etaplus(1:r,t,i)
        end do
    end do

end subroutine distsimnata

