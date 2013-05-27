subroutine predobs(ymiss,tvz, tvh, zt, ht, at, pt, p, n, m, theta, thetavar)

    implicit none

    integer, intent(in) :: p, m, n,tvz,tvh
    integer ::  t,i
    integer, intent(in), dimension(n,p) :: ymiss
    double precision, intent(in), dimension(p,m,(n-1)*tvz+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*tvh+1) :: ht
    double precision, intent(in), dimension(m,n) :: at
    double precision, intent(in), dimension(m,m,n) :: pt
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, intent(inout), dimension(p,p,n) :: thetavar
    double precision, dimension(p,m) :: pm
    double precision, dimension(p) :: help

    do t = 1, n
        if(sum(ymiss(t,1:p))>0) then
            call dgemv('n',p,m,1.0d0,zt(1:p,1:m,(t-1)*tvz+1),p,at(1:m,t),1,0.0d0,help,1)
            call dsymm('r','u',p,m,1.0d0,pt(1:m,1:m,t),m,zt(1:p,1:m,(t-1)*tvz+1),p,0.0d0,pm,p)
            call dgemm('n','t',p,p,m,1.0d0,pm,p,zt(1:p,1:m,(t-1)*tvz+1),p,0.0d0,thetavar(1:p,1:p,t),p)
            thetavar(1:p,1:p,t) = thetavar(1:p,1:p,t) + ht(1:p,1:p,(t-1)*tvh+1)
            do i=1, p
                if(ymiss(t,i).EQ.1) then
                    theta(t,i) = help(i)
                else
                    thetavar(i,1:p,t)=-1000.0d0
                    thetavar(1:p,i,t)=-1000.0d0
                end if
            end do

        end if
    end do

end subroutine predobs


subroutine signaltheta(tvz, zt, ahat, vt, p, n, m, theta, thetavar,d,states,m2)

    implicit none

    integer, intent(in) :: p, m, n,tvz,d,m2 !,tvh
    integer, intent(in), dimension(m2) :: states
    integer ::  t
    double precision, intent(in), dimension(p,m,(n-1)*tvz+1) :: zt
    !double precision, intent(in), dimension(p,p,(n-1)*tvh+1) :: ht
    double precision, intent(in), dimension(m,n) :: ahat
    double precision, intent(in), dimension(m,m,n) :: vt
    double precision, intent(inout), dimension(n,p) :: theta
    double precision, intent(inout), dimension(p,p,n) :: thetavar
    double precision, dimension(p,m2) :: pm

    do t = (d+1), n
        call dgemv('n',p,m2,1.0d0,zt(1:p,states,(t-1)*tvz+1),p,ahat(states,t),1,0.0d0,theta(t,1:p),1)
        call dsymm('r','u',p,m2,1.0d0,vt(states,states,t),m,zt(1:p,states,(t-1)*tvz+1),p,0.0d0,pm,p)
        call dgemm('n','t',p,p,m2,1.0d0,pm,p,zt(1:p,states,(t-1)*tvz+1),p,0.0d0,thetavar(1:p,1:p,t),p)
        !thetavar(1:p,1:p,t) = thetavar(1:p,1:p,t) + ht(1:p,1:p,(t-1)*tvh+1)
    end do

end subroutine signaltheta

subroutine signalthetang(tvz, zt, asim, w, n, m, thetahat, thetavar,states,m2,nsim)

    implicit none

    integer, intent(in) :: m, n,tvz,m2,nsim
    integer, intent(in), dimension(m2) :: states
    integer ::  t,i
    double precision, intent(in), dimension(1,m,(n-1)*tvz+1) :: zt
    double precision, intent(in), dimension(m,n,nsim) :: asim
    double precision, intent(in), dimension(nsim) :: w
    double precision, intent(inout), dimension(n) :: thetahat
    double precision, intent(inout), dimension(n) :: thetavar
    double precision :: thetasim

    double precision, external :: ddot

    do t = 1, n
        do i = 1, nsim
            thetasim = ddot(m,zt(1,states,(t-1)*tvz+1),1,asim(states,t,i),1)
            thetahat(t) = thetahat(t)  + thetasim*w(i)
            thetavar(t) = thetavar(t) + w(i)*thetasim**2
        end do
    end do
    thetavar = thetavar - thetahat**2




end subroutine signalthetang


