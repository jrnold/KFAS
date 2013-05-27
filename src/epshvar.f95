subroutine epshvar(ymiss, timevar, zt, tt, ft, kt, &
kinf, finf, d, m, n, tolf, ht,epshatvar)

    !Subroutine for computing the covariance matrix of smoothed disturbance vectors epsilon of linear gaussian state space model
    !Called by Fortran subroutine efloglik

    implicit none

    integer, intent(in) :: d, m, n
    integer :: t, i
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(1,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(1,1,n) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(1,n) :: ft
    double precision, intent(in), dimension(m,1,n) :: kt
    double precision, intent(in),dimension(m,1,d) ::  kinf
    double precision, intent(in), dimension(1,d) ::  finf
    double precision, intent(in) :: tolf
    double precision, intent(inout), dimension(1,1,n) :: epshatvar
    double precision, dimension(m,m) :: linf, lt,nrec,im,mm
    double precision, dimension(m) :: rhelp

    double precision, external :: ddot


    epshatvar(1,1,1:n) =  ht(1,1,1:n)


    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do


    nrec = 0.0d0

    do t = n, d+1, -1 !do until diffuse tts
        if(ymiss(t,1).EQ.0 .AND. ft(1,t) .GT. 0.0d0) then
            call dsymv('u',m,1.0d0/ft(1,t)**2,nrec,m,kt(1:m,1,t),1,0.0d0,rhelp,1)
            epshatvar(1,1,t) = ht(1,1,t)-(ht(1,1,t)**2)*(1.0d0/ft(1,t)+ddot(m,kt(1:m,1,t),1,rhelp,1))

            lt = im
            call dger(m,m,-1.0d0/ft(1,t),kt(1:m,1,t),1,zt(1,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = I -kz
            call dsymm('l','u',m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
            call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
            call dger(m,m,1.0d0/ft(1,t),zt(1,1:m,(t-1)*timevar(1)+1),1,zt(1,1:m,(t-1)*timevar(1)+1),1,nrec,m) ! n = n+z'z/f
        end if
        if(t.GT.1) then
            call dsymm('l','u',m,m,1.0d0,nrec,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,mm,m) !n*t
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,mm,m,0.0d0,nrec,m) !n_t,p = t'nt
        end if
  
    end do





    do t= d, 1, -1
        if(ymiss(t,1).EQ.0) then
            if(finf(1,t).GT. tolf) then
                call dsymv('u',m,1.0d0,nrec,m,kinf(1:m,1,t),1,0.0d0,rhelp,1)
                epshatvar(1,1,t) = ht(1,1,t)-(ht(1,1,t)**2)*ddot(m,kinf(1:m,1,t),1,rhelp,1)/finf(1,t)**2
                linf = im
                call dger(m,m,-1.0d0/finf(1,t),kinf(1:m,1,t),1,zt(1,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf
                call dsymm('l','u',m,m,1.0d0,nrec,m,linf,m,0.0d0,mm,m) !mm= nt0*linf
                call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm
            else
                if(ft(1,t).GT. 0.0d0) then
                    call dsymv('u',m,1.0d0,nrec,m,kt(1:m,1,t),1,0.0d0,rhelp,1)
                    epshatvar(1,1,t) =  ht(1,1,t)-(ht(1,1,t)**2.0d0)*(1.0d0/ft(1,t)+ddot(m,kt(1:m,1,t),1,rhelp,1)&
                    /ft(1,t)**2)
                    lt= im
                    call dger(m,m,-1.0d0/ft(1,t),kt(1:m,1,t),1,zt(1,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                    call dgemm('t','n',m,m,m,1.0d0,lt,m,nrec,m,0.0d0,mm,m) !mm =lt'*nt0
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,lt,m,0.0d0,nrec,m) !nt0 = lt'*nt0*lt
                    call dger(m,m,1.0d0/ft(1,t),zt(1,1:m,(t-1)*timevar(1)+1),1,zt(1,1:m,(t-1)*timevar(1)+1),1,nrec,m)  !nt0 = z'z/ft+lt'*nt0*lt
                end if
            end if
        end if

        if(t.GT.1) then
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t
        end if

    end do



end subroutine epshvar
