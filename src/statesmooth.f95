subroutine statesmooth(ymiss, timevar, zt, tt, p, n, m, d, j,&
at, pt, vt, ft, kt, pinf, kinf,finf, tolf, ahat, vvt)

    implicit none

    integer, intent(in) :: d, j, p, m, n
    integer :: t, i
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(inout), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,n+1) :: at
    double precision, intent(in), dimension(m,m,n+1) :: pt
    double precision, intent(in), dimension(p,n) ::  vt,ft
    double precision, intent(in), dimension(m,p,n) :: kt
    double precision, intent(in), dimension(m,m,d+1) ::  pinf
    double precision, intent(in),dimension(m,p,d) ::  kinf
    double precision, intent(in), dimension(p,d) ::  finf
    double precision, intent(in) :: tolf
    double precision, intent(inout), dimension(m,n) :: ahat
    double precision, intent(inout), dimension(m,m,n) :: vvt
    double precision, dimension(m,m) :: linf,lt,l0
    double precision, dimension(m,m) :: nrec,nrec1,nrec2,im,mm,mm2
    double precision, dimension(m) :: rrec,rrec1,rhelp
    double precision, external :: ddot

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do

    rrec = 0.0d0
    nrec = 0.0d0



    do t = n, d+1, -1 !do until diffuse starts
        do i = p, 1 , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then

                    lt = im
                    call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = I -kz !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                    rrec = rhelp + vt(i,t)/ft(i,t)*zt(i,1:m,(t-1)*timevar(1)+1)
                    call dsymm('l','u',m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
                    call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
                    call dger(m,m,(1.0d0/ft(i,t)),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m) ! n = n+z'z/f
                end if
            end if
        end do


        call dcopy(m,at(1:m,t),1,ahat(1:m,t),1) !ahat = at
        call dsymv('u',m,1.0d0,pt(1:m,1:m,t),m,rrec,1,1.0d0,ahat(1:m,t),1) !ahat = ahat+pt*r_t-1
        vvt(1:m,1:m,t) = pt(1:m,1:m,t)
        call dsymm('l','u',m,m,1.0d0,pt(1:m,1:m,t),m,nrec,m,0.0d0,mm,m) !pt*n_t-1
        call dsymm('r','u',m,m,-1.0d0,pt(1:m,1:m,t),m,mm,m,1.0d0,vvt(1:m,1:m,t),m) !pt*n_t-1*pt
        if(t.GT.1) then
            call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
            rrec = rhelp
            call dsymm('l','u',m,m,1.0d0,nrec,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,mm,m) !n*t
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,mm,m,0.0d0,nrec,m) !n_t,p = t'nt
        end if
    end do


    if(d.GT.0) then
        t=d


        do i = p, (j+1) , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                    lt = im
                    call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = i -kz !!!!!!!!!!
                    call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                    rrec=rhelp
                    call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1)
                    call dgemm('n','n',m,m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
                    call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
                    call dger(m,m,(1.0d0)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m) ! n = n+z'z/f
                end if
            end if
        end do

        rrec1 = 0.0d0
        nrec1 = 0.0d0
        nrec2 = 0.0d0

        do i = j, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t).GT.tolf) then
                    linf = im
                    call dger(m,m,-1.0d0/finf(i,t),kinf(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf = I- kinf*z/finf !!!!!!!!!!
                    rhelp = -kt(1:m,i,t)
                    call daxpy(m,ft(i,t)/finf(i,t),kinf(1:m,i,t),1,rhelp,1) !rhelp = -kt + ft/finf*kinf
                    l0=0.0d0
                    call dger(m,m,(1.0d0/finf(i,t)),rhelp,1,zt(i,1:m,(t-1)*timevar(1)+1),1,l0,m) !l0=  (-kt + ft/finf*kinf)*z/finf


                    call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                    call dcopy(m,rhelp,1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                    call daxpy(m,(vt(i,t)/finf(i,t)),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                    rrec = rhelp

                    call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
                    call dger(m,m,-1.0d0*ft(i,t)/(finf(i,t)**2.0d0),zt(i,1:m,(t-1)*timevar(1)+1)&
                    ,1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*ft/finf^2
                    call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*ft/finf^2 + l0'*nt0*l0
                    call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
                    call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
                    call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linf,m,0.0d0,mm,m) !mm = nt1'*linf
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
      
                    call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linf,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
                    call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
                    call dger(m,m,(1.0d0)/finf(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec1,m)
                    !nt1 = linf'nt1'linf + z'z/finf
                    call dsymm('l','u',m,m,1.0d0,nrec,m,linf,m,0.0d0,mm,m) !mm= nt0*linf
                    call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf
                    call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm
      
                else
                    if(ft(i,t).GT.0.0d0) then
                        lt= im
                        call dger(m,m,(-1.0d0)/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                        call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp
                        call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
                        call dgemv('t',m,m,1.0d0,lt,m,rrec1,1,0.0d0,rhelp,1)
                        rrec1=rhelp

                        call dgemm('t','n',m,m,m,1.0d0,lt,m,nrec,m,0.0d0,mm,m) !mm =lt'*nt0
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,lt,m,0.0d0,nrec,m) !nt0 = lt'*nt0*lt
                        call dger(m,m,(1.0d0)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m)  !nt0 = z'z/ft+lt'*nt0*lt
                        call dgemm('n','n',m,m,m,1.0d0,nrec1,m,lt,m,0.0d0,mm,m) !mm = nt1*lt
                        nrec1 = mm
                        call dgemm('n','n',m,m,m,1.0d0,nrec2,m,lt,m,0.0d0,mm,m) !mm = nt1*lt
                        nrec2 = mm
                    end if
                end if
            end if
        end do

        call dcopy(m,at(1:m,t),1,ahat(1:m,t),1) !ahat = at
        call dgemv('n',m,m,1.0d0,pt(1:m,1:m,t),m,rrec,1,1.0d0,ahat(1:m,t),1) !ahat = at + pt * rt0_t
        call dgemv('n',m,m,1.0d0,pinf(1:m,1:m,t),m,rrec1,1,1.0d0,ahat(1:m,t),1) !ahat = at + pt * rt0_t + pinf*rt1_t
        vvt(1:m,1:m,t) = pt(1:m,1:m,t)
        call dgemm('n','n',m,m,m,1.0d0,pt(1:m,1:m,t),m,nrec,m,0.0d0,mm,m) !mm = pt*nt0
        call dgemm('n','n',m,m,m,-1.0d0,mm,m,pt(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = pt - pt*nt0*pt
        call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nrec1,m,0.0d0,mm,m) !mm = pinf*nt1
        call dgemm('n','n',m,m,m,-1.0d0,mm,m,pt(1:m,1:m,t),m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pt
        vvt(1:m,1:m,t) = vvt(1:m,1:m,t) + mm2 + transpose(mm2) !vvt = pt - pt*nt0*pt  -pinf*nt1*pt - t(pinf*nt1*pt)
        call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nrec2,m,0.0d0,mm,m) !mm = pinf*nt2
        call dgemm('n','n',m,m,m,-1.0d0,mm,m,pinf(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = vvt - pinf*nt2*pinf
        if(t.GT.1) then  !LIS�TTY 6.11.11
            call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
            rrec = rhelp
            call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
            rrec1 = rhelp
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec2,m,0.0d0,mm,m) !mm =t'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec2,m) !nt2 = t'*nt2*t
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec1,m,0.0d0,mm,m) !mm =t'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec1,m) !nt2 = t'*nt2*t
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t
        end if

        do t=(d-1), 1, -1
            do i = p, 1, -1
                if(ymiss(t,i).EQ.0) then
                    if(finf(i,t).GT. tolf) then
                        linf = im !linf = I
                        call dger(m,m,-1.0d0/finf(i,t),kinf(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf
                        rhelp = -kt(1:m,i,t)
                        call daxpy(m,ft(i,t)/finf(i,t),kinf(1:m,i,t),1,rhelp,1)
                        l0=0.0d0
                        call dger(m,m,(1.0d0/finf(i,t)),rhelp,1,zt(i,1:m,(t-1)*timevar(1)+1),1,l0,m) !l0

                        call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                        call dcopy(m,rhelp,1,rrec1,1)
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                        call daxpy(m,vt(i,t)/finf(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec1,1)!rt1
                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                        rrec = rhelp
            
                        call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,linf,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
            
                        call dger(m,m,-1.0d0*ft(i,t)/(finf(i,t)**2.0d0),&
                        zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec2,m) !nt2 = linf'nt2'linf - z'z*ft/finf^2
  
                        call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
           
                        call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf - z'z*ft/finf^2 + l0'*nt0*l0
                        call dgemm('t','n',m,m,m,1.0d0,linf,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
                        call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linf,m,0.0d0,mm,m) !mm = nt1'*linf
                        call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
          
                        call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linf,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
                        call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
                        call dger(m,m,(1.0d0)/finf(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec1,m)
                        !nt1 = linf'nt1'linf + z'z/finf
                        call dsymm('l','u',m,m,1.0d0,nrec,m,linf,m,0.0d0,mm,m) !mm= nt0*linf
                        call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf
            
                        call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm
            
           
                    else
                        if(ft(i,t).GT.0.0d0) then !lis�tty 12.1.2012
                            lt= im
                            call dger(m,m,(-1.0d0)/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
                            rrec = rhelp
                            call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
                            call dgemv('t',m,m,1.0d0,lt,m,rrec1,1,0.0d0,rhelp,1)
                            rrec1=rhelp
               
                            call dgemm('t','n',m,m,m,1.0d0,lt,m,nrec,m,0.0d0,mm,m) !mm =lt'*nt0
                            call dgemm('n','n',m,m,m,1.0d0,mm,m,lt,m,0.0d0,nrec,m) !nt0 = lt'*nt0*lt
                            call dger(m,m,(1.0d0)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m)  !nt0 = z'z/ft+lt'*nt0*lt
                            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,lt,m,0.0d0,mm,m) !mm = nt1*lt
                            nrec1 = mm
                            call dgemm('n','n',m,m,m,1.0d0,nrec2,m,lt,m,0.0d0,mm,m) !mm = nt2*lt
                            nrec2 = mm
                        end if
                    end if
                end if
            end do
          
      

      
            call dcopy(m,at(1:m,t),1,ahat(1:m,t),1) !ahat = at
            call dgemv('n',m,m,1.0d0,pt(1:m,1:m,t),m,rrec,1,1.0d0,ahat(1:m,t),1) !ahat = at + pt * rt0_t
            call dgemv('n',m,m,1.0d0,pinf(1:m,1:m,t),m,rrec1,1,1.0d0,ahat(1:m,t),1) !ahat = at + pt * rt0_t + pinf*rt1_t
      
            vvt(1:m,1:m,t) = pt(1:m,1:m,t)
            call dgemm('n','n',m,m,m,1.0d0,pt(1:m,1:m,t),m,nrec,m,0.0d0,mm,m) !mm = pt*nt0
            call dgemm('n','n',m,m,m,-1.0d0,mm,m,pt(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = pt - pt*nt0*pt
            call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nrec1,m,0.0d0,mm,m) !mm = pinf*nt1
            call dgemm('n','n',m,m,m,-1.0d0,mm,m,pt(1:m,1:m,t),m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pt
            vvt(1:m,1:m,t) = vvt(1:m,1:m,t) + mm2 + transpose(mm2) !vvt = pt - pt*nt0*pt  -pinf*nt1*pt - t(pinf*nt1*pt)
            call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nrec2,m,0.0d0,mm,m) !mm = pinf*nt2
            call dgemm('n','n',m,m,m,-1.0d0,mm,m,pinf(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = vvt - pinf*nt2*pinf
      
            if(t.GT.1) then
                call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
                rrec = rhelp
                call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
                rrec1 = rhelp
                call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec2,m,0.0d0,mm,m) !mm =t'*nt2
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec2,m) !nt2 = t'*nt2*t
                call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec1,m,0.0d0,mm,m) !mm =t'*nt2
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec1,m) !nt2 = t'*nt2*t
                call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(3)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t
            end if

        end do
    end if


end subroutine statesmooth

