subroutine distsmooth(ymiss, timevar, zt, ht, tt, rtv, qt, &
p, n, m, r, d, j, vt, ft, kt, kinf, finf, tolf,&
epshat,epshatvar,etahat,etahatvar,aug)
    !Subroutine for disturbance smoothing
    !Called by R function disturbanceSmooth

    implicit none

    integer, intent(in) :: d, j, p, m, n,r,aug
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer :: t, i
    double precision, intent(in) :: tolf
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(p,n) ::  vt,ft
    double precision, intent(in), dimension(m,p,n) :: kt
    double precision, intent(in),dimension(m,p,d) ::  kinf
    double precision, intent(in), dimension(p,d) ::  finf
    double precision, intent(inout), dimension(p,n) :: epshat
    double precision, intent(inout), dimension(p,n) :: epshatvar
    double precision, intent(inout), dimension(r,n) :: etahat
    double precision, intent(inout), dimension(r,r,n) :: etahatvar
    double precision, dimension(m,m) :: linf,lt, nrec,im,mr,mr2,mm
    double precision, dimension(m) :: rrec,rhelp,help

    double precision, external :: ddot

    if(aug.EQ.0) then
        do i = 1, p
            do t = 1, n
                epshatvar(i,t) =  ht(i,i,(t-1)*timevar(2)+1)
            end do
        end do
    end if
    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do

    rrec = 0.0d0
    nrec = 0.0d0

    do t = n, d+1, -1 !do until diffuse tts

        call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(1:r,t),1)
        etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(5)+1)
        call dsymm('r','l',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','n',m,r,m,1.0d0,nrec,m,mr,m,0.0d0,mr2,m)
        call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)

        if(t.GT.0) then
            call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
            rrec = rhelp
            call dsymm('l','u',m,m,1.0d0,nrec,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,mm,m) !n*t
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,mm,m,0.0d0,nrec,m) !n_t,p = t'nt
        end if

        do i = p, 1 , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                if(aug.EQ.0) then
                    epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)/ft(i,t)*(vt(i,t)-ddot(m,kt(1:m,i,t),1,rrec,1))
                    call dsymv('u',m,1.0d0/ft(i,t)**2,nrec,m,kt(1:m,i,t),1,0.0d0,rhelp,1)
                    epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                    (1.0d0/ft(i,t)+ddot(m,kt(1:m,i,t),1,rhelp,1))
                end if
                end if
                lt = im
                call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = I -kz
                call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                rrec = rhelp + vt(i,t)/ft(i,t)*zt(i,1:m,(t-1)*timevar(1)+1)
                !call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1)
                call dsymm('l','u',m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
                call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
                call dger(m,m,1.0d0/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m) ! n = n+z'z/f
            
            end if
      
        end do



     ! epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)/ft(i,t)-ddot(m,kt(1:m,i,t),1,rrec,1))
 
      
    end do

    if(d.GT.0) then
        t=d



        call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(1:r,t),1)
        etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(5)+1)
        call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','n',m,r,m,1.0d0,nrec,m,mr,m,0.0d0,mr2,m)
        call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)

        if(t.GT.0) then
            call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
            rrec = rhelp
            call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t
        end if

        do i = p, (j+1) , -1
           if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                if(aug.EQ.0) then
                    epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)/ft(i,t)*(vt(i,t)-ddot(m,kt(1:m,i,t),1,rrec,1))
                    call dsymv('u',m,1.0d0/ft(i,t)**2,nrec,m,kt(1:m,i,t),1,0.0d0,rhelp,1)
                    epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                    (1.0d0/ft(i,t)+ddot(m,kt(1:m,i,t),1,rhelp,1))
                end if
                lt = im
                call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = i -kz
                call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                rrec=rhelp
                call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1)
                call dgemm('n','n',m,m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
                call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
                call dger(m,m,1.0d0/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m) ! n = n+z'z/f
            end if
            end if
        end do

        do i = j, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t).GT.tolf) then
                    if(aug.EQ.0) then
                        epshat(i,t) = -ht(i,i,(t-1)*timevar(2)+1)*ddot(m,kinf(1:m,i,t),1,rrec,1)/finf(i,t)
                        call dsymv('u',m,1.0d0,nrec,m,kinf(1:m,i,t),1,0.0d0,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                        ddot(m,kinf(1:m,i,t),1,rhelp,1)/finf(i,t)**2
                    end if
                    linf = im
                    call dger(m,m,-1.0d0/finf(i,t),kinf(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf
                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt
                    rrec = rhelp
                    call dsymm('l','u',m,m,1.0d0,nrec,m,linf,m,0.0d0,mm,m) !mm= nt*linf
                    call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec,m) !nt = linf'*mm
       


                else

                    if(ft(i,t).GT.0.0d0) then
                        if(aug.EQ.0) then
                            epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)/ft(i,t)-&
                            ddot(m,kt(1:m,i,t),1,rrec,1)/ft(i,t)) !ONKO OIKEASSA PAIKASSA?
                            call dsymv('u',m,1.0d0,nrec,m,kt(1:m,i,t),1,0.0d0,rhelp,1)
                            epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                            (1.0d0/ft(i,t)+ddot(m,kt(1:m,i,t),1,rhelp,1)/ft(i,t)**2)
                        end if
                        lt= im
                        call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                        call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp
                        call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
         
                        call dgemm('t','n',m,m,m,1.0d0,lt,m,nrec,m,0.0d0,mm,m) !mm =lt'*nt
                        call dgemm('n','n',m,m,m,1.0d0,mm,m,lt,m,0.0d0,nrec,m) !nt = lt'*nt*lt
                        call dger(m,m,1.0d0/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m)  !nt = z'z/ft+lt'*nt*lt


                    end if
                end if
            end if
        end do



        do t=(d-1), 1, -1

            call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,rrec,1,0.0d0,help,1)
            call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(1:r,t),1)
            etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(5)+1)
            call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
            call dgemm('n','n',m,r,m,1.0d0,nrec,m,mr,m,0.0d0,mr2,m)
            call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)

 if(t.GT.0) then
                call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
                rrec = rhelp
                call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
                call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t
            end if

            do i =p, 1, -1
                if(ymiss(t,i).EQ.0) then
                    if(finf(i,t).GT. tolf) then
                        if(aug.EQ.0) then
                            epshat(i,t) = -ht(i,i,(t-1)*timevar(2)+1)*ddot(m,kinf(1:m,i,t),1,rrec,1)/finf(i,t)
                            call dsymv('u',m,1.0d0,nrec,m,kinf(1:m,i,t),1,0.0d0,rhelp,1)
                            epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                            ddot(m,kinf(1:m,i,t),1,rhelp,1)/finf(i,t)**2
                        end if
                        linf = im
                        call dger(m,m,-1.0d0/finf(i,t),kinf(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf
                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt
                        rrec = rhelp
                        call dsymm('l','u',m,m,1.0d0,nrec,m,linf,m,0.0d0,mm,m) !mm= nt*linf
                        call dgemm('t','n',m,m,m,1.0d0,linf,m,mm,m,0.0d0,nrec,m) !nt = linf'*mm

           
      
                    else
                        if(ft(i,t).GT.0.0d0) then
                            if(aug.EQ.0) then
                                epshat(i,t) = ht(i,i,(t-1)*timevar(2)+1)*(vt(i,t)/ft(i,t)-&
                                ddot(m,kt(1:m,i,t),1,rrec,1)/ft(i,t))
                                call dsymv('u',m,1.0d0,nrec,m,kt(1:m,i,t),1,0.0d0,rhelp,1)
                                epshatvar(i,t) = ht(i,i,(t-1)*timevar(2)+1)-(ht(i,i,(t-1)*timevar(2)+1)**2)*&
                                (1.0d0/ft(i,t)+ddot(m,kt(1:m,i,t),1,rhelp,1)/ft(i,t)**2 )
                            end if
                            lt= im
                            call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
                            rrec = rhelp
                            call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0

                            call dgemm('t','n',m,m,m,1.0d0,lt,m,nrec,m,0.0d0,mm,m) !mm =lt'*nt
                            call dgemm('n','n',m,m,m,1.0d0,mm,m,lt,m,0.0d0,nrec,m) !nt = lt'*nt*lt
                            call dger(m,m,1.0d0/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m)  !nt = z'z/ft+lt'*nt*lt

            
                        end if
                    end if
                end if
            end do
          
      

        end do
    end if


end subroutine distsmooth
