subroutine dssima(ymiss, timevar, zt, tt, rtv, qt, vt, ft, kt, &
kinf, finf, d, j, p, m, n, r, tolf, etahat)


    implicit none

    integer, intent(in) :: d, j, p, m, n,r
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer :: t, i
    double precision, intent(in) :: tolf
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rtv
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(p,n) ::  vt,ft
    double precision, intent(in), dimension(m,p,n) :: kt
    double precision, intent(in),dimension(m,p,d) ::  kinf
    double precision, intent(in), dimension(p,d) ::  finf
    double precision, intent(inout), dimension(r,n) :: etahat
    double precision, dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
    double precision, dimension(m,m) :: linf,l0, lt,im
    double precision, dimension(m,d+1) :: rt0,rt1
    double precision, dimension(m) :: rrec,rrec1,rhelp,help

    double precision, external :: ddot

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do

    rrec = 0.0d0

    rt(1:m,n+1) = 0.0d0

    do t = n, d+1, -1 !do until diffuse tts
            call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,rt(1:m,t+1),1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(1:r,t),1)
        do i = p, 1 , -1
        if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                lt = im
                call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = I -kz
                call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                rrec = rhelp + vt(i,t)/ft(i,t)*zt(i,1:m,(t-1)*timevar(1)+1)
            end if
            end if
        end do

        call dcopy(m,rrec,1,rt(1:m,t),1) !r_t-1 = r_t,0
   
        if(t.GT.1) then
            call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
            rrec = rhelp
        end if

  

 
      
    end do

    if(d.GT.0) then
        t=d
        rt0(1:m,d+1)=rt(1:m,d+1)
        call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,rt0(1:m,t+1),1,0.0d0,help,1)
        call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(1:r,t),1)

        do i = p, (j+1) , -1
            if(ymiss(t,i).EQ.0) then
                if(ft(i,t) .GT. 0.0d0) then
                lt = im
                call dger(m,m,-1.0d0/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !l = i -kz
                call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                rrec=rhelp
                call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1)
                   end if
            end if
        end do

        rrec1 = 0.0d0
  

        do i = j, 1, -1
            if(ymiss(t,i).EQ.0) then
                if(finf(i,t).GT.tolf) then

         
                    linf = im
                    call dger(m,m,-1.0d0/finf(i,t),kinf(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf
                    rhelp = -kt(1:m,i,t)
                    call daxpy(m,ft(i,t)/finf(i,t),kinf(1:m,i,t),1,rhelp,1)
                    l0=0.0d0
                    call dger(m,m,(1.0d0/finf(i,t)),rhelp,1,zt(i,1:m,(t-1)*timevar(1)+1),1,l0,m) !l0
                    call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                    call dcopy(m,rhelp,1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                    call daxpy(m,(vt(i,t)/finf(i,t)),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec1,1)
                    call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                    rrec = rhelp
             
      
                else
        
                    if(ft(i,t).GT.0.0d0) then
                        lt= im
                        call dger(m,m,(-1.0d0)/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                        call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1)
                        rrec = rhelp
                        call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
                        call dgemv('t',m,m,1.0d0,lt,m,rrec1,1,0.0d0,rhelp,1)
                        rrec1=rhelp
                    end if
                end if
            end if
        end do
        rt0(1:m,t) = rrec
        rt1(1:m,t) = rrec1
  
   
        call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
        rrec = rhelp
        call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
        rrec1 = rhelp
   

   

        do t=(d-1), 1, -1
                    call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(4)+1),m,rt0(1:m,t+1),1,0.0d0,help,1)
            call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,help,1,0.0d0,etahat(1:r,t),1)

            do i = p, 1, -1
                if(ymiss(t,i).EQ.0) then
                    if(finf(i,t).GT. tolf) then
            
                        linf = im
                        call dger(m,m,-1.0d0/finf(i,t),kinf(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,linf,m) !linf
                        rhelp = -kt(1:m,i,t)
                        call daxpy(m,ft(i,t)/finf(i,t),kinf(1:m,i,t),1,rhelp,1)
                        l0=0.0d0
                        call dger(m,m,(1.0d0/finf(i,t)),rhelp,1,zt(i,1:m,(t-1)*timevar(1)+1),1,l0,m) !l0

                        call dgemv('t',m,m,1.0d0,linf,m,rrec1,1,0.0d0,rhelp,1) !rt1
                        call dcopy(m,rhelp,1,rrec1,1)
                        call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
                        call daxpy(m,vt(i,t)/finf(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec1,1)
                        call dgemv('t',m,m,1.0d0,linf,m,rrec,1,0.0d0,rhelp,1) !rt0
                        rrec = rhelp
        
                    else
                        if(ft(i,t).GT. 0.0d0) then
               
                            lt= im
                            call dger(m,m,(-1.0d0)/ft(i,t),kt(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(1)+1),1,lt,m) !lt = I -Kt*Z/Ft
                            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
                            rrec = rhelp
                            call daxpy(m,vt(i,t)/ft(i,t),zt(i,1:m,(t-1)*timevar(1)+1),1,rrec,1) !r0 = Z'vt/Ft - Lt'r0
                            call dgemv('t',m,m,1.0d0,lt,m,rrec1,1,0.0d0,rhelp,1)
                            rrec1=rhelp
               
                        end if

                    end if
                end if
            end do
          
      
            rt0(1:m,t) = rrec
            rt1(1:m,t) = rrec1
      
      
            if(t.GT.1) then
                call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
                rrec = rhelp
                call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(3)+1),m,rrec1,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
                rrec1 = rhelp
         
            end if

     
        end do
    end if


end subroutine dssima
