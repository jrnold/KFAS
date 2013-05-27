subroutine kfilter(yt, ymiss, timevar, zt, ht,tt, rt, qt, a1, p1, p1inf, p,n,m,r,d,j,&
at, pt, vt, ft,kt, pinf, finf, kinf, lik, tolf,rankp)
    !Subroutine for Kalman filtering of linear gaussian state space model
    !Called by R function kalmanFilter
    !Called by Fortran subroutines alphasim
    implicit none

    integer, intent(in) ::  p, m, r, n
    integer, intent(inout) :: d, j, rankp
    integer ::  t, i
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    double precision, intent(in), dimension(n,p) :: yt
    double precision, intent(in), dimension(p,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(p,p,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in) :: tolf
    double precision, intent(inout), dimension(m,n+1) :: at
    double precision, intent(inout), dimension(m,m,n+1) :: pt,pinf
    double precision, intent(inout), dimension(p,n) :: vt,ft,finf
    double precision, intent(inout), dimension(m,p,n) :: kt,kinf
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: arec
    double precision, dimension(m,m) :: prec, pirec,im,mm
    double precision, dimension(m,r) :: mr

    double precision, external :: ddot

    lik = 0.0d0

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do
    j=0
    d=0
    pinf(1:m,1:m,1)=p1inf(1:m,1:m)
    if(maxval(pinf(1:m,1:m,1)) .GT.  0.0d0) then

        pt(1:m,1:m,1) = p1
        prec = pt(1:m,1:m,1)
        pirec = pinf(1:m,1:m,1)
        at(1:m,1) = a1
        arec = a1
        diffuse: do while(d .LT. n)
            d = d+1
            do j=1, p
                if(ymiss(d,j) .EQ. 0) then
                    vt(j,d) = yt(d,j) - ddot(m,zt(j,1:m,(d-1)*timevar(1)+1),1,arec,1) !arec
                    call dsymv('u',m,1.0d0,prec,m,zt(j,1:m,(d-1)*timevar(1)+1),1,0.0d0,kt(1:m,j,d),1) ! kt_t,i = pt_t,i*t(z_t,i)
                    ft(j,d) = ddot(m,zt(j,1:m,(d-1)*timevar(1)+1),1,kt(1:m,j,d),1)  + ht(j,j,(d-1)*timevar(2)+1)
                    call dsymv('u',m,1.0d0,pirec,m,zt(j,1:m,(d-1)*timevar(1)+1),1,0.0d0,kinf(1:m,j,d),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                    finf(j,d) = ddot(m,zt(j,1:m,(d-1)*timevar(1)+1),1,kinf(1:m,j,d),1)! finf
                    if (finf(j,d) .GT. tolf) then
                        call daxpy(m,vt(j,d)/finf(j,d),kinf(1:m,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                        call dsyr('u',m,ft(j,d)/(finf(j,d)**2),kinf(1:m,j,d),1,prec,m) !prec = prec +  kinf*kinf'*ft/finf^2
                        call dsyr2('u',m,-1.0d0/finf(j,d),kt(1:m,j,d),1,kinf(1:m,j,d),1,prec,m) !prec = prec -(kt*kinf'+kinf*kt')/finf
                        call dsyr('u',m,-1.0d0/finf(j,d),kinf(1:m,j,d),1,pirec,m) !pirec = pirec -kinf*kinf'/finf
                        lik = lik - 0.5d0*log(finf(j,d))
                        rankp = rankp -1
                        do i = 1, m
                            if(pirec(i,i) .LT. tolf) then
                                pirec(i,1:m) = 0.0d0
                                pirec(1:m,i) = 0.0d0
                            end if
                        end do
                    else
                        finf(j,d) = 0.0d0
                        if(ft(j,d) .GT. 0.0d0) then
                            call daxpy(m,vt(j,d)/ft(j,d),kt(1:m,j,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                            call dsyr('u',m,(-1.0d0)/ft(j,d),kt(1:m,j,d),1,prec,m) !prec = prec -kt*kt'/ft
                            lik = lik - 0.5d0*(log(ft(j,d)) + vt(j,d)**2/ft(j,d))
                        end if
                    end if
                    if(rankp .EQ. 0) then
                        exit diffuse
                    end if
                end if
            end do

            call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(1:m,d+1),1)  !at(:,t+1) = matmul(tt,a_rec)
            call dcopy(m,at(1:m,d+1),1,arec,1) ! a_rec = at(:,t+1)
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,pt(1:m,1:m,d+1),m)

            if(r.GT.1) then
                call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(5)+1),r,rt(1:m,1:r,(d-1)*timevar(4)+1),m,0.0d0,mr,m)
                call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(d-1)*timevar(4)+1),m,1.0d0,pt(1:m,1:m,d+1),m)
            else
                call dger(m,m,qt(1,1,(d-1)*timevar(5)+1),rt(1:m,1,(d-1)*timevar(4)+1),1,rt(1:m,1,(d-1)*timevar(4)+1),1,&
                pt(1:m,1:m,d+1),m)
            end if
            prec = pt(1:m,1:m,d+1)
            call dsymm('r','u',m,m,1.0d0,pirec,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,pinf(1:m,1:m,d+1),m)
            pirec = pinf(1:m,1:m,d+1)
        end do diffuse       


            !non-diffuse filtering begins
        do i = j+1, p
            if(ymiss(d,i).EQ.0) then
                vt(i,d) = yt(d,i) - ddot(m,zt(i,1:m,(d-1)*timevar(1)+1),1,arec,1) !vt
                call dsymv('u',m,1.0d0,prec,m,zt(i,1:m,(d-1)*timevar(1)+1),1,0.0d0,kt(1:m,i,d),1) ! p symmetric!
                ft(i,d) = ddot(m,zt(i,1:m,(d-1)*timevar(1)+1),1,kt(1:m,i,d),1)  +  ht(i,i,(d-1)*timevar(2)+1)
                if (ft(i,d) .GT.  0.0d0) then !ft.NE.0
                    call daxpy(m,vt(i,d)/ft(i,d),kt(1:m,i,d),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    call dsyr('u',m,-1.0d0/ft(i,d),kt(1:m,i,d),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,t)
                    lik = lik - 0.5d0*(log(ft(i,d)) + vt(i,d)**2/ft(i,d))
                end if
            end if
        end do             
            
      
        call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(1:m,d+1),1)  !at(:,t+1) = matmul(tt,a_rec)
      
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,pt(1:m,1:m,d+1),m)
      
        if(r.GT.1) then
            call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(5)+1),r,rt(1:m,1:r,(d-1)*timevar(4)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(d-1)*timevar(4)+1),m,1.0d0,pt(1:m,1:m,d+1),m)
        else
            call dger(m,m,qt(1,1,(d-1)*timevar(5)+1),rt(1:m,1,(d-1)*timevar(4)+1),1,rt(1:m,1,(d-1)*timevar(4)+1),1,&
            pt(1:m,1:m,d+1),m)
        end if
      
        call dcopy(m,at(1:m,d+1),1,arec,1) ! a_rec =at(:,t+1)
        prec = pt(1:m,1:m,d+1)

    end if

    !Non-diffuse filtering continues from t=d+1, i=1


    if(d.EQ.0) then
        prec = p1
        arec = a1!   call dcopy(m,a1,1,arec,1)
        at(1:m,1) = a1 !call dcopy(m,a1,1,at(1:m,1),1) !at(:,1) = a1
        pt(1:m,1:m,1) = p1
    else
        if(d .EQ. n .AND. j .EQ. p+1) then
            j = p
        end if
    end if
    do t = d+1, n
        do i = 1, p
            if(ymiss(t,i).EQ.0) then
                vt(i,t) = yt(t,i) - ddot(m,zt(i,1:m,(t-1)*timevar(1)+1),1,arec,1) !variate vt
                call dsymv('u',m,1.0d0,prec,m,zt(i,1:m,(t-1)*timevar(1)+1),1,0.0d0,kt(1:m,i,t),1) ! p symmetric!
                ft(i,t) = ddot(m,zt(i,1:m,(t-1)*timevar(1)+1),1,kt(1:m,i,t),1) +  ht(i,i,(t-1)*timevar(2)+1)
                if (ft(i,t) .GT.  0.0d0) then !ft.NE.0
                    call daxpy(m,vt(i,t)/ft(i,t),kt(1:m,i,t),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                    call dsyr('u',m,-1.0d0/ft(i,t),kt(1:m,i,t),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
                    lik = lik - 0.5d0*(log(ft(i,t)) + vt(i,t)**2/ft(i,t))
                end if
            end if
        end do
   
        call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at(1:m,t+1),1)  !at(:,t+1) = matmul(tt,a_rec)
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,pt(1:m,1:m,t+1),m)
  
        if(r.GT.1) then
            call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,rt(1:m,1:r,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(t-1)*timevar(4)+1),m,1.0d0,pt(1:m,1:m,t+1),m)
        else
            call dger(m,m,qt(1,1,(t-1)*timevar(5)+1),rt(1:m,1,(t-1)*timevar(4)+1),1,rt(1:m,1,(t-1)*timevar(4)+1),&
            1,pt(1:m,1:m,t+1),m)
        end if
   
        call dcopy(m,at(1:m,t+1),1,arec,1) ! a_rec =at(:,t+1)
        prec = pt(1:m,1:m,t+1)
    end do



end subroutine kfilter
