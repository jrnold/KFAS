subroutine glogliku(yt, ymiss, timevar, zt, ht,tt, rt, qt, a1, p1, p1inf,&
m, r, n, lik, tolf,rankp)
    !Subroutine for computing the log-Likelihood of variate linear gaussian state space model
    !Called by R function gSSMlogLikU

    implicit none

    integer, intent(in) ::  m, r, n
    integer, intent(in), dimension(n,1) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: rankp
    integer ::  t, i, d
    double precision, intent(in), dimension(n,1) :: yt
    double precision, intent(in), dimension(1,m,(n-1)*timevar(1)+1) :: zt
    double precision, intent(in), dimension(1,1,(n-1)*timevar(2)+1) :: ht
    double precision, intent(in), dimension(m,m,(n-1)*timevar(3)+1) :: tt
    double precision, intent(in), dimension(m,r,(n-1)*timevar(4)+1) :: rt
    double precision, intent(in), dimension(r,r,(n-1)*timevar(5)+1) :: qt
    double precision, intent(in), dimension(m) :: a1
    double precision, intent(in), dimension(m,m) ::  p1,p1inf
    double precision, intent(in) :: tolf
    double precision, intent(inout) :: lik
    double precision, dimension(m) :: at,arec
    double precision  :: vt,ft,finf
    double precision, dimension(m,1) :: kt,kinf
    double precision, dimension(m,m) :: pt, pinf,mm,im,prec,pirec
    double precision, dimension(m,r) :: mr

    double precision, external :: ddot


    at=0.0d0
    pt=0.0d0
    vt=0.0d0
    ft=0.0d0
    kt=0.0d0
    pinf=0.0d0
    kinf=0.0d0
    finf=0.0d0
    lik = 0.0d0

    pinf(1:m,1:m)=p1inf

    im = 0.0d0
    do i = 1, m
        im(i,i) = 1.0d0
    end do
    d=0
    if(maxval(p1inf) .GT.  0.0d0) then

        pt(1:m,1:m) = p1
        prec = pt(1:m,1:m)
        pirec = pinf(1:m,1:m)
        at(1:m) = a1
        arec = a1
        diffuse: do while(d .LT. n)
            d = d+1
            if(ymiss(d,1).EQ.0) then

                vt = yt(d,1) - ddot(m,zt(1,1:m,(d-1)*timevar(1)+1),1,arec,1) !arec
                call dsymv('u',m,1.0d0,prec,m,zt(1,1:m,(d-1)*timevar(1)+1),1,0.0d0,kt(1:m,1),1) ! kt_t,i = pt_t,i*t(z_t,i)
                ft = ddot(m,zt(1,1:m,(d-1)*timevar(1)+1),1,kt(1:m,1),1)+ ht(1,1,(d-1)*timevar(2)+1)
        
                call dsymv('u',m,1.0d0,pirec,m,zt(1,1:m,(d-1)*timevar(1)+1),1,0.0d0,kinf(1:m,1),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
                finf = ddot(m,zt(1,1:m,(d-1)*timevar(1)+1),1,kinf(1:m,1),1)! finf
         
                ! call dsymv('u',m,1.0d0,pirec,m,zt(1,1:m,(d-1)*timevar+1),1,0.0d0,m1,1)

         
                if (finf .GT. tolf) then
                    call daxpy(m,vt/finf,kinf(1:m,1),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
                    call dsyr('u',m,ft/finf**2,kinf(1:m,1),1,prec,m) !prec = prec +  kinf*kinf'*ft/finf^2
                    call dsyr2('u',m,-1.0d0/finf,kt(1:m,1),1,kinf(1:m,1),1,prec,m) !prec = prec -(kt*kinf'+kinf*kt')/finf
                    !call dger(m,m,(-1.0d0/finf(j,d)),kinf(1:m,j,d),1,kinf(1:m,j,d),1,pirec,m)
                    call dsyr('u',m,(-1.0d0/finf),kinf(1:m,1),1,pirec,m) !pirec = pirec -kinf*kinf'/finf
                    lik = lik - 0.5d0*log(finf)
                    do i = 1, m
                            if(pirec(i,i) .LT. tolf) then
                                pirec(i,1:m) = 0.0d0
                                pirec(1:m,i) = 0.0d0
                            end if
                        end do
                    rankp = rankp -1
                else
                    if (ft .GT. 0.0d0) then
                        call daxpy(m,vt/ft,kt(1:m,1),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)/ft(i,t)
                        call dsyr('u',m,(-1.0d0)/ft,kt(1:m,1),1,prec,m) !prec = prec -kt*kt'/ft
                        lik = lik - 0.5d0*(log(ft) + vt**2/ft)
                    end if
                end if
                if(rankp .EQ. 0) then
                    exit diffuse
                end if

            end if
           
            call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(1:m),1)  !at(:,t+1) = matmul(tt,a_rec)
            call dcopy(m,at(1:m),1,arec,1) ! a_rec = at(:,t+1)
            call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,pt(1:m,1:m),m)
      
     
            call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(5)+1),r,rt(1:m,1:r,(d-1)*timevar(4)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(d-1)*timevar(4)+1),m,1.0d0,pt(1:m,1:m),m)
       
            prec = pt(1:m,1:m)
            call dsymm('r','u',m,m,1.0d0,pirec,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
            call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,pinf(1:m,1:m),m)
            pirec = pinf(1:m,1:m)
     
     
        end do diffuse

        !non-diffuse filtering begins
        prec = prec
      
        call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(3)+1),m,arec,1,0.0d0,at(1:m),1)  !at(:,t+1) = matmul(tt,a_rec)
  
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(3)+1),m,0.0d0,pt(1:m,1:m),m)
 
        call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(5)+1),r,rt(1:m,1:r,(d-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(d-1)*timevar(4)+1),m,1.0d0,pt(1:m,1:m),m)
    
   
        call dcopy(m,at(1:m),1,arec,1) ! a_rec =at(:,t+1)
        prec = pt(1:m,1:m)
  
    end if

   
    !Non-diffuse filtering continues from t=d+1, i=1


    if(d.EQ.0) then
        prec = p1
        arec = a1!   call dcopy(m,a1,1,arec,1)
       !!at(1:m) = a1 !call dcopy(m,a1,1,at(1:m),1) !at(:,1) = a1
       !!pt(1:m,1:m) = p1
    end if
    do t = d+1, n
        if(ymiss(t,1).EQ.0) then
            vt = yt(t,1) - ddot(m,zt(1,1:m,(t-1)*timevar(1)+1),1,arec,1) !variate vt
            call dsymv('u',m,1.0d0,prec,m,zt(1,1:m,(t-1)*timevar(1)+1),1,0.0d0,kt(1:m,1),1) ! p symmetric!
            ft = ddot(m,zt(1,1:m,(t-1)*timevar(1)+1),1,kt(1:m,1),1)+ht(1,1,(t-1)*timevar(2)+1)
            if (ft .GT.  0.0d0) then !ft.NE.0
                call daxpy(m,vt/ft,kt(1:m,1),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t)
                call dsyr('u',m,-1.0d0/ft,kt(1:m,1),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
                lik = lik - 0.5d0*(log(ft) + vt**2/ft)
            end if
        end if
   
        call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(3)+1),m,arec,1,0.0d0,at(1:m),1)  !at(:,t+1) = matmul(tt,a_rec)
        call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,mm,m)
        call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-1)*timevar(3)+1),m,0.0d0,pt(1:m,1:m),m)
  
  
     
        call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(5)+1),r,rt(1:m,1:r,(t-1)*timevar(4)+1),m,0.0d0,mr,m)
        call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(t-1)*timevar(4)+1),m,1.0d0,pt(1:m,1:m),m)
      
  
   
        call dcopy(m,at(1:m),1,arec,1) ! a_rec =at(:,t+1)
        prec = pt(1:m,1:m)
    end do



end subroutine glogliku


