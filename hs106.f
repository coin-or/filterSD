
christen this file hs106.f

      program hs106_driver
      
      implicit double precision (a-h, o-z)

      parameter (n=8,m=6,nm=n+m,mbar=5,maxa=25)
      parameter (mxws=30000, mxlws=5000)
      parameter (t=25.D-4,tm=-t)

      double precision x(nm),bl(nm),bu(nm),al(nm),ws(mxws),v(mbar)
      integer lws(0:mxlws-1)
      character cstype(m)

      common/defaultc/ainfty,ubd,mlp,mxf
      common/wsc/kk,ll,kkk,lll,mxws_,mxlws_
      common/refactorc/mc,mxmc
      common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt
      common/ngrc/mxgr

c  setting indices and pointers of la(0:*) and constant elements of a(*)
      data lws(0:25)
     *  /26,1,2,3,4,5,6,7,8,4,6,4,5,7,5,8,1,4,6,2,4,5,7,3,5,8/
      data lws(26:33)/1,9,11,14,16,19,23,26/
      data ws(2:26)/3*1.D0,5*0.D0,t,t,tm,t,t,
     *  -1.D-2,1.D-2,0.D0,833.33252D0,3*0.D0,125.D1,4*0.D0/

      data bl(1:n)/1.D2,2*1.D3,5*1.D1/
      data bu(1:n)/3*1.D4,5*1.D3/
      data x(1:n)/3*5.D3,2.D2,35.D1,15.D1,225.D0,425.D0/
      data al(1:n)/n*1.D-2/

c  solution and multipliers are
c     data x(1:n)/579.3167, 1359.943, 5110.071, 182.0174,
c    *  295.5985, 217.9799, 286.4162, 395.5979/
c     al = [1964.046;5210.645;5110.092;8.475914E-03;9.578792E-03;0.01]

      mxws_=mxws
      mxlws_=mxlws
      maxu=1
      maxiu=0
      maxla=maxa+m+3
      maxg=min(n,mbar)+1
      kmax=n

c  initializing second copy of a(*)
      do i=maxu+1,maxu+maxa
        ws(maxa+i)=ws(i)
      enddo

      do i=n+1,nm
        bl(i)=-ainfty
        bu(i)=0.D0
      enddo

c  user information passed to subroutines
      ws(1)=t

      nv=1
      v(1)=1.D0

      fmin=-ainfty
      maxit=60
      iprint=1
      nout=0
      mxmc=25
      mxgr=100
      rho=1.D2
      ubd=1.D5
      htol=1.D-6
      rgtol=1.D-4

c  call to check derivatives (now commented out since derivatives are correct)

c     call checkd(n,m,x,al,ws,lws,maxa,maxla,maxu,maxiu,
c    *  mxws,mxlws,1.D-8)
c     stop

      call filterSD(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv,
     *  maxa,maxla,maxu,maxiu,kmax,maxg,rho,htol,rgtol,maxit,iprint,
     *  nout,ifail)

      print 1,'ifail =',ifail
      print 1,'number of function and gradient calls =',nft,ngt
      print 4,'x =',(x(i),i=1,n)
      print 4,'al =',(al(i),i=1,nm)

      stop

    1 format(A,15I5)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
      end

      subroutine functions(n,m,x,f,c,user,iuser)
      implicit double precision (a-h,o-z)
      dimension x(*),c(*),user(*),iuser(*)
c     print *,'enter functions'
c objective function
      f=x(1)+x(2)+x(3)
c constraint functions
      t=user(1)
      c(1)=t*(x(4)+x(6))-1.D0
      c(2)=t*(x(5)+x(7)-x(4))-1.D0
      c(3)=1.D-2*(x(8)-x(5))-1.D0
      c(4)=1.D2*x(1)+833.33252D0*x(4)-x(1)*x(6)-83333.333D0
      c(5)=125.D1*(x(5)-x(4))+x(2)*(x(4)-x(7))
      c(6)=125.D4-25.D2*x(5)+x(3)*(x(5)-x(8))
      return
    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(20I4))
    4 format(A/(6E15.7))
      end

      subroutine gradients(n,m,x,a,user,iuser)
      implicit double precision(a-h,o-z)
      dimension x(*),a(*),user(*),iuser(*)
c     print *,'enter gradients'
c  grad f
c     do i=1,3
c       a(i)=1.D0
c     enddo
c     do i=4,n
c       a(i)=0.D0
c     enddo
c  Jacobian matrix
      a(16)=1.D2-x(6)
      a(18)=-x(1)
      a(19)=x(4)-x(7)
      a(20)=x(2)-125.D1
      a(22)=-x(2)
      a(23)=x(5)-x(8)
      a(24)=x(3)-25.D2
      a(25)=-x(3)
      return
    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(20I4))
    4 format(A/(6E15.7))
      end

