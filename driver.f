
      program cute_driver
      
      implicit double precision (a-h, o-z)

      parameter(nmax=12005,mmax=11997,nmmax=nmax+mmax,kmx=nmax,mbar=5)
c  this one for GASOIL
c     parameter(mxws=35000000,mxlws=1000000)
c  this one for PINENE
c     parameter(mxws=12000000,mxlws=1000000)
c  this one for MARINE and METHANOL
c     parameter(mxws=9000000,mxlws=1000000)
c  this one otherwise
      parameter(mxws=5000000,mxlws=1000000)
c  mxws >= 2*mxlws should hold
      dimension x(nmmax),al(nmmax),bl(nmmax),bu(nmmax),v(mbar),
     *  ws(mxws),lws(mxlws)
      character cstype(mmax)
      character*10 pname
      character ch

      common/epsc/eps,tol,emin
      common/defaultc/ainfty,ubd,mlp,mxf
      common/wsc/kk,ll,kkk,lll,mxws_,mxlws_
      common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt
      common/refactorc/mc,mxmc
      common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
      common/pnamec/pname
      common/ngrc/mxgr
      common/maxac/maxa

      mxws_=mxws
      mxlws_=mxlws
      maxa=mxlws

      call initialize(n,m,x,bl,bu,nmax,mmax,ws,lws,ws(mxlws+1))
      print *,pname

      v(1)=1.D0
      nv=1

      maxu=n
      maxiu=2*maxa
      maxla=maxa+m+3
      maxg=min(n,mbar)+1
      kmax=min(kmx,n)

      fmin=-ainfty
      maxit=999
      iprint=1
      nout=0
      mxmc=25
      mxgr=100
      mxgr=400
c     mxf=100
      mlp=50
      ubd=1.D1
      rho=1.D1
      htol=1.D-6
      rgtol=1.D-5
c     tol=1.D-8

c     do i=1,n
c       al(i)=1.D-2*x(i)
c       al(i)=1.D-2
c     enddo
c     call checkd(n,m,x,al,ws,lws,maxa,maxla,maxu,maxiu,
c    *  mxws,mxlws,1.D-8)
c     stop

   10 continue

      call filterSD(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv,
     *  maxa,maxla,maxu,maxiu,kmax,maxg,rho,htol,rgtol,maxit,iprint,
     *  nout,ifail)

      print 1,'ifail =',ifail
      if(ifail.eq.4.and.h.gt.ubd)then
        ubd=11.D-1*h
        goto10
      endif

c     print *,'cstype =',(cstype(i),i=1,m)
      print 1,'number of function and gradient calls =',nft,ngt
c     print 4,'x =',(x(i),i=1,n)
c     print 4,'al =',(al(i),i=1,n+m)

c     open(99,status='old',err=998)
c 997 continue
c     read(99,*,end=999)ch
c     goto997
c 998 continue
c     open(99)
c 999 continue
      nout1=99
      nout1=0
      if(ifail.eq.0)then
        if(abs(f).lt.1.D5.and.abs(f).ge.1.D0)then
          write(nout1,1111)pname,n,m,f,h,rgnorm,k,itn,nft,ngt
        else
          write(nout1,2222)pname,n,m,f,h,rgnorm,k,itn,nft,ngt
        endif
      elseif(ifail.eq.3)then
        write(nout1,3333)pname,n,m,hJt,h,rgnorm,k,itn,nft,ngt,ifail
      else
        write(nout1,3333)pname,n,m,f,h,rgnorm,k,itn,nft,ngt,ifail
      endif

      stop

    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))

 1111 format(A9,I4,I5,' ',G14.7,' ',E9.3,E10.3,2I4,2I6)
 2222 format(A9,I4,I5,' ',E14.6,' ',E9.3,E10.3,2I4,2I6)
 3333 format(A9,I4,I5,' ',E14.6,' ',E9.3,E10.3,2I4,2I6,' fail',I2)
      end
