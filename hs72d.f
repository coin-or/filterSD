
      program hs72d_driver
      
c program to drive the HS72 test problem, modified to give linear constraints,
c using dense storage format

      implicit double precision (a-h, o-z)

      parameter (maxa=12,n=4,m=2,nm=n+m,mlp=n,mxws=30000,mxlws=5000)

      dimension a(maxa),x(n),bl(nm),bu(nm),g(n),r(nm),
     * w(nm),e(nm),ls(nm),alp(mlp),lp(mlp),ws(mxws),lws(mxlws),v(n)

      character cws

      common/wsc/kk,ll,kkk,lll,mxws_,mxlws_
      common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
      common/mxm1c/mxm1

      data a/4*0.D0,4.D0,2.25D0,1.D0,0.25D0,0.16D0,0.36D0,2*0.64D0/

      parameter(ainfty=1.D20,tol=1.D-12)

      mxws_=mxws
      mxlws_=mxlws
      kk=0
      ll=0

c  set mxm1 (max size of non-trivial block of basis matrix (see denseL.f))
      mxm1=min(m+1,n)

c  set stride
      la=n

      do i=1,n
        x(i)=1.D0
        bl(i)=1.D0/((5-i)*1.D5)
        bu(i)=1.D3
      enddo
      bl(5)=-ainfty
      bl(6)=-ainfty
      bu(5)=4.01D-2
      bu(6)=1.0085D-2

      kmax=4
      maxg=5
      fmin=-ainfty
      rgtol=1.D-5
      mode=0
      mxgr=100
      iprint=1
      nout=0
      v(1)=1.D0
      nv=1

c     x(1)=0.5170432D-02
c     x(2)=0.5569570D-02
c     x(3)=0.5404878D-02
c     x(4)=0.5927444D-02

c     do i=1,n
c       g(i)=1.D-2
c     enddo
c     call checkg(n,x,g,r,w,ws,lws,ch,tol)
c     stop

      call glcpd(n,m,k,kmax,maxg,a,la,x,bl,bu,f,fmin,g,r,w,e,ls,alp,lp,
     * mlp,ipeq,ws,lws,ch,v,nv,rgtol,mode,ifail,mxgr,iprint,nout)

      write(nout,1)'total number of function and gradient calls =',
     *  nfn,ngr
      write(nout,4)'x =',(x(i),i=1,n)
      write(nout,1)'ifail,ipeq,k =',ifail,ipeq,k
      write(nout,4)'al =',(r(abs(ls(j))),j=1,n)

    1 format(A,6I5)
    4 format(A/(5E15.7))
      stop
      end

      subroutine funct(n,x,f,ws,lws,cws)
      implicit double precision (a-h,o-z)
      dimension x(*),ws(*),lws(*)
      character cws(*)
c     print *,'enter funct'
c     print 4,'x =',(x(i),i=1,n)
      f=1.D0+1.D0/x(1)+1.D0/x(2)+1.D0/x(3)+1.D0/x(4)
c     print *,'f =',f
      return
      end

      subroutine grad(n,x,g,ws,lws,cws)
      implicit double precision(a-h,o-z)
      dimension x(*),g(*),ws(*),lws(*)
      character cws(*)
      do i=1,n
        g(i)=-1.D0/x(i)**2
      enddo
      return
      end
