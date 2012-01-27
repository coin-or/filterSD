Christen this file user.f

      subroutine initialize(n,m,x,bl,bu,nmax,mmax,ws,lws,iuser)
      implicit double precision (a-h,o-z)
      dimension x(*),bl(*),bu(*),ws(*),lws(*),iuser(*)
      character*10 pname
      common/pnamec/pname
      common/maxac/maxa
c     print *,'enter initialize'
c  open SIF output file
      open(file='OUTSDIF.d',unit=7)
c  read SIF output file and set up the constraints
      call csetup(7,6,n,m,x,bl,bu,nmax,ws,ws,ws,
     *  bl(nmax+1),bu(nmax+1),mmax,.false.,.false.,.false.)
      do i=1,m
        bl(n+i)=bl(nmax+i)
        bu(n+i)=bu(nmax+i)
      enddo
c     print 4,'bl =',(bl(i),i=1,n+m)
c     print 4,'bu =',(bu(i),i=1,n+m)
      call unames(n,pname,ws)
c  call CUTE's coordinate format into sparse vector format ...
c  ... note nonzero column indices are in ascending order
c  ... followed by n zero column indices (objective function)
      call csgr(n,m,.false.,1,x,x,mxa,maxa,ws,iuser,lws(1))
c     print 1,'mxa =',mxa
c     print 3,'i =',(iuser(k),k=1,mxa)
c     print 3,'j =',(lws(k),k=1,mxa)
c     print 4,'a(ij) =',(ws(k),k=1,mxa)
      if(3*mxa+m+3.gt.maxa)then
        print *,'not enough space for CUTE Jacobian indices'
        print *,'mxa =',mxa,'   maxa =',maxa
        stop
      endif
      maxa=mxa
      lws(2*maxa+1)=maxa+1
      ip=3*maxa+3
      lws(ip-1)=1
      lws(ip)=n+1
      k=1
      do j=1,m
   10   continue
        if(lws(k).eq.j)then
          k=k+1
          goto10
        endif
        lws(ip+j)=k+n
      enddo
c     print 3,'pointers',(lws(ip+i),i=-1,m)
      do i=maxa-n+1,maxa
        lws(maxa+n+1+i)=iuser(i)
      enddo
      do i=1,maxa-n
        lws(2*maxa+n+1+i)=iuser(i)
      enddo
c     print 3,'indices',(lws(i),i=2*maxa+2,3*maxa+1)
    1 format(A,15I4)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
      return
      end

      subroutine functions(n,m,x,f,c,user,iuser)
      implicit double precision (a-h,o-z)
      dimension x(*),c(*),user(*),iuser(*)
c     print *,'enter functions'
c  call CUTE's objective and constraint function
c     print 4,'x =',(x(i),i=1,n)
      call cfn(n,m,x,f,m,c)
c     print *,'f =',f
c     print 4,'c =',(c(i),i=1,m)
    4 format(A/(5E15.7))
      return
      end

      subroutine gradients(n,m,x,a,user,iuser)
      implicit double precision(a-h,o-z)
      dimension x(*),a(*),user(*),iuser(*)
      common/maxac/maxa
c     print *,'enter gradients' 
c  call CUTE's sparse gradient/Jacobian evaluation function
      call csgr(n,m,.false.,1,x,x,mxa,maxa,a,iuser,iuser(maxa+1))
c     print 4,'a(ij) =',(a(k),k=1,mxa)
      do i=1,n
        user(i)=a(maxa-n+i)
      enddo
      do i=maxa-n,1,-1
        a(n+i)=a(i)
      enddo
      do i=1,n
        a(i)=user(i)
      enddo
c     print 4,'new a(ij) =',(a(k),k=1,mxa)
    4 format(A/(5E15.7))
      return
      end
