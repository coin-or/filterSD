
christen this file schurQR.f
cut here >>>>>>>>>>>>>>>>>>>

c  Copyright  (C) 2011 Roger Fletcher
c  Current version dated 17 January 2012

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c***************** sparse matrix routines for manipulating L *******************

c        **********************************************************
c        Basis matrix routines for LCP solvers with sparse matrices
c        **********************************************************

c  These routines form and update L-Implicit-U factors LPB=U of a matrix B
c  whose columns are the normal vectors of the active constraints. In this
c  method only the unit lower triangular matrix L and the diagonal of U (in
c  addition to the row permutation P) is stored. B is represented in block form

c    | I  A_2 |    where the last m1 columns (A_2 and A_1) come from the
c    | 0  A_1 |    general constraint normals (columns of the matrix A in bqpd)

c  and the remaining unit columns come from simple bounds. The matrix A must be
c  specified in sparse format and the user is referred to the file  sparseA.f.

c  The data structure used for L is that of a profile or skyline scheme, in
c  which the nontrivial rows of L are stored as dense row spikes. The use of
c  a Tarjan+spk1 ordering algorithm to control the length of these spikes has
c  proved quite effective.

c  In schurQR.f, the factors are updated by the Schur complement method with
c  QR factors. This is based on the block factorization
c
c            | B  V | = | L  0 | | U  V |
c            | E  0 |   | S  I | | 0  C |
c
c  where V are columns of the constraint normals [I A] that have been added to
c  the active set, and E has unit rows in which the unit element marks
c  columns of B that have been removed from B. C=-E.inv(B).V is the
c  Schur complement matrix, which is independent of how L and U are defined,
c  and its QR factors are stored. The current dimension of C is stored in the
c  parameter ms of   common/refactorc/mc,mxmc   and the user must set a
c  maximum permitted value of mc in mxmc (mxmc <= n). The current basis matrix
c  is refactored if mc would exceed mxmc, or if issues of numerical stability
c  arise. Typically mxmc=25 is suitable.


c  Workspace
c  *********
c  The user needs to supply storage for the row spikes in the LIU data
c  structure of L, Also storage for matrices in the Schur complement scheme
c  is required. The amount of storage required is unknown a-priori.
c  Storage for schurQR.f is situated at the end of the workspace arrays ws
c  and lws in bqpd. Allow as much space for ws as you can afford: the routine
c  will report if there is not enough. So far 10^6 locations has proved
c  adequate for problems of up to 5000 variables.

c  The user is also allowed to reserve storage in ws and lws, for use in the
c  user-supplied routine gdotx. This storage is situated at the start of the
c  arrays ws and lws. The user specifies the amount required by
c  setting the parameters kk and ll in the common block
c     common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  Storage required by the LCP solver is also required: the amount is set by
c  the LCP solver in kkk and lll. The user MUST set mxws and mxlws to be
c  the total amount of real and integer workspace for the arrays ws and lws.

c  Other information
c  *****************

c  The methodology behind the L-Implicit-U factors and the row spike storage
c  scheme for L is described in the references
c    Fletcher R., Dense Factors of Sparse Matrices, in "Approximation Theory
c    and Optimization. Tributes to M.J.D. Powell", (M.D. Buhmann and A. Iserles,
c    eds), Cambridge University Press (1997), pp. 145-166.
c  and
c    Fletcher R., Block Triangular Orderings and Factors for Sparse Matrices
c    in LP, in "Numerical analysis 1997" (D.F. Griffiths, D.J. Higham and
c    G.A. Watson, eds.), Pitman Research Notes in Mathematics 380, (1998),
c    Longman, Harlow, pp. 91-110.

c  The file contains routines for solving systems with B or its transpose
c  which might be of use in association with bqpd. These routines are
c  documented below.

c  Steepest edge coefficients e(i) are also updated in these routines

      subroutine start_up(n,nm,nmi,a,la,nk,e,ls,aa,ll,mode,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      common/wsc/kk,ll_,kkk,lll,mxws,mxlws
      common/epsc/eps,tol,emin
      common/schurc/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1,
     *  nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv,lv1,le,le1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/mc,mxmc
      mxmc=min(n,mxmc)
c  set storage map for sparse factors
      ns=n
      ns1=ns+1
      nt=ns+n
      nt1=nt+1
      nu=nt+n
      nu1=nu+1
      nx=nu+n
      nx1=nx+1
      nq=nx+n
      nq1=nq+1
      nr=nq+mxmc**2
      nr1=nr+1
      ny=nr+mxmc*(mxmc+1)/2
      ny1=ny+1
      nz=ny+mxmc+1
      nz1=nz+1
      np=nz+mxmc+1
      np1=np+1
      nprof=mxws-kk-kkk-np
c     print *,'nprof =',nprof
      if(nprof.le.0)then
        write(nout,*)'not enough real workspace in ws'
        write(nout,*)'you give mxws as',mxws
        write(nout,*)'mxws must be much greater than',mxws-nprof
        ifail=7
        return
      endif
      lc=n
      lc1=lc+1
      li=lc+n
      li1=li+1
      lm=li+nmi
      lm1=lm+1
      lp=lm+n
      lp1=lp+1
      lq=lp+n
      lq1=lq+1
      lr=lq+n
      lr1=lr+1
      ls_=lr+n
      ls1=ls_+1
      lt=ls_+n
      lt1=lt+1
      lv=lt+n
      lv1=lv+1
      le=lv+mxmc+1
      le1=le+1
      lleft=mxlws-ll_-lll-le-mxmc-1
      if(lleft.lt.0)then
        write(nout,*)'not enough integer workspace in lws'
        write(nout,*)'you give mxlws as',mxlws
        write(nout,*)'minimum value for mxlws is',mxlws-lleft
        ifail=7
        return
      endif
      m=nm-n
      mp=-1
      mq=-1
c     write(nout,*)'ls',(ls(ij),ij=1,nk)
      if(mode.eq.3)then
        if(nk.lt.n)then
c  reset ls from e
          do j=1,nk
            i=-ls(j)
            if(i.gt.0)e(i)=-e(i)
          enddo
          j=0
          nk=nmi
          do i=1,nmi
            if(e(i).ne.0.D0)then
              j=j+1
              if(e(i).gt.0.D0)then
                ls(j)=i
              else
                ls(j)=-i
                e(i)=-e(i)
              endif
            else
              ls(nk)=i
              nk=nk-1
            endif
          enddo
          if(j.ne.n)then
            write(nout,*)'malfunction in reset sequence in start_up'
            stop
          endif
        endif
c  reset lr, lc, li, m1 and m2 from ls
        do i=li+n+1,li+nm
          ll(i)=0
        enddo
        m1=n
        m2=0
        do j=1,n
          i=abs(ls(j))
          if(i.gt.n)then
            ll(lc+m1)=i
            ll(li+i)=m1
            m1=m1-1
          else
            m2=m2+1
            lii=ll(li+i)
            lrm2=ll(m2)
            call iexch(ll(lii),ll(m2))
            call iexch(ll(li+i),ll(li+lrm2))
          endif
        enddo
        m1=n-m1
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,ifail)
        if(ifail.ge.1)then
c         write(nout,*)'failure in re_order (1)'
          if(ifail.eq.7)return
          mode=2
          goto1
        endif
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,aa,ifail)
        if(ifail.eq.7)return
        call check_L(n,aa,ll(lp1),ifail)
        if(ifail.eq.1)then
          mode=2
          goto1
        endif
        call EBspace(n,ll(lp1),ll(lq1),ll(ls1),ll,aa(np1),
     *    neb,nprof,ifail)
        if(ifail.gt.0)return
        neb=np+neb
        neb1=neb+1
        mc=0
        do i=1,m2
          ll(lm+i)=ll(i)
        enddo
        do i=m2+1,n
          ll(lm+i)=ll(lc+i)
        enddo
        return
      endif
    1 continue
      if(emin.eq.0.D0)then
c  set a lower bound on e(i): setting emin=0.D0 will force emin to be recalculated: do this only if mode<3
        emin=1.D0
        do i=1,nmi-n
          emin=max(emin,ailen(n,a,la,i))
        enddo
        emin=1.D0/emin
      endif
      do i=1,n
        ll(i)=i
        ll(li+i)=i
        e(i)=1.D0
      enddo
      do i=n+1,nm
        ll(li+i)=0
        e(i)=0.D0
      enddo
      nu_=0
      if(mode.ne.0)then
c  shift designated bounds to end and order the resulting rows and columns
        do j=1,nk
          i=abs(ls(j))
          if(i.le.n)then
            nn=n-nu_
            nu_=nu_+1
            call iexch(ls(nu_),ls(j))
            ii=ll(li+i)
            ll(ii)=ll(nn)
            ll(li+ll(ii))=ii
            ll(nn)=i
            ll(li+i)=nn
          endif
        enddo
        call order(n,nu_,nk,la,ll,ls,ll(li1),ll(lp1),ll(lq1),ll(lr1),
     *    aa(np1),nprof,ifail)
        if(ifail.gt.0)return
      endif
      call factor(n,nm,nu_,nk,a,la,e,ls,aa(ns1),aa(nt1),aa(nu1),
     *  aa(nx1),ll,ll(lc1),ll(li1),ll(lm1),ll(lp1),ll(lq1),ll(lr1),
     *  ll(ls1),aa(np1),nprof,aa,ifail)
      call EBspace(n,ll(lp1),ll(lq1),ll(ls1),ll,aa(np1),
     *  neb,nprof,ifail)
      if(ifail.gt.0)return
      neb=np+neb
      neb1=neb+1
      mc=0
      do i=1,m2
        ll(lm+i)=ll(i)
      enddo
      do i=m2+1,n
        ll(lm+i)=ll(lc+i)
      enddo
      if(ifail.gt.0)return
    3 format(A/(15I5))
    4 format(A/(5E15.7))
c     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
c     emax=0.D0
c     do i=1,nm
c       if(e(i).gt.0.D0)then
c         call eptsol(n,a,la,i,a,aa(nq1),aa(nr1),aa(neb1),aa(ny1),
c    *    aa(ns1),aa(nu1),aa(nx1),aa,aa(np1),
c    *    ll,ll(lc1),ll(li1),ll(lv1),ll(le1),ll(lp1),ll(lq1),ei)
c         emax=max(emax,abs(ei-e(i)))
c       endif
c     enddo
c     if(emax.ge.tol)
c    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine refactor(n,nm,a,la,aa,ll,ifail)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),aa(*),ll(*)
      common/schurc/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1,
     *  nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv,lv1,le,le1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/noutc/nout
c     write(nout,*)'refactor'
      ifail=1
      return
      end

      subroutine pivot(p,q,n,nm,a,la,e,aa,ll,ifail,npv)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),aa(*),ll(*)
      common/noutc/nout
      common/iprintc/iprint
      common/schurc/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1,
     *  nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv,lv1,le,le1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/mxm1c/mxm1
      common/refactorc/mc,mxmc
      common/epsc/eps,tol,emin
      common/pqc/pc,qr,lmp
c     write(nout,*)'pivot: p,q =',p,q
      call updateSE(p,q,n,nm,a,la,e,aa(nq1),aa(nr1),aa(neb1),
     *  aa(ny1),aa(nz1),aa(ns1),aa(nt1),aa(nu1),aa(nx1),aa,aa(np1),ll,
     *  ll(lc1),ll(li1),ll(lp1),ll(lq1),ll(lm1),ll(lv1),ll(le1),ifail)
      if(ifail.eq.1)return
      if(mc.eq.mxmc.and.pc.eq.0.and.qr.eq.0)then
c  reset permutations and refactorize L
        mc1=mc+1
        ll(le+mc1)=p
        ll(lv+mc1)=q
c       print 3,'le =',(ll(le+i),i=1,mc1)
c       print 3,'lv =',(ll(lv+i),i=1,mc1)
        do i=1,mc1
          p=ll(le+i)
          q=ll(lv+i)
          ip=ll(li+p)
          if(p.gt.n)then
            m2=m2+1
            qq=ll(lc+m2)
            ll(lc+ip)=qq
            ll(li+qq)=ip
            ll(li+p)=0
          else
            ll(ip)=ll(m2)
            ll(li+ll(ip))=ip
            ll(m2)=p
            ll(li+p)=m2
          endif
          if(q.gt.n)then
            ll(lc+m2)=q
            ll(li+q)=m2
            m2=m2-1
          else
            iq=ll(li+q)
            ll(iq)=ll(m2)
            ll(li+ll(iq))=iq
            ll(m2)=q
            ll(li+q)=m2
          endif
        enddo
c       print 3,'lr =',(ll(i),i=1,n)
c       print 3,'lc =',(ll(lc+i),i=m2+1,n)
c       print 3,'li =',(ll(li+i),i=1,nm)
        m1=n-m2
c       call checkperms(n,ll,ll(lc1),ll(li1))
c       mp=-1
c       mq=-1
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,ifail)
        if(ifail.ge.1)then
c         print *,'no traversal in re_order (3)'
          ifail=11
          return
          stop
        endif
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,aa,ifail)
        if(ifail.eq.7)return
        call EBspace(n,ll(lp1),ll(lq1),ll(ls1),ll,aa(np1),
     *    neb,nprof,ifail)
        if(ifail.gt.0)return
        neb=np+neb
        neb1=neb+1
        mc=0
        do i=1,m2
          ll(lm+i)=ll(i)
        enddo
        do i=m2+1,n
          ll(lm+i)=ll(lc+i)
        enddo
      else
        call updateQR(p,q,n,a,la,aa(nq1),aa(nr1),aa(neb1),aa(nx1),
     *    aa(ny1),aa(nz1),ll,ll(lc1),ll(li1),ll(lv1),ll(le1),ll(lm1),
     *    ifail)
        if(ifail.gt.0)return
      endif
      npv=npv+1
      mp=-1
      mq=-1
c     call check_L(n,aa,ll(lp1),ifail)
c     print 4,'e =',(e(i),i=1,nm)
c     print 3,'lm =',(ll(lm+i),i=1,n)
      return
c  check Steepest Edge coefficients
      emax=0.D0
      do j=1,n
        i=ll(lm+j)
        call eptsol(n,a,la,i,a,aa(nq1),aa(nr1),aa(neb1),aa(ny1),
     *  aa(ns1),aa(nu1),aa(nx1),aa,aa(np1),
     *  ll,ll(lc1),ll(li1),ll(lv1),ll(le1),ll(lp1),ll(lq1),ei)
        emax=max(emax,abs(ei-e(i)))
c       if(abs(ei-e(i)).gt.tol)then
c         print *,'error in steepest edge coefficient =',i,ei,e(i)
c         print 4,'s =',(aa(ns+i),i=1,n)
c         if(abs(ei-e(i)).gt.1.D-6)stop
c       endif
      enddo
      if(emax.gt.tol)then
        print 2,'max error in steepest edge coefficients =',emax
c       if(emax.gt.1.D-2)stop
      endif
      return
    2 format(A,5E15.7)
    3 format(A/(15I5))
    4 format(A/(5E15.7))
    5 format((5E15.7))
      end

      subroutine fbsub(n,jmin,jmax,a,la,q,b,x,ls,aa,ll,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
    9 format(A/(15I5))
      dimension a(*),la(*),b(*),x(*),ls(*),aa(*),ll(*)

c  solves a system  B.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   jmin,jmax  now redundant 
c   a,la   specification of QP problem data (as for bqpd)
c   q   an integer which, if in the range 1:n+m, specifies that the rhs vector
c       b is to be column q of the matrix A of general constraint normals.
c       In this case the parameter b is not referenced by fbsub.
c       If q=0 then b is taken as the vector given in the parameter b.
c   b(n)  must be set to the r.h.s. vector b in natural order (but only if q=0)
c   x(n+m)  contains the solution x, set according to the index number of that
c           component (in the range 1:n for a simple bound and  n+1:n+m
c           for a general constraint)
c   ls(*)  now redundant. Previously ls was an index vector, listing the
c       components of x that are required. Now all the solution x is provided,
c       set as described above.
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   save   now redundant

      common/noutc/nout
      common/schurc/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1,
     *  nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv,lv1,le,le1
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'fbsub  q =',q
      if(q.eq.0)then
        do i=1,n
          aa(nt+ll(li+i))=b(i)
        enddo
      endif
      call  aqsol(n,a,la,q,aa(nq1),aa(nr1),aa(neb1),aa(nz1),aa(nt1),
     *  aa(nu1),aa(nx1),aa,aa(np1),ll,ll(lc1),ll(li1),ll(lv1),
     *  ll(le1),ll(lp1),ll(lq1))
      do j=1,n
        x(ll(lm+j))=aa(nt+j)
      enddo
c     print *,'x =',(x(i),i=1,18)
      return
      end

      subroutine ztg(n,k,rg,lv,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rg(*),lv(*),aa(*),ll(*)
      common/schurc/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1,
     *  nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv_,lv1,le,le1
c     print *,'aa =',(aa(nu+i),i=1,18)
      do j=1,k
        rg(j)=aa(nu+ll(li+lv(j)))
      enddo
      return
      end

      subroutine tfbsub(n,a,la,p,b,x,aa,ll,ep,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),aa(*),ll(*)

c  solves a system  Bt.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   a,la   specification of QP problem data (as for bqpd)
c   p    an integer which, if in the range 1:n+m, specifies that the rhs vector
c        b is a unit vector appropriate to the position of p in the current
c        ordering. In this case b is not referenced by tfbsub.
c   b(n+m)  If p=0, this must be set to the r.h.s. vector b. Only the components
c        of b need be set, according to the index number of each component (in
c        the range 1:n for a simple bound and n+1:n+m for a general constraint)
c   x(n)  contains the solution x (in natural ordering)
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   ep  if p>0, ep contains the L2 norm of the solution
c   save  now redundant

      common/noutc/nout
      common/schurc/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1,
     *  nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv,lv1,le,le1
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'tfbsub  p =',p
      call eptsol(n,a,la,p,b,aa(nq1),aa(nr1),aa(neb1),aa(ny1),
     *  aa(ns1),aa(nu1),aa(nx1),aa,aa(np1),ll,ll(lc1),ll(li1),
     *  ll(lv1),ll(le1),ll(lp1),ll(lq1),ep)
      do i=1,n
        x(ll(i))=aa(ns+i)
      enddo
c     print 4,'x =',(x(i),i=1,n)
    4 format(A/(5E15.7))
      return
      end

      subroutine newg
      common/factorc/m1,m2,mp,mq,lastr,irow
      mq=-1
      return
      end

c******** The following routines are internal to schurQR.f **************

      subroutine check_L(n,d,p,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension d(*),p(*)
      common/noutc/nout
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/epsc/eps,tol,emin
      write(nout,*)'check_L'
      ifail=1
      dmin=1.D37
      do k=nu+1,n
        dmin=min(dmin,abs(d(k)))
c       if(abs(d(k)).le.tol)return
      enddo
      write(nout,*)'dmin =',dmin
c     len=0
c     do i=1,n
c       len=len+p(i)
c     enddo
c     write(nout,*)m1*(m1+1)/2,len+m1
c     write(nout,*)'m1 =',m1,'   file length =',len,'   total =',len+m1
      ifail=0
      return
      end

      subroutine aqsol(n,a,la,q,Q_,R,EB,z,t,u,x,d,ws,
     *  lr,lc,li,lv,le,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),Q_(*),R(*),EB(*),z(*),t(*),u(*),x(*),
     *  d(*),ws(*),lr(*),lc(*),li(*),lv(*),le(*),pp(*),qq(*)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/mc,mxmc
      common/pqc/pc,qr,lmp
c     print *,'aqsol  q =',q
      if(q.gt.0)then
c       print *,'q,n,li(q),m2',q,n,li(q),m2
        if(q.le.n.and.li(q).le.m2.or.q.gt.n.and.li(q).gt.0)then
c  q is in the starting active set (and hence in row qr of E)
          do qr=1,mc
            if(q.eq.le(qr))goto10
          enddo
          print *,'malfunction: q not in E'
          stop
   10     continue
        else
c  q is a new column
          qr=0
        endif
c  form t=Bk^{-1}.aq, else form t=Bk^{-1}.t
c  scatter a_q into t
        liq=li(q)
        do i=1,n
          t(i)=0.D0
        enddo
        if(q.le.n)then
          t(liq)=1.D0
        else
          call iscatter(a,la,q-n,li,t,n)
        endif
      endif
c     print 4,'t=',(t(i),i=1,n)
      if(mc.gt.0)then
c  form u=E.B^{-1}.t and possibly z=-u
        if(q.eq.0.or.q.gt.n.and.qr.eq.0)then
          i1=1
          do i=1,mc
c           print 4,'EB =',(EB(j),j=i1,i1+m1-1)
            u(i)=scpr(0.D0,EB(i1),t(m2+1),m1)
            if(le(i).le.n)u(i)=u(i)+t(li(le(i)))
            z(i)=-u(i)
            i1=i1+m1
          enddo
        elseif(qr.gt.0)then
          do i=1,mc
            u(i)=0.D0
          enddo
          u(qr)=1.D0
        else
c         print 4,'EB =',(EB(j),j=1,m1*mc)
          liq=liq-m2
          do i=1,mc
            u(i)=EB(liq)
            z(i)=-u(i)
            liq=liq+m1
          enddo
        endif
c       print 4,'u=',(u(i),i=1,mc)
c  form x=C^{-1}.u
        call Qtprod(mc,mxmc,Q_,u,x)
        mm=mc*(3-mc)/2+(mc-1)*mxmc
        call rsol(mc,mm,mxmc,R,x)
c       print 4,'x=',(x(i),i=1,mc)
c  accumulate t=t+V.x
        do i=1,mc
          lvi=lv(i)
          if(lvi.le.n)then
            t(li(lvi))=t(li(lvi))+x(i)
          else
            call isaipy(x(i),a,la,lvi-n,t,n,lr,li)
          endif
        enddo
      endif
c     print 4,'t in natural order =',(t(li(i)),i=1,n)
c  finally t:=B^{-1}.t-E'.x
      call aqsol0(n,a,la,0,t,u,d,ws,lr,lc,li,pp,qq)
      do i=1,mc
        lei=li(le(i))
        t(lei)=t(lei)-x(i)
      enddo
c     print 4,'t in column order',(t(i),i=1,n)
      mq=q
      return
    3 format(A/(15I5))
    4 format(A/(5E15.7))
      end

      subroutine aqsol0(n,a,la,q,tn,xn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),tn(*),xn(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*)
      common/factorc/m1,m2,mp,mq,lastr,irow
      if(q.gt.0)then
        do i=1,n
          tn(i)=0.D0
        enddo
        if(q.le.n)then
          tn(li(q))=1.D0
        else
          call iscatter(a,la,q-n,li,tn,n)
        endif
      endif
c     print *,'tn =',(tn(i),i=1,n)
      do i=n,m2+1,-1
        ir=lr(i)
        pri=pp(ir)
        if(pri.eq.0)then
          xn(i)=tn(i)/d(i)
        else
          xn(i)=(scpr(tn(i),ws(qq(ir)+1),tn(i-pri),pri))/d(i)
        endif
        call isaipy(-xn(i),a,la,lc(i)-n,tn,n,lr,li)
      enddo
      do i=m2+1,n
        tn(i)=xn(i)
      enddo
c     print *,'tn =',(tn(i),i=1,n)
      return
      end

      subroutine eptsol(n,a,la,p,b,Q_,R,EB,y,s,u,x,d,ws,
     *  lr,lc,li,lv,le,pp,qq,ep)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),Q_(*),R(*),EB(*),y(*),s(*),u(*),x(*),
     *  d(*),ws(*),lr(*),lc(*),li(*),lv(*),le(*),pp(*),qq(*)
      common/epsc/eps,tol,emin
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/mc,mxmc
      common/pqc/pc,qr,lmp
c     print *,'eptsol  p =',p
c  column ordering is that defined by Bk = B + (V-B.E').E
c  row order is same as for B
      if(p.eq.0)then
c       print 3,'lr =',(lr(i),i=1,m2)
c       print 3,'lc =',(lc(i),i=m2+1,n)
c       print 3,'le =',(le(i),i=1,mc)
c       print 3,'lv =',(lv(i),i=1,mc)
        do i=1,mc
          x(i)=b(le(i))
          b(le(i))=0.D0
        enddo
        call eptsol0(n,a,la,0,b,s,d,ws,lr,lc,li,pp,qq)
        do i=1,mc
          b(le(i))=x(i)
          if(lv(i).le.n)then
            x(i)=s(li(lv(i)))-b(lv(i))
          else
            x(i)=aiscpri(n,a,la,lv(i)-n,s,-b(lv(i)),lr,li)
          endif
        enddo
      else
        if(p.le.n.and.li(p).le.m2.or.p.gt.n.and.li(p).gt.0)then
c  p is in the starting active set (set pc=0)
          pc=0
          lmp=li(p)
        else
c  p is in V  (pc indicates where p is in V)
          do pc=1,mc
            if(p.eq.lv(pc))goto10
          enddo
        print 1,'p,pc,li(p),m1,m2 =',p,pc,li(p),m1,m2
        print 3,'le =',(le(i),i=1,mc)
        print 3,'lv =',(lv(i),i=1,mc)
          print *,'malfunction: p not in V'
          stop
   10     continue
          lmp=li(le(pc))
        endif
c  get s=Bk^{-T}.ep
        if(pc.eq.0)then
          call eptsol0(n,a,la,p,a,s,d,ws,lr,lc,li,pp,qq)
c       print 4,'s0 ordered by lr',(s(i),i=1,n)
c       print 4,'s0 in natural order',(s(li(i)),i=1,n)
          m1mc=m1*mc
          do i=1,m1
            EB(m1mc+i)=s(m2+i)
          enddo
c         print 1,'eptsol: p =',p
c         print 4,'EB is',(s(i),i=m2+1,n)
c         print 4,'EB is',(EB(i),i=1,m1mc+m1)
c  form s=B^{-T}.ep and then y=-V'.s
          do i=1,mc
            if(lv(i).le.n)then
              x(i)=s(li(lv(i)))
            else
              x(i)=aiscpri(n,a,la,lv(i)-n,s,0.D0,lr,li)
            endif
            y(i)=-x(i)
          enddo
        else
c  this is pc>0: set s=0 and x=-e_pc
          do i=1,n
            s(i)=0.D0
          enddo
          do i=1,mc
            x(i)=0.D0
          enddo
          x(pc)=-1.D0
        endif
      endif
c     print 4,'x =',(x(i),i=1,mc)
      if(mc.gt.0)then
c  form u=C^{-T}.x and accumulate EB'.u into s
        call rtsol(mc,mm,mxmc,R,x)
        call Qprod(mc,mxmc,Q_,x,u)
        i1=1
        do i=1,mc
          call mysaxpy(u(i),EB(i1),s(m2+1),m1)
          if(le(i).le.n)s(li(le(i)))=u(i)
          i1=i1+m1
        enddo
      endif
c     print 4,'sk in natural order=',(s(li(i)),i=1,n)
c     print 4,'s =',(s(i),i=1,n)
      mp=p
      if(p.gt.0)ep=xlen(0.D0,s,n)
      return
    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(15I5))
    4 format(A/(5E15.7))
      end

      subroutine eptsol0(n,a,la,p,b,sn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),sn(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*)
      common/epsc/eps,tol,emin
      common/factorc/m1,m2,mp,mq,lastr,irow
      if(p.eq.0)then
        do i=1,m2
          sn(i)=b(lr(i))
        enddo
        do i=m2+1,n
          sn(i)=0.D0
        enddo
        do i=m2+1,n
          j=lc(i)
          sn(i)=-aiscpri(n,a,la,j-n,sn,-b(j),lr,li)/d(i)
          ir=lr(i)
          pri=pp(ir)
          if(pri.gt.0)call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
        enddo
      else
        do i=1,n
          sn(i)=0.D0
        enddo
        pr=li(p)
        if(p.le.n)then
          if(pr.gt.m2)goto1
          sn(pr)=1.D0
          do i=m2+1,n
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/d(i)
            ir=lr(i)
            pri=pp(ir)
            if(pri.gt.0)call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
          enddo
        else
          if(pr.le.m2)goto1
          do i=m2+1,n
            bi=0.D0
            if(i.eq.pr)bi=-1.D0
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,bi,lr,li)/d(i)
            ir=lr(i)
            pri=pp(ir)
            if(pri.gt.0)call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
          enddo
        endif
      endif
c     print *,'sn =',(sn(i),i=1,n)
      return
    1 continue
      print *,'malfunction detected in eptsol0: p =',p
      stop
      end

      subroutine order(n,nu,nc,la,lr,ls,li,p,q,r,ws,mxws,ifail)
      implicit integer (c-t)
      double precision ws
      dimension la(0:*),lr(*),ls(*),li(*),p(*),q(*),r(*),ws(*)
      common/noutc/nout
c     character star(1000,80)
c     write(nout,*)'order'
c  spk1 ordering on full matrix
      ifail=0
      if(nu.eq.n)return
c  set row and column counts and row-wise data structure
      nn=n-nu
      ii=mxws/nn
      do j=1,nn
        rowj=lr(j)
        p(rowj)=(j-1)*ii
        r(rowj)=0
      enddo
      do j=nn+1,n
        r(lr(j))=0
      enddo
    1 continue
      do i=nu+1,nc
        coli=abs(ls(i))
        li(coli)=0
        jp=la(0)+coli-n
        do j=la(jp),la(jp+1)-1
          rowj=la(j)
          if(li(rowj).le.nn)then
            li(coli)=li(coli)+1
            r(rowj)=r(rowj)+1
            ij=p(rowj)+r(rowj)
            if(ij.gt.mxws)then
              ij=mxws
              ifail=1
            endif
            ws(ij)=dble(coli)
          endif
        enddo
      enddo
c  check for no overlaps
      qrj=0
      do j=1,nn
        rowj=lr(j)
        if(p(rowj).lt.qrj)ifail=1
        qrj=p(rowj)+r(rowj)
        q(rowj)=qrj
        p(rowj)=p(rowj)+1
      enddo
      if(ifail.eq.1.or.qrj.gt.mxws)then
        qrj=0
        do j=1,nn
          rowj=lr(j)
          p(rowj)=qrj
          qrj=qrj+r(rowj)
          r(rowj)=0
        enddo
        if(qrj.gt.mxws)then
          write(nout,*)'not enough space for ws in order:  mxws =',mxws
          ifail=7
          return
        endif
        ifail=0
        goto1
      endif
      ifirstc=nu+1
      ifirstr=1
    2 continue
c  move zero-column-count columns to lhs and find minimum column count
      mcc=n
      do i=ifirstc,nc
        coli=abs(ls(i))
        if(li(coli).eq.0)then
          call iexch(ls(i),ls(ifirstc))
          li(coli)=ifirstr-1
          ifirstc=ifirstc+1
        else
          mcc=min(mcc,li(coli))
        endif
      enddo
c     write(nout,*)'ifirstc,ifirstr,mcc',ifirstc,ifirstr,mcc
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'ls =',(ls(i),i=nu+1,nc)
c     write(nout,*)'row counts =',(r(lr(j)),j=1,n)
c     write(nout,*)'column counts =',(li(abs(ls(i))),i=nu+1,nc)
      if(ifirstc.gt.nc)goto4
c  apply tie-break rule
      tie=0
      do i=ifirstc,nc
        coli=abs(ls(i))
        if(li(coli).eq.mcc)then
          ti=0
          jp=la(0)+coli-n
          do j=la(jp),la(jp+1)-1
            rowj=la(j)
            if(li(rowj).ge.ifirstr)ti=ti+r(rowj)
          enddo
          if(ti.gt.tie)then
            tie=ti
            mccc=coli
          endif
        endif
      enddo
c     write(nout,*)'tie,mccc',tie,mccc
c  permute rows of m-c-c column to top and update column counts
      jp=la(0)+mccc-n
      do j=la(jp),la(jp+1)-1
        rowj=la(j)
        jr=li(rowj)
        if(jr.lt.ifirstr)goto3
        if(jr.gt.nn)goto3
        lr(jr)=lr(ifirstr)
        li(lr(jr))=jr
        lr(ifirstr)=rowj
        li(rowj)=ifirstr
        ifirstr=ifirstr+1
        do i=p(rowj),q(rowj)
          coli=int(ws(i))
          li(coli)=li(coli)-1
        enddo
    3   continue
      enddo
      goto2
    4 continue
c  print star diagram
c     if(nc-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'spk1 ordering'
c     ij=li(abs(ls(nc)))
c     do i=1,ij
c       do j=1,nc-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,nc-nu
c       jp=la(0)+abs(ls(nu+j))-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=1,ij
c       write(nout,*)(star(i,j),j=1,nc-nu)
c     enddo
c     write(nout,*)'lr =',(lr(i),i=1,n)
c     write(nout,*)'ls =',(ls(i),i=nu+1,nc)
c     write(nout,*)'lower profile =',(li(abs(ls(i))),i=nu+1,nc)
      return
      end

      subroutine factor(n,nm,nu,nc,a,la,e,ls,sn,tn,un,xn,lr,lc,li,
     *  mao,p,q,r,s,ws,mxws,d,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      integer coli,r,s,rowi,rowp,tl,tu
      dimension a(*),la(0:*),e(*),ls(*),sn(*),tn(*),un(*),xn(*),
     *  lr(*),lc(*),li(*),mao(*),p(*),q(*),r(*),s(*),ws(*),d(*)
c     character star(1000,80)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/iprintc/iprint
      common/refactorc/mc,mxmc
      common/epsc/eps,tol,emin
      common/noutc/nout
      parameter (thresh=1.D-1)
c  factorize LPA=U when A is rectangular
c    p(row) stores the number of stored elements of a natural row
c    q(row) stores the base address in ws of a natural row
c    r(row) stores the previous row stored in ws (or 0 if the first row in ws)
c    s(row) stores the next row stored in ws (or 0 if the last row in ws)
c    li(n+*) stores the lower profile of the sparse matrix
c    irow stores the natural row number of the initial row stored in ws
c    lastr stores the natural row number of the previous row put into ws
c     write(nout,*)'factor'
      nup=0
      lastr=0
      irow=0
      do i=1,n
        p(i)=0
      enddo
      m1=0
      tl=1
      do ii=nu+1,nc
        coli=abs(ls(ii))
c       write(nout,*)'coli =',coli
        tu=li(coli)
        do i=1,n
          tn(i)=0.D0
        enddo
        call iscatter(a,la,coli-n,li,tn,n)
        do i=m1,1,-1
          rowi=lr(i)
          pri=p(rowi)
          if(pri.eq.0)then
            xn(i)=tn(i)/d(i)
          else
            xn(i)=(scpr(tn(i),ws(q(rowi)+1),tn(i-pri),pri))/d(i)
          endif
          call isaipy(-xn(i),a,la,lc(i)-n,tn,n,lr,li)
        enddo
        do i=1,m1
          tn(i)=xn(i)
        enddo
        m1p=m1+1
c       write(nout,*)'lr =',(lr(i),i=1,n)
c       write(nout,*)'tn =',(tn(i),i=1,tu)
c  threshold pivot selection
        call linf(tu-m1,tn(m1p),z,iz)
        if(z.le.tol)then
          li(coli)=0
          goto2
        endif
        zz=max(tol,z*thresh)
        do i=tl,tu
          q(lr(i))=m1p
        enddo
c       write(nout,*)'q =',(q(lr(i)),i=m1p,tu)
        iz=iz+m1
        if(iz.lt.tl)then
          z=0.D0
          qri=m1p
          do j=m1p,tu
            tnj=abs(tn(j))
            if(tnj.ge.zz)then
              qrj=q(lr(j))
              if(qrj.eq.qri)then
                if(tnj.gt.z)then
                  z=tnj
                  iz=j
                endif
              elseif(qrj.gt.qri)then
                z=tnj
                iz=j
                qri=qrj
              endif
            endif
          enddo
        endif
        tl=tu+1
c       write(nout,*)'zz,z,iz,m1,qri',zz,z,iz,m1,qri
        if(iz.gt.m1p)then
          call rexch(tn(m1p),tn(iz))
          call iexch(lr(m1p),lr(iz))
          li(lr(m1p))=m1p
          li(lr(iz))=iz
        endif
        rowp=lr(m1p)
c  reset q values
        qrp=q(rowp)
        do i=m1p+1,tu
          if(abs(tn(i)).gt.tol)then
            rowi=lr(i)
            if(qrp.lt.q(rowi))q(rowi)=qrp
          endif
        enddo
        tnp=tn(m1p)
        do i=1,n
          sn(i)=0.D0
        enddo
        sn(m1p)=1.D0
        do i=1,m1
          sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/d(i)
          rowi=lr(i)
          pri=p(rowi)
          if(pri.gt.0)call mysaxpy(sn(i),ws(q(rowi)+1),sn(i-pri),pri)
        enddo
c       write(nout,*)'sn =',(sn(i),i=1,m1)
c  update steepest edge coefficients
        ep=e(rowp)
        e(rowp)=0.D0
        eq=2.D0/ep
        do i=1,n
          un(i)=eq*sn(i)
        enddo
        do i=m1,1,-1
          rowi=lr(i)
          pri=p(rowi)
          if(pri.eq.0)then
            xn(i)=un(i)/d(i)
          else
            xn(i)=(scpr(un(i),ws(q(rowi)+1),un(i-pri),pri))/d(i)
          endif
          call isaipy(-xn(i),a,la,lc(i)-n,un,n,lr,li)
        enddo
        do i=1,m1
          un(i)=xn(i)
        enddo
c       write(nout,*)'un =',(un(i),i=1,n)
        eq=ep/tnp
        do i=1,nm
          if(e(i).gt.0.D0)then
            j=li(i)
            ei=e(i)
            wi=tn(j)*eq
            awi=abs(wi)
            if(ei.ge.awi)then
              wi=wi/ei
              e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-un(j)/ei))))
            else
              wi=ei/wi
              e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-un(j)/ei))))
            endif
          endif
        enddo
        e(coli)=max(emin,abs(eq))
        do j=qrp,m1
          if(abs(sn(j)).gt.tol)goto1
        enddo
        j=m1p
    1   continue
        pri=m1p-j
        if(pri.gt.0)then
          call newslot(rowp,pri,lastr,irow,p,q,r,s,ws,mxws,i,ifail)
          if(ifail.gt.0)return
          p(rowp)=pri
          i=q(rowp)
          do j=j,m1
            i=i+1
            ws(i)=sn(j)
          enddo
        endif
        m1=m1p
        ls(m1)=ls(ii)
        lc(m1)=coli
        li(coli)=m1
        d(m1)=tnp
    2   continue
      enddo
c  complete ls and reorder lr, lc and d
      do i=m1+1,n
        ls(i)=lr(i)
      enddo
      j=n
      do i=1,nm
        if(e(i).eq.0.D0)then
          j=j+1
          ls(j)=i
        endif
      enddo
      m2=n-m1
      do i=n,m2+1,-1
        lc(i)=lc(i-m2)
        li(lc(i))=i
        lr(i)=lr(i-m2)
        li(lr(i))=i
        d(i)=d(i-m2)
      enddo
      do i=1,m2
        lr(i)=ls(m1+i)
        li(lr(i))=i
      enddo
c  reset mao
      ilast=n
      ii=ilast
      do i=ilast,m2+1,-1
        mao(i)=ilast
        ii=min(ii,i-p(lr(i)))
        if(ii.eq.i)ilast=i-1
      enddo
c     write(nout,*)'PAQ factors:  m1 =',m1
c     write(nout,4)'d =',(d(ij),ij=m2+1,n)
c     do ij=m2+1,n
c       rowp=lr(ij)
c       if(p(rowp).ne.0)then
c         write(nout,*)'L(',rowp,')',
c    *      (ws(k),k=q(rowp)+1,q(rowp)+p(rowp))
c       endif
c     enddo
c  print star diagram
c     write(nout,*)'factored ordering:  m1 =',m1
c     if(m1.gt.80.or.n.gt.1000)stop
c     do i=1,n
c       do j=1,m1
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,m1
c       jp=la(0)+lc(m2+j)-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=m2+1,n
c       write(nout,*)(star(i,j),j=1,m1)
c     enddo
c     write(nout,9)'ls =',(ls(j),j=1,n)
c     write(nout,*)'s.e. coeffs =',(e(i),i=1,nm)
c     write(nout,9)'lr =',(lr(j),j=1,n)
c     write(nout,9)'lc =',(lc(j),j=m2+1,n)
c     write(nout,9)'li =',(li(j),j=1,nm)
c     write(nout,9)'mao =',(mao(j),j=m2+1,n)
c     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
    4 format(A/(5E15.7))
    9 format(A/(15I5))
      return
      end

      subroutine re_order(n,nm,a,la,point,lr,lc,li,mao,p,q,r,s,
     *  t,ws,mxws,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension a(*),la(*),point(0:*),lr(*),lc(*),li(*),mao(*),
     *  p(*),q(*),r(*),s(*),t(*),ws(*)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/noutc/nout
      logical backtrack
c     character star(1000,80)
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'initial ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo
c     write(nout,*)'re_order'
      if(nu.eq.n)then
        ifail=0
        return
      endif
      m=nm-n
c  transversal search
      do iq=nu+1,n
        backtrack=.false.
        istack=nu
        inode=iq
        nodec=lc(inode)
        nodec_n=nodec-n
        lap=point(nodec_n+1)-point(nodec_n)
c       write(nout,*)'column node =',nodec,'  look-ahead rows =',
c    *    (la(j),j=point(nodec_n),point(nodec_n)+lap-1)
c  look-ahead loop
    1   continue
          lap=lap-1
          nextr=la(point(nodec_n)+lap)
          inext=li(nextr)
          if(inext.ge.iq)goto4
          if(lap.gt.0)goto1
          li(nodec)=0
    2   continue
c  reassignment depth first search
        t(inode)=point(nodec_n+1)-point(nodec_n)
c       write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *    (la(j),j=point(nodec_n),point(nodec_n)+t(inode)-1)
    3   continue
c  examine successor nodes
        if(t(inode).eq.0)then
          if(istack.eq.nu)then
            ifail=1
c           ifail=iq
c           write(nout,*)'exit: ifail =',iq
c           print *,'lc(iq) =',lc(iq)
            return
          endif
          istack=istack-1
          backtrack=.true.
          if(istack.eq.nu)then
            inode=iq
          else
            inode=mao(istack)
          endif
c         write(nout,*)'backtrack to node at address =',inode
          nodec=lc(inode)
          nodec_n=nodec-n
c         write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *      (la(j),j=point(nodec_n),point(nodec_n)+t(inode)-1)
          goto3
        endif
        t(inode)=t(inode)-1
        nextr=la(point(nodec_n)+t(inode))
        inext=li(nextr)
        if(inext.le.nu)goto3
        if(t(inext).ge.0)goto3
c  extend depth first search
c       write(nout,*)'nextr,inext',nextr,inext
        inode=inext
c       write(nout,*)'put node address on stack'
        istack=istack+1
        mao(istack)=inode
c       write(nout,*)'stack =',(mao(j),j=nu+1,istack)
        nodec=lc(inode)
        nodec_n=nodec-n
        lap=li(nodec)
        if(lap.eq.0)goto2
c       write(nout,*)'column node =',nodec,'  look-ahead rows =',
c    *    (la(j),j=point(nodec_n),point(nodec_n)+lap-1)
        goto1
    4   continue
c       write(nout,*)'new assignment found in row',nextr
c       write(nout,*)'istack,inext,nextr',istack,inext,nextr
c       if(istack.gt.nu)write(nout,*)'stack =',(mao(j),j=nu+1,istack)
        li(nodec)=lap
c  perform row permutation
        lr(inext)=lr(iq)
        li(lr(inext))=inext
        inode=iq
        do i=nu+1,istack
          inext=mao(i)
          lr(inode)=lr(inext)
          li(lr(inode))=inode
          inode=inext
        enddo
        lr(inode)=nextr
        li(nextr)=inode
c       write(nout,*)'lr =',(lr(j),j=nu+1,n)
c       write(nout,*)'look-ahead lengths =',(li(lc(j)),j=nu+1,iq)
        t(iq)=-1
        if(backtrack.or.istack.gt.nu+1)then
          do i=nu+1,iq-1
            t(i)=-1
          enddo
        endif
        do i=1,n
          if(li(i).gt.n)then
            write(nout,*)'iq =',iq
            stop
          endif
        enddo
      enddo
c     write(nout,*)'transversal found'
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'transversal ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo

c  tarjan ordering
      do i=1,n
        q(i)=0
        r(i)=0
      enddo
c  reset li and pair off columns with rows
      do i=nu+1,n
        nodec=lc(i)
        li(nodec)=i
        t(lr(i))=nodec
        s(i)=0
      enddo
      do i=nu+1,n
        noder=lr(i)
        nodec=t(noder)
        lc(noder)=point(nodec-n+1)-point(nodec-n)
        li(nodec)=-1
      enddo
      ifath=nu
      istack=n+1
c  tarjan loop
   10 continue
        istack=istack-1
        inode=istack
        noder=lr(inode)
        if(lc(noder).eq.0)then
          write(nout,*)'malfunction: zero length'
          stop
        endif
        nodec=t(noder)
   11   continue
        li(nodec)=lc(noder)
        mao(inode)=istack
c       write(nout,*)'put new node',noder,' on stack'
c       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
c       write(nout,*)'ifath,istack =',ifath,istack
c       write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *    (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
   12   continue
          if(li(nodec).eq.0)then
c           write(nout,*)'backtrack to previous nodes'
   13       continue
              if(inode.eq.n)goto14
              inext=inode+1
              nextr=lr(inext)
              if(mao(inode).lt.mao(inext))goto14
              inode=inext
              noder=nextr
              nodec=t(noder)
              if(li(nodec).eq.0)goto13
c           write(nout,*)'stack =',(lr(j),j=istack,n)
c           write(nout,*)'lengths =',(li(t(lr(j))),j=istack,n)
c           write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *        (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
            goto12
          endif
c  examine successors of current node
          li(nodec)=li(nodec)-1
          nextr=la(point(nodec-n)+li(nodec))
          inext=li(nextr)
          if(inext.le.ifath)goto12
          q(nextr)=q(nextr)+1
          nextc=t(nextr)
c         write(nout,*)'nextc,nextr,inext',nextc,nextr,inext
          if(li(nextc).ge.0)then
            mx=mao(inext)
            if(mao(inode).ge.mx)goto12
            do j=istack,n
              if(mao(j).eq.mx)goto12
              mao(j)=mx
            enddo
            write(nout,*)'malfunction'
            stop
          endif
          nodec=nextc
          noder=nextr
          istack=istack-1
          inode=istack
          lr(inext)=lr(inode)
          li(lr(inext))=inext
          lr(inode)=noder
          li(noder)=inode
          goto11
   14   continue
c       write(nout,*)'strong component identified'
c       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
c       write(nout,*)'ifath,istack,inode =',ifath,istack,inode,n
c  shift forward strong component
        inext=istack-1
        ir=inode-inext
        do j=istack,inode
          mao(j)=lr(j)
        enddo
        do j=inext+ir,ifath+1+ir,-1
          lr(j)=lr(j-ir)
          li(lr(j))=j
        enddo
        mx=ifath+ir
        iq=inext-ifath
        ifath=ifath+1
        do j=ifath,mx
          lr(j)=mao(j+iq)
          li(lr(j))=j
          mao(j)=mx
        enddo
        istack=inode+1
        ifath=mx
c       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
c       write(nout,*)'ifath,istack =',ifath,istack
        if(istack.le.n)then
          inode=istack
          noder=lr(inode)
          nodec=t(noder)
          nodec_n=nodec-n
c         write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *      (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
          goto12
        endif
      if(ifath.lt.n)goto10
c  end of tarjan process
c  reset lc and li
      do i=nu+1,n
        lc(i)=t(lr(i))
        li(lc(i))=i
      enddo
c     write(nout,*)'mao =',(mao(j),j=nu+1,n)
c     write(nout,*)'q =',(q(j),j=1,n)
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c     write(nout,*)'li =',(li(j),j=1,n+m)
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'tarjan ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo
c  set up pointers for row-wise sparse structure
      p(1)=1
      do i=1,n-1
        p(i+1)=p(i)+q(i)
        q(i)=p(i)-1
      enddo
      if(p(n)+q(n).gt.mxws)then
        print *,'not enough space for ws in re_order'
        ifail=7
        return
      endif
      q(n)=p(n)-1
      i=nu+1
   20 continue
      if(i.eq.mao(i))then
        t(i)=i
      else
c  spk1 ordering on tarjan block
c  set row and column counts
        do inode=i,mao(i)
          nodec=lc(inode)
          do j=point(nodec-n),point(nodec-n+1)-1
            noder=la(j)
            if(li(noder).ge.i)then
              q(noder)=q(noder)+1
              ws(q(noder))=dble(nodec)
              s(inode)=s(inode)+1
            endif
          enddo
        enddo
c       print *,'r-c counts: i =',i,'   mao(i) =',mao(i)
c       print *,'q =',(q(j),j=i,mao(i))
c       print *,'s =',(s(j),j=i,mao(i))
c  find minimum-column-count column
        mcc=n
        do inode=i,mao(i)
          noder=lr(inode)
          r(noder)=q(noder)-p(noder)+1
          mcc=min(mcc,s(inode))
        enddo
c     write(nout,*)'i,mao(i),mcc',i,mao(i),mcc
c     write(nout,*)'p =',(p(lr(j)),j=i,mao(i))
c     write(nout,*)'q =',(q(lr(j)),j=i,mao(i))
c     write(nout,*)'r =',(r(lr(j)),j=i,mao(i))
c     write(nout,*)'s =',(s(j),j=i,mao(i))
c  check for fully dense block
        if(mcc.gt.mao(i)-i)then
          do inode=i,mao(i)
            t(inode)=mao(i)
          enddo
          goto22
        endif
c  determine spk1 ordering
        ifirstr=i
        ifirstc=i
   21   continue
c  apply tie-break rule
        tie=0
        do inode=ifirstc,mao(i)
          if(s(inode).eq.mcc)then
            nodec=lc(inode)-n
            ti=0
            do j=point(nodec),point(nodec+1)-1
              noder=la(j)
              if(li(noder).ge.ifirstr)ti=ti+r(noder)
            enddo
            if(ti.gt.tie)then
              tie=ti
              mccc=nodec
            endif
          endif
        enddo
c       write(nout,*)'tie,mccc',tie,mccc+n
c  permute rows of m-c-c column to top and update column counts
        do j=point(mccc),point(mccc+1)-1
          noder=la(j)
          ir=li(noder)
          if(ir.ge.ifirstr)then
            lr(ir)=lr(ifirstr)
            li(lr(ir))=ir
            lr(ifirstr)=noder
            li(noder)=ifirstr
            ifirstr=ifirstr+1
            do ir=p(noder),q(noder)
              inode=li(int(ws(ir)))
              s(inode)=s(inode)-1
            enddo
          endif
        enddo
c       write(nout,*)'s =',(s(ij),ij=i,mao(i))
c       write(nout,*)'lr =',(lr(ij),ij=i,mao(i))
c  move zero-column-count columns to lhs and find minimum column count
        mcc=n
        do inode=ifirstc,mao(i)
          if(s(inode).eq.0)then
            nodec=lc(inode)
            lc(inode)=lc(ifirstc)
            li(lc(inode))=inode
            lc(ifirstc)=nodec
            li(nodec)=ifirstc
            s(inode)=s(ifirstc)
            t(ifirstc)=ifirstr-1
            ifirstc=ifirstc+1
          else
            mcc=min(mcc,s(inode))
          endif
        enddo
c       write(nout,*)'lc =',(lc(ij),ij=i,mao(i))
c       write(nout,*)'ifirstc,mcc',ifirstc,mcc
        if(ifirstc.lt.mao(i))goto21
      endif
   22 continue
      i=mao(i)+1
      if(i.le.n)goto20
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'tarjan + spk1 ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo
c     write(nout,*)'lr =',(lr(j),j=nu+1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c     write(nout,*)'lower profile =',(t(j),j=nu+1,n)
      ifail=0
      return
      end

      subroutine re_factor(n,nm,a,la,lr,lc,li,mao,p,q,r,s,
     *  t,ws,mxws,d,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension a(*),la(0:*),lr(*),lc(*),li(*),mao(*),
     *  p(*),q(*),r(*),s(*),t(*),d(*),ws(*)
c     character star(1000,80)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/iprintc/iprint
      common/refactorc/mc,mxmc
      common/epsc/eps,tol,emin
      common/noutc/nout
      double precision thresh,tol
      parameter (thresh=1.D-1)
c  factorize LPA=U
c    p(row) stores the number of stored elements of a natural row
c    q(row) stores the base address in ws of a natural row
c    r(row) stores the previous row stored in ws (or 0 if the first row in ws)
c    s(row) stores the next row stored in ws (or 0 if the last row in ws)
c    t(*) stores the lower profile of the sparse matrix
c    irow stores the natural row number of the initial row stored in ws
c    lastr stores the natural row number of the previous row put into ws
c     write(nout,*)'re_factor'
      nup=0
      m=nm-n
      lastr=0
      irow=0
      do i=1,n
        p(i)=0
      enddo
      if(m1.eq.0)return
      i=nu+1
    1 continue
      if(i.eq.mao(i))then
        d(i)=aij(lr(i),lc(i)-n,a,la)
        if(d(i).eq.0.D0)d(i)=eps
c       write(nout,*)'row,col,d(i) =',lr(i),lc(i),d(i)
      else
c       write(nout,*)'lc =',(lc(j),j=i,mao(i))
        do inode=i,mao(i)-1
          nodec=lc(inode)-n
          im=inode-1
c  form L.a_q
          z=0.
c         write(nout,*)'inode,t(inode)',inode,t(inode)
          do j=inode,t(inode)
            rowj=lr(j)
            prj=p(rowj)
            if(prj.gt.0)then
              d(j)=aiscpri2(n,a,la,rowj,nodec,ws(q(rowj)+1),1.D0,im,
     *          prj,li)
            else
              d(j)=aij(rowj,nodec,a,la)
            endif
            z=max(z,abs(d(j)))
          enddo
c         write(nout,*)'d =',(d(ij),ij=inode,t(inode))
c  threshold pivot selection
          zz=z*thresh
          z=0.D0
          pri=n
          do j=inode,t(inode)
            dj=abs(d(j))
            if(dj.ge.zz)then
              prj=p(lr(j))
              if(prj.eq.pri)then
                if(dj.gt.z)then
                  z=dj
                  iz=j
                endif
              elseif(prj.lt.pri)then
                z=dj
                iz=j
                pri=prj
              endif
            endif
          enddo
c       write(nout,*)'zz,z,iz,pri',zz,z,iz,pri
          if(iz.gt.inode)then
c  pivot interchange
            call rexch(d(inode),d(iz))
            call iexch(lr(inode),lr(iz))
            li(lr(iz))=iz
            li(lr(inode))=inode
          endif
          if(d(inode).eq.0.D0)d(inode)=eps
c  update L
          qri=q(lr(inode))
          zz=-d(inode)
          do j=inode+1,t(inode)
            z=d(j)/zz
            rowj=lr(j)
            prj=p(rowj)
            qrj=q(rowj)
c  find space available in-situ in ws
            if(prj.eq.0)then
              len=0
            elseif(s(rowj).eq.0)then
              len=mxws-qrj
            else
              len=q(s(rowj))-qrj
            endif
            if(abs(z).le.tol)then
c  special case of a zero multiplier
              if(prj.eq.0)goto2
              len_=prj+1
              if(len_.gt.len)then
                call newslot(rowj,len_,lastr,irow,p,q,r,s,ws,mxws,qrj,
     *            ifail)
                if(ifail.gt.0)return
                qrj_=q(rowj)
                do k=1,prj
                  ws(qrj_+k)=ws(qrj+k)
                enddo
                ws(qrj_+len_)=z
              else
                ws(qrj+len_)=z
              endif
              p(rowj)=len_
              goto2
            endif
            len_=max(pri,prj)+1
            if(len_.gt.len.or.pri.gt.prj)then
c  create a new slot and use saxpyz ...
              call newslot(rowj,len_,lastr,irow,p,q,r,s,ws,mxws,qrj,
     *          ifail)
              if(ifail.gt.0)return
              qrj_=q(rowj)
              len=prj-pri
              if(len.ge.0)then
                do k=1,len
                  ws(qrj_+k)=ws(qrj+k)
                enddo
                len=len+1
                call saxpyz(z,ws(qri+1),ws(qrj+len),ws(qrj_+len),
     *            len_-len)
              else
                len=-len
                do k=1,len
                  ws(qrj_+k)=z*ws(qri+k)
                enddo
                len=len+1
                call saxpyz(z,ws(qri+len),ws(qrj+1),ws(qrj_+len),
     *            len_-len)
              endif
              ws(qrj_+len_)=z
            else
c  ... else saxpy in-situ
              if(pri.gt.0)
     *          call mysaxpy(z,ws(qri+1),ws(qrj+prj-pri+1),pri)
              ws(qrj+len_)=z
            endif
            p(rowj)=len_
c           do rj=1,n
c             if(p(rj).ne.0)then
c               write(nout,*)'storage for row',rj,'  p,q,r,s =',
c    *            p(rj),q(rj),r(rj),s(rj)
c             endif
c           enddo
    2       continue
          enddo
c         write(nout,*)'lr =',(lr(j),j=i,mao(i))
c         do j=i,mao(i)
c           rowj=lr(j)
c           if(p(rowj).ne.0)then
c             write(nout,*)'L(',rowj,')',
c    *          (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c           endif
c         enddo
        enddo
        inode=mao(i)
        noder=lr(inode)
        pri=p(noder)
        if(pri.gt.0)then
         d(inode)=aiscpri2(n,a,la,noder,lc(inode)-n,ws(q(noder)+1),
     *     1.D0,inode-1,pri,li)
        else
          d(inode)=aij(noder,lc(inode)-n,a,la)
        endif
        if(d(inode).eq.0.D0)d(inode)=eps
      endif
      i=mao(i)+1
      if(i.le.n)goto1
c     write(nout,*)'PAQ factors:  nu =',nu
c     write(nout,*)'column perm =',(lc(j),j=nu+1,n)
c     write(nout,*)'row perm =',(lr(j),j=nu+1,n)
c     write(nout,*)'d =',(d(ij),ij=nu+1,n)
c     do j=nu+1,n
c       rowj=lr(j)
c       if(p(rowj).ne.0)then
c         write(nout,*)'L(',rowj,')',
c    *      (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c       endif
c     enddo
c     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
c  print star diagram
c     if(m1.gt.80.or.n.gt.1000)stop
c     write(nout,*)'factored tarjan + spk1 ordering:  nu =',nu
c     do i=1,n
c       do j=1,m1
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,m1
c       jp=la(0)+lc(nu+j)-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,m1)
c     enddo
c     write(nout,*)'lr =',(lr(j),j=nu+1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
      mp=-1
      mq=-1
      ifail=0
      return
      end

      function aiscpri2(n,a,la,rowi,coli,ws,di,im,pri,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),ws(*),li(*)
      integer rowi,coli,rowj,pri
      aiscpri2=0.D0
      jp=la(0)+coli
      do j=la(jp),la(jp+1)-1
        rowj=la(j)
        if(rowj.eq.rowi)then
          aiscpri2=aiscpri2+di*a(j)
        else
          ir=li(rowj)-im
          if(ir.gt.0)goto1
          ir=ir+pri
          if(ir.gt.0)aiscpri2=aiscpri2+ws(ir)*a(j)
        endif
    1   continue
      enddo
      return
      end

      subroutine updateSE(p,q,n,nm,a,la,e,Q_,R,EB,y,z,s,t,u,x,
     *  d,ws,lr,lc,li,pp,qq,lm,lv,le,ifail)
      implicit double precision (a-h,o-z)
      integer p,q,pc,qr,pp,qq
      dimension a(*),la(0:*),e(*),Q_(*),R(*),EB(*),y(*),z(*),s(*),t(*),
     *  x(*),u(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*),lm(*),lv(*),le(*)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/mc,mxmc
      common/epsc/eps,tol,emin
      common/pqc/pc,qr,lmp
    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(15I5))
    4 format(A/(5E15.7))
    5 format((5E15.7))
c     print 1,'updateSE: p,q =',p,q
c     print 3,'lm =',(lm(i),i=1,n)
      if(p.ne.mp)then
        call eptsol(n,a,la,p,a,Q_,R,EB,y,s,u,x,d,ws,
     *    lr,lc,li,lv,le,pp,qq,ep)
c       print 4,'sk in natural ordering',(s(li(i)),i=1,n)
c       print 4,'sk ordered by lr',(s(i),i=1,n)
      else
        ep=e(p)
      endif
      if(q.ne.mq)then
        call aqsol(n,a,la,q,Q_,R,EB,z,t,u,x,d,ws,
     *    lr,lc,li,lv,le,pp,qq)
c       print 4,'tk =',(t(i),i=1,n)
      endif

c  update steepest edge coefficients
      tp=t(lmp)
      eq=2.D0/ep
      do i=1,n
        u(i)=eq*s(i)
      enddo
      call aqsol(n,a,la,0,Q_,R,EB,x,u,s,x,d,ws,
     *  lr,lc,li,lv,le,pp,qq)
c     print 4,'u in natural order =',(u(li(i)),i=1,n)
c     print 4,'u ordered by lr =',(u(i),i=1,n)
      eq=ep/tp
      do j=1,n
        i=lm(j)
        if(e(i).eq.0.D0)then
          print *,'malfunction: e(i)=0.D0, i =',i
          return
        endif
        ei=e(i)
        wi=t(j)*eq
        awi=abs(wi)
        if(ei.ge.awi)then
          wi=wi/ei
          e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-u(j)/ei))))
        else
          wi=ei/wi
          e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-u(j)/ei))))
        endif
      enddo
      e(p)=0.D0
      e(q)=max(emin,abs(eq))
c     print 4,'e =',(e(i),i=1,nm)
      return
      end

      subroutine updateQR(p,q,n,a,la,Q_,R,EB,x,y,z,
     *  lr,lc,li,lv,le,lm,ifail)
      implicit double precision (a-h,o-z)
      integer p,q,pc,qr
      dimension a(*),la(0:*),Q_(*),R(*),EB(*),x(*),y(*),z(*),
     *  lr(*),lc(*),li(*),lv(*),le(*),lm(*)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/mc,mxmc
      common/epsc/eps,tol,emin
      common/pqc/pc,qr,lmp
    1 format(A,15I5)
    2 format(A,6E15.7)
    3 format(A/(15I5))
    4 format(A/(5E15.7))
    5 format((5E15.7))
c     print 1,'updateQR:  p,q =',p,q,pc,qr
c     print *,'mc =',mc
c  x is only used for checking
      m1mc=m1*mc
c     print 4,'EB on entry to updateQR =',(EB(i),i=1,m1mc+m1)
      if(pc.eq.0)then
c  p is in the starting active set
        if(qr.eq.0)then
c  q is a new column: extend QR case
c         print *,'extend  ',p,q
          mc1=mc+1
c         print 4,'EB for c/s p =',(EB(i),i=m1mc+1,m1mc+m1)
          ii=m1mc-m2
          if(q.le.n)then
            if(li(q).gt.m2)sum=-EB(li(q)+ii)
          else
            sum=0.D0
            jp=la(0)+q-n
            do j=la(jp),la(jp+1)-1
              i=la(j)
              if(i.eq.p)then
                sum=sum-a(j)
              elseif(li(i).gt.m2)then
                sum=sum-EB(li(i)+ii)*a(j)
              endif
            enddo
          endif
          y(mc1)=sum
          if(mc.eq.0)then
            Q_(1)=1.D0
            R(1)=y(1)
c           print 2,'R(1) =',R(1)
            mc=1
          else
c  new row and column of C
c           print 2,'y =',(y(i),i=1,mc1)
c           print 2,'z =',(z(i),i=1,mc)
            ic=mc1
            mmx=mc*mxmc
            do i=1,mc
              Q_(ic)=0.D0
              Q_(mmx+i)=0.D0
              ic=ic+mxmc
            enddo
            Q_(ic)=1.D0 
            ii=1
            ic=1
            mmx1=mmx+1
            do i=1,mc
              call angle(R(ii),y(i),cos,sin)
              call rot(mc-i,R(ii+1),y(i+1),cos,sin)
              call rot(mc1,Q_(ic),Q_(mmx1),cos,sin)
              ii=ii+mxmc-i+1
              ic=ic+mxmc
            enddo
            z(mc1)=y(mc1)
c           print 2,'z =',(z(i),i=1,mc1)
            ii=mc1
            ic=1
            do i=1,mc1
              R(ii)=scpr(0.D0,Q_(ic),z,mc1)
              ii=ii+mxmc-i
              ic=ic+mxmc
            enddo
            mc=mc1
          endif
          lv(mc)=q
          le(mc)=p
          lm(lmp)=q
          goto10
        else
c  row replacement case: p is in the starting active set and
c  q is in the starting active set (and hence in row qr of E)
c         print *,'row     ',p,q
c         print 3,'lr =',(lr(i),i=1,n)
c         print 3,'lc =',(lc(i),i=m2+1,n)
c         print 3,'le =',(le(i),i=1,mc)
c         print 3,'lv =',(lv(i),i=1,mc)
c         print 3,'lm =',(lm(i),i=1,n)
          mc1=mc
          mc=mc-1
          j=mc1-qr
          if(mc.gt.0)then
            icp=mxmc*mc+mc1
            ic=icp
            call rexch(Q_(icp),Q_(icp-j))
            ii=mc1*(3-mc1)/2+(mc1-1)*mxmc
            z(mc1)=R(ii)
            do i=mc,1,-1
              ic=ic-mxmc
              call rexch(Q_(ic),Q_(ic-j))
              call angle(Q_(icp),Q_(ic),cos,sin)
              call rot(mc,Q_(icp-mc),Q_(ic-mc),cos,sin)
              ii=ii-mxmc+i-1
              z(i)=sin*R(ii)
              R(ii)=-cos*R(ii)
              call rot(mc1-i,z(i+1),R(ii+1),cos,sin)
            enddo
            ic=1
            icp=icp-mc
            do i=1,mc
              call angle(R(ii),y(i),cos,sin)
              call rot(mc1-i,R(ii+1),y(i+1),cos,sin)
              call rot(mc1,Q_(ic),Q_(icp),cos,sin)
              ii=ii+mxmc-i+1
              ic=ic+mxmc
            enddo
            R(ii)=y(mc1)
            le(qr)=le(mc1)
          else
            Q_(1)=1.D0
            R(1)=y(1)
          endif
          le(mc1)=p
c         j=m1*(mc1-qr)
          j=m1*j
          do i=m1mc-m1+1,m1mc
            EB(i-j)=EB(i)
            EB(i)=EB(i+m1)
          enddo
          mc=mc1
        endif
      else
c  p is not in the starting active set (and hence in column pc of V)
c  first remove column pc from V
        call ishift(lv(pc),mc-pc,1)
        ic=pc
        do i=1,pc-1
          call r_shift(R(ic),mc-pc,1)
          ic=ic+mxmc-i
        enddo
        mc1=mc
        mc=mc-1
        ii=ic+1
        iip=ii+mxmc-pc
        ic=(pc-1)*mxmc+1
        do i=pc,mc
          call angle(R(ii),R(iip),cos,sin)
          call rot(mc-i,R(ii+1),R(iip+1),cos,sin)
          call r_shift(R(ii-1),mc1-i,1)
          ii=iip+1
          iip=ii+mxmc-i-1
          icp=ic+mxmc
          call rot(mc1,Q_(ic),Q_(icp),cos,sin)
          ic=icp
        enddo
        if(qr.eq.0)then
c         print *,'column  ',p,q
c  q is a new column (column interchange case)
c         print 2,'z =',(z(i),i=1,mc1)
          ii=mc1
          ic=1
          do i=1,mc1
            R(ii)=scpr(0.D0,Q_(ic),z,mc1)
            ii=ii+mxmc-i
            ic=ic+mxmc
          enddo
          mc=mc1
          lv(mc)=q
        else
c  q is in the starting active set and hence in row qr of E (contract QR case)
c         print *,'contract',p,q
          icp=mxmc*mc+mc1
          ic=icp
          j=mc1-qr
          call rexch(Q_(icp),Q_(icp-j))
          ii=mc*(3-mc)/2+(mc-1)*mxmc
          do i=mc,1,-1
            ic=ic-mxmc
            call rexch(Q_(ic),Q_(ic-j))
            call angle(Q_(icp),Q_(ic),cos,sin)
            call rot(mc,Q_(icp-mc),Q_(ic-mc),cos,sin)
            y(i)=sin*R(ii)
            R(ii)=-cos*R(ii)
            call rot(mc-i,y(i+1),R(ii+1),cos,sin)
            ii=ii-mxmc+i-2
          enddo
          le(qr)=le(mc1)
          j=m1*(mc1-qr)
          do i=m1mc-m1+1,m1mc
            EB(i-j)=EB(i)
          enddo
        endif
      endif
c  reset lm (except when extending)
      do i=1,m2
        lm(i)=lr(i)
      enddo
      do i=m2+1,n
        lm(i)=lc(i)
      enddo
      do i=1,mc
        lm(li(le(i)))=lv(i)
      enddo
c     print 3,'le =',(le(i),i=1,mc)
c     print 3,'lv =',(lv(i),i=1,mc)
c     print 3,'lm =',(lm(i),i=1,n)
   10 continue
      ifail=0
      if(mc.eq.0)return
c     print 4,'Q transpose =',(Q_(j),j=1,mc)
c     do i=1,mc-1
c       print 5,(Q_(j),j=i*mxmc+1,i*mxmc+mc)
c     enddo
c     print 4,'R =',(R(j),j=1,mc)
c     ii=mxmc+1
c     do i=1,mc-1
c       print 5,(R(j),j=ii,ii+mc-i-1)
c       ii=ii+mxmc-i
c     enddo
      sum=R(mc*(3-mc)/2+(mc-1)*mxmc)
      if(sum-sum.ne.0.D0.or.sum.eq.0.D0)ifail=1
      return
c  check QR factors
c     print 4,'EB in check =',(EB(i),i=1,m1mc)
      do j=1,mc
        do i=1,n
          x(i)=0.D0
        enddo
        if(lv(j).le.n)then
          x(li(lv(j)))=1.D0
        else
          call iscatter(a,la,lv(j)-n,li,x,n)
        endif
        do i=1,mc
          z(i)=0.D0
        enddo
        i1=j
        do i=1,j
          call mysaxpy(R(i1),Q_((i-1)*mxmc+1),z,mc)
          i1=i1+mxmc-i
        enddo
c       print 4,'z =',(z(i),i=1,mc)
c       print 4,'x =',(x(i),i=1,n)
        emax=0.D0
        do i=1,mc
          cij=-scpr(0.D0,EB((i-1)*m1+1),x(m2+1),m1)
          if(le(i).le.n)cij=cij-x(li(le(i)))
          qrij=z(i)
          emax=max(emax,abs(cij-qrij))
c         if(abs(cij-qrij).gt.tol)then
c           print *,'error in C=QR:  i,j =',i,j
c           print *,'cij,qrij =',cij,qrij
c           if(abs(cij-qrij).gt.1.D-6)stop
c         endif
        enddo
      enddo
      if(emax.gt.tol)print *,'max error in C=QR =',emax
      if(emax.gt.tol)stop
      return
      end

      subroutine newslot(row,len,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *  ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      parameter (igap=10)
      dimension p(*),q(*),r(*),s(*),ws(*)
      common/noutc/nout
c     write(nout,*)'newslot: row =',row,'   len =',len
c     write(nout,*)'irow,lastr,mxws =',irow,lastr,mxws
      ifail=0
      if(lastr.eq.0)then
        if(mxws.lt.len)then
          write(nout,*)'insufficient space available for profile'
          ifail=7
        else
          irow=row
          q(row)=0
          r(row)=0
          s(row)=0
          lastr=row
        endif
        return
      endif
      igp=igap
    1 continue
      len_=len+igp
      thisr=lastr
    2 continue
      qrow=q(thisr)+p(thisr)
      nextr=s(thisr)
c     write(nout,*)'thisr,nextr,qrow,p(thisr),len_',
c    *  thisr,nextr,qrow,p(thisr),len_
      if(nextr.ne.0)then
        if(q(nextr).ge.qrow+len_)then
c  free slot after this row
          goto4
        else
          thisr=nextr
          if(thisr.ne.lastr)goto2
        endif
      else
        if(mxws-qrow.ge.len_)then
c  free slot at end of ws
          goto4
        elseif(q(irow).ge.len_)then
c  free slot at beginning of ws
          qrow=0
          thisr=0
          nextr=irow
          irow=row
          igp=0
          goto4
        endif
        thisr=irow
        if(thisr.ne.lastr)goto2
      endif
c  no free space: try minimum value of len
      if(igp.gt.0)then
        igp=0
        goto1
      endif
c  compress ws
      thisr=irow
      qrow=0
    3 continue
      call r_shift(ws(qrow+1),p(thisr),q(thisr)-qrow)
      q(thisr)=qrow
      qrow=qrow+p(thisr)
      if(s(thisr).ne.0)then
        thisr=s(thisr)
        goto3
      endif
      if(mxws.lt.qrow+len_)then
        write(nout,*)'insufficient space available for profile'
        write(nout,*)'mxws,qrow,len_',mxws,qrow,len_
        ifail=7
        return
      endif
c  insert at end of compressed file
      nextr=0
    4 continue
      qr_=q(row)
      q(row)=qrow+igp
      if(p(row).gt.0)then
        if(r(row).eq.thisr.or.s(row).eq.nextr)return
c  insert after row thisr and take out old row
        call erase(row,lastr,irow,r,s)
      endif
      lastr=row
      r(row)=thisr
      if(thisr.gt.0)s(thisr)=row
      s(row)=nextr
      if(nextr.gt.0)r(nextr)=row
      i=0
      return
      end

      subroutine erase(row,lastr,irow,r,s)
c  remove slot for row from the data file
      implicit integer (i-s)
      dimension r(*),s(*)
      common/noutc/nout
c     write(nout,*)'erase: row,irow,lastr =',row,irow,lastr
      if(r(row).eq.0)then
        if(s(row).eq.0)then
          irow=0
          lastr=0
          return
        endif
        irow=s(row)
        r(irow)=0
      elseif(s(row).eq.0)then
        s(r(row))=0
      else
        s(r(row))=s(row)
        r(s(row))=r(row)
      endif
      if(row.eq.lastr)lastr=irow
      return
      end

      subroutine trim_(rowi,pri,qri,q,ws)
c  trim leading zeros off slot for row i
      implicit double precision (a-h,s-z), integer (i-r)
      dimension q(*),ws(*)
      common/epsc/eps,tol,emin
    1 continue
      qri=qri+1
      pri=pri-1
      if(pri.eq.0)return
      if(abs(ws(qri+1)).le.tol)goto1
      q(rowi)=qri
      return
      end

      subroutine EBspace(n,p,q,s,lr,ws,neb,nprof,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension p(*),q(*),s(*),lr(*),ws(*)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/mc,mxmc
      common/noutc/nout
      ifail=0
c  find length and last entry in file
      len=0
      last=0
      do i=m2+1,n
        row=lr(i)
        if(p(row).gt.0)then
          len=len+p(row)
          last=max(last,q(row)+p(row))
        endif
      enddo
      neb=m1*(mxmc+1)
c     print *,'nprof =',nprof
c     print *,'m1,len,last,neb',m1,len,last,neb
      if(last+neb.le.nprof)then
        neb=last
        return
      elseif(len+neb.le.nprof)then
c  compress the file
        thisr=irow
        qrow=0
    1   continue
        call r_shift(ws(qrow+1),p(thisr),q(thisr)-qrow)
        q(thisr)=qrow
        qrow=qrow+p(thisr)
        if(s(thisr).ne.0)then
          thisr=s(thisr)
          goto1
        endif
        neb=len
        return
      endif
      write(nout,*)'not enough additional space for E.B^{-1} matrix'
      write(nout,*)'space required is at least len+m1*(mxmc+1)'
      write(nout,*)'space left =',nprof-len,'  len,m1 =',len,m1
      ifail=7
    9 format(A/(15I5))
      return
      end

      subroutine checkperms(n,lr,lc,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension lr(*),lc(*),li(*)
      common/factorc/m1,m2,mp,mq,lastr,irow
        do i=1,n
          if(lr(li(i)).ne.i)print *,'wrong perm 1'
          if(lr(li(i)).ne.i)stop
          if(li(lr(i)).ne.i)print *,'wrong perm 2'
          if(li(lr(i)).ne.i)stop
        enddo
        do i=m2+1,n
          if(li(lc(i)).ne.i)print *,'wrong perm 3'
          if(li(lc(i)).ne.i)stop
        enddo
      return
      end

      subroutine checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      implicit double precision (a-h,r-z), integer (i-q)
      integer r,s,rowj,thisr
      dimension a(*),la(*),lr(*),lc(*),li(*),p(*),q(*),r(*),s(*),ws(*),
     *  d(*)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/noutc/nout
      common/epsc/eps,tol,emin
c  check indexing
      do j=1,nu
        if(p(lr(j)).ne.0)then
          write(nout,*)'p(lr(j)).ne.0'
          goto11
        endif
      enddo
      np=0
      do i=nu+1,n
        if(p(lr(i)).gt.0)np=np+1
      enddo
      if(irow.gt.0)then
        if(r(irow).ne.0)then
          write(nout,*)'r(irow).ne.0'
          goto11
        endif
        thisr=irow
    1   continue
        if(p(thisr).le.0)then
          write(nout,*)'p(thisr).le.0'
          goto11
        endif
        np=np-1
        nextr=s(thisr)
        if(nextr.eq.0)then
          if(q(thisr)+p(thisr).gt.mxws)then
            write(nout,*)'q(thisr)+p(thisr).gt.mxws'
            goto11
          endif
        else
          if(r(nextr).ne.thisr)then
            write(nout,*)'r(nextr).ne.thisr'
            goto11
          endif
          if(nextr.ne.s(thisr))then
            write(nout,*)'nextr.ne.s(thisr)'
            goto11
          endif
          if(q(thisr)+p(thisr).gt.q(nextr))then
            write(nout,*)'q(thisr)+p(thisr).gt.q(nextr)'
            goto11
          endif
          thisr=nextr
          goto1
        endif
      endif
      if(np.ne.0)then
        write(nout,*)'np.ne.0'
        goto11
      endif
      last=0
      emax=0.D0
      length=0
      do inode=nu+1,n
        nodec=lc(inode)
c  form L.a_q
        rowj=lr(inode)
        prj=p(rowj)
        length=length+prj
        if(prj.lt.0)then
          write(nout,*)'prj.lt.0'
          goto11
        elseif(prj.eq.0)then
          e=abs(aij(rowj,nodec-n,a,la)-d(inode))
        else
          e=abs(d(inode)-aiscpri2(n,a,la,rowj,nodec-n,ws(q(rowj)+1),
     *      1.D0,inode-1,prj,li))
        endif
c       if(e.gt.tol)write(nout,*)'error =',e,
c    *    '  inode,nodec,rowj =',inode,nodec,rowj
        emax=max(emax,e)
        do j=inode+1,n
          rowj=lr(j)
          prj=p(rowj)
          if(prj.gt.0)then
            e=abs(aiscpri2(n,a,la,rowj,nodec-n,ws(q(rowj)+1),1.D0,j-1,
     *         prj,li))
          else
            e=abs(aij(rowj,nodec-n,a,la))
          endif
c         if(e.gt.tol)write(nout,*)'error =',e,
c    *      '  inode,nodec,j,rowj =',inode,nodec,j,rowj
          emax=max(emax,e)
        enddo
      enddo
      write(nout,*)'checkout:  m1 =',m1,'  file length =',length
      if(emax.gt.tol)write(nout,*)'error =',emax
      return
   11 continue
      write(nout,*)'thisr,nextr =',thisr,nextr
      write(nout,*)'i,p(i),q(i),r(i),s(i):  irow =',irow
      do i=1,n
        if(p(i).ne.0)write(nout,*)i,p(i),q(i),r(i),s(i)
      enddo
      stop
      end
