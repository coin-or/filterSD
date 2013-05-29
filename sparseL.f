christen this file sparseL.f
cut here >>>>>>>>>>>>>>>>>
c***************** sparse matrix routines for manipulating L *******************

c           ***************************************************
c           Basis matrix routines for bqpd with sparse matrices
c           ***************************************************

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
c  proved quite effective. The factors are updated by a variant of the
c  Fletcher-Matthews method, which has proved very reliable in practice.
c  However the B matrix is re-factored every 30 updates to control growth in
c  the total spike length.

c  Workspace
c  *********
c  The user needs to supply storage for the rows of L, although the amount
c  required is unknown a-priori.
c  sparse.f requires
c     5*n+nprof          locations of real workspace, and
c     9*n+m              locations of integer workspace
c  where nprof is the space required for storing the row spikes of the L matrix.
c  Storage for sparseL.f is situated at the end of the workspace arrays ws
c  and lws in bqpd.
c  Allow as much space for nprof as you can afford: the routine will report if
c  there is not enough. So far 10^6 locations has proved adequate for problems
c  of up to 5000 variables.

c  In addition the current version of bqpd.f requires
c     kmax*(kmax+9)/2+2*n+m   locations of real workspace in ws
c     kmax                    locations of integer workspace in lws
c  The user is also allowed to reserve storage in ws and lws, for use in the
c  user-supplied routine gdotx. This storage is situated at the start of the
c  arrays ws and lws. The user specifies the amount required by
c  setting the parameters kk and ll in the common block
c     common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  The user MUST also set mxws and mxlws to be (respectively) the total amount
c  of real and integer workspace for the arrays ws and lws.

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

c  Copyright, University of Dundee (R.Fletcher), January 1998
c  Current version dated 16/04/02

      subroutine start_up(n,nm,nmi,a,la,nk,e,ls,aa,ll,mode,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      common/wsc/kk,ll_,kkk,lll,mxws,mxlws
      common/epsc/eps,tol,emin
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/nup,nfreq
      nfreq=min(30,nfreq)
      nup=0
      ns=kk+kkk+5*n
      nt=ll_+lll+8*n+nmi
      nprof=mxws-ns
      if(nprof.le.0.or.nt.gt.mxlws)then
        write(nout,*)'not enough real (ws) or integer (lws) workspace'
        write(nout,*)'you give values for mxws and mxlws as',mxws,mxlws
        write(nout,*)'minimum values for mxws and mxlws are',ns,nt
        ifail=7
        return
      endif
    3 format(A/(20I5))
    4 format(A/(5E15.7))
c  set storage map for sparse factors
      ns=n
      ns1=ns+1
      nt=ns+n
      nt1=nt+1
      nu=nt+n
      nu1=nu+1
      nx=nu+n
      nx1=nx+1
      np=nx+n
      np1=np+1
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
      m=nm-n
      mp=-1
      mq=-1
c     write(nout,*)'ls',(ls(ij),ij=1,nk)
      if(mode.ge.3)then
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
        if(nk.eq.n)return
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
      do i=n+1,nmi
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
      call factor(n,nmi,nu_,nk,a,la,e,ls,aa(ns1),aa(nt1),aa(nu1),
     *  aa(nx1),ll,ll(lc1),ll(li1),ll(lm1),ll(lp1),ll(lq1),ll(lr1),
     *  ll(ls1),aa(np1),nprof,aa,ifail)
      if(ifail.gt.0)return
c     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
c     emax=0.D0
c     do i=1,nm
c       if(e(i).gt.0.D0)then
c         call eptsol(n,a,la,i,a,aa(ns1),aa(nt1),aa,aa(np1),
c    *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
c         ei=xlen(0.D0,aa(ns1),n)
c         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
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
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/noutc/nout
c     write(nout,*)'refactor'
      m=nm-n
      call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *  ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *  nprof,ifail)
      if(ifail.ge.1)then
c       write(nout,*)'failure in re_order (2)'
        return
      endif
      call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),ll(lm1),
     *  ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *  nprof,aa,ifail)
      if(ifail.eq.7)return
      call check_L(n,aa,ll(lp1),ifail)
      return
      end

      subroutine pivot(p,q,n,nm,a,la,e,aa,ll,ifail,info)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),aa(*),ll(*),info(*)
      common/noutc/nout
      common/iprintc/iprint
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/mxm1c/mxm1
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
c     write(nout,*)'pivot: p,q =',p,q
      ifail=0
      if(p.ne.mp)then
        call eptsol(n,a,la,p,a,aa(ns1),aa(nt1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        if(p.gt.n)then
          e(p)=xlen(0.D0,aa(ns1+m2),m1)
        else
          e(p)=xlen(1.D0,aa(ns1+m2),m1)
        endif
        epp=e(p)
        mp=p
      endif
      if(q.ne.mq)then
        call aqsol(n,a,la,q,a,aa(nt1),aa(nx1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        mq=q
      endif
c  update steepest edge coefficients
      tp=aa(nt+ll(li+p))
      if(tp.eq.0.D0)tp=eps
      ep=e(p)
      eq=2.D0/ep
c     do i=1,m2-1
c       aa(nu+i)=0.D0
c     enddo
c     do i=m2,n
      do i=1,n
        aa(nu+i)=eq*aa(ns+i)
      enddo
      call aqsol(n,a,la,-1,a,aa(nu1),aa(nx1),aa,aa(np1),
     *  ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
c     write(nout,*)'row perm',(ll(ij),ij=1,n)
c     write(nout,*)'column perm',(ll(lc+ij),ij=m2+1,n)
c     write(nout,*)'s =',(aa(ns+ij),ij=1,n)
c     write(nout,*)'t =',(aa(nt+ij),ij=1,n)
c     write(nout,*)'u =',(aa(nu+ij),ij=1,n)
      e(p)=0.D0
      eq=ep/tp
      do i=1,nm
        if(e(i).gt.0.D0)then
          j=ll(li+i)
          ei=e(i)
          wi=aa(nt+j)*eq
          awi=abs(wi)
          if(ei.ge.awi)then
            wi=wi/ei
            e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-aa(nu+j)/ei))))
          else
            wi=ei/wi
            e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-aa(nu+j)/ei))))
          endif
        endif
      enddo
      e(q)=max(emin,abs(eq))
      info(1)=info(1)+1
      if(nup.ge.nfreq)then
c     if(nup.ge.30)then
c  refactorize L
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
        m1=n-m2
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,ifail)
        if(ifail.ge.1)then
c         write(nout,*)'failure in re_order (3)'
          return
        endif
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,aa,ifail)
      else
c  update L
        call update_L(p,q,n,nm,a,la,ll,ll(lc1),ll(li1),ll(lm1),ll(lp1),
     *    ll(lq1),ll(lr1),ll(ls1),aa(np1),nprof,aa,aa(ns1),ifail)
      endif
      if(ifail.eq.7)return
      mp=-1
      mq=-1
      call check_L(n,aa,ll(lp1),ifail)
c     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
c     emax=0.D0
c     do i=1,nm
c       if(e(i).gt.0.D0)then
c         call eptsol(n,a,la,i,a,aa(ns1),aa(nt1),aa,aa(np1),
c    *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
c         ei=xlen(0.D0,aa(ns1),n)
c         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
c         emax=max(emax,abs(ei-e(i)))
c       endif
c     enddo
c     if(emax.ge.tol)
c    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine fbsub(n,jmin,jmax,a,la,q,b,x,ls,aa,ll,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),ls(*),aa(*),ll(*)

c  solves a system  B.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   jmin,jmax  (see description of ls below)
c   a,la   specification of QP problem data (as for bqpd)
c   q   an integer which, if in the range 1:n+m, specifies that the rhs vector
c       b is to be column q of the matrix A of general constraint normals.
c       In this case the parameter b is not referenced by fbsub.
c       If q=0 then b is taken as the vector given in the parameter b.
c   b(n)  must be set to the r.h.s. vector b (but only if q=0)
c   x(n+m)  contains the required part of the solution x, set according to the
c       index number of that component (in the range 1:n for a simple bound and
c       n+1:n+m for a general constraint)
c   ls(*)  an index vector, listing the components of x that are required.
c       Only the absolute value of the elements of ls are used (this allows
c       the possibility of using of the contents of the ls parameter of bqpd).
c       Elements of x in the range abs(ls(j)), j=jmin:jmax are set by fbsub.
c       These contortions allow bqpd to be independent of the basis matrix code.
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   save   indicates if fbsub is to save its copy of the solution for possible
c       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'fbsub  q =',q
      if(save)then
        if(q.ne.mq)then
          call aqsol(n,a,la,q,b,aa(nt1),aa(nx1),aa,aa(np1),
     *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
          mq=q
        endif
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nt+ll(li+i))
        enddo
      else
        call aqsol(n,a,la,q,b,aa(nu1),aa(nx1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nu+ll(li+i))
        enddo
      endif
      return
      end

      subroutine ztg(n,k,rg,lv,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rg(*),lv(*),aa(*),ll(*)
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
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
c   ep  if p.ne.0 and save is true, ep contains the l_2 length of x on exit
c   save  indicates if tfbsub is to save its copy of the solution for possible
c       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'tfbsub  p =',p
      if(save)then
        if(p.ne.mp)then
          call eptsol(n,a,la,p,b,aa(ns1),aa(nt1),aa,aa(np1),
     *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
          mp=p
        endif
        do i=1,n
          x(ll(i))=aa(ns+i)
        enddo
        if(p.gt.n)then
          ep=xlen(0.D0,aa(ns1+m2),m1)
        elseif(p.gt.0)then
          ep=xlen(1.D0,aa(ns1+m2),m1)
        endif
      else
        call eptsol(n,a,la,p,b,aa(nu1),aa(nt1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        do i=1,n
          x(ll(i))=aa(nu+i)
        enddo
      endif
c     write(nout,*)'x =',(x(i),i=1,n)
      return
      end

      subroutine newg
      common/factorc/m1,m2,mp,mq,lastr,irow
      mq=-1
      return
      end

c******** The following routines are internal to sparseL.f **************

      subroutine check_L(n,d,p,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension d(*),p(*)
      common/noutc/nout
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/epsc/eps,tol,emin
c     write(nout,*)'check_L'
      ifail=1
c     dmin=1.D37
      do k=nu+1,n
c       dmin=min(dmin,abs(d(k)))
        if(abs(d(k)).le.tol)return
      enddo
c     write(nout,*)'dmin =',dmin
c     len=0
c     do i=1,n
c       len=len+p(i)
c     enddo
c     write(nout,*)m1*(m1+1)/2,len+m1
c     write(nout,*)'m1 =',m1,'   file length =',len,'   total =',len+m1
      ifail=0
      return
      end

      subroutine aqsol(n,a,la,q,b,tn,xn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),tn(*),xn(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*)
      common/noutc/nout
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'aqsol  q =',q
      if(q.gt.0)then
        do i=1,n
          tn(i)=0.D0
        enddo
        if(q.le.n)then
          tn(li(q))=1.D0
        else
          call iscatter(a,la,q-n,li,tn,n)
        endif
      elseif(q.eq.0)then
        do i=1,n
          tn(li(i))=b(i)
        enddo
      endif
c     write(nout,*)'tn =',(tn(i),i=1,n)
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
c     write(nout,*)'tn =',(tn(i),i=1,n)
      return
      end

      subroutine eptsol(n,a,la,p,b,sn,tn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),sn(*),tn(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*)
      common/noutc/nout
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'eptsol  p =',p
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
c     write(nout,*)'sn =',(sn(i),i=1,n)
      return
    1 continue
      write(nout,*)'malfunction detected in eptsol: p =',p
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
      common/refactorc/nup,nfreq
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
c     write(nout,*)'d =',(d(ij),ij=m2+1,n)
c     do j=m2+1,n
c       rowp=lr(j)
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
c     write(nout,*)'ls =',(ls(j),j=1,n)
c     write(nout,*)'s.e. coeffs =',(e(i),i=1,nm)
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'lc =',(lc(j),j=m2+1,n)
c     write(nout,*)'mao =',(mao(j),j=m2+1,n)
c     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
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
      common/refactorc/nup,nfreq
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

      subroutine update_L(pp,qq,n,nm,a,la,lr,lc,li,mao,p,q,r,s,
     *  ws,mxws,d,sn,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),lr(*),lc(*),li(*),mao(*),
     *  p(*),q(*),r(*),s(*),ws(*),d(*),sn(*)
c     character star(1000,80)
      double precision l11,l21
      integer r,s,rowim,rowi,rowj,rrj
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/refactorc/nup,nfreq
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/noutc/nout
      parameter (thresh=1.D-1,growth=1.D1)
c     write(nout,*)'update_L:  p,q =',pp,qq
      nup=nup+1
      if(qq.gt.n)then
        ilast=nu
        jp=la(0)+qq-n
        do j=la(jp),la(jp+1)-1
          ip=li(la(j))
          if(ip.gt.nu)ilast=max(ilast,mao(ip))
        enddo
        qqq=qq
      else
c  row flma procedure to remove row qq (includes qq amongst the unit vectors)
        iq=li(qq)
        if(iq.le.nu)goto99
        ilast=mao(iq)
        l11=1.D0
        u11=d(iq)
        ss=-sn(iq)
        nu=nu+1
        do i=iq,nu+1,-1
          lr(i)=lr(i-1)
          li(lr(i))=i
          sn(i)=sn(i-1)
          d(i)=d(i-1)
        enddo
        lr(nu)=qq
        li(qq)=nu
c  update mao
        do j=iq-1,nu,-1
          if(mao(j).lt.ilast)goto5
        enddo
        j=nu-1
    5   continue
        do j=j,nu,-1
          mao(j+1)=mao(j)+1
        enddo
        prq=p(qq)
        if(prq.gt.0)qrq=q(qq)
        do i=iq+1,ilast
          im=i-1
          rowi=lr(i)
          pri=p(rowi)
          u22=d(i)
          if(prq.gt.0)then
            u12=aiscpri2(n,a,la,qq,lc(i)-n,ws(qrq+1),l11,im,prq,li)
          else
            u12=l11*aij(qq,lc(i)-n,a,la)
          endif
          if(abs(u12).le.tol)u12=0.D0
          if(pri.gt.0)then
            qri=q(rowi)
            is=im-iq
            ii=pri-is
            if(ii.le.0)then
              l21=0.
            else
              l21=ws(qri+ii)
              if(abs(l21).le.tol)l21=0.D0
              if(ii.eq.1)then
                call trim_(rowi,pri,qri,q,ws)
                if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
                if(s(rowi).eq.0)then
                  qr_=mxws
                else
                  qr_=q(s(rowi))
                endif
                if(qri+pri.ge.qr_)then
                  call r_shift(ws(qri),pri,1)
                  qri=qri-1
                  q(rowi)=qri
                endif
              else
                pri=pri-1
                call r_shift(ws(qri+ii),is,1)
              endif
            endif
            p(rowi)=pri
          else
            l21=0.D0
          endif
          rr=-l21/l11
          del=rr*u12+u22
          test=abs(rr)*max(abs(u11),abs(u22))
c         write(nout,*)'l11,l21,u11,u12,u22,del,test',
c    *      l11,l21,u11,u12,u22,del,test
          is=pri-prq
          if(is.lt.0)test=test*growth
          if(u12.eq.0.D0.and.is.gt.0)test=test*thresh
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'rowq,prq,qrq =',qq,prq,qrq
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowq =',(ws(qrq+ij),ij=1,prq)
c           write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
          if(abs(del).le.test)then
c  no-perm operation for row flma
c           write(nout,*)'no-perm operation for row flma'
            if(is.gt.0)then
              pr_=prq
              prq=pri+1
              call newslot(qq,prq,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
              if(ifail.gt.0)return
              qrq=q(qq)
              qri=q(rowi)
              call r_shift(ws(qrq+1),pri,qri-qrq)
              call mysaxpy(rr,ws(qr_+1),ws(qri+is+1),pr_)
            else
              if(prq.eq.0)then
                call erase(rowi,lastr,irow,r,s)
                p(rowi)=0
                call newslot(qq,1,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
                if(ifail.gt.0)return
                prq=1
                qrq=q(qq)
              else
                is=-is
                do j=1,is
                  ws(qrq+j)=rr*ws(qrq+j)
                enddo
                if(pri.gt.0)then
                  call saxpyx(rr,ws(qrq+is+1),ws(qri+1),pri)
                else
                  call newslot(rowi,1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *              ifail)
                  if(ifail.gt.0)return
                  qri=q(rowi)
                  qrq=q(qq)
                endif
                if(abs(ws(qrq+1)).le.tol)call trim_(qq,prq,qrq,q,ws)
c  rename qq as rowi and vice-versa
                if(qri.lt.qrq)then
                  if(s(rowi).eq.qq)then
                    r(qq)=r(rowi)
                    r(rowi)=qq
                    s(rowi)=s(qq)
                    s(qq)=rowi
                  else
                    call iexch(r(qq),r(rowi))
                    call iexch(s(qq),s(rowi))
                    r(s(qq))=qq
                    s(r(rowi))=rowi
                  endif
                  if(r(qq).gt.0)then
                    s(r(qq))=qq
                  else
                    irow=qq
                  endif
                  if(s(rowi).gt.0)r(s(rowi))=rowi
                else
                  if(s(qq).eq.rowi)then
                    r(rowi)=r(qq)
                    r(qq)=rowi
                    s(qq)=s(rowi)
                    s(rowi)=qq
                  else
                    call iexch(r(rowi),r(qq))
                    call iexch(s(rowi),s(qq))
                    r(s(rowi))=rowi
                    s(r(qq))=qq
                  endif 
                  if(r(rowi).gt.0)then
                    s(r(rowi))=rowi
                  else
                    irow=rowi
                  endif
                  if(s(qq).gt.0)r(s(qq))=qq
                endif
                call iexch(pri,prq)
                call iexch(qri,qrq)
                call iexch(q(rowi),q(qq))
                if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
                prq=prq+1
              endif
            endif
            p(rowi)=pri
            p(qq)=prq
            ws(qrq+prq)=1.D0
            d(i)=rr*u11
            u11=u22
            l11=l21
          else
c  perm operation for row flma
c           write(nout,*)'perm operation for row flma'
            if(rr.ne.0.D0)then
              if(is.ge.0)then
                if(prq.gt.0)then
                  call mysaxpy(rr,ws(qrq+1),ws(qri+is+1),prq)
                  if(abs(ws(qri+1)).le.tol)call trim_(rowi,pri,qri,q,ws)
                  if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
                endif
                is=pri-prq
              else
                pr_=pri
                pri=prq
                call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)
                if(ifail.gt.0)return
                qrq=q(qq)
                qri=q(rowi)
                is=-is
                do j=1,is
                  ws(qri+j)=rr*ws(qrq+j)
                enddo
                call saxpyz(rr,ws(qrq+is+1),ws(qr_+1),ws(qri+is+1),pr_)
                is=0
              endif
            endif
            p(rowi)=pri
            if(u12.ne.0.D0)then
              u12=-u12/del
              if(is.gt.0)then
                pr_=prq
                prq=pri+1
                call newslot(qq,prq,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)
                if(ifail.gt.0)return
                qrq=q(qq)
                qri=q(rowi)
                do j=1,is
                  ws(qrq+j)=u12*ws(qri+j)
                enddo
                call saxpyz(u12,ws(qri+is+1),ws(qr_+1),ws(qrq+is+1),pr_)
                ws(qrq+prq)=u12
                goto7
              else
                if(pri.gt.0)then
                  is=-is
                  call mysaxpy(u12,ws(qri+1),ws(qrq+is+1),pri)
                  if(abs(ws(qrq+1)).le.tol)then
                    call trim_(qq,prq,qrq,q,ws)
                    if(prq.eq.0)call erase(qq,lastr,irow,r,s)
                    p(qq)=prq
                  endif
                endif
              endif
            endif
            if(prq.gt.0.or.u12.ne.0.D0)then
              if(prq.eq.0)then
                len=0
              elseif(s(qq).eq.0)then
                len=mxws-qrq
              else
                len=q(s(qq))-qrq
              endif
              if(len.eq.prq)then
                call newslot(qq,prq+1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)     
                if(ifail.gt.0)return
                qrq=q(qq)
                qri=q(rowi)
                call r_shift(ws(qrq+1),prq,qr_-qrq)
              endif
              prq=prq+1
              ws(qrq+prq)=u12
            endif
   7        continue
            p(rowi)=pri
            p(qq)=prq
            d(i)=del
            u11=u11*u22/del
            call iexch(lc(i),lc(im))
          endif
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'rowq,prq,qrq =',qq,prq,qrq
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowq* =',(ws(qrq+ij),ij=1,prq)
c           write(nout,*)'rowi* =',(ws(qri+ij),ij=1,pri)
        enddo
        if(prq.gt.0)then
c         write(nout,*)'ss,l11,ilast,n,prq',ss,l11,ilast,n,prq
c         write(nout,*)'sn =',(sn(ij),ij=nu+1,n)
          call mysaxpy(ss/l11,ws(qrq+1),sn(ilast-prq+1),prq)
          call erase(qq,lastr,irow,r,s)
          p(qq)=0
        endif
        qqq=lc(ilast)
        do i=ilast,nu+1,-1
          lc(i)=lc(i-1)
          li(lc(i))=i
        enddo
c       if(pp.le.n)then
c         ip=li(pp)
c         write(nout,*)'check sn'
c         do i=nu+1,ilast
c           nodec=lc(i)
c           u12=aiscpri2(n,a,la,pp,lc(i)-n,sn(nu+1),1.D0,ilast,
c             ilast-nu,li)
c           if(abs(u12).gt.tol)write(nout,*)'error,nodec =',u12,nodec
c         enddo
c       endif
c       write(nout,*)'intermediate PAQ factors:  new q =',qqq
c       write(nout,*)'lr =',(lr(j),j=nu+1,n)
c       write(nout,*)'lc =',(lc(j),j=nu+1,n)
c       write(nout,*)'d =',(d(ij),ij=nu+1,n)
c       do j=nu+1,n
c         rowj=lr(j)
c         if(p(rowj).ne.0)then
c           write(nout,*)'L(',rowj,')',
c    *        (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c         endif
c       enddo
c       call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      endif
      ip=li(pp)
      if(pp.gt.n)then
        li(pp)=0
        if(pp.eq.qqq)goto30
        if(ip.le.nu)goto99
        iout=ip
        rowim=lr(ip)
        prim=p(rowim)
        if(prim.gt.0)qrim=q(rowim)
      else
        if(ip.gt.nu.or.p(pp).gt.0)goto99
        lr(ip)=lr(nu)
        li(lr(ip))=ip
c  check for growth in sn
c       write(nout,*)'sn =',(sn(i),i=nu+1,n)
        iout=ilast
        i=nu+1
        if(i.gt.ilast)goto13
   11   continue
          do j=i,mao(i)
            if(abs(sn(j)).gt.growth)then
              iout=i-1
              goto13
            endif
          enddo
          i=mao(i)+1
          if(i.le.ilast)goto11
   13   continue
        do j=nu+1,iout
          if(abs(sn(j)).gt.tol)goto14
        enddo
        j=iout+1
   14   continue
        rowim=pp
        prim=iout-j+1
        if(prim.gt.0)then
          call newslot(pp,prim,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
          if(ifail.gt.0)return
          p(pp)=prim
          qrim=q(pp)
          ii=qrim
          do j=j,iout
            ii=ii+1
            ws(ii)=sn(j)
          enddo
        endif
        do i=nu,iout-1
          lr(i)=lr(i+1)
          li(lr(i))=i
          lc(i)=lc(i+1)
          li(lc(i))=i
          d(i)=d(i+1)
        enddo
        lr(iout)=pp
        li(pp)=iout
c       write(nout,*)'lr =',(lr(ij),ij=nu,iout)
c       write(nout,*)'lc =',(lc(ij),ij=nu,iout-1)
c       if(prim.gt.0)write(nout,*)'L(',pp,') =',(ws(qrim+j),j=1,prim)
        nu=nu-1
      endif
c     write(nout,*)'iout,ilast,rowim,prim =',iout,ilast,rowim,prim
c  column flma operations to restore L to triangular form
      iswap=0
      do i=iout+1,ilast
        im=i-1
        lc(im)=lc(i)
        li(lc(im))=im
        rowi=lr(i)
        pri=p(rowi)
c       if(pri.gt.0)write(nout,*)'L(',rowi,') =',(ws(q(rowi)+j),j=1,pri)
        u22=d(i)
        if(prim.gt.0)then
          u12=aiscpri2(n,a,la,rowim,lc(i)-n,ws(qrim+1),1.D0,im-1,prim,
     *      li)
          if(abs(u12).le.tol)u12=0.D0
        else
          u12=aij(rowim,lc(i)-n,a,la)
        endif
        if(pri.gt.0)then
c         write(nout,*)'pri,iswap',pri,iswap
          qri=q(rowi)
          ii=pri-iswap
          if(ii.le.0)then
            l21=0.D0
          else
            l21=ws(qri+ii)
            if(abs(l21).le.tol)l21=0.D0
            if(ii.eq.1)then
              call trim_(rowi,pri,qri,q,ws)
              if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
              if(s(rowi).eq.0)then
                qr_=mxws
              else
                qr_=q(s(rowi))
              endif
              if(qri+pri.ge.qr_)then
                call r_shift(ws(qri),pri,1)
                qri=qri-1
                q(rowi)=qri
              endif
            else
              pri=pri-1
              call r_shift(ws(qri+ii),iswap,1)
            endif
            p(rowi)=pri
c           write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
          endif
        else
          l21=0.D0
        endif
        del=u22-l21*u12
        test=abs(u12)*max(1.D0,abs(l21))
c       write(nout,*)'l21,u12,u22,del,test',l21,u12,u22,del,test
        is=pri-prim
        if(is.gt.0)test=growth*test
        if(l21.eq.0.D0.and.is.lt.0)test=thresh*test
c         write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
c         write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c         write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c         do j=1,n
c           if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c         enddo
c         write(nout,*)'rowim =',(ws(qrim+ij),ij=1,prim)
c         write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
        if(abs(del).le.test)then
c  no-perm operation for column flma
c         write(nout,*)'no-perm operation for column flma'
          rr=-u22/u12
          l21=l21+rr
          if(abs(l21).le.tol)l21=0.D0
          if(is.ge.0)then
            if(prim.gt.0)then
              call mysaxpy(rr,ws(qrim+1),ws(qri+is+1),prim)
              if(abs(ws(qri+1)).le.tol)call trim_(rowi,pri,qri,q,ws)
              if(pri.eq.0)then
                call erase(rowi,lastr,irow,r,s)
                p(rowi)=0
              endif
            endif
            if(pri.gt.0.or.l21.ne.0.D0)then
              if(pri.eq.0)then
                len=0
              elseif(s(rowi).eq.0)then
                len=mxws-qri
              else
                len=q(s(rowi))-qri
              endif
              if(len.eq.pri)then
                call newslot(rowi,pri+1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)
                if(ifail.gt.0)return
                qrim=q(rowim)
                qri=q(rowi)
                call r_shift(ws(qri+1),pri,qr_-qri)
              endif
              pri=pri+1
              ws(qri+pri)=l21
            endif
          else
            pr_=pri
            pri=prim+1
            call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)      
            if(ifail.gt.0)return
            qrim=q(rowim)
            qri=q(rowi)
            is=-is
            do j=1,is
              ws(qri+j)=rr*ws(qrim+j)
            enddo
            call saxpyz(rr,ws(qrim+is+1),ws(qr_+1),ws(qri+is+1),pr_)
            ws(qri+pri)=l21
          endif
c           write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
c           write(nout,*)'rowi* =',(ws(q(rowi)+ij),ij=1,p(rowi))
          p(rowi)=pri
          rowim=rowi
          prim=pri
          qrim=qri
          d(im)=u12
c  perform accumulated cyclic permutation in subsequent rows
          if(iswap.gt.0)then
            do j=i+1,ilast
              rowj=lr(j)
              prj=p(rowj)
              is=prj-j+i
              if(is.gt.0)then
                qrj=q(rowj)
                if(is.gt.iswap)then
                  ii=is-iswap
                  l21=ws(qrj+ii)
                  call r_shift(ws(qrj+ii),iswap,1)
                  ws(qrj+is)=l21
                  if(abs(ws(qrj+1)).le.tol)call trim_(rowj,prj,qrj,q,ws)
                  if(prj.eq.0)call erase(rowj,lastr,irow,r,s)
                else
                  prj=prj+1
                  rrj=r(rowj)
                  if(rrj.eq.0)then
                    len=qrj
                  else
                    len=qrj-q(rrj)-p(rrj)
                  endif
                  if(len.gt.0)then
                    call r_shift(ws(qrj),is,1)
                    ws(qrj+is)=0.D0
                    qrj=qrj-1
                    q(rowj)=qrj
                  else
                    call newslot(rowj,prj,lastr,irow,p,q,r,s,ws,mxws,
     *                qr_,ifail) 
                    if(ifail.gt.0)return
                    qrj=q(rowj)
                    qrim=q(rowim)
                    call r_shift(ws(qrj+1),is,qr_-qrj)
                    ws(qrj+is+1)=0.D0
                    call r_shift(ws(qrj+is+2),j-i,qr_-qrj-1)
                  endif
                endif
                p(rowj)=prj
c               write(nout,*)'L(',rowj,')* =',(ws(qrj+ij),ij=1,prj)
              endif
            enddo
          endif
          iswap=0
        else
c  perm operation for column flma
c         write(nout,*)'perm operation for column flma'
          rr=-l21
          if(rr.ne.0.D0)then
            if(is.ge.0)then
              if(prim.gt.0)then
                call mysaxpy(rr,ws(qrim+1),ws(qri+is+1),prim)
                if(abs(ws(qri+1)).le.tol)call trim_(rowi,pri,qri,q,ws)
                if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
              endif
              is=pri-prim
            else
              pr_=pri
              pri=prim
              call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *          ifail)
              if(ifail.gt.0)return
              qrim=q(rowim)
              qri=q(rowi)
              is=-is
              do j=1,is
                ws(qri+j)=rr*ws(qrim+j)
              enddo
              call saxpyz(rr,ws(qrim+is+1),ws(qr_+1),ws(qri+is+1),pr_)
              is=0
            endif
          endif
          p(rowi)=pri
          if(u12.ne.0.D0)then
            u12=-u12/del
            if(is.gt.0)then
              pr_=prim
              prim=pri+1
              call newslot(rowim,prim,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *          ifail)
              if(ifail.gt.0)return
              qrim=q(rowim)
              qri=q(rowi)
              do j=1,is
                ws(qrim+j)=u12*ws(qri+j)
              enddo
              call saxpyz(u12,ws(qri+is+1),ws(qr_+1),ws(qrim+is+1),pr_)
              ws(qrim+prim)=u12
              goto27
            else
              if(pri.gt.0)then
                is=-is
                call mysaxpy(u12,ws(qri+1),ws(qrim+is+1),pri)
                if(abs(ws(qrim+1)).le.tol)then
                  call trim_(rowim,prim,qrim,q,ws)
                  if(prim.eq.0)call erase(rowim,lastr,irow,r,s)
                  p(rowim)=prim
                endif
              endif
            endif
          endif
          if(prim.gt.0.or.u12.ne.0.D0)then
            if(prim.eq.0)then
              len=0
            elseif(s(rowim).eq.0)then
              len=mxws-qrim
            else
              len=q(s(rowim))-qrim
            endif
            if(len.eq.prim)then
              call newslot(rowim,prim+1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *          ifail)
              if(ifail.gt.0)return
              qrim=q(rowim)
              qri=q(rowi)
              call r_shift(ws(qrim+1),prim,qr_-qrim)
            endif
            prim=prim+1
            ws(qrim+prim)=u12
          endif
   27     continue
          p(rowim)=prim
          p(rowi)=pri
c           write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
c           write(nout,*)'rowi* =',(ws(q(rowi)+ij),ij=1,p(rowi))
          d(im)=del
          call iexch(lr(i),lr(i-1))
          call iexch(li(lr(i)),li(lr(i-1)))
          iswap=iswap+1
        endif
      enddo
      lc(ilast)=qqq
      li(qqq)=ilast
c     write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
c     write(nout,*)'ilast,prim,qrim',ilast,prim,qrim
      if(prim.gt.0)then
       d(ilast)=aiscpri2(n,a,la,rowim,qqq-n,ws(qrim+1),1.D0,ilast-1,
     *    prim,li)
      else
        d(ilast)=aij(rowim,qqq-n,a,la)
      endif
c  reset mao
      iout=ilast
      do i=ilast,nu+1,-1
        mao(i)=ilast
        iout=min(iout,i-p(lr(i)))
        if(iout.eq.i)ilast=i-1
      enddo
   30 continue
      m1=n-nu
c     write(nout,*)'PAQ factors:  nu =',nu
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
c     write(nout,*)'updated ordering:  nu =',nu
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
c     write(nout,*)'mao =',(mao(j),j=nu+1,n)
      return
   99 continue
      write(nout,*)'malfunction in update_L:  p,q =',pp,qq
      stop
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
