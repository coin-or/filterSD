christen this file    sparseA.f

c  Copyright (C) 1996 Roger Fletcher

c  Current version dated 9 December 2010

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c  ******************************************
c  Specification of A in sparse matrix format
c  ******************************************

c  The matrix A contains gradients of the linear terms in the objective
c  function (column 0) and the general constraints (columns 1:m).
c  No explicit reference to simple bound constraints is required in A.
c  The information is set in the parameters a(*) (double precision real) and
c  la(*) (integer).

c  In this sparse format, these vectors have dimension  a(1:maxa)  and
c   la(0:maxla-1),  where maxa is at least nnza (the number of nonzero elements
c  in A), and maxla is at least  nnza+m+3.  la(0) and the last m+2 elementss
c  in la are pointers.

c  The vectors a(.) and la(.) must be set as follows:

c  a(j) and la(j) for j=1,nnza are set to the values and row indices (resp.)
c  of all the nonzero elements of A. Entries for each column are grouped
c  together in increasing column order. Within each column group, it is
c  not necessary to have the row indices in increasing order.

c  la(0) is a pointer which points to the start of the pointer information in
c  la. la(0) must be set to nnza+1 (or a larger value if it is desired to
c  allow for future increases to nnza).

c  The last m+2 elements of la(.) contain pointers to the first elements in
c  each of the column groupings. Thus la(la(0)+i)) for i=0,m is set to the
c  location in a(.) containing the first nonzero element for column i of A.
c  Also la(la(0)+m+1)) is set to nnza+1 (the first unused location in a(.)).

c  Copyright, University of Dundee (R.Fletcher), June 1996
c  Current version dated 31/01/07

      subroutine saipy(s,a,la,i,y,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*)
c  saxpy with column i of A
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        y(ir)=y(ir)+s*a(j)
      enddo
    1 format(A,15I2)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))

      return
      end

      subroutine msaipy(s,a,la,i,y,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*)
c  saxpy with modulus of column i of A
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        y(ir)=y(ir)+s*abs(a(j))
      enddo
      return
      end

c     subroutine daipy(s,a,la,i,y,n)
c     DOUBLE PRECISION a(*),y(*),d
c     dimension la(0:*)
c     if(s.eq.0.D0)return
c     d=dble(s)
c     jp=la(0)+i
c     do j=la(jp),la(jp+1)-1
c       ir=la(j)
c       y(ir)=y(ir)+d*dble(a(j))
c     enddo
c     return
c     end

      subroutine isaipy(s,a,la,i,y,n,lr,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        y(ir)=y(ir)+s*a(j)
      enddo
      return
      end

c the old isaipy was what might be called isaipy2

      subroutine isaipy1(s,a,la,i,y,n,lr,li,m1)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A_1
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        if(ir.le.m1)y(ir)=y(ir)+s*a(j)
      enddo
      return
      end

c     subroutine ssaipy(s,a,la,i,y,n)
c     implicit double precision (a-h,o-z)
c     dimension a(*),la(0:*),y(*)
c  saxpy with squares of column i of A
c     if(s.eq.0.D0)return
c     jp=la(0)+i
c     do j=la(jp),la(jp+1)-1
c       ir=la(j)
c       y(ir)=y(ir)+s*(a(j))**2
c     enddo
c     return
c     end

      function aiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*)
c  scalar product with column i of A
      aiscpr=b
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        aiscpr=aiscpr+x(ir)*a(j)
      enddo
      return
      end

      function daiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*)
      DOUBLE PRECISION daiscpr
      daiscpr=dble(b)
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        daiscpr=daiscpr+dble(x(ir))*dble(a(j))
      enddo
      return
      end

      function aiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A
      aiscpri=b
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        aiscpri=aiscpri+x(ir)*a(j)
      enddo
      return
      end

      function daiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),lr(*),li(*)
      DOUBLE PRECISION daiscpri
      daiscpri=dble(b)
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        daiscpri=daiscpri+dble(x(ir))*dble(a(j))
      enddo
      return
      end

c the old aiscpri was what might be called aiscpri2

      function aiscpri1(n,a,la,i,x,b,lr,li,m1)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A_1
      aiscpri1=b
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        if(ir.le.m1)aiscpri1=aiscpri1+x(ir)*a(j)
      enddo
      return
      end

      function ailen(n,a,la,i)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)
c  L2 length of column i of A
      ailen=0.D0
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ailen=ailen+a(j)**2
      enddo
      ailen=sqrt(ailen)
      return
      end

      subroutine iscatter(a,la,i,li,an,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),li(*),an(*)
c  indirect scatter into previously zeroed vector an
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        an(li(la(j)))=a(j)
      enddo
      return
      end

      subroutine iunscatter(a,la,i,li,an,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),li(*),an(*)
c  undo effect of iscatter
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        an(li(la(j)))=0.D0
      enddo
      return
      end

      function aij(i,j,a,la)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)
c  get element A(i,j)
      jp=la(0)+j
      do ij=la(jp),la(jp+1)-1
        ir=la(ij)
        if(ir.eq.i)then
          aij=a(ij)
          return
        endif
      enddo
      aij=0.D0
      return
      end

      subroutine setaij(aij,i,j,a,la)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)
c  set element A(i,j)
      jp=la(0)+j
      do jj=la(jp+1)-1,la(jp),-1
        ir=la(jj)
        if(ir.eq.i)then
          a(jj)=aij
          return
        endif
      enddo
      if(aij.eq.0.D0)return
      print *,'malfunction: no slot for A(i,j)'
      stop
      end

      subroutine cscale(n,m,a,la,x,bl,bu,s,menu,ifail)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),bl(*),bu(*),s(*)

c     Constraint scaling procedure for use prior to calling bqpd when using
c     sparseA.f

c     Parameters are set as for bqpd, except for s, menu and ifail

c     The user must set the parameter menu to control how the
c     x-variables are scaled (or equivalently how constraints i = 1:n
c     are scaled), as follows

c     menu = 1 indicates that a unit scaling applies to the x-variables

c     menu = 2 the user provides estimates s(i)>0 of the magnitude of
c              x(i) for i = 1:n. In this case the elements  x(i), bl(i), bu(i)
c              are divided by s(i) for i = 1:n.

c     In all cases, cscale goes on to scale the general constraints, in
c     such a way that the normal vector of each nontrivial constraint in
c     the scaled problem has an l_2 norm of unity. This scaling is also
c     applied to the right hand sides  bl(i), bu(i) for i = n+1:n+m.
c     The scaled data overwrites the original data.

c     cscale also scales the constant vector of the quadratic function,
c     which is found in a(1:n). However if a non-unit x-variable scaling
c     is used, it is necessary for the user to scale the Hessian matrix
c     G appropriately. This can be done by passing the x-variable scale
c     factors s(i) i = 1:n into the subroutine gdotx using the
c     parameter ws, and multiplying G(i,j) by s(i)*s(j) (possibly
c     implicitly).

c     cscale sets ifail = 1 to indicate that some s(i)< = 0,
c             and ifail = 2 to indicate an incorrect setting of menu.
c       Otherwise ifail = 0.

      integer pjp

      ifail=2
      if(menu.lt.1.or.menu.gt.2)return
      pjp=la(0)
c     z=1.D0/log(2.D0)
      if(menu.eq.1)then
        do j=1,n
          s(j)=1.D0
        enddo
      else
        ifail=1
        do j=1,n
          if(s(j).le.0.D0)return
        enddo
c       if(menu.eq.2)then
c         do j=1,n
c           s(j)=2.D0**nint(log(s(j))*z)
c         enddo
c       endif
        do j=1,n
          if(s(j).ne.1.D0)then
            x(j)=x(j)/s(j)
            bl(j)=bl(j)/s(j)
            bu(j)=bu(j)/s(j)
          endif
        enddo
        do j=1,la(pjp+1)-1
          a(j)=a(j)*s(la(j))
        enddo
      endif
      do i=1,m
        t=0.D0
        do j=la(pjp+i),la(pjp+i+1)-1
          a(j)=s(la(j))*a(j)
          t=t+a(j)**2
        enddo
        t=sqrt(t)
        if(t.eq.0.D0)then
          s(n+i)=1.D0
        else
c         t=2.D0**nint(log(t)*z)
          s(n+i)=t
          do j=la(pjp+i),la(pjp+i+1)-1
            a(j)=a(j)/t 
          enddo
          bl(n+i)=bl(n+i)/t
          bu(n+i)=bu(n+i)/t
        endif
      enddo
      ifail=0
      return
      end

      subroutine modify(n,m,sigma,s,a,la,maxa,iws)
      implicit double precision (a-h,o-z)
      dimension s(*),a(*),la(0:*),iws(*)

c  Modifies the sparse data structure to add an extra variable and duplicate
c  the general constraints, to enable scaled L-infinity QPs to be solved.
c  Scale factors given in s(1:m) and the coefficient of the objective function in sigma
c  For unscaled problems set s=ones and sigma=1.
c  Needs m+1 locations of integer workspace in iws(*)

      n1=n+1
      n1m=n1+m
      m1=m+1
      la0=la(0)
      nextra=la(la0+m1)-la(la0+1)+m1+m
      ij=la(la0+m1)+nextra
c     print 1,'la0,nextra,ij',la0,nextra,ij
      if(ij-1.gt.maxa)then
        print *,'not enough space:  reset maxa to at least ',ij-1
        stop
      endif
      do i=1,m1
        iws(i)=la(la0+i)
      enddo
      la0=ij
      la(la0+m1+m)=ij
c  set lower bounds
      do i=m,1,-1
        ij=ij-1
c       a(ij)=-1.D0
        a(ij)=-s(i)
        la(ij)=n1
        do j=iws(i+1)-1,iws(i),-1
          ij=ij-1
          a(ij)=-a(j)
          la(ij)=la(j)
        enddo
        la(la0+m+i)=ij
      enddo
c  set upper bounds
      do i=m,1,-1
        ij=ij-1
c       a(ij)=-1.D0
        a(ij)=-s(i)
        la(ij)=n1
        do j=iws(i+1)-1,iws(i),-1
          ij=ij-1
          a(ij)=a(j)
          la(ij)=la(j)
        enddo
        la(la0+i)=ij
      enddo
      ij=ij-1
      la(ij)=n1
      a(ij)=sigma
      la(la0)=1
      la(0)=la0
c     print 3,'pointers =',(la(i),i=la0,la0+m1+m)
c     print 4,'a =',(a(i),i=1,la(la0+m1+m)-1)
c     print 3,'la =',(la(i),i=1,la(la0+m1+m)-1)
    1 format(A,15I4)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
      return
      end

      subroutine restore(n,m,a,la)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)

c  restores the changes made by subroutine modify

      la0=la(0)
      do i=1,m
        do j=la(la0+i),la(la0+i+1)-2
          la(j-i)=la(j)
          a(j-i)=a(j)
        enddo
        la(la0+i)=la(la0+i)-i
      enddo
      la(la0+m+1)=la(la0+m+1)-m-1
c     print 3,'pointers =',(la(i),i=la0,la0+m+1)
c     print 4,'a =',(a(i),i=1,la(la0+m+1)-1)
c     print 3,'la =',(la(i),i=1,la(la0+m+1)-1)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
      return
      end

      subroutine extend_la(n,m,la,lax)
      implicit double precision (a-h,o-z)
      dimension la(0:*),lax(0:*)

c  Modifies the sparse data structure to add an extra variable and duplicate
c  the general constraints, to enable scaled L-infinity QPs to be solved.
c  The gradient vector is assumed to have a single non-zero entry (n+1).

      n1=n+1
      la0=la(0)
      lax0=2*(la(la0+m+1)-la(la0+1)+m+1)
      lax(0)=lax0
      lax(1)=n1
      lax(lax0)=1
      ijx=2
      jp=lax0+1
      do k=1,2
        do j=1,m
          lax(jp)=ijx
          do ij=la(la0+j),la(la0+j+1)-1
            lax(ijx)=la(ij)
            ijx=ijx+1
          enddo
          lax(ijx)=n1
          ijx=ijx+1
          jp=jp+1
        enddo
      enddo
      lax(jp)=ijx
c     print 3,'lax pointers =',(lax(i),i=lax0,lax0+m+m+1)
c     print 3,'lax =',(lax(i),i=1,lax(lax0+m+m+1)-1)
    3 format(A/(20I4))
      return
      end

      subroutine extend_a(n,m,a,la,ax,lax,s,sigma)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),ax(*),lax(0:*),s(*)

c  Extends the sparse data values to enable scaled L-infinity QPs to be solved.
c  Scale factors given in s(1:m) and the coefficient of the objective function in sigma
c  For unscaled problems set s=ones and sigma=1.

      n1=n+1
      la0=la(0)
      ax(1)=sigma
      ijx=2
c  set lower bounds
      do j=1,m
        do ij=la(la0+j),la(la0+j+1)-1
          ax(ijx)=-a(ij)
          ijx=ijx+1
        enddo
        ax(ijx)=-s(j)
        ijx=ijx+1
      enddo
c  set upper bounds
      do j=1,m
        do ij=la(la0+j),la(la0+j+1)-1
          ax(ijx)=a(ij)
          ijx=ijx+1
        enddo
        ax(ijx)=-s(j)
        ijx=ijx+1
      enddo
c     print 4,'ax =',(ax(i),i=1,lax(lax(0)+m+m+1)-1)
    4 format(A/(5E15.7))
      return
      end
