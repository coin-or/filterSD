
      subroutine checkg(n,x,h,a,b,user,iuser,cuser,tol)
      implicit double precision (a-h,o-z)
      dimension x(*),h(*),a(*),b(*),user(*),iuser(*)
      character cuser(*)

c  Copyright (C) 2010 Roger Fletcher

c  Current version 11 January 2011

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c  gradient checking subroutine for use with the glcpd code

c  Parameters
c  **********
c  n    set the number of variables
c  x(n) set any vector of variables in x(i), i=1,n
c  h(n) set a vector of differencing increments in h(i), i=1,n
c        (h(i) is the difference interval for x(i) and must be nonzero,
c        say 0.01 times the typical magnitude of x(i))
c  a(n)  double precision real workspce
c  b(n)  double precision real workspce
c  user(*)  double precision real user storage, passed through to funct and grad
c  iuser(*) integer user storage, passed through to funct and grad
c  cuser(*) character user storage, passed through to funct and grad
c  tol  tolerace (>0) for reporting inconsistencies in derivatives (eg 1.D-12)

c  Usage
c  *****
c  The user must write subroutines 'funct' and 'grad' as for glcpd
c  Write a driver program for your problem,  but replace the call of glcpd
c  by a call of checkg (having set differencing increments in h).
c     The program will report any inconsistencies in the gradients.
c  If the difference quotient estimate lies between the gradient component
c  at x and x+[0,...,0,h(i),0,...,0]' then the derivative is assumed to be
c  correct. Small errors in this comparison may be ignored. If no errors are
c  reported then the call of glcpd may be restored.

      print *,'entering checkg'

      call funct(n,x,fx,user,iuser,cuser)
      call grad(n,x,a,user,iuser,cuser)
      do i=1,n
        xi=x(i)
        x(i)=x(i)+h(i)
        call funct(n,x,fxh,user,iuser,cuser)
        call grad(n,x,b,user,iuser,cuser)
        dfi=(fxh-fx)/h(i)
        if((dfi.ge.a(i)-tol.and.dfi.le.b(i)+tol).or.
     *    (dfi.ge.b(i)-tol.and.dfi.le.a(i)+tol))goto10
        print 1,'derivative inconsistency in variable',i
        print 2,'deriv at x, diff quotient, deriv at x+h =',
     *        a(i),dfi,b(i)
        stop
   10   continue
        x(i)=xi
      enddo
      print *,'exiting checkg'
    1 format(A,15I5)
    2 format(A,6E15.7)
    4 format(A/(6E15.7))
      return
      end
