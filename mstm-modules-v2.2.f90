!
!  numerical constants
!
!
!  last revised: 15 January 2011
!
      module numconstants
      implicit none
      integer :: print_intermediate_results
      integer, allocatable :: monen(:)
      integer, private :: nmax=0
      real(8) :: pi
      real(8), allocatable :: bcof(:,:),fnr(:),vwh_coef(:,:,:,:)
      real(8), allocatable :: vcc_const(:,:,:),fnm1_const(:,:),fn_const(:,:),fnp1_const(:,:)
      data pi/3.141592653589793/

      contains

         subroutine init(notd)
         implicit none
         integer :: notd,l,n,ierr,nbc,m,mm1,mp1,np1,nm1,nn1,mn
         real(8) :: fnorm1,fnorm2
!
!  bcof(n,l)=((n+l)!/(n!l!))^(1/2)
!
         if(notd.le.nmax) return
         nmax=max(nmax,notd)
         nbc=6*notd+6
         if(allocated(fnr)) deallocate(monen,fnr,bcof)
         allocate (monen(0:2*notd),bcof(0:nbc,0:nbc),fnr(0:2*nbc),stat=ierr)
!         write(*,'('' nmax, bcof status:'',2i5)') nmax,ierr
         do n=0,2*notd
            monen(n)=(-1)**n
         enddo
         fnr(0)=0.d0
         do n=1,2*nbc
            fnr(n)=dsqrt(dble(n))
         enddo
         bcof(0,0)=1.d0
         do n=0,nbc-1
            do l=n+1,nbc
               bcof(n,l)=fnr(n+l)*bcof(n,l-1)/fnr(l)
               bcof(l,n)=bcof(n,l)
            enddo
            bcof(n+1,n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n,n)/fnr(n+1)/fnr(n+1)
         enddo
         if(allocated(vwh_coef)) deallocate(vwh_coef)
         allocate(vwh_coef(-notd:notd,1:notd,-1:1,-1:1))
!
!  constants used for calculation of svwf functions.
!
         do n=1,notd
            nn1=n*(n+1)
            np1=n+1
            nm1=n-1
            fnorm1=-.5d0/fnr(n+n+1)/fnr(n)/fnr(n+1)
            fnorm2=-.5d0*fnr(n+n+1)/fnr(n)/fnr(n+1)
            m=-n
            mp1=m+1
            mm1=m-1
            vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
            vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
            vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
            vwh_coef(m,n,-1,-1)=0.d0
            vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
            vwh_coef(m,n, 0,-1)=0.d0
            vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
            vwh_coef(m,n,-1, 0)=-0.d0
            vwh_coef(m,n, 0, 0)=-fnorm2*m
            do m=-n+1,-1
               mp1=m+1
               mm1=m-1
               vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
               vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
               vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
               vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
               vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
               vwh_coef(m,n, 0,-1)=fnorm1*np1*fnr(n+m)*fnr(n-m)
               vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
               vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
               vwh_coef(m,n, 0, 0)=-fnorm2*m
            enddo
            do m=0,n-1
               mp1=m+1
               mm1=m-1
               vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
               vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
               vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
               vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
               vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
               vwh_coef(m,n, 0,-1)=fnorm1*np1*fnr(n+m)*fnr(n-m)
               vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
               vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
               vwh_coef(m,n, 0, 0)=-fnorm2*m
            enddo
            m=n
            mp1=m+1
            mm1=m-1
            vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
            vwh_coef(m,n, 1,-1)=0.d0
            vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
            vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
            vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
            vwh_coef(m,n, 0,-1)=0.d0
            vwh_coef(m,n, 1, 0)=-0.d0
            vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
            vwh_coef(m,n, 0, 0)=-fnorm2*m
         enddo
         end subroutine init

      end module numconstants
!
!  special function for the multiple sphere problem
!
      module specialfuncs
      implicit none
      contains

         subroutine timewrite(iunit,char1,time)
         use intrinsics
         implicit none
         integer :: iunit
         real(8) :: time,time2
         character(*) :: char1
         if(time.gt.3600.d0) then
            time2=time/3600.d0
            write(iunit,'(a,f9.3,'' hours'')') char1,time2
         elseif(time.gt.60.d0) then
            time2=time/60.d0
            write(iunit,'(a,f9.2,'' min'')') char1,time2
         else
            write(iunit,'(a,f9.2,'' sec'')') char1,time
         endif
         call flush(iunit)
         end subroutine timewrite
!
!  ricatti-bessel function psi(n), real argument
!
         subroutine ricbessel(n,ds,eps,nmax,psi)
         implicit none
         integer :: n,nmax,ns,i
         real(8) :: ds,dns,sn,psi(0:n),psit,ds2,sum,eps,err
         if(int(ds).lt.n) then
            ns=nint(ds+4.*(ds**.3333d0)+17)
            ns=max(n+10,ns)
            dns=0.d0
            do i=ns-1,n,-1
               sn=dble(i+1)/ds
               dns=sn-1.d0/(dns+sn)
            enddo
            psi(n)=dns
            psi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
            do i=n-2,1,-1
               sn=dble(i+1)/ds
               psi(i)=sn-1.d0/(psi(i+1)+sn)
            enddo
            psit=dsin(ds)
            psi(0)=psit
            ds2=ds*ds
            sum=psit*psit/ds2
            do i=1,n
               psit=psit/(dble(i)/ds+psi(i))
               sum=sum+dble(i+i+1)*psit*psit/ds2
               err=dabs(1.d0-sum)
               psi(i)=psit
               if(err.lt.eps) then
                  nmax=i
                  return
               endif
            enddo
            nmax=n
         else
            psi(0)=dsin(ds)
            psi(1)=psi(0)/ds-dcos(ds)
            do i=1,n-1
               sn=dble(i+i+1)/ds
               psi(i+1)=sn*psi(i)-psi(i-1)
            enddo
            nmax=n
         endif
         end subroutine ricbessel
!
!  ricatti-hankel function xi(n), real argument
!
!
!  last revised: 15 January 2011
!
         subroutine richankel(n,ds,xi)
         implicit none
         integer :: n,i,ns
         real(8) :: ds,dns,sn,chi0,chi1,chi2,psi,psi0,psi1
         complex(8) :: xi(0:n)
         if(int(ds).lt.n) then
            ns=nint(ds+4.*(ds**.3333)+17)
            ns=max(n+10,ns)
            dns=0.d0
            do i=ns-1,n,-1
               sn=dble(i+1)/ds
               dns=sn-1.d0/(dns+sn)
            enddo
            xi(n)=dns
            xi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
            do i=n-2,1,-1
               sn=dble(i+1)/ds
               xi(i)=sn-1.d0/(xi(i+1)+sn)
            enddo
            chi0=-dcos(ds)
            psi=dsin(ds)
            chi1=chi0/ds-psi
            xi(0)=dcmplx(psi,chi0)
            do i=1,n
               chi2=dble(i+i+1)/ds*chi1-chi0
               psi=psi/(dble(i)/ds+xi(i))
               xi(i)=dcmplx(psi,chi1)
               chi0=chi1
               chi1=chi2
            enddo
            return
         else
            chi0=-dcos(ds)
            psi0=dsin(ds)
            chi1=chi0/ds-psi0
            psi1=psi0/ds+chi0
            xi(0)=dcmplx(psi0,chi0)
            xi(1)=dcmplx(psi1,chi1)
            do i=1,n-1
               sn=dble(i+i+1)/ds
               xi(i+1)=sn*xi(i)-xi(i-1)
            enddo
            return
         endif
         end subroutine richankel
!
!  ricatti-bessel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!
         subroutine cricbessel(n,ds,psi)
         implicit none
         integer :: n,i
         complex(8) :: ds,psi(0:n),chi(0:n)
         call cspherebessel(n,ds,psi,chi)
         do i=0,n
            psi(i)=psi(i)*ds
         enddo
         return
         end subroutine cricbessel
!
!  ricatti-hankel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!  7 october 2011: forces upwards recurrence for real argument ds
!
         subroutine crichankel(n,ds,xi)
         implicit none
         integer :: n,i,i1
         complex(8) :: ds,psi(0:n),chi(0:n),xi(0:n),ci
         data ci/(0.d0,1.d0)/
         xi(0)=-ci*cdexp(ci*ds)
         xi(1)=-cdexp(ci*ds)*(ci+ds)/ds
         if(dimag(ds).eq.0.d0) then
            do i=1,n-1
               i1=i+1
               xi(i1)=dble(i+i1)/ds*xi(i)-xi(i-1)
            enddo
            return
         endif
         if(cdabs(xi(0)).lt.1.d-10) then
            do i=1,n-1
               i1=i+1
               xi(i1)=dble(i+i1)/ds*xi(i)-xi(i-1)
            enddo
            return
         else
            call cspherebessel(n,ds,psi,chi)
            do i=1,n-1
               i1=i+1
               xi(i1)=(psi(i1)+ci*chi(i1))*ds
            enddo
            return
         endif
         end subroutine crichankel
!
!     ==========================================================
!     Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!              for a complex argument
!     Input :  z --- Complex argument
!              n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!     Output:  CSJ(n) --- jn(z)
!              CSY(n) --- yn(z)
!              NM --- Highest order computed
!     Routines called:
!              MSTA1 and MSTA2 for computing the starting
!              point for backward recurrence
!     ==========================================================
!
!    obtained from, and copywrited by, Jian-Ming Jin
!    http://jin.ece.uiuc.edu/
!
!
!  last revised: 15 January 2011
!
         subroutine cspherebessel(n,z,csj,csy)
         implicit none
         integer :: n,nm,k,m
         real(8) :: a0
         complex(8) :: z,csj(0:n),csy(0:n),csa,csb,cs,cf0,cf1,cf
         a0=cdabs(z)
         nm=n
         if (a0.lt.1.0d-60) then
            csj=(0.d0,0.d0)
            csy=(-1.d300,0.d0)
            csy(0)=(1.d0,0.d0)
            return
         endif
         csj=(0.d0,0.d0)
         csj(0)=cdsin(z)/z
         csj(1)=(csj(0)-cdcos(z))/z
         if (n.ge.2) then
            csa=csj(0)
            csb=csj(1)
            m=msta1(a0,200)
            if (m.lt.n) then
               nm=m
            else
               m=msta2(a0,n,15)
            endif
            cf0=0.0d0
            cf1=1.0d0-100
            do k=m,0,-1
               cf=(2.0d0*k+3.0d0)*cf1/z-cf0
               if (k.le.nm) csj(k)=cf
               cf0=cf1
               cf1=cf
            enddo
            if (cdabs(csa).gt.cdabs(csb)) cs=csa/cf
            if (cdabs(csa).le.cdabs(csb)) cs=csb/cf0
            do k=0,min(nm,n)
               csj(k)=cs*csj(k)
            enddo
         endif
         csy=(1.d200,0.d0)
         csy(0)=-cdcos(z)/z
         csy(1)=(csy(0)-cdsin(z))/z
         do k=2,min(nm,n)
            if (cdabs(csj(k-1)).gt.cdabs(csj(k-2))) then
               csy(k)=(csj(k)*csy(k-1)-1.0d0/(z*z))/csj(k-1)
            else
               csy(k)=(csj(k)*csy(k-2)-(2.0d0*k-1.0d0)/z**3)/csj(k-2)
            endif
         enddo
         end subroutine cspherebessel
!
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that the magnitude of
!              Jn(x) at that point is about 10^(-MP)
!     Input :  x     --- Argument of Jn(x)
!              MP    --- Value of magnitude
!     Output:  MSTA1 --- Starting point
!     ===================================================
!
!
!  last revised: 15 January 2011
!
         integer function msta1(x,mp)
         implicit none
         integer :: mp,n0,n1,it,nn
         real(8) :: x, a0,f1,f,f0
         a0=dabs(x)
         n0=int(1.1*a0)+1
         f0=envj(n0,a0)-mp
         n1=n0+5
         f1=envj(n1,a0)-mp
         do it=1,20
            nn=n1-(n1-n0)/(1.0d0-f0/f1)
            f=envj(nn,a0)-mp
            if(abs(nn-n1).lt.1) exit
            n0=n1
            f0=f1
            n1=nn
            f1=f
         enddo
         msta1=nn
         end function msta1
!
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that all Jn(x) has MP
!              significant digits
!     Input :  x  --- Argument of Jn(x)
!              n  --- Order of Jn(x)
!              MP --- Significant digit
!     Output:  MSTA2 --- Starting point
!     ===================================================
!
!
!  last revised: 15 January 2011
!
         integer function msta2(x,n,mp)
         implicit none
         integer :: n,mp,n0,n1,it,nn
         real(8) :: x,a0,hmp,ejn,obj,f0,f1,f
         a0=dabs(x)
         hmp=0.5d0*dble(mp)
         ejn=envj(n,a0)
         if (ejn.le.hmp) then
            obj=mp
            n0=int(1.1*a0)
         else
            obj=hmp+ejn
            n0=n
         endif
         f0=envj(n0,a0)-obj
         n1=n0+5
         f1=envj(n1,a0)-obj
         do it=1,20
            nn=n1-(n1-n0)/(1.0d0-f0/f1)
            f=envj(nn,a0)-obj
            if (abs(nn-n1).lt.1) exit
            n0=n1
            f0=f1
            n1=nn
            f1=f
         enddo
         msta2=nn+10
         end function msta2

         real(8) function envj(n,x)
         implicit none
         integer :: n
         real(8) :: x
         n=max(1,abs(n))
         envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
         end function envj
!
!    vector coupling coefficients vc(w) = C(m,n|k,l|m+k,w), w = |n-l|,... n+l
!    uses downwards and upwards recurrence
!
!
!  last revised: 15 January 2011
!
         subroutine vcfunc(m,n,k,l,vcn)
         use numconstants
         implicit none
         integer :: m,n,k,l,wmax,wmin,w,mk
         real(8) :: vcn(0:n+l),t1,t2,t3,vcmax,vctest,rat
         vcn=0.d0
         wmax=n+l
         wmin=max(abs(n-l),abs(m+k))
         vcn(wmax)=bcof(n+m,l+k)*bcof(n-m,l-k)/bcof(n+n,l+l)
         if(wmin.eq.wmax) return
         vcn(wmax-1)=vcn(wmax)*(l*m-k*n)*fnr(2*(l+n)-1)/fnr(l)/fnr(n)&
        &  /fnr(n+l+m+k)/fnr(n+l-m-k)
         if(wmin.eq.wmax-1) return
         mk=m+k
         vcmax=abs(vcn(wmax))+abs(vcn(wmax-1))
!
!  a downwards recurrence is used initially
!
         do w=wmax,wmin+2,-1
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)&
        &     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))&
        &    /dble(2*w*(w-1))
            t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1)&
        &     *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3)&
        &     *fnr(2*w-1))
            vcn(w-2)=(t2*vcn(w-1)-vcn(w)/t1)/t3
            if(mod(wmax-w,2).eq.1) then
               vctest=abs(vcn(w-2))+abs(vcn(w-1))
               vcmax=max(vcmax,vctest)
               rat=vctest/vcmax
!
!  if/when the coefficients start to decrease in magnitude, an upwards recurrence takes over
!
               if(rat.lt.0.01d0) exit
            endif
         enddo
         if(w-2.gt.wmin) then
            wmax=w-3
            call vcfuncuprec(m,n,k,l,wmax,vcn)
         endif
         end subroutine vcfunc
!
!  upwards VC coefficient recurrence
!
!
!  last revised: 15 January 2011
!
         subroutine vcfuncuprec(m,n,k,l,wmax,vcn)
         use numconstants
         implicit none
         integer :: m,n,k,l,wmax,wmin,w,mk,nl,m1,n1,l1,k1,w1,w2
         real(8) :: vcn(0:n+l),t1,t2,t3,vc1
         mk=abs(m+k)
         nl=abs(n-l)
         if(nl.ge.mk) then
            w=nl
            if(n.ge.l) then
               m1=m
               n1=n
               l1=l
               k1=k
            else
               m1=k
               n1=l
               k1=m
               l1=n
            endif
            vc1=(-1)**(k1+l1)*bcof(l1+k1,w-m1-k1) &
               *bcof(l1-k1,w+m1+k1)/bcof(l1+l1,w+w+1)
         else
            w=mk
            if(m+k.ge.0) then
               vc1=(-1)**(n+m)*bcof(n-l+w,l-k)*bcof(l-n+w,n-m) &
                  /bcof(w+w+1,n+l-w)
            else
               vc1=(-1)**(l+k)*bcof(n-l+w,l+k)*bcof(l-n+w,n+m) &
                 /bcof(w+w+1,n+l-w)
            endif
         endif
         w1=w
         vcn(w)=vc1
         w=w1+1
         mk=m+k
         w2=min(wmax,n+l)
         if(w2.gt.w1) then
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk) &
              *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            if(w1.eq.0) then
               t2=.5*dble(m-k)
            else
               t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1)) &
                 /dble(2*w*(w-1))
            endif
            vcn(w)=t1*t2*vcn(w1)
         endif
         do w=w1+2,w2
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk) &
              *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1)) &
             /dble(2*w*(w-1))
            t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1) &
              *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3) &
              *fnr(2*w-1))
            vcn(w)=t1*(t2*vcn(w-1)-t3*vcn(w-2))
         enddo
         end subroutine vcfuncuprec
!
!  Normalized associated legendre functions
!
!
!  last revised: 15 January 2011
!
         subroutine normalizedlegendre(cbe,mmax,nmax,dc)
         use numconstants
         implicit none
         integer :: nmax,mmax,m,n,np1,nm1,im
         real(8) :: dc(-mmax:mmax,0:nmax),cbe,sbe
         sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
         dc=0.d0
         do m=0,mmax
            dc(m,m)=(-1)**m*(0.5d0*sbe)**m*bcof(m,m)
            if(m.eq.nmax) exit
            dc(m,m+1)=fnr(m+m+1)*cbe*dc(m,m)
            do n=m+1,nmax-1
               dc(m,n+1)=(-fnr(n-m)*fnr(n+m)*dc(m,n-1)+dble(n+n+1)*cbe*dc(m,n)) &
                         /(fnr(n+1-m)*fnr(n+1+m))
            enddo
         enddo
         do m=1,mmax
            im=(-1)**m
            do n=m,nmax
               dc(-m,n)=im*dc(m,n)
            enddo
         enddo
         end subroutine normalizedlegendre
!
!  Generalized spherical functions
!
!  dc(m,n*(n+1)+k)=(-1)^(m + k)((n - k)!(n + k)!/(n - m)!/(n + m)!)^(1/2)
!  ((1 + x)/2)^((m + k)/2)((1 - x)/2)^((k - m)/2)JacobiP[n - k, k - m, k + m, x]
!
!  for |m| <= kmax, n=0,1,...nmax, |k| <= n
!
!
!  last revised: 15 January 2011
!
         subroutine rotcoef(cbe,kmax,nmax,dc)
         use numconstants
         implicit none
         integer :: kmax,nmax,k,m,in,n,knmax,nn1,kn,im,m1
         real(8) :: cbe,sbe,dc(-kmax:kmax,0:nmax*(nmax+2)),cbe2,sbe2,dk0(-nmax-1:nmax+1),&
                    dk01(-nmax-1:nmax+1),sben,dkt,fmn,dkm0,dkm1,dkn1
         sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         in=1
         dk0(0)=1.d0
         sben=1.d0
         dc(0,0)=1.d0
         dk01(0)=0.
         do n=1,nmax
            knmax=min(n,kmax)
            nn1=n*(n+1)
            in=-in
            sben=sben*sbe/2.d0
            dk0(n)=in*sben*bcof(n,n)
            dk0(-n)=in*dk0(n)
            dk01(n)=0.
            dk01(-n)=0.
            dc(0,nn1+n)=dk0(n)
            dc(0,nn1-n)=dk0(-n)
            do k=-n+1,n-1
               kn=nn1+k
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)&
                     /(fnr(n+k)*fnr(n-k))
               dc(0,kn)=dk0(k)
            enddo
            im=1
            do m=1,knmax
               im=-im
               fmn=1.d0/fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.
               do k=-n,n
                  kn=nn1+k
                  dkm1=dkm0
                  dkm0=dc(m1,kn)
                  if(k.eq.n) then
                     dkn1=0.
                  else
                     dkn1=dc(m1,kn+1)
                  endif
                  dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                          -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1  &
                          -dble(k)*sbe*dc(m1,kn))*fmn
                  dc(-m,nn1-k)=dc(m,kn)*(-1)**(k)*im
               enddo
            enddo
         enddo
         end subroutine rotcoef

         subroutine rotcoefvecarg(narg,cbe,kmax,nmax,dc)
         use numconstants
         implicit none
         integer :: kmax,nmax,k,m,in,n,knmax,nn1,kn,im,m1,narg
         real(8) :: cbe(narg),sbe(narg),dc(-kmax:kmax,0:nmax*(nmax+2),narg), &
                    cbe2(narg),sbe2(narg),dk0(-nmax-1:nmax+1,narg),&
                    dk01(-nmax-1:nmax+1,narg),sben(narg),dkt(narg), &
                    fmn,dkm0(narg),dkm1(narg),dkn1(narg)
         sbe=sqrt((1.d0+cbe)*(1.d0-cbe))
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         in=1
         dk0(0,:)=1.d0
         sben=1.d0
         dc(0,0,:)=1.d0
         dk01(0,:)=0.
         do n=1,nmax
            knmax=min(n,kmax)
            nn1=n*(n+1)
            in=-in
            sben=sben*sbe/2.d0
            dk0(n,:)=in*sben(:)*bcof(n,n)
            dk0(-n,:)=in*dk0(n,:)
            dk01(n,:)=0.
            dk01(-n,:)=0.
            dc(0,nn1+n,:)=dk0(n,:)
            dc(0,nn1-n,:)=dk0(-n,:)
            do k=-n+1,n-1
               kn=nn1+k
               dkt(:)=dk01(k,:)
               dk01(k,:)=dk0(k,:)
               dk0(k,:)=(cbe(:)*dble(n+n-1)*dk01(k,:)-fnr(n-k-1)*fnr(n+k-1)*dkt(:)) &
                     /(fnr(n+k)*fnr(n-k))
               dc(0,kn,:)=dk0(k,:)
            enddo
            im=1
            do m=1,knmax
               im=-im
               fmn=1.d0/fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.
               do k=-n,n
                  kn=nn1+k
                  dkm1=dkm0
                  dkm0(:)=dc(m1,kn,:)
                  if(k.eq.n) then
                     dkn1=0.
                  else
                     dkn1(:)=dc(m1,kn+1,:)
                  endif
                  dc(m,kn,:)=(fnr(n+k)*fnr(n-k+1)*cbe2(:)*dkm1(:) &
                          -fnr(n-k)*fnr(n+k+1)*sbe2(:)*dkn1(:)  &
                          -dble(k)*sbe(:)*dc(m1,kn,:))*fmn
                  dc(-m,nn1-k,:)=dc(m,kn,:)*(-1)**(k)*im
               enddo
            enddo
         enddo
         end subroutine rotcoefvecarg
!
!  tau are the vector spherical harmonic functions, normalized
!
!
!  last revised: 15 January 2011
!
         subroutine taufunc(cb,nmax,tau)
         use numconstants
         implicit none
         integer :: nmax,n,m,p,nn1,mn
         real(8) :: drot(-1:1,0:nmax*(nmax+2)),tau(0:nmax+1,nmax,2),cb,fnm
         call rotcoef(cb,1,nmax,drot)
         do n=1,nmax
            nn1=n*(n+1)
            fnm=sqrt(dble(n+n+1)/2.d0)/4.d0
            do m=-n,-1
               mn=nn1+m
               tau(n+1,-m,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(n+1,-m,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
            do m=0,n
               mn=nn1+m
               tau(m,n,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(m,n,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
         enddo
         end subroutine taufunc
!
! vector spherical harmonic function
! november 2011
!

         subroutine pifunc(cb,ephi,nmax,ndim,pivec)
         use numconstants
         implicit none
         integer :: nmax,n,m,p,nn1,mn,ndim
         real(8) :: drot(-1:1,0:nmax*(nmax+2)),tau(2),cb,fnm
         complex(8) :: pivec(0:ndim+1,ndim,2),ephi,ephim(-nmax:nmax),cin
         call rotcoef(cb,1,nmax,drot)
         ephim(0)=1.d0
         do m=1,nmax
            ephim(m)=ephi*ephim(m-1)
            ephim(-m)=dconjg(ephim(m))
         enddo
         do n=1,nmax
            cin=(0.d0,-1.d0)**(n+1)
            nn1=n*(n+1)
            fnm=sqrt(dble(n+n+1)/2.d0)/4.d0
            do m=-n,-1
               mn=nn1+m
               tau(1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(2)=-fnm*(drot(-1,mn)+drot(1,mn))
               pivec(n+1,-m,1)=cin*tau(1)*ephim(m)
               pivec(n+1,-m,2)=cin*tau(2)*ephim(m)
            enddo
            do m=0,n
               mn=nn1+m
               tau(1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(2)=-fnm*(drot(-1,mn)+drot(1,mn))
               pivec(m,n,1)=cin*tau(1)*ephim(m)
               pivec(m,n,2)=cin*tau(2)*ephim(m)
            enddo
         enddo
         end subroutine pifunc
!
!  regular vswf expansion coefficients for a plane wave.
!  alpha, beta: incident azimuth and polar angles.
!
!
!  last revised: 15 January 2011
!
         subroutine planewavecoef(alpha,beta,nodr,pmnp0)
         use numconstants
         implicit none
         integer :: nodr,m,n,p,k,ierr
         real(8) :: alpha,beta,cb,sb,ca,sa
         real(8), allocatable :: tau(:,:,:)
         complex(8) :: ealpha,ci,cin
         complex(8), allocatable :: ealpham(:)
         complex(8) :: pmnp0(0:nodr+1,nodr,2,2)
         data ci/(0.d0,1.d0)/
         call init(nodr)
         allocate(ealpham(-nodr:nodr))
         allocate(tau(0:nodr+1,nodr,2))
         cb=cos(beta)
         sb=sqrt((1.d0-cb)*(1.d0+cb))
         ca=cos(alpha)
         sa=sin(alpha)
         ealpha=dcmplx(ca,sa)
         call taufunc(cb,nodr,tau)
         call ephicoef(ealpha,nodr,ealpham)
         do n=1,nodr
            cin=4.d0*ci**(n+1)
            do p=1,2
               do m=-n,-1
                  pmnp0(n+1,-m,p,1)=-cin*tau(n+1,-m,p)*ealpham(-m)
                  pmnp0(n+1,-m,p,2)=ci*cin*tau(n+1,-m,3-p)*ealpham(-m)
               enddo
               do m=0,n
                  pmnp0(m,n,p,1)=-cin*tau(m,n,p)*ealpham(-m)
                  pmnp0(m,n,p,2)=ci*cin*tau(m,n,3-p)*ealpham(-m)
               enddo
            enddo
         enddo
         deallocate(ealpham,tau)
         end subroutine planewavecoef
!
!  regular vswf expansion coefficients for a gaussian beam, localized approximation.
!  cbeam = 1/(k omega)
!
!
!  last revised: 15 January 2011
!
         subroutine gaussianbeamcoef(alpha,beta,cbeam,nodr,pmnp0)
         use numconstants
         implicit none
         integer :: nodr,m,n,p,k,ierr
         real(8) :: alpha,beta,cbeam,gbn
         complex(8) :: pmnp0(0:nodr+1,nodr,2,2)
         call planewavecoef(alpha,beta,nodr,pmnp0)
         do n=1,nodr
            gbn=dexp(-((dble(n)+.5d0)*cbeam)**2.)
            do p=1,2
               do k=1,2
                  do m=-n,-1
                     pmnp0(n+1,-m,p,k)=pmnp0(n+1,-m,p,k)*gbn
                  enddo
                  do m=0,n
                     pmnp0(m,n,p,k)=pmnp0(m,n,p,k)*gbn
                  enddo
               enddo
            enddo
         enddo
         end subroutine gaussianbeamcoef
!
!  plane wave expansion coefficients at sphere origins.  uses a phase shift.
!
!
!  last revised: 15 January 2011
!
         subroutine sphereplanewavecoef(nsphere,neqns,nodr,nodrmax,alpha,beta,rpos,pmnp)
         implicit none
         integer :: m,n,p,nsphere,i,l,nodr(nsphere),nblk,nboff,nodrmax,neqns,k
         real(8) :: alpha,beta,cb,sb,ca,sa,rpos(3,nsphere)
         complex(8) :: ci,phasefac, pmnp(neqns,2)
         complex(8) :: pmnp0(0:nodrmax+1,nodrmax,2,2)
         data ci/(0.d0,1.d0)/
         call planewavecoef(alpha,beta,nodrmax,pmnp0)
         cb=cos(beta)
         sb=sqrt((1.d0-cb)*(1.d0+cb))
         ca=cos(alpha)
         sa=sin(alpha)
         l=0
         do i=1,nsphere
            phasefac=cdexp(ci*((ca*rpos(1,i)+sa*rpos(2,i))*sb+rpos(3,i)*cb))
            do p=1,2
               do n=1,nodr(i)
                  do m=0,nodr(i)+1
                     l=l+1
                     do k=1,2
                        pmnp(l,k)=phasefac*pmnp0(m,n,p,k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         end subroutine sphereplanewavecoef
!
! this computes the normalized translation coefficients for an
! axial translation of positive distance r.  For itype=1 or 3, the translation
! uses the spherical Bessel or Hankel functions as a basis function,
! respectively.    They are related to the coefficients appearing in
! M&M JOSA 96 by
!
! J^{ij}_{mnp mlq} = (E_{ml}/E_{mn})^(1/2) ac(s,n,l*(l+1)+m)
!
! where
!
!   E_{mn} = n(n+1)(n+m)!/((2n+1)(n-m)!)
!   s=mod(p+q,2)+1 (i.e., s=1 for the A coefficient, =2 for the B
!   coefficient)
!
!  The calculation procedure is based on the derivation
!  of the addition theorem for vector harmonics, appearing in
!  Fuller and Mackowski, proc. Light Scattering by Nonspherical
!  Particles, NASA/GISS Sept. 1998.
!
!  revised: 10 october 2011: used F90 vector arithmetic and precalculation
!           of various constants.
!
         subroutine axialtrancoef(itype,r,ri,nmax,lmax,ac)
         use numconstants
         implicit none
         integer :: itype,nmax,lmax,n,l,m,p,w,n21,ll1,nlmin,lblk,wmin,wmax,ml
         integer :: iadd,nlmax
         integer, save :: nlmax0
         real(8) :: r
         complex(8) :: ri,ci,z,xi(0:nmax+lmax)
         complex(8) :: ac(nmax,lmax*(lmax+3)/2,2)
         data ci,nlmax0/(0.d0,1.d0),0/
         nlmax=max(nmax,lmax)
         if(nlmax.gt.nlmax0) then
            nlmax0=nlmax
            call axialtrancoefinit(nlmax)
         endif
         if(r.eq.0.d0) then
            ac=(0.d0,0.d0)
            if(itype.ne.1) return
            do m=0,min(nmax,lmax)
               do n=max(1,m),min(nmax,lmax)
                  iadd=atcadd(m,n,lmax)
                  ac(n,iadd,l)=1.
               enddo
            enddo
            return
         endif
         z=r*ri
         if(itype.eq.1) then
            call cricbessel(nmax+lmax,z,xi)
         else
            call crichankel(nmax+lmax,z,xi)
         endif
         xi=xi/z
         do n=1,nmax
            do l=1,lmax
               wmin=abs(n-l)
               wmax=n+l
               do m=0,min(n,l)
                  iadd=atcadd(m,l,lmax)
                  ml=l*(l+1)/2+m
                  ac(n,iadd,1)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2))
                  ac(n,iadd,2)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2))
               enddo
            enddo
         enddo
         end subroutine axialtrancoef
!
!  axial translation coefficients calculated by the diamond recurrence formula
!  new: 10 october 2011
!
         subroutine axialtrancoefrecurrence(itype,r,ri,nmax,lmax,ac)
         use numconstants
         implicit none
         integer :: itype,nmax,lmax,n,l,m,p,q,w,n21,ll1,nlmin,lblk, &
                    wmin,wmax,ml,m1,np1,nm1,iaddp1,iaddm1,lm1,lp1
         integer :: iadd,nlmax,iadd0,iadd1
         integer, save :: nlmax0
         real(8) :: r,fnp1,fn,fnm1,flp1,fl,flm1
         complex(8) :: ri,ci,z,xi(0:nmax+lmax)
         complex(8) :: ac(nmax,lmax*(lmax+3)/2,2)
         data ci,nlmax0/(0.d0,1.d0),0/
         nlmax=max(nmax,lmax)
         nlmin=min(nmax,lmax)
         if(nlmax.gt.nlmax0) then
            nlmax0=nlmax
            call axialtrancoefinit(nlmax)
         endif

         if(r.eq.0.d0) then
            ac=(0.d0,0.d0)
            if(itype.ne.1) return
            do m=0,nlmin
               m1=max(1,m)
               do n=m1,nlmin
                  iadd=atcadd(m,n,lmax)
                  ac(n,iadd,l)=1.
               enddo
            enddo
            return
         endif
         z=r*ri
         if(itype.eq.1) then
            call cricbessel(nmax+lmax,z,xi)
         else
            call crichankel(nmax+lmax,z,xi)
         endif
         xi=xi/z

         lm1=lmax-1
         do m=0,nlmin
            m1=max(1,abs(m))
            lp1=m1+1
            iadd0=atcadd(m,m1,lmax)
            iadd1=atcadd(m,lmax,lmax)
            iaddp1=iadd0+1
            iaddm1=iadd1-1

            iadd=iadd0-1
            n=m1
            do l=m1,lmax
               wmin=abs(n-l)
               wmax=n+l
               iadd=iadd+1
               ml=l*(l+1)/2+m
               ac(n,iadd,1)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2))
               ac(n,iadd,2)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2))
            enddo
            l=lmax
            iadd=iadd1
            ml=l*(l+1)/2+m
            do n=m1+1,nmax
               wmin=abs(n-l)
               wmax=n+l
               ac(n,iadd,1)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2))
               ac(n,iadd,2)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2))
            enddo
            if(m1.eq.nlmin) cycle

            do n=m1,nmax-1
               np1=n+1
               nm1=n-1
               do p=1,2
                  q=3-p
                  ac(np1,iadd0:iaddm1,p)= &
                    - ac(n,iaddp1:iadd1,p)*fnp1_const(m,m1:lm1) &
                    + (fn_const(m,m1:lm1)-fn_const(m,n))*ci*ac(n,iadd0:iaddm1,q)
                  ac(np1,iaddp1:iaddm1,p)=ac(np1,iaddp1:iaddm1,p) &
                    + ac(n,iadd0:iadd1-2,p)*fnm1_const(m,lp1:lm1)
                  if(n.gt.m1) then
                     ac(np1,iadd0:iaddm1,p)=ac(np1,iadd0:iaddm1,p) &
                       + ac(nm1,iadd0:iaddm1,p)*fnm1_const(m,n)
                  endif
                  ac(np1,iadd0:iaddm1,p)=ac(np1,iadd0:iaddm1,p)/fnp1_const(m,n)
               enddo
            enddo
         enddo
         end subroutine axialtrancoefrecurrence
!
!  constants for translation coefficient calculation
!
         subroutine axialtrancoefinit(nmax)
         use numconstants
         implicit none
         integer :: nmax,m,n,l,w,n21,ml,ll1,wmin,wmax,nlmin,lp1,lm1
         real(8) :: c1,c2,vc1(0:2*nmax),vc2(0:2*nmax),alnw
         complex(8) :: ci,inlw
         data ci/(0.d0,1.d0)/
         if(allocated(vcc_const)) deallocate(vcc_const,fnm1_const,fn_const,fnp1_const)
         allocate(vcc_const(nmax,nmax*(nmax+1)/2+nmax,0:2*nmax),fnm1_const(0:nmax,nmax), &
                  fn_const(0:nmax,nmax),fnp1_const(0:nmax,nmax))
         do n=1,nmax
            n21=n+n+1
            do l=1,nmax
               c1=fnr(n21)*fnr(l+l+1)
               ll1=l*(l+1)/2
               call vcfunc(-1,n,1,l,vc2)
               wmin=abs(n-l)
               wmax=n+l
               nlmin=min(l,n)
               do m=0,nlmin
                  ml=ll1+m
                  c2=-c1*(-1)**m
                  call vcfunc(-m,n,m,l,vc1)
                  do w=wmin,wmax
                     inlw=ci**(n-l+w)
                     vcc_const(n,ml,w)=c2*vc1(w)*vc2(w)*(dble(inlw)+dimag(inlw))
                  enddo
               enddo
            enddo
         enddo
         fnm1_const=0.
         fn_const=0.
         fnp1_const=0.
         do m=0,nmax
            do l=max(1,m),nmax
               lp1=l+1
               lm1=l-1
               fnm1_const(m,l)=fnr(lm1)*fnr(lp1)*fnr(l-m)*fnr(l+m)/fnr(lm1+l)/fnr(l+lp1)/dble(l)
               fn_const(m,l)=dble(m)/dble(l)/dble(lp1)
               fnp1_const(m,l)=fnr(l)*fnr(l+2)*fnr(lp1-m)*fnr(lp1+m)/fnr(l+lp1)/fnr(l+l+3)/dble(lp1)
            enddo
         enddo
         end subroutine axialtrancoefinit
!
!  test to determine convergence of regular vswf addition theorem for max. order lmax
!  and translation distance r w/ refractive index ri.
!
!
!  last revised: 15 January 2011
!
         subroutine tranordertest(r,ri,lmax,eps,nmax)
         use numconstants
         implicit none
         integer :: itype,nmax,lmax,n,l,m,p,w,n21,ll1,nlmin,lblk,wmin,wmax
         integer, parameter :: nlim=200
         integer :: iadd
         real(8) :: r,alnw,sum,eps
         real(8) :: vc1(0:nlim+lmax)
         complex(8) :: ri,ci,z,a,b,c
         complex(8) :: xi(0:nlim+lmax)
         data ci/(0.d0,1.d0)/
         if(r.eq.0.d0) then
            nmax=lmax
            return
         endif
         z=r*ri
         sum=0.d0
         do n=1,nlim
            call init(n+lmax)
            call cricbessel(n+lmax,z,xi)
            do l=0,n+lmax
               xi(l)=xi(l)/z*ci**l
            enddo
            n21=n+n+1
            l=lmax
            c=fnr(n21)*fnr(l+l+1)*ci**(n-l)
            call vcfunc(-1,n,1,l,vc1)
            wmin=abs(n-l)
            wmax=n+l
            m=1
            a=0.
            b=0.
            do w=wmin,wmax
               alnw=vc1(w)*vc1(w)
               if(mod(n+l+w,2).eq.0) then
                  a=a+alnw*xi(w)
               else
                  b=b+alnw*xi(w)
               endif
            enddo
            a=c*a
            b=c*b
            sum=sum+a*conjg(a)+b*conjg(b)
            if(abs(1.d0-sum).lt.eps) exit
         enddo
         nmax=min(n,nlim)
         nmax=max(nmax,lmax)
         end subroutine tranordertest
!
!  address for axial translation coefficient
!
!
!  last revised: 15 January 2011
!
         integer function atcadd(m,n,ntot)
         implicit none
         integer :: m,n,ntot
         atcadd=n-ntot+(max(1,m)*(1+2*ntot-max(1,m)))/2+ntot*min(1,m)
         end function atcadd
!
!   gentrancoef: calculates the vwh translation coefficients for
!   a general translation from one origin to another
!
!   input: itype: integer, =1, regular, =3, outgoing type harmonics
!          xptran: real, dim 3 vector: x,y,z components of translation, in units
!                   of 1/k
!          ri: complex, refractive index of medium
!          nrow0,nrow1,ncol0,ncol1: integer, starting and stopping row and column order
!          iaddrow0,iaddcol0: address offset for row and column order (see below)
!   output: ac(p,mn,kl): complex translation matrix.  calculated for mode p=1,2 (A or B type),
!           order n=nrow0,nrow1, degree m=-n,n
!           order l=ncol0,ncol1, degree k=-n,n
!           address is given by
!           mn=m+n*(n+1)-(nrow0-1)*(nrow0+1)+iaddrow0
!           kl=k+l*(l+1)-(ncol0-1)*(ncol0+1)+iaddcol0
!           that is, if iaddrow0=0 the address is mn=1 for n=nrow0 and m=-n.
!
!
!  last revised: 15 January 2011
!
         subroutine gentrancoef(itype,xptran,ri,nrow0,nrow1,ncol0,ncol1, &
                               iaddrow0,iaddcol0,ac)
         use numconstants
         implicit none
         integer :: itype,nrow0,nrow1,ncol0,ncol1,iaddrow0,iaddcol0,kmax
         integer :: ntot,nblkr0,nblkr1,nblkc0,nblkc1
         integer :: v,vw,w,wmax,wmin,n,l,m,k,p,nn1,ll1,mn,kl,m1m
         real(8) :: vc1(0:nrow1+ncol1),vc2(0:nrow1+ncol1),&
                    xptran(3),r,ct,ct0
         real(8) :: drot(0:0,0:(nrow1+ncol1)*(nrow1+ncol1+2))
         complex(8) :: ri,ci,ephi,ac(2,nrow1*(nrow1+2)-(nrow0-1)*(nrow0+1)-iaddrow0,&
                       ncol1*(ncol1+2)-(ncol0-1)*(ncol0+1)-iaddcol0),&
                       z,c,a,b
         complex(8) :: ephim(-(nrow1+ncol1):nrow1+ncol1),jnc(0:nrow1+ncol1)
         data ci/(0.d0,1.d0)/
         call cartosphere(xptran,r,ct,ephi)
         ntot=nrow1+ncol1
         nblkr0=(nrow0-1)*(nrow0+1)
         nblkr1=nrow1*(nrow1+2)
         nblkc0=(ncol0-1)*(ncol0+1)
         nblkc1=ncol1*(ncol1+2)
         if(r.eq.0.d0) then
            do n=nblkr0+1,nblkr1
               mn=n-nblkr0+iaddrow0
               do l=nblkc0+1,nblkc1
                  kl=l-nblkc0+iaddcol0
                  do p=1,2
                     ac(p,mn,kl)=0.d0
                  enddo
               enddo
               if(n.gt.nblkc0.and.n.le.nblkc1.and.itype.eq.1) then
                  ac(1,mn,n-nblkc0+iaddcol0)=1.d0
               endif
            enddo
            return
         endif
         kmax=0
         ct0=ct
         call rotcoef(ct0,kmax,ntot,drot)
         call ephicoef(ephi,ntot,ephim)
         z=ri*r
         if(itype.eq.1) then
            call cricbessel(ntot,z,jnc)
         else
            call crichankel(ntot,z,jnc)
         endif
         do n=0,ntot
            c=ci**n
            jnc(n)=c*jnc(n)/z
         enddo
         do l=ncol0,ncol1
            ll1=l*(l+1)
            do n=nrow0,nrow1
               nn1=n*(n+1)
               wmax=n+l
               call vcfunc(-1,n,1,l,vc2)
               c=-ci**(n-l)*fnr(n+n+1)*fnr(l+l+1)
               do k=-l,l
                  kl=ll1+k-nblkc0+iaddcol0
                  do m=-n,n
                     m1m=(-1)**m
                     mn=nn1+m-nblkr0+iaddrow0
                     v=k-m
                     call vcfunc(-m,n,k,l,vc1)
                     a=0.
                     b=0.
                     wmin=max(abs(v),abs(n-l))
                     do w=wmax,wmin,-1
                        vw=w*(w+1)+v
                        if(mod(wmax-w,2).eq.0) then
                           a=a+vc1(w)*vc2(w)*jnc(w)*drot(0,vw)
                        else
                           b=b+vc1(w)*vc2(w)*jnc(w)*drot(0,vw)
                        endif
                     enddo
                     ac(1,mn,kl)=a*c*m1m*ephim(v)
                     ac(2,mn,kl)=b*c*m1m*ephim(v)
                  enddo
               enddo
            enddo
         enddo
         return
         end subroutine gentrancoef
!
! cartosphere takes the cartesian point (x,y,z) = xp(1), xp(2), xp(3)
! and converts to polar form: r: radius, ct: cos(theta), ep = exp(i phi)
!
!
!  last revised: 15 January 2011
!
         subroutine cartosphere(xp,r,ct,ep)
         implicit none
         real(8) :: xp(3),r,ct
         complex(8) :: ep
         r=xp(1)*xp(1)+xp(2)*xp(2)+xp(3)*xp(3)
         if(r.eq.0.d0) then
            ct=1.d0
            ep=(1.d0,0.d0)
            return
         endif
         r=sqrt(r)
         ct=xp(3)/r
         if(xp(1).eq.0.d0.and.xp(2).eq.0.d0) then
            ep=(1.d0,0.d0)
         else
            ep=dcmplx(xp(1),xp(2))/sqrt(xp(1)*xp(1)+xp(2)*xp(2))
         endif
         return
         end subroutine cartosphere
!
! ephicoef returns the complex array epm(m) = exp(i m phi) for
! m=-nodr,nodr.   ep =exp(i phi), and epm is dimensioned epm(-nd:nd)
!
!
!  last revised: 15 January 2011
!
         subroutine ephicoef(ep,nodr,epm)
         implicit none
         integer :: nodr,m
         complex(8) :: ep,epm(-nodr:nodr)
         epm(0)=(1.d0,0.d0)
         do m=1,nodr
            epm(m)=ep*epm(m-1)
            epm(-m)=dconjg(epm(m))
         enddo
         return
         end subroutine ephicoef
!
!  test to determine max order of vswf expansion of a plane wave at distance r
!
!
!  last revised: 15 January 2011
!
         subroutine planewavetruncationorder(r,eps,nodr)
         implicit none
         integer :: nodr,n1,n
         real(8) :: r,eps,err
         real(8), allocatable :: jn(:)
         complex(8) :: sum, ci,eir
         data ci/(0.d0,1.d0)/
         n1=max(10,int(3.*r+1))
         allocate(jn(0:n1))
         call ricbessel(n1,r,-1.d0,n1,jn)
         jn(0:n1)=jn(0:n1)/r
         eir=cdexp(-ci*r)
         sum=jn(0)*eir
         do n=1,n1
            sum=sum+ci**n*dble(n+n+1)*jn(n)*eir
            err=cdabs(1.d0-sum)
            if(err.lt.eps) then
               nodr=n
               deallocate(jn)
               return
            endif
         enddo
         nodr=n1
         deallocate(jn)
         end subroutine planewavetruncationorder
!
!  calculates the cartesian components of the vswf at position rpos, in ref. index ri.
!
!
!  original: 15 January 2011
!  revised: 23 February 2011: multiplied by root 2
!
         subroutine vwhcalc(rpos,ri,nodr,itype,vwh)
         use numconstants
         implicit none
         integer :: nodr,itype,m,n,p,nodrp1,nodrm1,nn1,mn,np1,nm1,mp1,mm1,ndim, &
                    nblkp
         integer, save :: nodrmax
         real(8) ::  rpos(3),r,ct,fnorm1,fnorm2
         real(8) pmn(0:0,0:(nodr+1)*(nodr+3))
         complex(8) :: ci,vwh(3,2,1:*),ri,ephi,a,b,a1,b1,z1,a2,b2,z2
         complex(8)  :: a1vec(-nodr:nodr), &
                       b1vec(-nodr:nodr),z1vec(-nodr:nodr),a2vec(-nodr:nodr), &
                       b2vec(-nodr:nodr),z2vec(-nodr:nodr)
         complex(8) :: umn(-nodr-2:nodr+2,0:nodr+1), hn(0:nodr+1), ephim(-nodr-1:nodr+1)
         data ci,nodrmax/(0.d0,1.d0),0/
         if(nodr.gt.nodrmax) then
            nodrmax=nodr
            call init(nodr+2)
         endif
         call cartosphere(rpos,r,ct,ephi)
         if(r.le.1.d-4) then
            vwh(:,:,1:nodr*(nodr+1))=(0.d0,0.d0)
            if(itype.eq.3) return
            vwh(1,1,1)=.5d0*fnr(2)/fnr(3)
            vwh(2,1,1)=-.5d0*ci*fnr(2)/fnr(3)
            vwh(3,1,2)=1.d0*fnr(2)/fnr(6)
            vwh(1,1,3)=-.5d0*fnr(2)/fnr(3)
            vwh(2,1,3)=-.5d0*ci*fnr(2)/fnr(3)
            return
         endif
         nodrp1=nodr+1
         nodrm1=nodr-1
         a=ri*r
         if(itype.eq.1) then
            call cricbessel(nodrp1,a,hn)
         else
            call crichankel(nodrp1,a,hn)
         endif
         hn(0:nodrp1)=hn(0:nodrp1)/a
         call rotcoef(ct,0,nodrp1,pmn)
         call ephicoef(ephi,nodrp1,ephim)
         umn=0.d0
         umn(0,0)=hn(0)*fnr(2)
         do n=1,nodrp1
            nn1=n*(n+1)
            umn(-n:n,n)=fnr(2)*pmn(0,nn1-n:nn1+n)*ephim(-n:n)*hn(n)
            umn(-n-1,n)=0.d0
            umn(n+1,n)=0.d0
         enddo
         do n=1,nodr
            nn1=n*(n+1)
            np1=n+1
            nm1=n-1
            a1vec(-n:n)=vwh_coef(-n:n,n,1,1)*umn(-nm1:np1,np1) &
               +vwh_coef(-n:n,n,1,-1)*umn(-nm1:np1,nm1)
            b1vec(-n:n)=vwh_coef(-n:n,n,-1,1)*umn(-np1:nm1,np1) &
               +vwh_coef(-n:n,n,-1,-1)*umn(-np1:nm1,nm1)
            z1vec(-n:n)=vwh_coef(-n:n,n,0,1)*umn(-n:n,np1) &
               +vwh_coef(-n:n,n,0,-1)*umn(-n:n,nm1)
            a2vec(-n:n)=vwh_coef(-n:n,n,1,0)*umn(-nm1:np1,n)
            b2vec(-n:n)=vwh_coef(-n:n,n,-1,0)*umn(-np1:nm1,n)
            z2vec(-n:n)=vwh_coef(-n:n,n,0,0)*umn(-n:n,n)
            vwh(1,1,nn1-n:nn1+n)=-0.5d0*(a1vec(-n:n)+b1vec(-n:n))
            vwh(2,1,nn1-n:nn1+n)=-0.5d0*ci*(-a1vec(-n:n)+b1vec(-n:n))
            vwh(3,1,nn1-n:nn1+n)=-z1vec(-n:n)
            vwh(1,2,nn1-n:nn1+n)=-0.5d0*ci*(a2vec(-n:n)+b2vec(-n:n))
            vwh(2,2,nn1-n:nn1+n)=-0.5d0*(a2vec(-n:n)-b2vec(-n:n))
            vwh(3,2,nn1-n:nn1+n)=-ci*z2vec(-n:n)
         enddo
         end subroutine vwhcalc
!
!  svwf calculation for an axial translation
!
!
!  original: 15 January 2011
!  revised: 23 February 2011: multiplied by root 2
!
         subroutine vwhaxialcalc(rpos,ri,nodr,itype,vwh)
         use numconstants
         implicit none
         integer :: nodr,itype,m,n,p,nodrp1,nodrm1,nn1,mn,np1,nm1,mp1,mm1,ndim, &
                    nblkp
         integer, save :: nodrmax
         real(8) ::  rpos(3),r,ct
         real(8) pmn(-2:2,0:nodr+1)
         complex(8) :: ci,vwh(3,2,2,1:nodr),ri,ephi,a,b,a1,b1,z1,a2,b2,z2
         complex(8) :: umn(-2:2,0:nodr+1), hn(0:nodr+1), ephim(-2:2)
         data ci,nodrmax/(0.d0,1.d0),0/
         if(nodr.gt.nodrmax) then
            nodrmax=nodr
            call init(nodr+2)
         endif
         call cartosphere(rpos,r,ct,ephi)
         if(r.le.1.d-4) then
            vwh(:,:,:,1:nodr)=(0.d0,0.d0)
            if(itype.eq.3) return
            vwh(1,1,1,1)=.5d0*fnr(2)/fnr(3)
            vwh(2,1,1,1)=-.5d0*ci*fnr(2)/fnr(3)
            vwh(1,1,2,1)=-.5d0*fnr(2)/fnr(3)
            vwh(2,1,2,1)=-.5d0*ci*fnr(2)/fnr(3)
            return
         endif
         nodrp1=nodr+1
         nodrm1=nodr-1
         a=ri*r
         if(itype.eq.1) then
            call cricbessel(nodrp1,a,hn)
         else
            call crichankel(nodrp1,a,hn)
         endif
         hn(0:nodrp1)=hn(0:nodrp1)/a
         call normalizedlegendre(ct,2,nodrp1,pmn)
         call ephicoef(ephi,2,ephim)
         umn(-2:2,0:nodrp1)=0.d0
         umn(0,0)=hn(0)*fnr(2)
         do n=1,nodrp1
            p=min(n,2)
            do m=-p,p
               umn(m,n)=fnr(2)*pmn(m,n)*ephim(m)*hn(n)
            enddo
         enddo
         vwh(:,:,:,1:nodr)=0.d0
         do n=1,nodr
            np1=n+1
            nm1=n-1
            m=-1
            mp1=m+1
            mm1=m-1
            a1=vwh_coef(m,n,1,1)*umn(mp1,np1) &
               +vwh_coef(m,n,1,-1)*umn(mp1,nm1)
            b1=vwh_coef(m,n,-1,1)*umn(mm1,np1) &
               +vwh_coef(m,n,-1,-1)*umn(mm1,nm1)
            z1=vwh_coef(m,n,0,1)*umn(m,np1) &
               +vwh_coef(m,n,0,-1)*umn(m,nm1)
            a2=vwh_coef(m,n,1,0)*umn(mp1,n)
            b2=vwh_coef(m,n,-1,0)*umn(mm1,n)
            z2=vwh_coef(m,n,0,0)*umn(m,n)
            vwh(1,1,1,n)=-0.5d0*(a1+b1)
            vwh(2,1,1,n)=-0.5d0*ci*(-a1+b1)
            vwh(3,1,1,n)=-z1
            vwh(1,2,1,n)=-0.5d0*ci*(a2+b2)
            vwh(2,2,1,n)=-0.5d0*(a2-b2)
            vwh(3,2,1,n)=-ci*z2
            m=1
            mp1=m+1
            mm1=m-1
            a1=vwh_coef(m,n,1,1)*umn(mp1,np1) &
               +vwh_coef(m,n,1,-1)*umn(mp1,nm1)
            b1=vwh_coef(m,n,-1,1)*umn(mm1,np1) &
               +vwh_coef(m,n,-1,-1)*umn(mm1,nm1)
            z1=vwh_coef(m,n,0,1)*umn(m,np1) &
               +vwh_coef(m,n,0,-1)*umn(m,nm1)
            a2=vwh_coef(m,n,1,0)*umn(mp1,n)
            b2=vwh_coef(m,n,-1,0)*umn(mm1,n)
            z2=vwh_coef(m,n,0,0)*umn(m,n)
            vwh(1,1,2,n)=-0.5d0*(a1+b1)
            vwh(2,1,2,n)=-0.5d0*ci*(-a1+b1)
            vwh(3,1,2,n)=-z1
            vwh(1,2,2,n)=-0.5d0*ci*(a2+b2)
            vwh(2,2,2,n)=-0.5d0*(a2-b2)
            vwh(3,2,2,n)=-ci*z2
         enddo
         return
         end subroutine vwhaxialcalc

      end module specialfuncs
!
!  module mpidata
!
!
!  last revised: 15 January 2011
!
      module mpidata
      implicit none
      integer :: group_comm,root_group_comm,base_rank,group_rank,root_group_rank, &
                 base_group,number_groups,proc_per_group,number_proc
      integer, allocatable :: mpi_sphere_index(:), mpi_sphere_number(:)

      contains
!
! allocates the processors into groups
!
!  last revised: 15 January 2011: original
!  20 April 2011: fixedorran=0 now looks for 2 groups.
!  10 october 2011: option for not storing matrices.  If fixorran=0, 2 groups, else
!     nproc groups
!  november 2011: near and far field translation differentiation
!
         subroutine mpisetup(nsphere,nodr,rpos,fixorran,maxmbperproc,istore, &
                    nfdistance,fftranpresent,iunit)
         use mpidefs
         use intrinsics
         use specialfuncs
         implicit none
         integer :: nsphere,numprocs,ierr,i,iunit,nodr(nsphere),fixorran, &
                    nodrmax,nodrmin,temp_comm,newgroup,j,rank,maxmbperproc, &
                    istore,nfspheres,fftranpresent,ffspheres
         integer, allocatable :: grouplist1(:),grouplist2(:)
         real(8) :: memrow(nsphere),memtot,maxmemproc,memperproc
         real(8) :: fp,sum,rpos(3,*),nfdistance,rij(3),r,avenfspheres,rmax, &
                    nfdistancei,aveffspheres
         maxmemproc=maxmbperproc*1.d6
         call ms_mpi(mpi_command='size',mpi_size=numprocs)
         call ms_mpi(mpi_command='rank',mpi_rank=rank)
         call ms_mpi(mpi_command='group',mpi_group=base_group)
         base_rank=rank
         number_proc=numprocs
         memrow=0.d0
         memtot=0.d0
!
!  compute the memory storage requirements
!
         avenfspheres=0.
         aveffspheres=0.
         rmax=0.
         if(1.eq.1) then
            do i=1,nsphere
               nfspheres=0
               do j=1,nsphere
                  rij(:)=rpos(:,i)-rpos(:,j)
                  if(j.ne.i) then
                     if(nfdistance.lt.0.) then
                        nfdistancei=(.5*dble(nodr(i)+nodr(j)))**2.
                     else
                        nfdistancei=nfdistance
                     endif
                     r=sqrt(dot_product(rij,rij))
                     rmax=max(rmax,r)
                     if(r.le.nfdistancei) then
                        nfspheres=nfspheres+1
                        nodrmax=max(nodr(j),nodr(i))
                        nodrmin=min(nodr(j),nodr(i))
                        memrow(i)=memrow(i)+(2*nodrmin+1)*(1+nodrmax*(nodrmax+2))*8.d0
                        memrow(i)=memrow(i)+nodr(i)*nodr(j)*(nodr(j)+3)*16.d0
                        memrow(i)=memrow(i)+(2*nodrmax+1)*16.d0
                     endif
                     if(r.gt.nfdistancei.and.istore.eq.2) then
                        memrow(i)=memrow(i)+2*nodrmax*(nodrmax+2)*16.d0
                     endif
                  endif
               enddo
               ffspheres=nsphere-1-nfspheres
               avenfspheres=avenfspheres+nfspheres
               aveffspheres=aveffspheres+ffspheres
               memtot=memtot+memrow(i)
            enddo
            if(aveffspheres.eq.0) then
               fftranpresent=0
            else
               fftranpresent=1
            endif
            proc_per_group=ceiling(memtot/maxmemproc)
            proc_per_group=min(proc_per_group,numprocs)
            proc_per_group=max(proc_per_group,1)
            do
               if(mod(numprocs,proc_per_group).eq.0) exit
               if(proc_per_group.eq.numprocs) exit
               proc_per_group=proc_per_group+1
            enddo
         endif
         avenfspheres=avenfspheres/dble(nsphere)
         if(rank.eq.0) then
            write(iunit,'('' average near field translations per sphere:'', f10.1)') avenfspheres
            call flush(iunit)
         endif
!
! no-store option
!
         if(istore.eq.0) then
            if(fixorran.eq.0) then
               proc_per_group=max(floor(dble(numprocs)/2.),1)
            else
               proc_per_group=1
            endif
            memrow=1.d0
            memtot=dble(nsphere)
         else
!
! only one or two groups for fixed orientation
!
            if(fixorran.eq.0) proc_per_group=max(floor(dble(numprocs)/2.),proc_per_group)
         endif
         number_groups=numprocs/proc_per_group
         if(allocated(mpi_sphere_index)) deallocate(mpi_sphere_index)
         if(allocated(mpi_sphere_number)) deallocate(mpi_sphere_number)
         allocate(mpi_sphere_index(0:proc_per_group-1),mpi_sphere_number(0:proc_per_group-1), &
                  grouplist1(proc_per_group),grouplist2(number_groups))
         memperproc=memtot/dble(proc_per_group)
!
!  associate the spheres with the processors in a group
!
         mpi_sphere_index(0)=0
         do j=1,proc_per_group-1
            memtot=0.d0
            do i=1,nsphere
               memtot=memtot+memrow(i)
               if(memtot.gt.dble(j)*memperproc) then
                  mpi_sphere_index(j)=i-1
                  exit
               endif
            enddo
         enddo
         do i=0,proc_per_group-2
            mpi_sphere_number(i)=mpi_sphere_index(i+1)-mpi_sphere_index(i)
         enddo
         mpi_sphere_number(proc_per_group-1)=nsphere-mpi_sphere_index(proc_per_group-1)
!
!  assign the sphere-based groups
!
         do i=0,number_groups-1
            do j=0,proc_per_group-1
               grouplist1(j+1)=i*proc_per_group+j
            enddo
            call ms_mpi(mpi_command='incl',mpi_group=base_group,mpi_size=proc_per_group, &
                        mpi_new_group_list=grouplist1,mpi_new_group=newgroup)
            call ms_mpi(mpi_command='create',mpi_group=newgroup,mpi_new_comm=temp_comm)
            if(rank.ge.grouplist1(1).and.rank.le.grouplist1(proc_per_group)) then
               group_comm=temp_comm
            endif
            grouplist2(i+1)=i*proc_per_group
         enddo
!
!  make a group associated with the rank 0 members of the sphere groups
!
         call ms_mpi(mpi_command='incl',mpi_group=base_group,mpi_size=number_groups, &
                     mpi_new_group_list=grouplist2,mpi_new_group=newgroup)
         call ms_mpi(mpi_command='create',mpi_group=newgroup,mpi_new_comm=temp_comm)
         group_rank=mod(rank,proc_per_group)
         root_group_rank=floor(dble(rank)/dble(proc_per_group))
         if(group_rank.eq.0) root_group_comm=temp_comm
         if(rank.eq.0) then
            if(istore.ge.1) then
               write(iunit,'('' number of processors, number groups, mb mem/processor:'',2i5,f9.3)') &
                 numprocs,number_groups,memperproc*1.d-6
               if(memperproc.gt.maxmemproc) then
                  write(iunit,'('' warning: set maximum memory/processor is exceeded!!'')')
               endif
            else
               write(iunit,'('' number of processors, number groups:'',2i5)') &
                 numprocs,number_groups
            endif
            call flush(iunit)
         endif
         deallocate(grouplist1,grouplist2)
         end subroutine mpisetup

      end module mpidata

!
! module spheredata: used to 1) input sphere data, 2) dimension sphere data
! arrays, and 3) provide common access to the data in other subroutines.
!
!
!  last revised: 15 January 2011
!
!  30 March 2011: added optical activity
!
      module spheredata
      use specialfuncs
      use mpidata
      implicit none
      integer, private :: numberspheres,numberiterations,fixedorrandom,numbertheta, &
               calcnf,nfplane,calctmatrix,runprintunit,calcamn,maxmemperproc, &
               trackiterations,nfoutdata,normalizesm,storetranmat,niterstep, &
               fftranpresent
      real(8), private :: lengthscalefactor,realriscalefactor,imriscalefactor,epsmie, &
                 epstran,epssoln,phideg,thetamindeg,thetamaxdeg,alphadeg, &
                 betadeg,epstcon,nfplanepos,nfplanevert(2,2),deltax,gammadeg,epspw, &
                 cgaussbeam,gaussbeamfocus(3),realchiralfactor,imchiralfactor,nfdistance
      character(30), private :: positionfile,outputfile,nfoutputfile,tmatrixfile,printfile, &
                                amnfile
      real(8), private :: xspmax,xvsp
      real(8), private, allocatable :: rpos(:,:),xsp(:)
      complex(8), private, allocatable :: ri(:,:)
      data numberiterations,fixedorrandom,numbertheta/2000,0,181/
      data calcamn,trackiterations,niterstep/1,1,20/
      data lengthscalefactor,realriscalefactor,imriscalefactor,epsmie, &
                 epstran,epssoln,phideg,thetamindeg,thetamaxdeg,alphadeg, &
                 betadeg,epstcon/1.d0,1.d0,1.d0,1.d-4,1.d-6,1.d-10,0.d0,0.d0, &
                 180.d0,0.d0,0.d0,1.d-6/
      data realchiralfactor,imchiralfactor/0.d0,0.d0/
      data normalizesm,storetranmat,nfdistance/0,1,-1.d0/
      data maxmemperproc/1500/
      data cgaussbeam/0.d0/
      data gaussbeamfocus/0.d0,0.d0,0.d0/
      data calcnf,calctmatrix,nfoutdata/0,1,1/
      data runprintunit/6/
      data positionfile,outputfile,tmatrixfile,printfile/'at_bottom','test.dat','tmatrix-temp.dat',' '/
      data nfoutputfile/'nf-temp.dat'/
      data amnfile/'amn-temp.dat'/

      contains
!
!  Find the number of data points in input unit iunit, and reposition the unit to the
!  point after record containing parmid
!
!
!  last revised: 15 January 2011
!
         subroutine numberinrecord(iunit,parmid,numrec)
         implicit none
         integer :: numrec,iunit
         character*1 :: a
         character*35 :: parmid
         character*10 :: rec
         numrec=0
         do
            read(iunit,"(a)",advance="no",err=100,eor=100) a
            if(a.ne.' '.and.a.ne.',') then
!
! start of a number
!
               numrec=numrec+1
!
! look for the delimeter
!
               do
                  read(iunit,"(a)",advance="no",err=100,eor=100) a
                  if(a.eq.' '.or.a.eq.',') exit
               enddo
            endif
         enddo
100      if(parmid.eq.'rewind') then
            rewind(iunit)
         else
            backspace(iunit)
            backspace(iunit)
            backspace(iunit)
            do
               read(iunit,'(a10)') rec
               if(rec.eq.parmid(1:10)) exit
            enddo
         endif
         end subroutine numberinrecord
!
!  inputdata:  reads parameters from inputfile
!              reads sphere data from position file
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fix output file initialization.
!  30 March 2011: added optical activity
!
!
         subroutine inputdata(inputfile,printdata)
         integer :: imax,i,j,ierr,iunit,numrec,nsphere,printdata
         real(8) :: rmax,rtoi,rposmean(3),rireal,riimag,dtemp,betareal,betaimag, &
                    rij,xij(3),rijmax
         real(8), allocatable :: sdat(:)
         complex(8) :: ribulk,beta
         character*35 :: parmid
         character*30 :: inputfile
!
!  cycle through parameter input operations
!
         open(1,file=inputfile)
         do
            read(1,'(a)',end=10) parmid
            parmid=parmid(:index(parmid,' '))
            if(parmid.eq.'number_spheres') then
               read(1,*) numberspheres
               cycle
            endif
            if(parmid.eq.'sphere_position_file') then
               read(1,'(a)') positionfile
               positionfile=positionfile(:index(positionfile,' '))
               cycle
            endif
            if(parmid.eq.'output_file') then
               read(1,'(a)') outputfile
               outputfile=outputfile(:index(outputfile,' '))
               cycle
            endif
            if(parmid.eq.'run_print_file') then
               read(1,'(a)') printfile
               printfile=printfile(:index(printfile,' '))
               if(printdata.eq.1) then
                  if((printfile.eq.' '.or.printfile.eq.'console')) then
                     printfile=' '
                     runprintunit=6
                  else
                     runprintunit=4
                     open(runprintunit,file=printfile)
                  endif
               else
                  runprintunit=6
               endif
               cycle
            endif
            if(parmid.eq.'length_scale_factor') then
               read(1,*) lengthscalefactor
               cycle
            endif
            if(parmid.eq.'real_ref_index_scale_factor') then
               read(1,*) realriscalefactor
               cycle
            endif
            if(parmid.eq.'imag_ref_index_scale_factor') then
               read(1,*) imriscalefactor
               cycle
            endif
            if(parmid.eq.'real_chiral_factor') then
               read(1,*) realchiralfactor
               cycle
            endif
            if(parmid.eq.'imag_chiral_factor') then
               read(1,*) imchiralfactor
               cycle
            endif
            if(parmid.eq.'mie_epsilon') then
               read(1,*) epsmie
               cycle
            endif
            if(parmid.eq.'translation_epsilon') then
               read(1,*) epstran
               cycle
            endif
            if(parmid.eq.'solution_epsilon') then
               read(1,*) epssoln
               cycle
            endif
            if(parmid.eq.'max_number_iterations') then
               read(1,*) numberiterations
               cycle
            endif
            if(parmid.eq.'max_memory_per_processor') then
               read(1,*) maxmemperproc
               cycle
            endif
            if(parmid.eq.'store_translation_matrix') then
               read(1,*) storetranmat
               cycle
            endif
            if(parmid.eq.'near_field_distance') then
               read(1,*) nfdistance
               cycle
            endif
            if(parmid.eq.'iterations_per_correction') then
               read(1,*) niterstep
               cycle
            endif
            if(parmid.eq.'fixed_or_random_orientation') then
               read(1,*) fixedorrandom
               cycle
            endif
            if(parmid.eq.'scattering_plane_angle_deg') then
               read(1,*) phideg
               cycle
            endif
            if(parmid.eq.'min_scattering_angle_deg') then
               read(1,*) thetamindeg
               cycle
            endif
            if(parmid.eq.'max_scattering_angle_deg') then
               read(1,*) thetamaxdeg
               cycle
            endif
            if(parmid.eq.'number_scattering_angles') then
               read(1,*) numbertheta
               cycle
            endif
            if(parmid.eq.'normalize_scattering_matrix') then
               read(1,*) normalizesm
               cycle
            endif
            if(parmid.eq.'incident_azimuth_angle_deg') then
               read(1,*) alphadeg
               cycle
            endif
            if(parmid.eq.'incident_polar_angle_deg') then
               read(1,*) betadeg
               cycle
            endif
            if(parmid.eq.'calculate_scattering_coefficients') then
               read(1,*) calcamn
               cycle
            endif
            if(parmid.eq.'scattering_coefficient_file') then
               read(1,'(a)') amnfile
               if(amnfile.eq.' ') then
                  amnfile='amn-temp.dat'
               else
                  amnfile=amnfile(:index(amnfile,' '))
               endif
               cycle
            endif
            if(parmid.eq.'track_iterations') then
               read(1,*) trackiterations
               cycle
            endif
            if(parmid.eq.'calculate_near_field') then
               read(1,*) calcnf
               cycle
            endif
            if(parmid.eq.'near_field_plane_coord') then
               read(1,*) nfplane
               cycle
            endif
            if(parmid.eq.'near_field_plane_position') then
               read(1,*) nfplanepos
               cycle
            endif
            if(parmid.eq.'near_field_plane_vertices') then
               read(1,*) nfplanevert
               cycle
            endif
            if(parmid.eq.'spacial_step_size') then
               read(1,*) deltax
               cycle
            endif
            if(parmid.eq.'polarization_angle_deg') then
               read(1,*) gammadeg
               cycle
            endif
            if(parmid.eq.'near_field_output_file') then
               read(1,'(a)') nfoutputfile
               if(nfoutputfile.eq.' ') then
                  nfoutputfile='nf-temp.dat'
               else
                  nfoutputfile=nfoutputfile(:index(nfoutputfile,' '))
               endif
               cycle
            endif
            if(parmid.eq.'near_field_output_data') then
               read(1,*) nfoutdata
               cycle
            endif
            if(parmid.eq.'plane_wave_epsilon') then
               read(1,*) epspw
               cycle
            endif
            if(parmid.eq.'gaussian_beam_constant') then
               read(1,*) cgaussbeam
               cycle
            endif
            if(parmid.eq.'gaussian_beam_focal_point') then
               read(1,*) gaussbeamfocus
               cycle
            endif
            if(parmid.eq.'t_matrix_convergence_epsilon') then
               read(1,*) epstcon
               cycle
            endif
            if(parmid.eq.'calculate_t_matrix') then
               read(1,*) calctmatrix
               cycle
            endif
            if(parmid.eq.'t_matrix_file') then
               read(1,'(a)') tmatrixfile
               if(tmatrixfile.eq.' ') then
                  tmatrixfile='tmatrix-temp.dat'
               else
                  tmatrixfile=tmatrixfile(:index(tmatrixfile,' '))
               endif
               cycle
            endif
            if(parmid.eq.'sphere_sizes_and_positions') exit
            if(parmid.eq.'end_of_options') exit
            write(*,'('' warning: unknown parameter ID:'',a35)') parmid
         enddo
!
!  end of parameter input options.   Input of sphere data follows
!
10       write(runprintunit,'('' input file is '',a30)') inputfile
         if(positionfile.ne.'at_bottom'.and.positionfile.ne.' ') then
            close(1)
            open(1,file=positionfile)
            parmid='rewind'
         endif
!
!  find number of records in position file
!
         call numberinrecord(1,parmid,numrec)
         if(printdata.eq.1) write(runprintunit,'('' position data has '',i3,'' records'')') numrec
         nsphere=numberspheres
         iunit=1
         allocate(sdat(numrec))
         allocate(xsp(0:nsphere),rpos(3,0:nsphere),ri(2,0:nsphere),stat=ierr)
         xvsp=0.d0
         do i=1,nsphere
            read(iunit,*,end=20) sdat
            xsp(i)=sdat(1)*lengthscalefactor
            rpos(1:3,i)=sdat(2:4)*lengthscalefactor
            if(numrec.gt.4) then
               rireal=sdat(5)*realriscalefactor
               riimag=sdat(6)*imriscalefactor
            else
               rireal=realriscalefactor
               riimag=imriscalefactor
            endif
            if(numrec.gt.6) then
               betareal=sdat(7)*realchiralfactor
               betaimag=sdat(8)*imchiralfactor
            else
               betareal=realchiralfactor
               betaimag=imchiralfactor
            endif
            ribulk=dcmplx(rireal,riimag)
            beta=dcmplx(betareal,betaimag)
            if(beta.eq.(0.d0,0.d0)) then
               ri(1,i)=ribulk
               ri(2,i)=ribulk
            else
               ri(1,i)=ribulk/(1.d0-beta*ribulk)
               ri(2,i)=ribulk/(1.d0+beta*ribulk)
            endif
            xvsp=xvsp+xsp(i)**3.d0
         enddo
20       nsphere=min(nsphere,i-1)
         close(iunit)
         deallocate(sdat)
         if(nsphere.ne.numberspheres.and.printdata.eq.1) then
            write(runprintunit,'('' warning: insufficient position points in file.'')')
            write(runprintunit,'('' number of spheres truncated to:'',i5)') nsphere
         endif
!
! check for overlapping spheres, and find maximum translation
!
         rijmax=0.
         do i=1,nsphere
            do j=i+1,nsphere
               xij=rpos(:,i)-rpos(:,j)
               rij=sqrt(dot_product(xij,xij))
               rijmax=max(rijmax,rij)
               if(rij/(xsp(i)+xsp(j)).lt..999d0) then
                  write(runprintunit,'('' warning: spheres '',i4,'' and '',i4 '' overlap. '',&
                  & '' scaled distance:'' f8.4)') i,j,rij/(xsp(i)+xsp(j))
               endif
            enddo
         enddo
         if(rijmax.gt.nfdistance) then
            fftranpresent=1
         else
            fftranpresent=0
         endif
         numberspheres=nsphere
         xvsp=xvsp**(1.d0/3.d0)
         gaussbeamfocus=gaussbeamfocus*lengthscalefactor
         if(nsphere.eq.1) then
            rposmean=rpos(:,1)
            rpos(:,1)=0.d0
            xspmax=xsp(1)
         else
            rposmean=0.d0
            do i=1,nsphere
               rposmean=rposmean+rpos(:,i)
            enddo
            rposmean=rposmean/dble(nsphere)
            rmax=0.d0
!
!  the target origin is defined as the GB focal point.
!
            do i=1,nsphere
!               rpos(1:3,i)=rpos(1:3,i)-rposmean(1:3)
               rpos(1:3,i)=rpos(1:3,i)-gaussbeamfocus(1:3)
               rtoi=dot_product(rpos(:,i),rpos(:,i))
               if(rtoi.gt.rmax) then
                  rmax=rtoi
                  imax=i
               endif
            enddo
            xspmax=sqrt(rmax)+xsp(imax)
         endif
!
!  xsp(0) is the circumscribing sphere size parameter
!
         xsp(0)=xspmax
         ri(1,0)=(1.d0,0.d0)
         ri(2,0)=(1.d0,0.d0)
         rpos(:,0)=0.d0
!
!  write run data to run file and output file
!
         if(printdata.eq.1) then
            call writerundata(runprintunit)
            call flush(runprintunit)
            open(1,file=outputfile,status='replace',action='write')
            call writerundata(1)
            close(1)
         endif
         end subroutine inputdata
!
!  writes run data to output unit iunit
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine writerundata(iunit)
         implicit none
         integer :: iunit,i
         character*1 :: lf
         if(iunit.ne.1) then
            lf = ' '
         else
            lf = '/'
         endif
         write(iunit,'('' number of spheres, volume size parameter:'' '//lf//',i5,e13.5)') &
                       numberspheres,xvsp
         write(iunit,'('' position file:'' '//lf//',a)') positionfile
         write(iunit,'('' output file:'' '//lf//',a)') outputfile
         write(iunit,'('' length, ref. indx. scale factors:'' '//lf//',3f8.3)') lengthscalefactor, &
                      realriscalefactor,imriscalefactor
         write(iunit,'('' chiral factors:'' '//lf//',2e13.5)')  &
                      realchiralfactor,imchiralfactor
         write(iunit,'('' thetamin, thetamax, num. theta:'' '//lf//',2f9.1,i5)') &
                      thetamindeg,thetamaxdeg,numbertheta
         write(iunit,'('' epsmie, epssoln, max number iterations:'' '//lf//',2e12.4,i5)') epsmie, &
                      epssoln, numberiterations
         if(fftranpresent.eq.1) then
            write(iunit,'('' far field kr, iterations/correction:'' '//lf//',e12.4,i5)') &
                 nfdistance,niterstep
         else
            write(iunit,'('' all translations computed exactly'' '//lf//')')
         endif
         if(cgaussbeam.ne.0.d0) then
            write(iunit,'('' gaussian incident beam: 1/width:'' '//lf//',f9.4,)') cgaussbeam
            write(iunit,'('' beam focal point:'' '//lf//',3f9.3,)') gaussbeamfocus
         else
            write(iunit,'('' plane wave incidence'')')
         endif
         if(fixedorrandom.eq.0) then
            write(iunit,'('' fixed orientation calculations'')')
            write(iunit,'('' scattering plane, incident alpha, beta:'' '//lf//',3f9.2)') &
                     phideg,alphadeg,betadeg
            write(iunit,'('' common expansion epsilon:'' '//lf//',e12.4)') epstran
            if(calcamn.eq.0) then
               write(iunit,'('' scattering coefficients read from file '' '//lf//',a)') amnfile
            else
               write(iunit,'('' scattering coefficients calculated, stored in file '' '//lf//',a)') amnfile
            endif
            if(calcnf.eq.1) then
               write(iunit,'('' near field calculated, stored in file '' '//lf//',a)') nfoutputfile
               write(iunit,'('' near field data output option: '' '//lf//',i4)') nfoutdata
               write(iunit,'('' near field plane, position: '' '//lf//', i4,f9.3)') nfplane, nfplanepos
               write(iunit,'('' near field plane vertices: '' '//lf//',4f9.3)') nfplanevert
               write(iunit,'('' spacial step size:'' '//lf//',f9.4)') deltax
               write(iunit,'('' polarization angle, deg.:'' '//lf//',f9.2)') gammadeg
               write(iunit,'('' plane wave epsilon:'' '//lf//',e13.5)') epspw
            endif
         else
            write(iunit,'('' random orientation calculations'')')
            if(calctmatrix.eq.0) then
               write(iunit,'('' t matrix read from file '' '//lf//',a)') tmatrixfile
            elseif(calctmatrix.eq.1) then
               write(iunit,'('' t matrix calculated, stored in file '' '//lf//',a)') tmatrixfile
               write(iunit,'('' t matrix convergence epsilon:'' '//lf//',e12.4)') epstcon
            else
               write(iunit,'('' t matrix calculated from end of file '' '//lf//',a)') tmatrixfile
               write(iunit,'('' t matrix convergence epsilon:'' '//lf//',e12.4)') epstcon
            endif
         endif
         end subroutine writerundata
!
!  getspheredata: retrieves sphere data
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine getspheredata(number_spheres, sphere_size_parameters, sphere_positions, &
              sphere_refractive_indices, volume_size_parameter)
         implicit none
         integer, optional :: number_spheres
         real(8), optional :: sphere_size_parameters(numberspheres), &
                              sphere_positions(3,numberspheres), volume_size_parameter
         complex(8), optional :: sphere_refractive_indices(2,numberspheres)
         if (present(number_spheres)) number_spheres=numberspheres
         if (present(sphere_size_parameters)) sphere_size_parameters(1:numberspheres)=xsp(1:numberspheres)
         if (present(sphere_positions)) sphere_positions(:,1:numberspheres)=rpos(:,1:numberspheres)
         if (present(sphere_refractive_indices)) &
               sphere_refractive_indices(:,1:numberspheres)=ri(:,1:numberspheres)
         if (present(volume_size_parameter)) volume_size_parameter=xvsp
         end subroutine getspheredata

         subroutine getspheredataone(sphere,sphere_size_parameter, sphere_position, &
              sphere_refractive_index)
         implicit none
         integer :: sphere
         real(8), optional :: sphere_size_parameter,sphere_position(3)
         complex(8), optional :: sphere_refractive_index(2)
         if (present(sphere_size_parameter)) sphere_size_parameter=xsp(sphere)
         if (present(sphere_position)) sphere_position(:)=rpos(:,sphere)
         if (present(sphere_refractive_index)) &
               sphere_refractive_index(:)=ri(:,sphere)
         end subroutine getspheredataone
!
!  setspheredata: sets sphere data
!
         subroutine setspheredata(number_spheres, sphere_size_parameters, sphere_positions, &
              sphere_refractive_indices, volume_size_parameter)
         implicit none
         integer :: i
         integer, optional :: number_spheres
         real(8), optional :: sphere_size_parameters(*), &
                              sphere_positions(3,*), volume_size_parameter
         complex(8), optional :: sphere_refractive_indices(2,*)
         if (present(number_spheres)) then
            numberspheres=number_spheres
            if(allocated(xsp)) deallocate(xsp,rpos,ri)
            allocate(xsp(0:numberspheres),rpos(3,0:numberspheres),ri(2,0:numberspheres))
         endif
         if (present(sphere_size_parameters))    xsp(1:numberspheres)       =sphere_size_parameters(1:numberspheres)
         if (present(sphere_positions))          rpos(:,1:numberspheres)    =sphere_positions(:,1:numberspheres)
         if (present(sphere_refractive_indices)) ri(:,1:numberspheres)      =sphere_refractive_indices(:,1:numberspheres)
         if (present(volume_size_parameter))     xvsp                       =volume_size_parameter
         end subroutine setspheredata
!
!  getrunparameters: retrieves run parameters read from input file
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine getrunparameters(number_spheres,sphere_position_file,output_file, &
                       length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       max_number_iterations,fixed_or_random_orientation,scattering_plane_angle_deg, &
                       min_scattering_angle_deg,max_scattering_angle_deg,number_scattering_angles, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,calculate_near_field, &
                       near_field_plane_coord,near_field_plane_position,near_field_plane_vertices, &
                       spacial_step_size,polarization_angle_deg,near_field_output_file, &
                       plane_wave_epsilon,t_matrix_convergence_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point,calculate_t_matrix,t_matrix_file,run_print_file, &
                       run_print_unit,calculate_scattering_coefficients,scattering_coefficient_file, &
                       max_memory_per_processor,track_iterations,near_field_output_data, &
                       real_chiral_factor,imag_chiral_factor,normalize_scattering_matrix, &
                       store_translation_matrix,near_field_distance, &
                       iterations_per_correction)
         implicit none
         integer, optional :: number_spheres,max_number_iterations,fixed_or_random_orientation, &
                              number_scattering_angles,calculate_near_field,near_field_plane_coord, &
                              calculate_t_matrix,run_print_unit,calculate_scattering_coefficients, &
                              max_memory_per_processor,track_iterations,near_field_output_data, &
                              normalize_scattering_matrix,store_translation_matrix, &
                              iterations_per_correction
         real(8), optional :: length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       scattering_plane_angle_deg, &
                       min_scattering_angle_deg,max_scattering_angle_deg, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,t_matrix_convergence_epsilon, &
                       near_field_plane_position,near_field_plane_vertices(2,2),spacial_step_size, &
                       polarization_angle_deg,plane_wave_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point(3),real_chiral_factor,imag_chiral_factor, &
                       near_field_distance
         character*30, optional :: sphere_position_file,output_file,near_field_output_file, &
                                   t_matrix_file,run_print_file,scattering_coefficient_file
         if(present(number_spheres))                    number_spheres                    =numberspheres
         if(present(sphere_position_file))              sphere_position_file              =positionfile
         if(present(output_file))                       output_file                       =outputfile
         if(present(length_scale_factor))               length_scale_factor               =lengthscalefactor
         if(present(real_ref_index_scale_factor))       real_ref_index_scale_factor       =realriscalefactor
         if(present(imag_ref_index_scale_factor))       imag_ref_index_scale_factor       =imriscalefactor
         if(present(mie_epsilon))                       mie_epsilon                       =epsmie
         if(present(translation_epsilon))               translation_epsilon               =epstran
         if(present(solution_epsilon))                  solution_epsilon                  =epssoln
         if(present(max_number_iterations))             max_number_iterations             =numberiterations
         if(present(track_iterations))                  track_iterations                  =trackiterations
         if(present(max_memory_per_processor))          max_memory_per_processor          =maxmemperproc
         if(present(fixed_or_random_orientation))       fixed_or_random_orientation       =fixedorrandom
         if(present(scattering_plane_angle_deg))        scattering_plane_angle_deg        =phideg
         if(present(min_scattering_angle_deg))          min_scattering_angle_deg          =thetamindeg
         if(present(max_scattering_angle_deg))          max_scattering_angle_deg          =thetamaxdeg
         if(present(number_scattering_angles))          number_scattering_angles          =numbertheta
         if(present(normalize_scattering_matrix))       normalize_scattering_matrix       =normalizesm
         if(present(incident_azimuth_angle_deg))        incident_azimuth_angle_deg        =alphadeg
         if(present(incident_polar_angle_deg))          incident_polar_angle_deg          =betadeg
         if(present(t_matrix_convergence_epsilon))      t_matrix_convergence_epsilon      =epstcon
         if(present(calculate_near_field))              calculate_near_field              =calcnf
         if(present(near_field_plane_coord))            near_field_plane_coord            =nfplane
         if(present(near_field_plane_position))         near_field_plane_position         =nfplanepos
         if(present(near_field_plane_vertices))         near_field_plane_vertices         =nfplanevert
         if(present(spacial_step_size))                 spacial_step_size                 =deltax
         if(present(polarization_angle_deg))            polarization_angle_deg            =gammadeg
         if(present(near_field_output_file))            near_field_output_file            =nfoutputfile
         if(present(near_field_output_data))            near_field_output_data            =nfoutdata
         if(present(plane_wave_epsilon))                plane_wave_epsilon                =epspw
         if(present(gaussian_beam_constant))            gaussian_beam_constant            =cgaussbeam
         if(present(gaussian_beam_focal_point))         gaussian_beam_focal_point         =gaussbeamfocus
         if(present(t_matrix_file))                     t_matrix_file                     =tmatrixfile
         if(present(calculate_t_matrix))                calculate_t_matrix                =calctmatrix
         if(present(run_print_file))                    run_print_file                    =printfile
         if(present(run_print_unit))                    run_print_unit                    =runprintunit
         if(present(calculate_scattering_coefficients)) calculate_scattering_coefficients =calcamn
         if(present(scattering_coefficient_file))       scattering_coefficient_file       =amnfile
         if(present(real_chiral_factor))                real_chiral_factor                =realchiralfactor
         if(present(imag_chiral_factor))                imag_chiral_factor                =imchiralfactor
         if(present(store_translation_matrix))          store_translation_matrix          =storetranmat
         if(present(near_field_distance))               near_field_distance               =nfdistance
         if(present(iterations_per_correction))         iterations_per_correction         =niterstep
         end subroutine getrunparameters
!
!  set run parameters: set run parameters
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine setrunparameters(number_spheres,sphere_position_file,output_file, &
                       length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       max_number_iterations,fixed_or_random_orientation,scattering_plane_angle_deg, &
                       min_scattering_angle_deg,max_scattering_angle_deg,number_scattering_angles, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,calculate_near_field, &
                       near_field_plane_coord,near_field_plane_position,near_field_plane_vertices, &
                       spacial_step_size,polarization_angle_deg,near_field_output_file, &
                       plane_wave_epsilon,t_matrix_convergence_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point,calculate_t_matrix,t_matrix_file,run_print_file, &
                       run_print_unit,calculate_scattering_coefficients,scattering_coefficient_file, &
                       max_memory_per_processor,track_iterations,near_field_output_data, &
                       real_chiral_factor,imag_chiral_factor,store_translation_matrix, &
                       near_field_distance,iterations_per_correction)
         implicit none
         integer, optional :: number_spheres,max_number_iterations,fixed_or_random_orientation, &
                              number_scattering_angles,calculate_near_field,near_field_plane_coord, &
                              calculate_t_matrix,run_print_unit,calculate_scattering_coefficients, &
                              max_memory_per_processor,track_iterations,near_field_output_data, &
                              store_translation_matrix,iterations_per_correction
         real(8), optional :: length_scale_factor,real_ref_index_scale_factor, &
                       imag_ref_index_scale_factor,mie_epsilon,translation_epsilon,solution_epsilon, &
                       scattering_plane_angle_deg,near_field_distance,&
                       min_scattering_angle_deg,max_scattering_angle_deg, &
                       incident_azimuth_angle_deg,incident_polar_angle_deg,t_matrix_convergence_epsilon, &
                       near_field_plane_position,near_field_plane_vertices(2,2),spacial_step_size, &
                       polarization_angle_deg,plane_wave_epsilon,gaussian_beam_constant, &
                       gaussian_beam_focal_point(3),real_chiral_factor,imag_chiral_factor
         character*30, optional :: sphere_position_file,output_file,near_field_output_file, &
                                   t_matrix_file,run_print_file,scattering_coefficient_file

         if(present(number_spheres))                     numberspheres        =number_spheres
         if(present(sphere_position_file))               positionfile         =sphere_position_file
         if(present(output_file))                        outputfile           =output_file
         if(present(length_scale_factor))                lengthscalefactor    =length_scale_factor
         if(present(real_ref_index_scale_factor))        realriscalefactor    =real_ref_index_scale_factor
         if(present(imag_ref_index_scale_factor))        imriscalefactor      =imag_ref_index_scale_factor
         if(present(mie_epsilon))                        epsmie               =mie_epsilon
         if(present(translation_epsilon))                epstran              =translation_epsilon
         if(present(solution_epsilon))                   epssoln              =solution_epsilon
         if(present(max_number_iterations))              numberiterations     =max_number_iterations
         if(present(track_iterations))                   trackiterations      =track_iterations
         if(present(max_memory_per_processor))           maxmemperproc        =max_memory_per_processor
         if(present(fixed_or_random_orientation))        fixedorrandom        =fixed_or_random_orientation
         if(present(scattering_plane_angle_deg))         phideg               =scattering_plane_angle_deg
         if(present(min_scattering_angle_deg))           thetamindeg          =min_scattering_angle_deg
         if(present(max_scattering_angle_deg))           thetamaxdeg          =max_scattering_angle_deg
         if(present(number_scattering_angles))           numbertheta          =number_scattering_angles
         if(present(incident_azimuth_angle_deg))         alphadeg             =incident_azimuth_angle_deg
         if(present(incident_polar_angle_deg))           betadeg              =incident_polar_angle_deg
         if(present(t_matrix_convergence_epsilon))       epstcon              =t_matrix_convergence_epsilon
         if(present(calculate_near_field))               calcnf               =calculate_near_field
         if(present(near_field_plane_coord))             nfplane              =near_field_plane_coord
         if(present(near_field_plane_position))          nfplanepos           =near_field_plane_position
         if(present(near_field_plane_vertices))          nfplanevert          =near_field_plane_vertices
         if(present(spacial_step_size))                  deltax               =spacial_step_size
         if(present(polarization_angle_deg))             gammadeg             =polarization_angle_deg
         if(present(near_field_output_file))             nfoutputfile         =near_field_output_file
         if(present(near_field_output_data))             nfoutdata            =near_field_output_data
         if(present(plane_wave_epsilon))                 epspw                =plane_wave_epsilon
         if(present(gaussian_beam_constant))             cgaussbeam           =gaussian_beam_constant
         if(present(gaussian_beam_focal_point))          gaussbeamfocus       =gaussian_beam_focal_point
         if(present(t_matrix_file))                      tmatrixfile          =t_matrix_file
         if(present(calculate_t_matrix))                 calctmatrix          =calculate_t_matrix
         if(present(run_print_file))                     printfile            =run_print_file
         if(present(run_print_unit))                     runprintunit         =run_print_unit
         if(present(calculate_scattering_coefficients))  calcamn              =calculate_scattering_coefficients
         if(present(scattering_coefficient_file))        amnfile              =scattering_coefficient_file
         if(present(real_chiral_factor))                 realchiralfactor     =real_chiral_factor
         if(present(imag_chiral_factor))                 imchiralfactor       =imag_chiral_factor
         if(present(store_translation_matrix))           storetranmat         =store_translation_matrix
         if(present(near_field_distance))                nfdistance           =near_field_distance
         if(present(iterations_per_correction))          niterstep            =iterations_per_correction
         end subroutine setrunparameters

      end module spheredata
!
! module miecoefdata: used to 1) calculate single sphere mie coefficient values,
! 2) store values in an allocated array, 3) provide common access to values, and
! 4) perform multiplication of coefficient values with vectors containing VWH scattering
! coefficients.
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
      module miecoefdata
      implicit none
      integer, private :: numeqns,maxorder
      integer, allocatable, private :: nodr(:),nodroffset(:),nblk(:),nblkoffset(:)
      real(8), allocatable, private :: qextmie(:),qabsmie(:)
      complex(8), allocatable, private :: anmie(:,:,:),cnmie(:,:,:)
      interface getmiedata
            module procedure getmiedataall, getmiedataone
      end interface getmiedata

      contains
!
!  calculation of the max order of sphere expansions and storage of mie coefficients
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine miecoefcalc(nsphere,xsp,ri,qeps)
         implicit none
         integer :: n,nodrn,nsphere,nodrtot,ierr,nblktot
         real(8) :: qext,qabs,qsca,qeps,xsp(nsphere)
         complex(8) :: ri(2,nsphere)
         complex(8), allocatable :: anp(:,:,:),cnp(:,:,:)
         if(allocated(nodr)) deallocate(nodr,nodroffset,nblk, &
                      nblkoffset,qextmie,qabsmie)
         allocate(nodr(nsphere),nodroffset(nsphere+1), &
                  nblk(nsphere),nblkoffset(nsphere+1),  &
                  qextmie(nsphere),qabsmie(nsphere),stat=ierr)
         if(ierr.ne.0) then
            write(*,'('' bad allocation in nodr: stat:'',i4)') ierr
         endif
         nodrtot=0
         nblktot=0
         maxorder=0
!
!  calculate the order limits and efficiencies
!
         do n=1,nsphere
            call mieoa(xsp(n),ri(1,n),nodrn,qeps,qext,qsca)
            nodroffset(n)=nodrtot
            nblkoffset(n)=nblktot
            nodr(n)=nodrn
            maxorder=max(maxorder,nodrn)
            nblk(n)=nodrn*(nodrn+2)*2
            nodrtot=nodrtot+nodrn
            nblktot=nblktot+nblk(n)
            qextmie(n)=qext
            qabsmie(n)=qext-qsca
         enddo
         nodroffset(nsphere+1)=nodrtot
         nblkoffset(nsphere+1)=nblktot
         numeqns=nblktot
!
! calculate the mie coefficients, and store in memory
!
         if(allocated(anmie)) deallocate(anmie,cnmie)
         allocate(anmie(2,2,nodrtot),cnmie(2,2,nodrtot),stat=ierr)
         if(ierr.ne.0) then
            write(*,'('' bad allocation in anmie: stat:'',i4)') ierr
         endif
         do n=1,nsphere
            if(abs(ri(1,n)-ri(2,n)).eq.0) then
               allocate(anp(2,1,nodr(n)),cnp(2,1,nodr(n)))
               call mieregular(xsp(n),ri(1,n),nodrn,qeps,qext,qsca,anp_mie=anp,cnp_mie=cnp)
               anmie(1,1,nodroffset(n)+1:nodroffset(n+1))=anp(1,1,1:nodr(n))
               anmie(2,2,nodroffset(n)+1:nodroffset(n+1))=anp(2,1,1:nodr(n))
               anmie(1,2,nodroffset(n)+1:nodroffset(n+1))=0.d0
               anmie(2,1,nodroffset(n)+1:nodroffset(n+1))=0.d0
               cnmie(1,1,nodroffset(n)+1:nodroffset(n+1))=cnp(1,1,1:nodr(n))
               cnmie(2,2,nodroffset(n)+1:nodroffset(n+1))=cnp(2,1,1:nodr(n))
               cnmie(1,2,nodroffset(n)+1:nodroffset(n+1))=0.d0
               cnmie(2,1,nodroffset(n)+1:nodroffset(n+1))=0.d0
               deallocate(anp,cnp)
            else
               allocate(anp(2,2,nodr(n)),cnp(2,2,nodr(n)))
               call mieoa(xsp(n),ri(1,n),nodrn,qeps,qext,qsca,anp_mie=anp,cnp_mie=cnp)
               anmie(1:2,1:2,nodroffset(n)+1:nodroffset(n+1))=anp(1:2,1:2,1:nodr(n))
               cnmie(1:2,1:2,nodroffset(n)+1:nodroffset(n+1))=cnp(1:2,1:2,1:nodr(n))
               deallocate(anp,cnp)
            endif
         enddo
         end subroutine miecoefcalc
!
!  retrieve the array of mie data
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine getmiedataall(sphere_order, sphere_block, &
                     sphere_order_offset, sphere_block_offset, sphere_qext, &
                     sphere_qabs, sphere_mie_coefficients, sphere_int_mie_coefficients, &
                     number_equations, max_order)
         use spheredata
         implicit none
         integer, optional :: sphere_order(:), sphere_block(:), sphere_order_offset(:), &
                     sphere_block_offset(:),number_equations, max_order
         integer :: i,nsphere
         real(8), optional :: sphere_qext(:), sphere_qabs(:)
         complex(8), optional :: sphere_mie_coefficients(:,:,:,:), &
                       sphere_int_mie_coefficients(:,:,:,:)
         call getspheredata(number_spheres=nsphere)
         if(present(sphere_order)) sphere_order=nodr
         if(present(sphere_block)) sphere_block=nblk
         if(present(sphere_order_offset)) sphere_order_offset=nodroffset
         if(present(sphere_block_offset)) sphere_block_offset=nblkoffset
         if(present(sphere_qext)) sphere_qext=qextmie
         if(present(sphere_qabs)) sphere_qabs=qabsmie
         if(present(number_equations)) number_equations=numeqns
         if(present(max_order)) max_order=maxorder
         if(present(sphere_mie_coefficients)) then
            do i=1,nsphere
               sphere_mie_coefficients(1:2,1:2,1:nodr(i),i)  &
                =anmie(1:2,1:2,nodroffset(i)+1:nodroffset(i+1))
            enddo
         endif
         if(present(sphere_int_mie_coefficients)) then
            do i=1,nsphere
               sphere_int_mie_coefficients(1:2,1:2,1:nodr(i),i)  &
                =cnmie(1:2,1:2,nodroffset(i)+1:nodroffset(i+1))
            enddo
         endif
         end subroutine getmiedataall
!
!  retrieve mie data for a single sphere
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine getmiedataone(which_sphere, sphere_order, sphere_block, &
                     sphere_order_offset, sphere_block_offset, sphere_qext, &
                     sphere_qabs, sphere_mie_coefficients, sphere_int_mie_coefficients, &
                     number_equations, max_order)
         use spheredata
         implicit none
         integer, optional :: sphere_order, sphere_block, sphere_order_offset, &
                     sphere_block_offset, number_equations, max_order
         integer :: which_sphere
         integer :: i,nsphere
         real(8), optional :: sphere_qext, sphere_qabs
         complex(8), optional :: sphere_mie_coefficients(:,:,:), sphere_int_mie_coefficients(:,:,:)
         i=which_sphere
         if(present(sphere_order)) sphere_order=nodr(i)
         if(present(sphere_block)) sphere_block=nblk(i)
         if(present(sphere_order_offset)) sphere_order_offset=nodroffset(i)
         if(present(sphere_block_offset)) sphere_block_offset=nblkoffset(i)
         if(present(sphere_qext)) sphere_qext=qextmie(i)
         if(present(sphere_qabs)) sphere_qabs=qabsmie(i)
         if(present(number_equations)) number_equations=numeqns
         if(present(max_order)) max_order=maxorder
         if(present(sphere_mie_coefficients)) &
            sphere_mie_coefficients(1:2,1:2,1:nodr(i))  &
                =anmie(1:2,1:2,nodroffset(i)+1:nodroffset(i+1))
         if(present(sphere_int_mie_coefficients)) &
            sphere_int_mie_coefficients(1:2,1:2,1:nodr(i))  &
                =cnmie(1:2,1:2,nodroffset(i)+1:nodroffset(i+1))
         end subroutine getmiedataone
!
!  retrieve mie coefficients for sphere n
!  30 March 2011: added optical activity
!
         function miecoef(n)
         implicit none
         integer :: n
         complex(8), dimension(2,2,nodr(n)) :: miecoef
         miecoef=anmie(1:2,1:2,nodroffset(n)+1:nodroffset(n+1))
         end function miecoef

         function internalmiecoef(n)
         implicit none
         integer :: n
         complex(8), dimension(2,2,nodr(n)) :: internalmiecoef
         internalmiecoef=cnmie(1:2,1:2,nodroffset(n)+1:nodroffset(n+1))
         end function internalmiecoef
!
!  multiples the solution vector cx by mie coefficients and returns in y
!  i1: starting sphere, i2: ending sphere
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine miecoeffmult(i1,i2,neqns,cx,cy)
         implicit none
         integer :: i1,i2,neqns,i,n,p,nodrvec(3),nodri,nblki,noffi,icon,q
         complex(8) :: cx(neqns),cy(neqns)
         complex(8), allocatable :: cxt(:,:,:),an1(:,:,:),cxtt(:,:,:)

         do i=i1,i2
            nodri=nodr(i)
            nblki=nblk(i)
            noffi=nblkoffset(i)
            allocate(cxt(0:nodri+1,nodri,2),cxtt(0:nodri+1,nodri,2),an1(2,2,nodri))
            cxt=reshape(cx(noffi+1:noffi+nblki),(/nodri+2,nodri,2/))
            cxtt=0.d0
            an1=miecoef(i)
            do n=1,nodri
               do p=1,2
                  cxtt(n+1,n:1:-1,p)=an1(p,1,n)*cxt(n+1,n:1:-1,1)+an1(p,2,n)*cxt(n+1,n:1:-1,2)
                  cxtt(0:n,n,p)=an1(p,1,n)*cxt(0:n,n,1)+an1(p,2,n)*cxt(0:n,n,2)
               enddo
            enddo
            cy(noffi+1:noffi+nblki)=reshape(cxtt(0:nodri+1,1:nodri,1:2),(/nblki/))
            deallocate(cxt,cxtt,an1)
         enddo
         end subroutine miecoeffmult

         subroutine internalmiecoeffmult(i1,i2,neqns,cx,cy)
         implicit none
         integer :: i1,i2,neqns,i,n,p,nodrvec(3),nodri,nblki,noffi,icon,q
         complex(8) :: cx(neqns),cy(neqns)
         complex(8), allocatable :: cxt(:,:,:),an1(:,:,:),cxtt(:,:,:)

         do i=i1,i2
            nodri=nodr(i)
            nblki=nblk(i)
            noffi=nblkoffset(i)
            allocate(cxt(0:nodri+1,nodri,2),cxtt(0:nodri+1,nodri,2),an1(2,2,nodri))
            cxt=reshape(cx(noffi+1:noffi+nblki),(/nodri+2,nodri,2/))
            cxtt=0.d0
            an1=internalmiecoef(n)
            do n=1,nodri
               do p=1,2
                  cxtt(n+1,n:1:-1,p)=an1(p,1,n)*cxt(n+1,n:1:-1,1)+an1(p,2,n)*cxt(n+1,n:1:-1,2)
                  cxtt(0:n,n,p)=an1(p,1,n)*cxt(0:n,n,1)+an1(p,2,n)*cxt(0:n,n,2)
               enddo
            enddo
            cy(noffi+1:noffi+nblki)=reshape(cxtt(0:nodri+1,1:nodri,1:2),(/nblki/))
            deallocate(cxt,cxtt,an1)
         enddo
         end subroutine internalmiecoeffmult
!
! single-sphere lorenz/mie coefficients
!
!
!  last revised: 15 January 2011
!
         subroutine mieregular(x,ri,nstop,qeps,qext,qsca,anp_mie,cnp_mie)
         use specialfuncs
         implicit none
         integer :: nstop,n,iancalc
         real(8) :: x,qeps,qext,qsca,prn,prp,qext1,err
         complex(8), optional :: anp_mie(2,*), cnp_mie(2,*)
         complex(8) :: ri,y,pcp,xip,da,db,na,nb,an1,an2,cn1,cn2
         complex(8), allocatable :: pc(:),xi(:)
!
!  modified LM criterion
!
         if(qeps.gt.0.) nstop=nint(x+4.*x**(1./3.))+15
!
!  user-set order limit
!
         if(qeps.lt.0) nstop=-qeps
!
!  basic calculations follow
!
         allocate(pc(0:nstop),xi(0:nstop))
         y=x*ri
         call cricbessel(nstop,y,pc)
         call richankel(nstop,x,xi)
         qsca=0.0
         qext=0.0
         do n=1,nstop
            prn=dble(xi(n))
            pcp=pc(n-1)-n*pc(n)/y
            xip=xi(n-1)-n*xi(n)/x
            prp=dble(xip)
            da=ri*xip*pc(n)-xi(n)*pcp
            db=ri*xi(n)*pcp-xip*pc(n)
            na=ri*prp*pc(n)-prn*pcp
            nb=ri*prn*pcp-prp*pc(n)
            an1=-na/da
            an2=-nb/db
            cn1=-dcmplx(0.d0,1.d0)*ri/na
            cn2=dcmplx(0.d0,1.d0)*ri/nb
            if(present(anp_mie)) then
               anp_mie(1,n)=an1
               anp_mie(2,n)=an2
            endif
            if(present(cnp_mie)) then
               cnp_mie(1,n)=cn1
               cnp_mie(2,n)=cn2
            endif
            qsca=qsca+(n+n+1)*(cdabs(an1)*cdabs(an1) &
                     +cdabs(an2)*cdabs(an2))
            qext1=-(n+n+1)*dble(an1+an2)
            qext=qext+qext1
            err=abs(qext1)/abs(qext)
            if(err.lt.qeps.or.n.eq.nstop) exit
         enddo
         nstop=n
         qsca=2./x/x*qsca
         qext=2./x/x*qext
         deallocate(pc,xi)
         end subroutine mieregular
!
! optically active lorenz/mie coefficients
! 30 March 2011
!
         subroutine mieoa(x,ri,nstop,qeps,qext,qsca,anp_mie,cnp_mie)
         use specialfuncs
         implicit none
         integer :: nstop
         real(8) :: x,qeps,qext,qsca,fn1,err
         complex(8) :: ri(2)
         complex(8), optional :: anp_mie(2,2,*),cnp_mie(2,2,*)
         integer :: n,i,p,q
         real(8) :: psi,psip,qext1
         complex (8) :: xri(2),xip,psicp,psic,wn(2),vn(2),an(2),bn(2), &
                        den,xi,anct(2,2),cnct(2,2),ri0,ci
         complex(8), allocatable :: psicn(:,:),xin(:)
         data ci/(0.d0,1.d0)/

         ri0=2.d0/(1.d0/ri(1)+1.d0/ri(2))
         if(qeps.ge.0.) then
            nstop=nint(x+4.*x**(1./3.))+5.
         else
            nstop=-qeps
         endif
         allocate(psicn(0:nstop+1,2),xin(0:nstop+1))
         do i=1,2
            xri(i)=x*ri(i)
            call cricbessel(nstop+1,xri(i),psicn(0,i))
         enddo
         call richankel(nstop+1,x,xin)
         qsca=0.0
         qext=0.0
         do n=1,nstop
            do i=1,2
               psic=psicn(n,i)
               psicp=psicn(n-1,i)-dble(n)*psic/xri(i)
               xi=xin(n)
               xip=xin(n-1)-dble(n)*xi/x
               psi=dble(xi)
               psip=dble(xip)
               wn(i)=ri0*psic*xip-xi*psicp
               vn(i)=psic*xip-ri0*xi*psicp
               an(i)=ri0*psic*psip-psi*psicp
               bn(i)=psic*psip-ri0*psi*psicp
            enddo
            den=wn(1)*vn(2)+wn(2)*vn(1)
            anct(1,1)=-(vn(1)*an(2)+vn(2)*an(1))/den
            anct(2,2)=-(wn(1)*bn(2)+wn(2)*bn(1))/den
            anct(1,2)=(wn(1)*an(2)-wn(2)*an(1))/den
            anct(2,1)=anct(1,2)
            den=an(1)*bn(2)+an(2)*bn(1)
            cnct(1,1)=-ci*ri(1)*bn(2)/den
            cnct(1,2)=-ci*ri(1)*an(2)/den
            cnct(2,1)=ri(2)*ri0*bn(1)/den
            cnct(2,2)=-ri(2)*ri0*an(1)/den
            if(present(anp_mie)) then
               do p=1,2
                  do q=1,2
                     anp_mie(p,q,n)=anct(p,q)
                     cnp_mie(p,q,n)=cnct(p,q)
                  enddo
               enddo
            endif
            qext1=0.d0
            fn1=n+n+1
            do p=1,2
               do q=1,2
                  qsca=qsca+fn1*cdabs(anct(p,q))*cdabs(anct(p,q))
               enddo
               qext1=qext1-fn1*dble(anct(p,p))
            enddo
            qext=qext+qext1
            err=abs(qext1)/abs(qext)
            if(err.lt.qeps.or.n.eq.nstop) exit
         enddo
         nstop=min(n,nstop)
         qsca=2./x/x*qsca
         qext=2./x/x*qext
         return
         end subroutine mieoa

      end module miecoefdata
!
!  module translation contains subroutines for VSWF translation and rotation
!
!
!  last revised: 15 January 2011
!
      module translation
      implicit none
      integer, private :: stored_max_order,store_tran_mat
      integer, allocatable, private :: nsizerot(:,:),nsizetran(:,:),nsizeephi(:,:), &
                 noffrot(:,:),nofftran(:,:),noffephi(:,:)
      real(8), private :: near_field_distance
      real(8), allocatable, private :: sphere_position(:,:)
      real(8), target, allocatable, private :: rotmatstore(:)
      complex(8), target, allocatable, private :: tranmatstore(:), ephimatstore(:)
      complex(8), allocatable, private :: rvec_temp(:,:),tvec_temp(:,:),c_temp(:,:,:), &
                  ct_temp(:,:,:),rvec2_temp(:,:),tvec2_temp(:,:),c2_temp(:,:,:), &
                                ct2_temp(:,:,:)

      contains
!
!  rotation of expansion coefficients amn by euler angles alpha,beta,gamma
!  idir=1: forward rotation, idir=-1, reverse rotation.
!
!
!  last revised: 15 January 2011
!
         subroutine rotvec(alpha,beta,gamma,nmax,mmax,amn,idir)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nmax,mmax,idir,k,n,m,in,kmax,kn,ka,na,p,im,m1
         real(8) :: dc(-nmax-1:nmax+1,-nmax-1:nmax+1),dk0(-nmax-1:nmax+1), &
                    dk01(-nmax-1:nmax+1),sbe,cbe,sbe2,cbe2,sben,dkt, &
                    fmn,dkm0,dkm1,alpha,beta,gamma
         complex(8) :: ealpha,amn(0:nmax+1,nmax,2),ealpham(-nmax:nmax), &
                       amnt(2,-nmax:nmax),a,b,ci,egamma,egammam(-nmax:nmax)
         data ci/(0.d0,1.d0)/
         call init(nmax)
         dc=0.d0
         dk01=0.d0
         dk0=0.d0
         ealpha=cdexp(ci*alpha)
         egamma=cdexp(ci*gamma)
         cbe=cos(beta)
         sbe=sqrt((1.d0+cbe)*(1.d0-cbe))
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         call ephicoef(ealpha,nmax,ealpham)
         call ephicoef(egamma,nmax,egammam)
         in=1
         dk0(0)=1.d0
         sben=1.d0
         dk01(0)=0.d0
         do n=1,nmax
            kmax=min(n,mmax)
            do k=-kmax,kmax
               if(k.le.-1) then
                  ka=n+1
                  na=-k
               else
                  ka=k
                  na=n
               endif
               if(idir.eq.1) then
                  amnt(1,k)=amn(ka,na,1)*ealpham(k)
                  amnt(2,k)=amn(ka,na,2)*ealpham(k)
               else
                  amnt(1,-k)=amn(ka,na,1)*egammam(k)
                  amnt(2,-k)=amn(ka,na,2)*egammam(k)
               endif
            enddo
            in=-in
            sben=sben*sbe/2.d0
            dk0(n)=in*sben*bcof(n,n)
            dk0(-n)=in*dk0(n)
            dk01(n)=0.d0
            dk01(-n)=0.d0
            dc(0,n)=dk0(n)
            dc(0,-n)=dk0(-n)
            do k=-n+1,n-1
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt) &
                     /(fnr(n+k)*fnr(n-k))
               dc(0,k)=dk0(k)
            enddo
            im=1
            do m=1,kmax
               im=-im
               fmn=1./fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.
               do k=-n,n
                  dkm1=dkm0
                  dkm0=dc(m1,k)
                  dc(m,k)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                         -fnr(n-k)*fnr(n+k+1)*sbe2*dc(m1,k+1) &
                         -k*sbe*dc(m1,k))*fmn
                  dc(-m,-k)=dc(m,k)*(-1)**(k)*im
               enddo
            enddo
            do m=-n,n
               if(m.le.-1) then
                  ka=n+1
                  na=-m
               else
                  ka=m
                  na=n
               endif
               a=0.
               b=0.
               do k=-kmax,kmax
                  a=a+dc(-k,-m)*amnt(1,k)
                  b=b+dc(-k,-m)*amnt(2,k)
               enddo
               if(idir.eq.1) then
                  amn(ka,na,1)=a*egammam(m)
                  amn(ka,na,2)=b*egammam(m)
               else
                  amn(ka,na,1)=a*ealpham(m)
                  amn(ka,na,2)=b*ealpham(m)
               endif
            enddo
         enddo
         end subroutine rotvec
!
!  sets up the stored translation matrices for mpi
!
!
!  last revised: 15 January 2011
!  november 2011: added near and far field translation
!
         subroutine mpirottranmtrxsetup(nsphere,nodr,rpos,ri,istore,nfdistance,&
                    runprintunit)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         implicit none
         integer :: nsphere,nodr(nsphere),i,j,nodrmax,nodrmin,n,ntotrot,ntottran,ntotephi, &
                    ierr,n1,n2,nt,rank,nsrank,runprintunit,isendok,tag,sendrank,numprocs,brank, &
                    nsend,istore
         real(8) :: rpos(3,nsphere),xij(3),r,ct,memused(1),memusedmax(1),memusedmin(1), &
                    nfdistance,nfdistancei
         real(8), allocatable :: rotmat(:,:)
         complex(8) :: ri,ephi
         complex(8), allocatable :: tranmat(:,:,:),ephimat(:),pivec(:,:,:)
         data isendok,tag/0,1/
         numprocs=proc_per_group
         rank=group_rank
         brank=base_rank
         nsrank=mpi_sphere_number(rank)
         nodrmax=maxval(nodr)
         call init(nodrmax)
         store_tran_mat=istore
         near_field_distance=nfdistance
         if(allocated(sphere_position)) deallocate(sphere_position)
         allocate(sphere_position(3,nsphere))
         sphere_position=rpos
         if(istore.eq.0) then
            return
         endif
         if(allocated(nsizerot)) deallocate(nsizerot,nsizetran,nsizeephi,noffrot,nofftran,noffephi)
         allocate(nsizerot(nsphere,nsphere),nsizetran(nsphere,nsphere),nsizeephi(nsphere,nsphere), &
                  noffrot(nsphere,nsphere),nofftran(nsphere,nsphere),noffephi(nsphere,nsphere))
         if(allocated(rvec_temp)) deallocate(rvec_temp,tvec_temp,c_temp,ct_temp, &
                     rvec2_temp,tvec2_temp,c2_temp,ct2_temp)
         allocate(rvec_temp(-nodrmax:nodrmax,2),tvec_temp(nodrmax,2), &
                  c_temp(-nodrmax:nodrmax,nodrmax,2),ct_temp(nodrmax,2,2), &
                       rvec2_temp(-nodrmax:nodrmax,2),tvec2_temp(nodrmax,2), &
                       c2_temp(-nodrmax:nodrmax,nodrmax,2),ct2_temp(nodrmax,2,2))
         stored_max_order=nodrmax
!
!  determine the memory requirements
!
         ntotrot=0
         ntottran=0
         ntotephi=0
         do i=mpi_sphere_index(rank)+1,mpi_sphere_index(rank)+nsrank
            do j=1,nsphere
               xij(:)=rpos(:,i)-rpos(:,j)
               if(j.ne.i) then
                  if(nfdistance.lt.0.) then
                     nfdistancei=(.5*dble(nodr(i)+nodr(j)))**2.
                  else
                     nfdistancei=nfdistance
                  endif
                  r=sqrt(dot_product(xij,xij))
                  if(r.le.nfdistancei) then
                     nodrmax=max(nodr(j),nodr(i))
                     nodrmin=min(nodr(j),nodr(i))
                     noffrot(i,j)=ntotrot
                     nofftran(i,j)=ntottran
                     noffephi(i,j)=ntotephi
                     nsizerot(i,j)=(2*nodrmin+1)*(1+nodrmax*(nodrmax+2))
                     nsizetran(i,j)=nodr(i)*nodr(j)*(nodr(j)+3)
                     nsizeephi(i,j)=2*nodrmax+1
                     ntotrot=ntotrot+nsizerot(i,j)
                     ntottran=ntottran+nsizetran(i,j)
                     ntotephi=ntotephi+nsizeephi(i,j)
                  endif
                  if(r.gt.nfdistancei.and.istore.eq.2) then
                     nodrmax=max(nodr(j),nodr(i))
                     nofftran(i,j)=ntottran
                     nsizetran(i,j)=2*nodrmax*(nodrmax+2)
                     ntottran=ntottran+nsizetran(i,j)
                  endif
               endif
            enddo
         enddo
         memused(1)=dble(8*ntotrot+16*(ntottran+ntotephi))*1.d-6
         nsend=1
         call ms_mpi(mpi_command='reduce',mpi_send_buf_dp=memused,mpi_recv_buf_dp=memusedmax,&
                     mpi_number=1,mpi_rank=0,mpi_operation=ms_mpi_max)
         call ms_mpi(mpi_command='reduce',mpi_send_buf_dp=memused,mpi_recv_buf_dp=memusedmin,&
                     mpi_number=1,mpi_rank=0,mpi_operation=ms_mpi_min)
         call ms_mpi(mpi_command='barrier')
         if(brank.eq.0) then
            write(runprintunit,'('' maximum translation matrix storage:'',f9.4,'' MB'')') memusedmax
            write(runprintunit,'('' minimum translation matrix storage:'',f9.4,'' MB'')') memusedmin
            call flush(runprintunit)
         endif
!
!  calculate the matrices and store in memory
!
         if(allocated(rotmatstore)) deallocate(rotmatstore,tranmatstore,ephimatstore)
         allocate(rotmatstore(ntotrot),stat=ierr)
         allocate(tranmatstore(ntottran),stat=ierr)
         allocate(ephimatstore(ntotephi),stat=ierr)
         do i=mpi_sphere_index(rank)+1,mpi_sphere_index(rank)+nsrank
            do j=1,nsphere
               if(j.ne.i) then
                  nodrmax=max(nodr(j),nodr(i))
                  nodrmin=min(nodr(j),nodr(i))
                  xij=rpos(:,i)-rpos(:,j)
                  call cartosphere(xij,r,ct,ephi)
                  if(nfdistance.lt.0.) then
                     nfdistancei=(.5*dble(nodr(i)+nodr(j)))**2.
                  else
                     nfdistancei=nfdistance
                  endif
                  if(r.le.nfdistancei) then
!
!  rotation matrix
!
                     n1=noffrot(i,j)+1
                     nt=nsizerot(i,j)
                     n2=n1+nt-1
                     allocate(rotmat(-nodrmin:nodrmin,0:nodrmax*(nodrmax+2)))
                     call rotcoef(ct,nodrmin,nodrmax,rotmat)
                     rotmatstore(n1:n2)=reshape(rotmat,(/nt/))
                     deallocate(rotmat)
!
!  axial translation matrix
!
                     n1=nofftran(i,j)+1
                     nt=nsizetran(i,j)
                     n2=n1+nt-1
                     allocate(tranmat(nodr(i),nodr(j)*(nodr(j)+3)/2,2))
                     call axialtrancoef(3,r,ri,nodr(i),nodr(j),tranmat)
                     tranmatstore(n1:n2)=reshape(tranmat,(/nt/))
                     deallocate(tranmat)
!
!  ephi matrix
!
                     n1=noffephi(i,j)+1
                     nt=nsizeephi(i,j)
                     n2=n1+nt-1
                     allocate(ephimat(-nodrmax:nodrmax))
                     call ephicoef(ephi,nodrmax,ephimat)
                     ephimatstore(n1:n2)=ephimat(-nodrmax:nodrmax)
                     deallocate(ephimat)
!
!  ff translation matrix storage
!
                  elseif(istore.eq.2) then
                     n1=nofftran(i,j)+1
                     nt=nsizetran(i,j)
                     n2=n1+nt-1
                     nodrmax=max(nodr(j),nodr(i))
                     allocate(pivec(0:nodrmax+1,nodrmax,2))
                     call pifunc(ct,ephi,nodrmax,nodrmax,pivec)
                     tranmatstore(n1:n2)=reshape(pivec,(/nt/))
                     deallocate(pivec)
                  endif
               endif
            enddo
         enddo
         end subroutine mpirottranmtrxsetup
!
!  clear the stored translation matrices
!
!
!  last revised: 15 January 2011
!
         subroutine rottranmtrxclear()
         implicit none
         if(allocated(rotmatstore)) deallocate(rotmatstore,tranmatstore,ephimatstore)
         if(allocated(sphere_position)) deallocate(sphere_position)
         end subroutine rottranmtrxclear
!
!  translation coefficient vector cx by xij in medium with ri by rotation-translation
!  itype: 1 or 3
!  icalc: =1, calculate matrices; = 0, use stored matrix
!  idir: =1, translation of xij, =-1, -xij (reverse)
!  itran=1, A(i-j) a(j), = -1, a(j) A(i-j)
!
!
!  last revised: 15 January 2011
!
         subroutine rottran(cx,cy,xij,ri,nodrx,nodry,itype,icalc,idir,itran)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,itype,icalc,idir,itran,nmax,nmin,n,m,p,nblk
         real(8) :: xij(3),r,ct
         complex(8) :: ri,ephi,cx(0:nodrx+1,nodrx,2),cy(0:nodry+1,nodry,2)
         real(8), allocatable, save :: rotmat(:,:)
         complex(8), allocatable, save  :: ephimat(:), tranmat(:,:,:)
         if(icalc.eq.1) then
            nmax=max(nodrx,nodry)
            nmin=min(nodrx,nodry)
            call cartosphere(xij,r,ct,ephi)
            if(r.lt.1.d-4) then
               do p=1,2
                  do n=1,nmin
                     do m=0,nmin+1
                        cy(m,n,p)=cy(m,n,p)+cx(m,n,p)
                     enddo
                  enddo
               enddo
               return
            endif
            if(allocated(ephimat)) deallocate(rotmat,ephimat,tranmat)
            if(nmax.gt.stored_max_order) then
               if(allocated(rvec_temp)) deallocate(rvec_temp,tvec_temp, &
                       c_temp,ct_temp,rvec2_temp,tvec2_temp, &
                       c2_temp,ct2_temp)
               allocate(rvec_temp(-nmax:nmax,2),tvec_temp(nmax,2), &
                       c_temp(-nmax:nmax,nmax,2),ct_temp(nmax,2,2), &
                       rvec2_temp(-nmax:nmax,2),tvec2_temp(nmax,2), &
                       c2_temp(-nmax:nmax,nmax,2),ct2_temp(nmax,2,2))
               stored_max_order=nmax
            endif
            nblk=(nodrx*(nodrx+3))/2
            allocate(rotmat(-nmin:nmin,0:nmax*(nmax+2)))
            allocate(ephimat(-nmax:nmax))
            allocate(tranmat(1:nodry,1:nblk,1:2))
            call rotcoef(ct,nmin,nmax,rotmat)
!            call axialtrancoef(itype,r,ri,nodry,nodrx,tranmat)
            call axialtrancoefrecurrence(itype,r,ri,nodry,nodrx,tranmat)
            call ephicoef(ephi,nmax,ephimat)
         endif
         call rottranmtrx(cx,cy,idir,itran,nodrx,nodry,ephimat,rotmat,tranmat)
         return
         end subroutine rottran
!
!  far field formula for outgoing SVWF translation
!  October 2011
!
         subroutine farfieldtranslation(cx,cy,xij,ri,nodrx,nodry,icase, &
                    stored_pivec_matrix)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,itype,icalc,icase,nmax,nmin,n,m,p,nblk,im
         real(8) :: xij(3),r,ct,xijt(3)
         complex(8) :: ri,ephi,cx(0:nodrx+1,nodrx,2),cy(0:nodry+1,nodry,2), &
                       cxt(0:nodrx+1,nodrx,2),cyt(0:nodry+1,nodry,2), &
                       sumx(2),c1,pivec(0:max(nodrx,nodry)+1,max(nodrx,nodry),2)
         complex(8), optional :: stored_pivec_matrix(0:max(nodrx,nodry)+1,max(nodrx,nodry),2)

         call cartosphere(xij,r,ct,ephi)
         nmax=max(nodrx,nodry)
         if(present(stored_pivec_matrix)) then
            pivec=stored_pivec_matrix
         else
            call pifunc(ct,ephi,nmax,nmax,pivec)
         endif
         if(icase.eq.1) then
            sumx(1)=sum(pivec(0:nodrx+1,1:nodrx,1:2)*cx(0:nodrx+1,1:nodrx,1:2))
            sumx(2)=sum(pivec(0:nodrx+1,1:nodrx,2:1:-1)*cx(0:nodrx+1,1:nodrx,1:2))
            sumx=sumx*cdexp((0.d0,1.d0)*ri*r)/((0.d0,1.d0)*ri*r)*8.d0
            cyt(0:nodry+1,1:nodry,1) = conjg(pivec(0:nodry+1,1:nodry,1))*sumx(1) &
                  +conjg(pivec(0:nodry+1,1:nodry,2))*sumx(2)
            cyt(0:nodry+1,1:nodry,2) = conjg(pivec(0:nodry+1,1:nodry,2))*sumx(1) &
                  +conjg(pivec(0:nodry+1,1:nodry,1))*sumx(2)
         else
            do n=1,nodrx
               do p=1,2
                  im=(-1)**(n+p)
                  cxt(n+1,1:n,p)=im*cx(n+1,1:n,p)
                  cxt(0:n,n,p)=im*cx(0:n,n,p)
               enddo
            enddo
            sumx(1)=sum(conjg(pivec(0:nodrx+1,1:nodrx,1:2))*cxt(0:nodrx+1,1:nodrx,1:2))
            sumx(2)=sum(conjg(pivec(0:nodrx+1,1:nodrx,2:1:-1))*cxt(0:nodrx+1,1:nodrx,1:2))
            sumx=sumx*cdexp((0.d0,1.d0)*ri*r)/((0.d0,1.d0)*ri*r)*8.d0
            cyt(0:nodry+1,1:nodry,1) = pivec(0:nodry+1,1:nodry,1)*sumx(1) &
                  +pivec(0:nodry+1,1:nodry,2)*sumx(2)
            cyt(0:nodry+1,1:nodry,2) = pivec(0:nodry+1,1:nodry,2)*sumx(1) &
                  +pivec(0:nodry+1,1:nodry,1)*sumx(2)
            do n=1,nodry
               do p=1,2
                  im=(-1)**(n+p)
                  cyt(n+1,1:n,p)=im*cyt(n+1,1:n,p)
                  cyt(0:n,n,p)=im*cyt(0:n,n,p)
               enddo
            enddo
         endif
         cy=cy+cyt
         end subroutine farfieldtranslation
!
! far field translation: normal and transpose, for bcgm solution
! october 2011
!
         subroutine farfieldtranslationtwovec(cx1,cx2,cy1,cy2,xij,ri,nodrx,nodry, &
                    stored_pivec_matrix)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,itype,icalc,icase,nmax,nmin,n,m,p,nblk,im
         real(8) :: xij(3),r,ct,xijt(3)
         complex(8) :: ri,ephi,cx1(0:nodrx+1,nodrx,2),cy1(0:nodry+1,nodry,2), &
                       cx2(0:nodrx+1,nodrx,2),cy2(0:nodry+1,nodry,2), &
                       cxt(0:nodrx+1,nodrx,2),cyt1(0:nodry+1,nodry,2), &
                       cyt2(0:nodry+1,nodry,2), &
                       sumx(2),c1,phasefunc, &
                       pivec(0:max(nodrx,nodry)+1,max(nodrx,nodry),2)
         complex(8), optional :: stored_pivec_matrix(0:max(nodrx,nodry)+1,max(nodrx,nodry),2)

         call cartosphere(xij,r,ct,ephi)
         nmax=max(nodrx,nodry)
         if(present(stored_pivec_matrix)) then
            pivec=stored_pivec_matrix
         else
            call pifunc(ct,ephi,nmax,nmax,pivec)
         endif
         phasefunc=cdexp((0.d0,1.d0)*ri*r)/((0.d0,1.d0)*ri*r)*8.d0
         sumx(1)=sum(pivec(0:nodrx+1,1:nodrx,1:2)*cx1(0:nodrx+1,1:nodrx,1:2))
         sumx(2)=sum(pivec(0:nodrx+1,1:nodrx,2:1:-1)*cx1(0:nodrx+1,1:nodrx,1:2))
         sumx=sumx*phasefunc
         cyt1(0:nodry+1,1:nodry,1) = conjg(pivec(0:nodry+1,1:nodry,1))*sumx(1) &
               +conjg(pivec(0:nodry+1,1:nodry,2))*sumx(2)
         cyt1(0:nodry+1,1:nodry,2) = conjg(pivec(0:nodry+1,1:nodry,2))*sumx(1) &
               +conjg(pivec(0:nodry+1,1:nodry,1))*sumx(2)
         do n=1,nodrx
            do p=1,2
               im=(-1)**(n+p)
               cxt(n+1,1:n,p)=im*cx2(n+1,1:n,p)
               cxt(0:n,n,p)=im*cx2(0:n,n,p)
            enddo
         enddo
         sumx(1)=sum(conjg(pivec(0:nodrx+1,1:nodrx,1:2))*cxt(0:nodrx+1,1:nodrx,1:2))
         sumx(2)=sum(conjg(pivec(0:nodrx+1,1:nodrx,2:1:-1))*cxt(0:nodrx+1,1:nodrx,1:2))
         sumx=sumx*phasefunc
         cyt2(0:nodry+1,1:nodry,1) = pivec(0:nodry+1,1:nodry,1)*sumx(1) &
               +pivec(0:nodry+1,1:nodry,2)*sumx(2)
         cyt2(0:nodry+1,1:nodry,2) = pivec(0:nodry+1,1:nodry,2)*sumx(1) &
               +pivec(0:nodry+1,1:nodry,1)*sumx(2)
         do n=1,nodry
            do p=1,2
               im=(-1)**(n+p)
               cyt2(n+1,1:n,p)=im*cyt2(n+1,1:n,p)
               cyt2(0:n,n,p)=im*cyt2(0:n,n,p)
            enddo
         enddo
         cy1=cy1+cyt1
         cy2=cy2+cyt2
         end subroutine farfieldtranslationtwovec
!
! correction term for hybrid bcgm solution: difference between exact and
! ff translation field
! november 2011
!
         subroutine fftranslationerror(cx,cy,jx,iy,nodrx,nodry)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,idir,itran,iy,jx,istore
         integer :: nr1,nr2,nt1,nt2,ne1,ne2
         real(8) :: xj(3),xi(3),xij(3),rij,nfdist
         complex(8) :: cx(0:nodrx+1,nodrx,2),cy(0:nodry+1,nodry,2), &
                       cyt(0:nodry+1,nodry,2)
         xj(:)=sphere_position(:,jx)
         xi(:)=sphere_position(:,iy)
         xij=xi-xj
         rij=sqrt(dot_product(xij,xij))
         if(near_field_distance.lt.0.) then
            nfdist=(.5*(nodrx+nodry))**2.
         else
            nfdist=near_field_distance
         endif
         if(rij.gt.nfdist) then
            cyt=0.d0
            call farfieldtranslation(cx,cyt,xij,(1.d0,0.d0),nodrx,nodry,1)
            cyt=-cyt
            call rottran(cx,cyt,xij,(1.d0,0.d0),nodrx,nodry,3,1,1,1)
            cy=cy+cyt
         endif
         end subroutine fftranslationerror
!
!  translation via stored or calculated matrices (replaces rottranstoredmatrix)
!
!  12 October 2011.
!     if rij> near_field_distance, the far field formula is
!     applied.
!
         subroutine rottranjtoi(cx,cy,jx,iy,nodrx,nodry,idir,itran)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,idir,itran,iy,jx,istore
         integer :: nr1,nr2,nt1,nt2,ne1,ne2
         real(8) :: xj(3),xi(3),xij(3),rij,nfdist
         complex(8) :: cx(0:nodrx+1,nodrx,2),cy(0:nodry+1,nodry,2)
         xj(:)=sphere_position(:,jx)
         xi(:)=sphere_position(:,iy)
         xij=xi-xj
         rij=sqrt(dot_product(xij,xij))
         if(near_field_distance.lt.0.) then
            nfdist=(.5*(nodrx+nodry))**2.
         else
            nfdist=near_field_distance
         endif
         if(rij.gt.nfdist) then
            if(store_tran_mat.eq.2) then
               nt1=nofftran(iy,jx)+1
               nt2=nt1+nsizetran(iy,jx)-1
               call farfieldtranslation(cx,cy,xij,(1.d0,0.d0),nodrx,nodry,itran, &
                    stored_pivec_matrix=tranmatstore(nt1:nt2))
            else
               call farfieldtranslation(cx,cy,xij,(1.d0,0.d0),nodrx,nodry,itran)
            endif
         else
            if(store_tran_mat.eq.0) then
               call rottran(cx,cy,xij,(1.d0,0.d0),nodrx,nodry,3,1,idir,itran)
            else
               nr1=noffrot(iy,jx)+1
               nr2=nr1+nsizerot(iy,jx)-1
               nt1=nofftran(iy,jx)+1
               nt2=nt1+nsizetran(iy,jx)-1
               ne1=noffephi(iy,jx)+1
               ne2=ne1+nsizeephi(iy,jx)-1
               call rottranmtrx(cx,cy,idir,itran,nodrx,nodry,ephimatstore(ne1:ne2), &
                                rotmatstore(nr1:nr2),tranmatstore(nt1:nt2))
            endif
         endif
         end subroutine rottranjtoi
!
! normal and transpose translation, for bcgm
! november 2011
!
         subroutine rottrantwojtoi(cx1,cx2,cy1,cy2,jx,iy,nodrx,nodry)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,idir,itran,iy,jx,istore
         integer :: nr1,nr2,nt1,nt2,ne1,ne2
         real(8) :: xj(3),xi(3),xij(3),rij,nfdist
         complex(8) :: cx1(0:nodrx+1,nodrx,2),cy1(0:nodry+1,nodry,2), &
                       cx2(0:nodrx+1,nodrx,2),cy2(0:nodry+1,nodry,2)
         xj(:)=sphere_position(:,jx)
         xi(:)=sphere_position(:,iy)
         xij=xi-xj
         rij=sqrt(dot_product(xij,xij))
         if(near_field_distance.lt.0.) then
            nfdist=(.5*(nodrx+nodry))**2.
         else
            nfdist=near_field_distance
         endif
         if(rij.gt.nfdist) then
            if(store_tran_mat.eq.2) then
               nt1=nofftran(iy,jx)+1
               nt2=nt1+nsizetran(iy,jx)-1
               call farfieldtranslationtwovec(cx1,cx2,cy1,cy2,xij,(1.d0,0.d0),nodrx,nodry, &
               stored_pivec_matrix=tranmatstore(nt1:nt2))
            else
               call farfieldtranslationtwovec(cx1,cx2,cy1,cy2,xij,(1.d0,0.d0),nodrx,nodry)
            endif
         else
            if(store_tran_mat.eq.0) then
               call rottran(cx1,cy1,xij,(1.d0,0.d0),nodrx,nodry,3,1,1,1)
               call rottran(cx2,cy2,xij,(1.d0,0.d0),nodrx,nodry,3,0,-1,-1)
            else
               nr1=noffrot(iy,jx)+1
               nr2=nr1+nsizerot(iy,jx)-1
               nt1=nofftran(iy,jx)+1
               nt2=nt1+nsizetran(iy,jx)-1
               ne1=noffephi(iy,jx)+1
               ne2=ne1+nsizeephi(iy,jx)-1
               call rottranmtrxtwovec(cx1,cx2,cy1,cy2,nodrx,nodry,ephimatstore(ne1:ne2), &
                                rotmatstore(nr1:nr2),tranmatstore(nt1:nt2))
            endif
         endif
         end subroutine rottrantwojtoi
!
!  the vectorized rotation-translation-rotation operation
!
!
!  last revised: 15 January 2011
!
         subroutine rottranmtrx(cx,cy,idir,itran,nodrx,nodry,ephimat,rotmat,tranmat)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,itype,icalc,idir,itran,nmax,nmin
         integer :: m,n,k,l,nn1,nn2,ll1,mn,kl,m1,p,n1,addr(2)
         real(8), target :: rotmat(-min(nodrx,nodry):min(nodrx,nodry), &
                            0:max(nodrx,nodry)*(max(nodrx,nodry)+2))
         real(8), pointer :: rmat(:,:)
         complex(8) :: cx(0:nodrx+1,nodrx,2),cy(0:nodry+1,nodry,2), &
                       ephimat(-max(nodrx,nodry):max(nodrx,nodry))
         complex(8), target :: tranmat(nodry,nodrx*(nodrx+3)/2,2)
         complex(8), pointer :: tmat1(:,:),tmat2(:,:)
         c_temp=(0.d0,0.d0)
         nmin=min(nodrx,nodry)
!
!  rotation to origin of target
!
         do n=1,nodrx
            nn1=n*(n+1)-n
            nn2=nn1+(2*n+1)-1
            n1=min(n,nodry)
            rmat=>rotmat(-n1:n1,nn1:nn2)
            do p=1,2
               rvec_temp(-n:-1,p)=cx(n+1,n:1:-1,p)
               rvec_temp(0:n,p)=cx(0:n,n,p)
               if(itran.eq.1) then
                  rvec_temp(-n:n,p)=rvec_temp(-n:n,p)*ephimat(-n:n)
               else
                  rvec_temp(-n:n,p)=rvec_temp(-n:n,p)*conjg(ephimat(-n:n))
               endif
            enddo
            c_temp(-n1:n1,n,1:2)=matmul(rmat,rvec_temp(-n:n,1:2))
         enddo
!
!  axial translation to target
!
         do m=0,nmin
            m1=max(1,m)
            nn1=atcadd(m,m1,nodrx)
            nn2=atcadd(m,nodrx,nodrx)
            tmat1=>tranmat(m1:nodry,nn1:nn2,1)
            tmat2=>tranmat(m1:nodry,nn1:nn2,2)
            tvec_temp(m1:nodrx,1)=idir*c_temp(m,m1:nodrx,1)
            tvec_temp(m1:nodrx,2)=c_temp(m,m1:nodrx,2)
            if(itran*idir.eq.-1) then
               tvec_temp(m1:nodrx,1)=tvec_temp(m1:nodrx,1)*monen(m1:nodrx)
               tvec_temp(m1:nodrx,2)=tvec_temp(m1:nodrx,2)*monen(m1:nodrx)
            endif
            ct_temp=(0.d0,0.d0)
            ct_temp(m1:nodry,1,1:2)=matmul(tmat1,tvec_temp(m1:nodrx,1:2))
            ct_temp(m1:nodry,2,1:2)=matmul(tmat2,tvec_temp(m1:nodrx,1:2))
            c_temp(m,m1:nodry,1)=idir*(ct_temp(m1:nodry,1,1)+ct_temp(m1:nodry,2,2))
            c_temp(m,m1:nodry,2)=ct_temp(m1:nodry,2,1)+ct_temp(m1:nodry,1,2)
            if(itran*idir.eq.-1) then
               c_temp(m,m1:nodry,1)=c_temp(m,m1:nodry,1)*monen(m1:nodry)
               c_temp(m,m1:nodry,2)=c_temp(m,m1:nodry,2)*monen(m1:nodry)
            endif
            if(m.gt.0) then
               tvec_temp(m1:nodrx,1)=idir*c_temp(-m,m1:nodrx,1)
               tvec_temp(m1:nodrx,2)=c_temp(-m,m1:nodrx,2)
               if(itran*idir.eq.-1) then
                  tvec_temp(m1:nodrx,1)=tvec_temp(m1:nodrx,1)*monen(m1:nodrx)
                  tvec_temp(m1:nodrx,2)=tvec_temp(m1:nodrx,2)*monen(m1:nodrx)
               endif
               ct_temp=(0.d0,0.d0)
               ct_temp(m1:nodry,1,1:2)=matmul(tmat1,tvec_temp(m1:nodrx,1:2))
               ct_temp(m1:nodry,2,1:2)=matmul(tmat2,tvec_temp(m1:nodrx,1:2))
               c_temp(-m,m1:nodry,1)=idir*(ct_temp(m1:nodry,1,1)-ct_temp(m1:nodry,2,2))
               c_temp(-m,m1:nodry,2)=-ct_temp(m1:nodry,2,1)+ct_temp(m1:nodry,1,2)
               if(itran*idir.eq.-1) then
                  c_temp(-m,m1:nodry,1)=c_temp(-m,m1:nodry,1)*monen(m1:nodry)
                  c_temp(-m,m1:nodry,2)=c_temp(-m,m1:nodry,2)*monen(m1:nodry)
               endif
            endif
         enddo
!
!  rotation back to original frame
!
         do n=1,nodry
            rvec_temp=(0.d0,0.d0)
            m1=min(n,nmin)
            nn1=n*(n+1)-n
            nn2=n*(n+1)+n
            rmat=>rotmat(-m1:m1,nn1:nn2)
            do p=1,2
               if(itran.eq.1) then
                  rvec_temp(-n:n,p)=matmul(c_temp(-m1:m1,n,p),rmat)*conjg(ephimat(-n:n))
               else
                  rvec_temp(-n:n,p)=matmul(c_temp(-m1:m1,n,p),rmat)*ephimat(-n:n)
               endif
               cy(n+1,n:1:-1,p)=cy(n+1,n:1:-1,p)+rvec_temp(-n:-1,p)
               cy(0:n,n,p)=cy(0:n,n,p)+rvec_temp(0:n,p)
            enddo
         enddo
         end subroutine rottranmtrx
!
! two vector rotation: normal and transpose
! november 2011
!
         subroutine rottranmtrxtwovec(cx1,cx2,cy1,cy2,nodrx,nodry, &
          ephimat,rotmat,tranmat)
         use numconstants
         use specialfuncs
         implicit none
         integer :: nodrx,nodry,itype,icalc,idir,itran,nmax,nmin
         integer :: m,n,k,l,nn1,nn2,ll1,mn,kl,m1,p,n1,addr(2)
         real(8), target :: rotmat(-min(nodrx,nodry):min(nodrx,nodry), &
                            0:max(nodrx,nodry)*(max(nodrx,nodry)+2))
         real(8), pointer :: rmat(:,:)
         complex(8) :: cx1(0:nodrx+1,nodrx,2),cy1(0:nodry+1,nodry,2), &
                       cx2(0:nodrx+1,nodrx,2),cy2(0:nodry+1,nodry,2), &
                       ephimat(-max(nodrx,nodry):max(nodrx,nodry))
         complex(8), target :: tranmat(nodry,nodrx*(nodrx+3)/2,2)
         complex(8), pointer :: tmat1(:,:),tmat2(:,:)
         c_temp=(0.d0,0.d0)
         nmin=min(nodrx,nodry)
!
!  rotation to origin of target
!
         do n=1,nodrx
            nn1=n*(n+1)-n
            nn2=nn1+(2*n+1)-1
            n1=min(n,nodry)
            rmat=>rotmat(-n1:n1,nn1:nn2)
            do p=1,2
               rvec_temp(-n:-1,p)=cx1(n+1,n:1:-1,p)
               rvec_temp(0:n,p)=cx1(0:n,n,p)
               rvec2_temp(-n:-1,p)=cx2(n+1,n:1:-1,p)
               rvec2_temp(0:n,p)=cx2(0:n,n,p)
               rvec_temp(-n:n,p)=rvec_temp(-n:n,p)*ephimat(-n:n)
               rvec2_temp(-n:n,p)=rvec2_temp(-n:n,p)*conjg(ephimat(-n:n))
            enddo
            c_temp(-n1:n1,n,1:2)=matmul(rmat,rvec_temp(-n:n,1:2))
            c2_temp(-n1:n1,n,1:2)=matmul(rmat,rvec2_temp(-n:n,1:2))
         enddo
!
!  axial translation to target
!
         do m=0,nmin
            m1=max(1,m)
            nn1=atcadd(m,m1,nodrx)
            nn2=atcadd(m,nodrx,nodrx)
            tmat1=>tranmat(m1:nodry,nn1:nn2,1)
            tmat2=>tranmat(m1:nodry,nn1:nn2,2)
            tvec_temp(m1:nodrx,1)=c_temp(m,m1:nodrx,1)
            tvec_temp(m1:nodrx,2)=c_temp(m,m1:nodrx,2)
            tvec2_temp(m1:nodrx,1)=-c2_temp(m,m1:nodrx,1)
            tvec2_temp(m1:nodrx,2)=c2_temp(m,m1:nodrx,2)
            ct_temp=(0.d0,0.d0)
            ct2_temp=(0.d0,0.d0)
            ct_temp(m1:nodry,1,1:2)=matmul(tmat1,tvec_temp(m1:nodrx,1:2))
            ct_temp(m1:nodry,2,1:2)=matmul(tmat2,tvec_temp(m1:nodrx,1:2))
            ct2_temp(m1:nodry,1,1:2)=matmul(tmat1,tvec2_temp(m1:nodrx,1:2))
            ct2_temp(m1:nodry,2,1:2)=matmul(tmat2,tvec2_temp(m1:nodrx,1:2))
            c_temp(m,m1:nodry,1)=(ct_temp(m1:nodry,1,1)+ct_temp(m1:nodry,2,2))
            c_temp(m,m1:nodry,2)=ct_temp(m1:nodry,2,1)+ct_temp(m1:nodry,1,2)
            c2_temp(m,m1:nodry,1)=-(ct2_temp(m1:nodry,1,1)+ct2_temp(m1:nodry,2,2))
            c2_temp(m,m1:nodry,2)=ct2_temp(m1:nodry,2,1)+ct2_temp(m1:nodry,1,2)
            if(m.gt.0) then
               tvec_temp(m1:nodrx,1)=c_temp(-m,m1:nodrx,1)
               tvec_temp(m1:nodrx,2)=c_temp(-m,m1:nodrx,2)
               tvec2_temp(m1:nodrx,1)=-c2_temp(-m,m1:nodrx,1)
               tvec2_temp(m1:nodrx,2)=c2_temp(-m,m1:nodrx,2)
               ct_temp=(0.d0,0.d0)
               ct2_temp=(0.d0,0.d0)
               ct_temp(m1:nodry,1,1:2)=matmul(tmat1,tvec_temp(m1:nodrx,1:2))
               ct_temp(m1:nodry,2,1:2)=matmul(tmat2,tvec_temp(m1:nodrx,1:2))
               ct2_temp(m1:nodry,1,1:2)=matmul(tmat1,tvec2_temp(m1:nodrx,1:2))
               ct2_temp(m1:nodry,2,1:2)=matmul(tmat2,tvec2_temp(m1:nodrx,1:2))
               c_temp(-m,m1:nodry,1)=(ct_temp(m1:nodry,1,1)-ct_temp(m1:nodry,2,2))
               c_temp(-m,m1:nodry,2)=-ct_temp(m1:nodry,2,1)+ct_temp(m1:nodry,1,2)
               c2_temp(-m,m1:nodry,1)=-(ct2_temp(m1:nodry,1,1)-ct2_temp(m1:nodry,2,2))
               c2_temp(-m,m1:nodry,2)=-ct2_temp(m1:nodry,2,1)+ct2_temp(m1:nodry,1,2)
            endif
         enddo
!
!  rotation back to original frame
!
         do n=1,nodry
            rvec_temp=(0.d0,0.d0)
            rvec2_temp=(0.d0,0.d0)
            m1=min(n,nmin)
            nn1=n*(n+1)-n
            nn2=n*(n+1)+n
            rmat=>rotmat(-m1:m1,nn1:nn2)
            do p=1,2
               rvec_temp(-n:n,p)=matmul(c_temp(-m1:m1,n,p),rmat)*conjg(ephimat(-n:n))
               rvec2_temp(-n:n,p)=matmul(c2_temp(-m1:m1,n,p),rmat)*ephimat(-n:n)
               cy1(n+1,n:1:-1,p)=cy1(n+1,n:1:-1,p)+rvec_temp(-n:-1,p)
               cy1(0:n,n,p)=cy1(0:n,n,p)+rvec_temp(0:n,p)
               cy2(n+1,n:1:-1,p)=cy2(n+1,n:1:-1,p)+rvec2_temp(-n:-1,p)
               cy2(0:n,n,p)=cy2(0:n,n,p)+rvec2_temp(0:n,p)
            enddo
         enddo
         end subroutine rottranmtrxtwovec
!
!  GB coefficients for sphere-centered expansions, obtained via translation
!
!  last revised: 15 January 2011
!
         subroutine spheregaussianbeamcoef(nsphere,neqns,nodr,alpha,beta,cbeam, &
                    rpos,rbeam,epstran,pmnp)
         use specialfuncs
         implicit none
         integer :: m,n,p,nsphere,i,l,nodr(nsphere),nblk,noff,nodrgb,neqns,k
         real(8) :: alpha,beta,cb,sb,ca,sa,rpos(3,nsphere),rmax,rbeam(3),xib(3),rib, &
                    cbeam,epstran
         complex(8) :: pmnp(neqns,2)
         complex(8), allocatable :: pmnp0(:,:,:,:)
         nodrgb=0
         rmax=0.d0
         do i=1,nsphere
            xib(:)=rpos(:,i)-rbeam(:)
            rib=sqrt(dot_product(xib,xib))
            rmax=max(rmax,rib)
            call tranordertest(rib,(1.d0,0.d0),nodr(i),epstran,n)
            nodrgb=max(n,nodrgb)
         enddo
         allocate(pmnp0(0:nodrgb+1,nodrgb,2,2))
         call gaussianbeamcoef(alpha,beta,cbeam,nodrgb,pmnp0)
         pmnp=0.d0
         noff=0
         do i=1,nsphere
            nblk=2*nodr(i)*(nodr(i)+2)
            xib(:)=rpos(:,i)-rbeam(:)
            do k=1,2
               call rottran(pmnp0(0:nodrgb+1,1:nodrgb,1:2,k),pmnp(noff+1:noff+nblk,k),xib, &
                    (1.d0,0.d0),nodrgb,nodr(i),1,1,1,1)
            enddo
            noff=noff+nblk
         enddo
         deallocate(pmnp0)
         end subroutine spheregaussianbeamcoef

      end module translation
!
! scatprops module: various subroutines for calculation of observables from the solution
!
!
!  last revised: 15 January 2011
!
      module scatprops
      implicit none
      contains
!
!  determination of maximum orders for target--based expansions
!
!
!  last revised: 15 January 2011
!
         subroutine tranorders(nsphere,nodr,rpos,eps,ntran,nodrt)
         use numconstants
         use specialfuncs
         use translation
         implicit none
         integer :: nsphere,nodr(nsphere),nodrt,ntran(nsphere),i
         real(8) :: rpos(3,nsphere),r,eps
         nodrt=0
         do i=1,nsphere
            r=sqrt(dot_product(rpos(:,i),rpos(:,i)))
            call tranordertest(r,(1.d0,0.d0),nodr(i),eps,ntran(i))
            if(print_intermediate_results.eq.1) &
               write(*,'('' i, nodr, ntran:'',3i7)') i,nodr(i),ntran(i)
            nodrt=max(nodrt,ntran(i))
         enddo
         end subroutine tranorders
!
!  translation of sphere-based expansions to common target origin
!
!
!  last revised: 15 January 2011
!
         subroutine amncommonorigin(neqns,nsphere,nodr,ntran,nodrt,rpos,amnp,amnp0)
         use specialfuncs
         use translation
         implicit none
         integer :: neqns,nsphere,nodr(nsphere),nodrt,i,m,n,p,nblk,ntran(nsphere),noff
         real(8) :: rpos(3,nsphere),r,eps,xij(3)
         complex(8) :: amnp(neqns),amnp0(0:nodrt+1,nodrt,2)
         complex(8), allocatable :: amnpt(:,:,:)
         amnp0=(0.d0,0.d0)
         noff=0
         do i=1,nsphere
            allocate(amnpt(0:ntran(i)+1,ntran(i),2))
            amnpt=(0.d0,0.d0)
            nblk=nodr(i)*(nodr(i)+2)*2
            xij=-rpos(:,i)
            call rottran(amnp(noff+1:noff+nblk),amnpt,xij,(1.d0,0.d0), &
                 nodr(i),ntran(i),1,1,1,1)
            do p=1,2
               do n=1,ntran(i)
                  do m=0,ntran(i)+1
                     amnp0(m,n,p)=amnp0(m,n,p)+amnpt(m,n,p)
                  enddo
               enddo
            enddo
            deallocate(amnpt)
            noff=noff+nblk
         enddo
         end subroutine amncommonorigin
!
!  sphereqeff computes the efficiency factors for the sphere, given an1: mie coefficients,
!  anp: scattering coefficients, pnp: incident field coefficients.
!
! This subroutine is specific to the OA model for the sphere.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: polarized and cross-polarized efficiency calculation
!  30 March 2011: added optical activity
!
         subroutine sphereqeff(nsphere,neqns,nodr,nodrmax,xsp,anp1,anp2,&
                    pnp1,pnp2,qext,qabs,qsca)
         use miecoefdata
         use spheredata
         implicit none
         integer :: nsphere,m,n,p,i,nodr(nsphere),nblk,noff,neqns,nodrmax
         real(8) :: xsp(nsphere),qext(nsphere),qabs(nsphere),qsca(nsphere), &
                    qe,qa,qs
         complex(8) :: anp1(neqns),pnp1(neqns),anp2(neqns),pnp2(neqns)
         complex(8) :: anmie(2,2,nodrmax)
         qext=0.d0
         qabs=0.d0
         qsca=0.d0
         noff=0
         do i=1,nsphere
            nblk=nodr(i)*(nodr(i)+2)*2
            call getmiedata(which_sphere=i,sphere_mie_coefficients=anmie)
            call qeffcalc(nodr(i),anp1(noff+1:noff+nblk),anp2(noff+1:noff+nblk), &
                    pnp1(noff+1:noff+nblk),pnp2(noff+1:noff+nblk),anmie,qe,qa,qs)
            noff=noff+nblk
            qext(i)=2.d0*qe/xsp(i)/xsp(i)
            qabs(i)=2.d0*qa/xsp(i)/xsp(i)
            qsca(i)=2.d0*qs/xsp(i)/xsp(i)
         enddo
         end subroutine sphereqeff
!
!  calculation of sphere efficiency factors for scattered and incident field
!  coefficient anp1, pnp1, anp2, pnp2  and mie coefficients anmie
!
!  original: 15 January 2011
!  revised: 21 February 2011: polarized and cross-polarized efficiency calculation
!  30 March 2011: added optical activity
!
         subroutine qeffcalc(nodr,anp1,anp2,pnp1,pnp2,anmie,qe,qa,qs)
         implicit none
         integer :: nodr,m,n,p,q
         real(8) :: qe,qa,qs,babs,aninv(2,2)
         complex(8) :: anp1(0:nodr+1,nodr,2),pnp1(0:nodr+1,nodr,2), &
                       anp2(0:nodr+1,nodr,2),pnp2(0:nodr+1,nodr,2),anmie(2,2,nodr), &
                       a
         qe=0.d0
         qa=0.d0
         qs=0.d0
         do n=1,nodr
            a=anmie(1,1,n)*anmie(2,2,n)-anmie(1,2,n)*anmie(1,2,n)
            do p=1,2
               do q=1,2
                  aninv(p,q)=(-1)**(p+q)*anmie(3-p,3-q,n)/a
               enddo
               aninv(p,p)=aninv(p,p)+1.d0
            enddo
            do p=1,2
!               babs=-(1.d0/anmie(p,n)+1.d0)
               do m=-n,-1
                  qe=qe-(anp1(n+1,-m,p)*conjg(pnp2(n+1,-m,p)) &
                       + anp2(n+1,-m,p)*conjg(pnp1(n+1,-m,p)))*.5d0
                  qs=qs+anp1(n+1,-m,p)*conjg(anp2(n+1,-m,p))
                  do q=1,2
                     qa=qa-conjg(anp1(n+1,-m,p))*aninv(p,q)*anp2(n+1,-m,q)
                  enddo
               enddo
               do m=0,n
                  qe=qe-(anp1(m,n,p)*conjg(pnp2(m,n,p)) &
                       +anp2(m,n,p)*conjg(pnp1(m,n,p)))*.5d0
                  qs=qs+anp1(m,n,p)*conjg(anp2(m,n,p))
                  do q=1,2
                     qa=qa-conjg(anp1(m,n,p))*aninv(p,q)*anp2(m,n,q)
                  enddo
               enddo
            enddo
         enddo
         end subroutine qeffcalc
!
!  scattering amplitude sa and matrix sm calculation
!
!  original: 15 January 2011
!  revised: 21 February 2011: S11 normalization changed
!
         subroutine scatteringmatrix(amn0,nodrt,xv,ct,phi,sa,sm)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nodrt,m,n,p,m1,n1,i,j
         real(8) :: xv,ct,phi,sm(4,4),tau(0:nodrt+1,nodrt,2),cphi,sphi,qsca
         complex(8) :: amn0(0:nodrt+1,nodrt,2,2),sa(4),ephi,ephim(-nodrt:nodrt), &
                       ci,cin,a,b,sp(4,4)
         data ci/(0.d0,1.d0)/


         call taufunc(ct,nodrt,tau)
         cphi=cos(phi)
         sphi=sin(phi)
         ephi=dcmplx(cphi,sphi)
         call ephicoef(ephi,nodrt,ephim)
         sa=(0.d0,0.d0)
         qsca=0.d0
         do n=1,nodrt
            cin=(-ci)**n
            do m=-n,n
               if(m.le.-1) then
                  m1=n+1
                  n1=-m
               else
                  m1=m
                  n1=n
               endif
               do p=1,2
                  qsca=qsca+amn0(m1,n1,p,1)*dconjg(amn0(m1,n1,p,1)) &
                           + amn0(m1,n1,p,2)*dconjg(amn0(m1,n1,p,2))
                  a=amn0(m1,n1,p,1)*cphi+amn0(m1,n1,p,2)*sphi
                  b=amn0(m1,n1,p,1)*sphi-amn0(m1,n1,p,2)*cphi
                  sa(1)=sa(1)+cin*tau(m1,n1,3-p)*b*ephim(m)
                  sa(2)=sa(2)+ci*cin*tau(m1,n1,p)*a*ephim(m)
                  sa(3)=sa(3)+ci*cin*tau(m1,n1,p)*b*ephim(m)
                  sa(4)=sa(4)+cin*tau(m1,n1,3-p)*a*ephim(m)
               enddo
            enddo
         enddo
         qsca=qsca*2.d0
         do i=1,4
            do j=1,4
               sp(i,j)=sa(i)*dconjg(sa(j))*16.d0/qsca
            enddo
         enddo
         sm(1,1)=sp(1,1)+sp(2,2)+sp(3,3)+sp(4,4)
         sm(1,2)=-sp(1,1)+sp(2,2)-sp(3,3)+sp(4,4)
         sm(2,1)=-sp(1,1)+sp(2,2)+sp(3,3)-sp(4,4)
         sm(2,2)=sp(1,1)+sp(2,2)-sp(3,3)-sp(4,4)
         sm(3,3)=2.*(sp(1,2)+sp(3,4))
         sm(3,4)=-2.*dimag(sp(1,2)+sp(3,4))
         sm(4,3)=2.*dimag(sp(1,2)-sp(3,4))
         sm(4,4)=2.*(sp(1,2)-sp(3,4))
         sm(1,3)=2.*(sp(2,3)+sp(1,4))
         sm(3,1)=2.*(sp(2,4)+sp(1,3))
         sm(1,4)=2.*dimag(sp(2,3)-sp(1,4))
         sm(4,1)=-2.*dimag(sp(2,4)+sp(1,3))
         sm(2,3)=2.*(sp(2,3)-sp(1,4))
         sm(3,2)=2.*(sp(2,4)-sp(1,3))
         sm(2,4)=2.*dimag(sp(2,3)+sp(1,4))
         sm(4,2)=-2.*dimag(sp(2,4)-sp(1,3))
!         do i=1,4
!            do j=1,4
!               if(i.ne.1.or.j.ne.1) then
!                  sm(i,j)=sm(i,j)/sm(1,1)
!               endif
!            enddo
!         enddo
         end subroutine scatteringmatrix
!   c                                                                               c
!   c  subroutine scatexp(amn0,nodrt,nodrg,gmn) computes the expansion coefficients c
!   c  for the spherical harmonic expansion of the scattering phase function from   c
!   c  the scattering coefficients amn0.  For a complete expansion, the max. order  c
!   c  of the phase function expansion (nodrg) will be 2*nodrt, where nodrt is      c
!   c  the max. order of the scattered field expansion.   In this code nodrg is     c
!   c  typically set to 1, so that the subroutine returns the first moments         c
!   c  of the phase function; gmn(1) and gmn(2).                                    c
!   c                                                                               c
!   c  The expansion coefficients are normalized so that gmn(0)=1                   c
!   c                                                                               c
!   c  gmn(1)/3 is the asymmetry parameter.                                         c
!   c                                                                               c
         subroutine s11expansion(amn0,nodrt,mmax,nodrg,gmn)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nodrt,m,n,p,ma,na,mmax,nodrg,w,w1,w2,u,uw,ww1, &
                    l1,l2,ka,la,k,l,q,ik
         real(8) :: vc1(0:nodrt*2+1),vc2(0:nodrt*2+1),g0
         complex(8) :: amn0(0:nodrt+1,nodrt,2,2),gmn(0:nodrg*(nodrg+3)/2), &
                       a(2,2),c,c2
         gmn=(0.d0,0.d0)
         do n=1,nodrt
            l1=max(1,n-nodrg)
            l2=min(nodrt,n+nodrg)
            do l=l1,l2
               c=sqrt(dble((n+n+1)*(l+l+1)))*dcmplx(0.d0,1.d0)**(l-n)
               w2=min(n+l,nodrg)
               call vcfunc(-1,l,1,n,vc2)
               do m=-n,n
                  if(m.le.-1) then
                     ma=n+1
                     na=-m
                  else
                     ma=m
                     na=n
                  endif
                  do k=-l,min(l,m)
                     if(k.le.-1) then
                        ka=l+1
                        la=-k
                     else
                        ka=k
                        la=l
                     endif
                     u=m-k
                     if(u.le.mmax) then
                        ik=(-1)**k
                        c2=ik*c
                        do p=1,2
                           do q=1,2
                              a(p,q)=c2*(amn0(ma,na,p,1)*conjg(amn0(ka,la,q,1)) &
                                    +amn0(ma,na,p,2)*conjg(amn0(ka,la,q,2)))
                           enddo
                        enddo
                        w1=max(abs(n-l),abs(u))
                        w2=min(n+l,nodrg)
                        call vcfunc(-k,l,m,n,vc1)
                        do w=w1,w2
                           uw=(w*(w+1))/2+u
                           do p=1,2
                              if(mod(n+l+w,2).eq.0) then
                                 q=p
                              else
                                 q=3-p
                              endif
                              gmn(uw)=gmn(uw)-vc1(w)*vc2(w)*a(p,q)
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
         g0=dble(gmn(0))
         gmn(0)=1.d0
         do w=1,nodrg
            ww1=(w*(w+1))/2
            gmn(ww1)=dcmplx(dble(gmn(ww1)),0.d0)/g0
            do u=1,min(mmax,w)
               uw=ww1+u
               gmn(uw)=(-1)**u*2.d0*gmn(uw)/g0
            enddo
         enddo
         end subroutine s11expansion
!
!  calculate azimuth--averaged scattering matrix from expansion, for cos(theta) = ct
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
         subroutine fosmcalc(ntot,s00,s02,sp22,sm22,ct,sm)
         use numconstants
         use specialfuncs
         integer :: ntot,w,i,j,ww1
         real(8) :: s00(4,4,0:ntot*2),s02(4,4,0:ntot*2),sp22(4,4,0:ntot*2),sm22(4,4,0:ntot*2), &
                    sm(4,4),dc(-2:2,0:2*ntot*(2*ntot+2)),ct
         call rotcoef(ct,2,2*ntot,dc)
         sm=0.d0
         do w=0,2*ntot
            ww1=w*(w+1)
            sm(:,:)=sm(:,:)+s00(:,:,w)*dc(0,ww1)+s02(:,:,w)*dc(0,ww1+2) &
                   +sp22(:,:,w)*dc(2,ww1+2)+sm22(:,:,w)*dc(-2,ww1+2)
         enddo
         sm=sm/s00(1,1,0)
!         do i=1,4
!            do j=1,4
!               if(i.ne.1.or.j.ne.1) then
!                  sm(i,j)=sm(i,j)/sm(1,1)
!               endif
!            enddo
!         enddo
         end subroutine fosmcalc
!
!  determine the generalized spherical function expansion for the azimuth-averaged scattering matrix
!  corresponding to the target-based scattering field expansion of amnp.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fixed flush call.
!
         subroutine fosmexpansion(ntot,amnp,s00,s02,sp22,sm22)
         use mpidefs
         use mpidata
         use specialfuncs
         use numconstants
         use spheredata
         integer :: ntot,n,p,m,l,wmin,wmax,m1m,q,m1mq,m1mnpl,w,m1w,fe,fo,i,j,wtot
         integer :: rank,numprocs,nl,nsend,runprintunit
         integer, allocatable :: nlindex(:),nlnum(:)
         real(8) :: s00(4,4,0:ntot*2),s02(4,4,0:ntot*2),sp22(4,4,0:ntot*2),sm22(4,4,0:ntot*2), &
                    cm1p1(0:ntot*2),cm1m1(0:ntot*2),cmmpm(0:ntot*2),cmmm2pm(0:ntot*2), &
                    cmmp2pm(0:ntot*2),sum,nlperproc
         complex(8) :: amnp(0:ntot+1,ntot,2,2),a1(-ntot-2:ntot+2,ntot,2),a2(-ntot-2:ntot+2,ntot,2), &
                       ci,fnl,a1122,a2112,a1p2,a1m2
         data ci/(0.d0,1.d0)/
         call init(2*ntot)
         call getrunparameters(run_print_unit=runprintunit)
         call ms_mpi(mpi_command='rank',mpi_rank=rank)
         call ms_mpi(mpi_command='size',mpi_size=numprocs)
         allocate(nlindex(0:numprocs-1),nlnum(0:numprocs-1))
         nlperproc=dble(ntot*ntot)/dble(numprocs)
         sum=0.
         do i=0,numprocs-1
            nlindex(i)=floor(sum)
            sum=sum+nlperproc
         enddo
         do i=0,numprocs-2
            nlnum(i)=nlindex(i+1)-nlindex(i)
         enddo
         nlnum(numprocs-1)=ntot*ntot-nlindex(numprocs-1)
         if(rank.eq.0) then
            write(runprintunit,'('' SM calc, orders per processor:'',f10.4)') nlperproc
            call flush(runprintunit)
         endif
         a1=(0.d0,0.d0)
         a2=(0.d0,0.d0)
         s00=0.d0
         s02=0.d0
         sp22=0.d0
         sm22=0.d0
         wtot=ntot+ntot
         do n=1,ntot
            do p=1,2
               do m=-n,-1
                  a1(m,n,p)=amnp(n+1,-m,p,1)
                  a2(m,n,p)=amnp(n+1,-m,p,2)
               enddo
               do m=0,n
                  a1(m,n,p)=amnp(m,n,p,1)
                  a2(m,n,p)=amnp(m,n,p,2)
               enddo
            enddo
         enddo
         do nl=nlindex(rank)+1,nlindex(rank)+nlnum(rank)
            n=floor((nl-1)/dble(ntot))+1
            l=mod(nl-1,ntot)+1
            wmin=abs(n-l)
            wmax=n+l
            fnl=sqrt(dble((n+n+1)*(l+l+1)))*ci**(l-n)
            call vcfunc(-1,n,1,l,cm1p1)
            call vcfunc(-1,n,-1,l,cm1m1)
            do m=-min(n,l+2),min(n,l+2)
               m1m=(-1)**m
               if(abs(m).le.l) then
                  call vcfunc(-m,n,m,l,cmmpm)
               else
                  cmmpm=0.d0
               endif
               if(abs(-2+m).le.l) then
                  call vcfunc(-m,n,-2+m,l,cmmm2pm)
               else
                  cmmm2pm=0.d0
               endif
               if(abs(2+m).le.l) then
                  call vcfunc(-m,n,2+m,l,cmmp2pm)
               else
                  cmmp2pm=0.d0
               endif
               do p=1,2
                  do q=1,2
                     m1mq=(-1)**(m+q)
                     m1mnpl=(-1)**(m+n+p+l)
                     a1122=(a1(m,n,p)*conjg(a1(m,l,q)) + a2(m,n,p)*conjg(a2(m,l,q)))
                     a2112=(a2(m,n,p)*conjg(a1(m,l,q)) - a1(m,n,p)*conjg(a2(m,l,q)))
                     a1p2=(a1(m,n,p)+ci*a2(m,n,p))*conjg(a1(m-2,l,q)-ci*a2(m-2,l,q))
                     a1m2=(a1(m,n,p)-ci*a2(m,n,p))*conjg(a1(m+2,l,q)+ci*a2(m+2,l,q))
                     do w=wmin,wmax
                        m1w=(-1)**w
                        if(mod(n+l+w+p+q,2).eq.0) then
                           fe=1
                           fo=0
                        else
                           fe=0
                           fo=1
                        endif
                        s00(1,1,w) = s00(1,1,w)-(m1m*fe*fnl*a1122*cm1p1(w)*cmmpm(w))/2.
                        s00(3,2,w) = s00(3,2,w)+ (ci/2.*m1m*fnl*fo*a1122*cm1p1(w)*cmmpm(w))
                        s00(4,2,w) = s00(4,2,w)+ dimag(-ci/2.*m1m*fnl*fo*a1122*cm1p1(w)*cmmpm(w))
                        s00(1,4,w) = s00(1,4,w)+ dimag(m1m*fe*fnl*(-a2112)*cm1p1(w)*cmmpm(w))/2.
                        s00(2,3,w) = s00(2,3,w)+ (m1m*fe*fnl*(-a2112)*cm1p1(w)*cmmpm(w))/2.
                        s00(4,3,w) = s00(4,3,w)+ dimag(ci/2.*m1m*fnl*fo*a2112*cm1p1(w)*cmmpm(w))
                        s00(4,4,w) = s00(4,4,w)+ (ci/2.*m1m*fnl*fo*a2112*cm1p1(w)*cmmpm(w))

                        if(w.lt.2) cycle

                        s02(2,1,w) = s02(2,1,w)-(m1mq*a1122*fe*fnl*cm1m1(w)*cmmpm(w))/2.
                        s02(3,1,w) = s02(3,1,w)+ (-ci/2.*m1mq*a1122*fnl*fo*cm1m1(w)*cmmpm(w))
                        s02(4,1,w) = s02(4,1,w)+ dimag(ci/2.*m1mq*a1122*fnl*fo*cm1m1(w)*cmmpm(w))
                        s02(1,3,w) = s02(1,3,w)-(m1mq*a2112*fe*fnl*cm1m1(w)*cmmpm(w))/2.
                        s02(2,4,w) = s02(2,4,w)-dimag(m1mq*a2112*fe*fnl*cm1m1(w)*cmmpm(w))/2.
                        s02(3,3,w) = s02(3,3,w)+ (ci/2.*m1mq*a2112*fnl*fo*cm1m1(w)*cmmpm(w))
                        s02(3,4,w) = s02(3,4,w)+ dimag(-ci/2.*m1mq*a2112*fnl*fo*cm1m1(w)*cmmpm(w))

                        s02(1,2,w) = s02(1,2,w)-(m1m*a1p2*fe*fnl*cm1p1(w)*cmmm2pm(w))/4.
                        s02(1,3,w) = s02(1,3,w)+ (-ci/4.*m1m*a1p2*fe*fnl*cm1p1(w)*cmmm2pm(w))
                        s02(2,4,w) = s02(2,4,w)+ dimag(-ci/4.*m1m*a1p2*fe*fnl*cm1p1(w)*cmmm2pm(w))
                        s02(3,1,w) = s02(3,1,w)+ (ci/4.*m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))
                        s02(4,1,w) = s02(4,1,w)+ dimag(-ci/4.*m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))
                        s02(4,3,w) = s02(4,3,w)+ dimag(m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))/4.
                        s02(4,4,w) = s02(4,4,w)+ (m1m*a1p2*fnl*fo*cm1p1(w)*cmmm2pm(w))/4.

                        sm22(1,4,w) = sm22(1,4,w)+ dimag(-ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sm22(2,2,w) = sm22(2,2,w)-(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sm22(2,3,w) = sm22(2,3,w)+ (-ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sm22(3,2,w) = sm22(3,2,w)+ (ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sm22(3,3,w) = sm22(3,3,w)-(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sm22(3,4,w) = sm22(3,4,w)+ dimag(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sm22(4,2,w) = sm22(4,2,w)+ dimag(-ci/8.*m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))

                        sp22(1,4,w) = sp22(1,4,w)+ dimag(-ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sp22(2,2,w) = sp22(2,2,w)-(m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sp22(2,3,w) = sp22(2,3,w)+ (-ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sp22(3,2,w) = sp22(3,2,w)+ (-ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))
                        sp22(3,3,w) = sp22(3,3,w)+ (m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sp22(3,4,w) = sp22(3,4,w)-dimag(m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                        sp22(4,2,w) = sp22(4,2,w)+ dimag(ci/8.*m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))

                        s02(1,2,w) = s02(1,2,w)-(m1m*a1m2*fe*fnl*cm1p1(w)*cmmp2pm(w))/4.
                        s02(1,3,w) = s02(1,3,w)+ (ci/4.*m1m*a1m2*fe*fnl*cm1p1(w)*cmmp2pm(w))
                        s02(2,4,w) = s02(2,4,w)+ dimag(ci/4.*m1m*a1m2*fe*fnl*cm1p1(w)*cmmp2pm(w))
                        s02(3,1,w) = s02(3,1,w)+ (ci/4.*m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))
                        s02(4,1,w) = s02(4,1,w)+ dimag(-ci/4.*m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))
                        s02(4,3,w) = s02(4,3,w)-dimag(m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))/4.
                        s02(4,4,w) = s02(4,4,w)-(m1m*a1m2*fnl*fo*cm1p1(w)*cmmp2pm(w))/4.

                        sm22(1,4,w) = sm22(1,4,w)+ dimag(ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sm22(2,2,w) = sm22(2,2,w)-(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sm22(2,3,w) = sm22(2,3,w)+ (ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sm22(3,2,w) = sm22(3,2,w)+ (-ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sm22(3,3,w) = sm22(3,3,w)-(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sm22(3,4,w) = sm22(3,4,w)+ dimag(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sm22(4,2,w) = sm22(4,2,w)+ dimag(ci/8.*m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))

                        sp22(1,4,w) = sp22(1,4,w)+ dimag(ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sp22(2,2,w) = sp22(2,2,w)-(m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sp22(2,3,w) = sp22(2,3,w)+ (ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sp22(3,2,w) = sp22(3,2,w)+ (ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                        sp22(3,3,w) = sp22(3,3,w)+ (m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sp22(3,4,w) = sp22(3,4,w)-dimag(m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        sp22(4,2,w) = sp22(4,2,w)+ dimag(-ci/8.*m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         call ms_mpi(mpi_command='barrier')
         nsend=4*4*(2*ntot+1)
         call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dp=s00,&
              mpi_number=nsend,mpi_operation=ms_mpi_sum)
         call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dp=s02,&
              mpi_number=nsend,mpi_operation=ms_mpi_sum)
         call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dp=sp22,&
              mpi_number=nsend,mpi_operation=ms_mpi_sum)
         call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dp=sm22,&
              mpi_number=nsend,mpi_operation=ms_mpi_sum)
!
!  a patch
!
         do i=3,4
            do j=1,i
               s00(j,i,0:wtot)=-s00(j,i,0:wtot)
               s02(j,i,0:wtot)=-s02(j,i,0:wtot)
               sm22(j,i,0:wtot)=-sm22(j,i,0:wtot)
               sp22(j,i,0:wtot)=-sp22(j,i,0:wtot)
            enddo
         enddo
         deallocate(nlindex,nlnum)
         end subroutine fosmexpansion
!
!  compute the coefficients for the GSF expansion of the random orientation
!  scattering matrix.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
         subroutine ranorientscatmatrix(xv,nsphere,nodr,nodrw,cbeam,tmatrixfile,&
                    sm,qext,qabs,qsca)
         use mpidefs
         use mpidata
         use intrinsics
         use specialfuncs
         use spheredata
         use numconstants
         implicit none
         integer :: nodr,nodrw,nodr2,m,n,p,k,l,q,s,t,v,u,w,nblk,kl,mn,nn1,tn, &
                    lmax,ll1,tvl,ku,k1,ns,ik,ik1,m1,nu,n1s,n1e,nu1,p1,n1max, &
                    in,n1,i,lt,kt,qt,nt,mt,ikm,klm,mnm,nodrt,nsphere, &
                    rank,iunit,numprocs
         real(8) :: sm(4,4,0:nodrw),fl,vc(0:4*nodr+2),xv,fl2,fc1,fc2,fc3,fc4, &
                    cbeam,gbn,qext(nsphere),qabs(nsphere),qsca(nsphere),qel, &
                    qal,qsl,fc(4),time1,time2,qsca0
         complex(8) :: ci,cin,a
         complex(8) :: aw(0:2,-1:1,0:nodrw),bw(0:2,-1:1,0:nodrw),cw(0:nodrw), &
                       dw(0:nodrw),pp(nodr,2,2), &
                       bm(2,nodr*(nodr+2),2),am(2,nodr+1,2),fm(3,nodr,2,nodr,2)
         complex(8), allocatable :: dm(:,:,:,:,:,:)
         complex(4), allocatable :: tc(:,:,:,:)
         integer :: nblkw,wv,sizedm,ierr,sizetm,nsend
         integer, allocatable :: windex(:),vindex(:),wvindex(:),wvnum(:)
         real(8) :: wvperproc,sum
         character*30 :: tmatrixfile
         data ci/(0.d0,1.d0)/
         call ms_mpi(mpi_command='rank',mpi_rank=rank)
         call ms_mpi(mpi_command='size',mpi_size=numprocs)
         call getrunparameters(run_print_unit=iunit)
         if(rank.eq.0) time1=mytime()
!
!  read the T matrix from the file
!
         if(rank.eq.0) then
            open(3,file=tmatrixfile)
            read(3,*) nodrt
         endif
         nodrt=nodr
         nblk=nodr*(nodr+2)
         sizetm=4*nblk*nblk
         allocate(tc(2,nblk,2,nblk))
         tc=(0.,0.)
         if(rank.eq.0) then
            qext=0.d0
            qabs=0.d0
            qsca=0.d0
            do l=1,nodr
               gbn=dexp(-((dble(l)+.5d0)*cbeam)**2.)
               do k=-l,l
                  kl=l*(l+1)+k
                  klm=l*(l+1)-k
                  do q=1,2
                     read(3,*) lt,kt,qt
                     do n=1,l
                        do m=-n,n
                           mn=n*(n+1)+m
                           mnm=n*(n+1)-m
                           read(3,*) nt,mt,fc
                           tc(1,mn,q,kl)=cmplx(fc(1),fc(2))
                           tc(2,mn,q,kl)=cmplx(fc(3),fc(4))
                           if(n.lt.l) then
                              ikm=(-1)**(m+k)
                              do p=1,2
                                 tc(q,klm,p,mnm)=tc(p,mn,q,kl)*ikm
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               enddo
               do i=1,nsphere
                  read(3,*) n,qel,qal,qsl
                  qext(i)=qext(i)+qel*gbn*gbn
                  qabs(i)=qabs(i)+qal*gbn*gbn
                  qsca(i)=qsca(i)+qsl*gbn*gbn
               enddo
            enddo
            close(3)
         endif
!
!  send to the other processors
!
         if(numprocs.gt.1) then
            call ms_mpi(mpi_command='barrier')
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qext,mpi_number=nsphere,mpi_rank=0)
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qabs,mpi_number=nsphere,mpi_rank=0)
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qsca,mpi_number=nsphere,mpi_rank=0)
            call ms_mpi(mpi_command='bcast',mpi_send_buf_c=tc,mpi_number=sizetm,mpi_rank=0)
         endif
         allocate(dm(-nodr-1:nodr+1,3,nodr,2,nodr,2))
         if(rank.eq.0) then
            time2=mytime()-time1
            call timewrite(iunit,' t matrix read time:',time2)
            time1=mytime()
         endif
         nodr2=nodr+nodr
         nblk=nodr*(nodr+2)
         dm=(0.d0,0.d0)
         sizedm=size(dm)
         call init(nodr2)
!
!  compute the GB modified T matrix
!
         do n=1,nodr
            gbn=dexp(-((dble(n)+.5d0)*cbeam)**2.)
            cin=ci**(n+1)
            pp(n,1,1) =-.5d0*cin*fnr(n+n+1)*gbn
            pp(n,2,1) =-pp(n,1,1)
            pp(n,1,2)=-pp(n,1,1)
            pp(n,2,2)=pp(n,2,1)
         enddo
         do n=1,nodr
            nn1=n*(n+1)
            do m=-n,n
               mn=nn1+m
               do p=1,2
                  do l=1,nodr
                     do k=-l,l
                        kl=l*(l+1)+k
                        a=tc(p,mn,1,kl)
                        tc(p,mn,1,kl)=tc(p,mn,1,kl)*pp(l,1,1)&
                           +tc(p,mn,2,kl)*pp(l,1,2)
                        tc(p,mn,2,kl)=a*pp(l,2,1)+tc(p,mn,2,kl)*pp(l,2,2)
                     enddo
                  enddo
               enddo
            enddo
         enddo
!
!  determine the distribution of work load among the processors
!
         nblkw=nodr2*(nodr2+2)+1
         allocate(windex(nblkw),vindex(nblkw),wvindex(0:numprocs-1),wvnum(0:numprocs-1))
         w=0
         do n=0,nodr2
            do m=-n,n
               w=w+1
               windex(w)=n
               vindex(w)=m
            enddo
         enddo
         wvperproc=dble(nblkw)/dble(numprocs)
         sum=0.
         do i=0,numprocs-1
            wvindex(i)=floor(sum)
            sum=sum+wvperproc
         enddo
         do i=0,numprocs-2
            wvnum(i)=wvindex(i+1)-wvindex(i)
         enddo
         wvnum(numprocs-1)=nblkw-wvindex(numprocs-1)
         if(rank.eq.0) then
            write(iunit,'('' d matrix calculation, order+degree per proc.:'',f9.2)') &
                wvperproc
            call flush(iunit)
         endif
!
!  the big loop
!
         do wv=wvindex(rank)+1,wvindex(rank)+wvnum(rank)
            w=windex(wv)
            v=vindex(wv)
            bm=(0.d0,0.d0)
            do n=1,nodr
               nn1=n*(n+1)
               do l=max(1,abs(w-n)),min(nodr,w+n)
                  am(1,l,1)=0.d0
                  am(1,l,2)=0.d0
                  am(2,l,1)=0.d0
                  am(2,l,2)=0.d0
               enddo
               do t=-n,n
                  tn=nn1+t
                  lmax=min(nodr,w+n)
                  call vcfunc(v,w,-t,n,vc)
                  do l=max(1,abs(v-t),abs(n-w)),lmax
                     ll1=l*(l+1)
                     tvl=ll1+t-v
                     do k=1,2
                        do p=1,2
                           am(k,l,p)=am(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                        enddo
                     enddo
                  enddo
               enddo
               do m=-n,n
                  mn=nn1+m
                  do k=1,2
                     u=m-(-3+2*k)
                     if(abs(u).le.w) then
                        lmax=min(nodr,w+n)
                        call vcfunc(-u,w,m,n,vc)
                        do l=max(1,abs(w-n)),lmax
                           fl=-(-1)**l*vc(l)/dble(l+l+1)
                           do p=1,2
                              bm(k,mn,p)=bm(k,mn,p)+am(k,l,p)*fl
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
            do u=-min(w,nodr+1),min(w,nodr+1)
               do ku=1,3
                  if(ku.eq.1) then
                     k=-1
                     k1=-1
                  elseif(ku.eq.2) then
                     k=1
                     k1=1
                  else
                     k=1
                     k1=-1
                  endif
                  m=u+k
                  ns=max(1,abs(m))
                  ik=(k+1)/2+1
                  ik1=(k1+1)/2+1
                  m1=u+k1
                  do n=ns,nodr
                     nu=n*(n+1)+m
                     n1s=max(1,abs(m1),n-nodrw)
                     n1e=min(nodr,n+nodrw)
                     do n1=n1s,n1e
                        cin=ci**(n-n1)
                        nu1=n1*(n1+1)+m1
                        fl=-fnr(n+n+1)*fnr(n1+n1+1)*dble(w+w+1)
                        do p=1,2
                           do p1=1,2
                              a=bm(ik,nu,p)*cin*fl*conjg(bm(ik1,nu1,p1))
                              dm(u,ku,n,p,n1,p1)=dm(u,ku,n,p,n1,p1)+a
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         deallocate(tc)
         call ms_mpi(mpi_command='barrier')
         call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dc=dm,&
              mpi_number=sizedm,mpi_operation=ms_mpi_sum)
         if(rank.eq.0) then
            time2=mytime()-time1
            call timewrite(iunit,' d matrix time:',time2)
            time1=mytime()
         endif
!
!  compute the expansion coefficients
!
         aw=0.d0
         bw=0.d0
         cw=0.d0
         dw=0.d0
         do w=0,nodrw
            do n=1,nodr
               n1s=max(1,abs(n-w))
               n1e=min(nodr,n+w)
               do n1=n1s,n1e
                  do k=1,3
                     do p=1,2
                        do p1=1,2
                           fm(k,n,p,n1,p1)=0.
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do u=-nodr-1,nodr+1
               do k=-1,1,2
                  m=u+k
                  ik=(k+1)/2+1
                  ns=max(1,abs(m))
                  do n=ns,nodr
                     n1max=min(w+n,nodr)
                     call vcfunc(m,n,0,w,vc)
                     do n1=ns,nodr
                        if((n+n1.lt.w).or.(abs(n-n1).gt.w)) cycle
                        fl=-(-1)**n*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                        do p=1,2
                           do p1=1,2
                              fm(ik,n,p,n1,p1)=fm(ik,n,p,n1,p1)+dm(u,ik,n,p,n1,p1)*fl
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               if(w.lt.2) cycle
               m=u+1
               m1=u-1
               ns=max(1,abs(m))
               n1s=max(1,abs(m1))
               do n=ns,nodr
                  n1max=min(w+n,nodr)
                  call vcfunc(m,n,-2,w,vc)
                  do n1=n1s,nodr
                     if((n+n1.lt.w).or.(abs(n-n1).gt.w)) cycle
                     fl=-(-1)**n*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     do p=1,2
                        do p1=1,2
                           fm(3,n,p,n1,p1)=fm(3,n,p,n1,p1)+dm(u,3,n,p,n1,p1)*fl
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do n=1,nodr
               n1s=max(1,abs(n-w))
               n1e=min(nodr,n+w)
               in=(-1)**n
               n1max=min(w+n,nodr)
               call vcfunc(1,n,0,w,vc)
               do n1=n1s,n1e
                  fl=2.d0*in*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                  i=mod(n+n1-w,2)+1
                  do p=1,2
                     p1=(2-i)*p+(i-1)*(3-p)
                     do k=-1,1,2
                        ik=(k+1)/2+1
                        aw(0,k,w)=aw(0,k,w)+fm(ik,n,p,n1,p1)*fl
                        bw(0,k,w)=bw(0,k,w)+fm(ik,n,p,n1,3-p1)*fl
                     enddo
                     bw(2,0,w)=bw(2,0,w)+fm(3,n,p,n1,3-p1)*fl
                     aw(2,0,w)=aw(2,0,w)+fm(3,n,p,n1,p1)*fl
                  enddo
               enddo
               if(w.lt.2) cycle
               call vcfunc(1,n,-2,w,vc)
               do n1=n1s,n1e
                  fl=2.d0*in*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                  i=mod(n+n1-w,2)+1
                  do p=1,2
                     p1=(2-i)*p+(i-1)*(3-p)
                     do k=-1,1,2
                        ik=(k+1)/2+1
                        aw(2,k,w)=aw(2,k,w)+fm(ik,n,p,n1,p1)*fl*(-1)**p1
                        bw(2,k,w)=bw(2,k,w)+fm(ik,n,p,n1,3-p1)*fl*(-1)**(3-p1)
                     enddo
                  enddo
                  fl2=2.*(-1)**(n1+w)*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                  do p=1,2
                     do p1=1,2
                        cw(w)=cw(w)+fm(3,n,p,n1,p1)*fl*(-1)**p1
                        dw(w)=dw(w)+fm(3,n,p,n1,p1)*fl2*(-1)**p
                     enddo
                  enddo
               enddo
            enddo
         enddo
         do w=0,nodrw
            do k=-1,1
               do i=0,2
                  aw(i,k,w)=aw(i,k,w)*2./xv/xv
                  bw(i,k,w)=bw(i,k,w)*2./xv/xv
               enddo
            enddo
            cw(w)=cw(w)*2./xv/xv
            dw(w)=dw(w)*2./xv/xv
         enddo
         do n=0,nodrw
            sm(1,1,n)=aw(0,-1,n)+aw(0,1,n)
            sm(1,2,n)=aw(2,-1,n)+aw(2,1,n)
            sm(1,3,n)=2.d0*dimag(aw(2,0,n))
            sm(1,4,n)=aw(0,1,n)-aw(0,-1,n)
            sm(2,2,n)=dw(n)
            sm(2,3,n)=dimag(dw(n))
            sm(2,4,n)=aw(2,1,n)-aw(2,-1,n)
            sm(3,2,n)=dimag(cw(n))
            sm(3,3,n)=cw(n)
            sm(3,4,n)=dimag(bw(2,-1,n)-bw(2,1,n))
            sm(4,4,n)=bw(0,1,n)-bw(0,-1,n)
         enddo
!
!  normalization
!
         qsca0=sm(1,1,0)
         do n=0,nodrw
            sm(1,1,n)=sm(1,1,n)/qsca0
            sm(1,2,n)=sm(1,2,n)/qsca0
            sm(1,3,n)=sm(1,3,n)/qsca0
            sm(1,4,n)=sm(1,4,n)/qsca0
            sm(2,2,n)=sm(2,2,n)/qsca0
            sm(2,3,n)=sm(2,3,n)/qsca0
            sm(2,4,n)=sm(2,4,n)/qsca0
            sm(3,2,n)=sm(3,2,n)/qsca0
            sm(3,3,n)=sm(3,3,n)/qsca0
            sm(3,4,n)=sm(3,4,n)/qsca0
            sm(4,4,n)=sm(4,4,n)/qsca0
         enddo
         call ms_mpi(mpi_command='barrier')
         if(rank.eq.0) then
            time2=mytime()-time1
            call timewrite(iunit,' scat matrix coef time:',time2)
         endif
         deallocate(windex,vindex,wvindex,wvnum,dm)
         end subroutine ranorientscatmatrix
!
!  calculation of the RO scattering matrix from the GSF expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
         subroutine ranorienscatmatrixcalc(nt,tmin,tmax,iscale,smc,nodrexp,sm)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nt,iscale,nodrexp,i,j,k,n,nn0,nnp2,nnm2
         real(8) :: tmin,tmax,smc(4,4,0:nodrexp),sm(4,4,nt),dc(-2:2,0:nodrexp*(nodrexp+2)), &
                    ct,th,qsca
         do i=1,nt
            if(nt.eq.1) then
               th=tmin
            else
               th=tmin+(tmax-tmin)*dble(i-1)/dble(nt-1)
            endif
            ct=cos(th*pi/180.d0)
!
!     dc is the normalized generalized spherical function
!     dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
!     where D^k_{mn} is defined in M&M JOSA 96
!
            call rotcoef(ct,2,nodrexp,dc)
            do j=1,4
               do k=j,4
                  sm(j,k,i)=0.d0
               enddo
            enddo
            do n=0,nodrexp
               nn0=n*(n+1)
               nnp2=nn0+2
               nnm2=nn0-2
               sm(1,1,i)=sm(1,1,i)+dc(0,nn0)*smc(1,1,n)
               sm(1,4,i)=sm(1,4,i)+dc(0,nn0)*smc(1,4,n)
               sm(4,4,i)=sm(4,4,i)+dc(0,nn0)*smc(4,4,n)
               if(n.ge.2) then
                  sm(1,2,i)=sm(1,2,i)+dc(2,nn0)*smc(1,2,n)
                  sm(2,4,i)=sm(2,4,i)+dc(2,nn0)*smc(2,4,n)
                  sm(3,4,i)=sm(3,4,i)+dc(2,nn0)*smc(3,4,n)
                  sm(1,3,i)=sm(1,3,i)+dc(2,nn0)*smc(1,3,n)
                  sm(2,2,i)=sm(2,2,i)+dc(2,nnm2)*smc(2,2,n)+dc(2,nnp2)*smc(3,3,n)
                  sm(2,3,i)=sm(2,3,i)+dc(2,nnp2)*smc(2,3,n)+dc(2,nnp2)*smc(3,2,n)
                  sm(3,3,i)=sm(3,3,i)-dc(2,nnm2)*smc(2,2,n)+dc(2,nnp2)*smc(3,3,n)
               endif
            enddo
!
!  discontiued scaling option: now done in main program
!
!            if(iscale.eq.1) then
!               do j=1,4
!                  do k=j,4
!                     if(j.ne.1.or.k.ne.1) then
!                        sm(j,k,i)=sm(j,k,i)/sm(1,1,i)
!                     endif
!                  enddo
!               enddo
!            endif
!
!    here are the VV and HH differential cross sections
!
!            gvv=.25*(sm(1,1)+sm(2,2)-2.*sm(1,2))
!            ghh=.25*(sm(1,1)+sm(2,2)+2.*sm(1,2))
!
         enddo
         return
         end subroutine ranorienscatmatrixcalc
      end module scatprops
!
! module nearfield contains local data and subroutines for near field calculation
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
      module nearfield
      implicit none
      integer, private :: axialinc,ndimpw,nodrpwmax
      integer, private, allocatable :: nblk_nf(:),noff_nf(:)
      real(8), private :: rplotmax
      complex(8), allocatable, private :: amnp_nf(:),cmnp_nf(:),pmnp_nf(:), &
               amn3mp_nf(:),cmn3mp_nf(:),pmn3mp_nf(:)
      contains
!
!  nearfieldspherepart calculates the field at point xg generated by the spheres
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine nearfieldspherepart(xg,nsphere,xsp,rpos,ri,&
                    nodr,neqns,insphere,efield,hfield)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nsphere,nodr(nsphere),neqns,i,insphere,nblki,n
         real(8) :: xg(3),xsp(nsphere),rpos(3,nsphere),x(3),r
         complex(8) :: ri(2,nsphere),vwh(3,neqns),efield(3),hfield(3),ri0,cn1,cn2
         complex(8), allocatable :: vwhleft(:,:,:),vwhright(:,:,:)

!
!  find if the point is inside a sphere
!
         insphere=0
         do i=1,nsphere
            x=xg(:)-rpos(:,i)
            r=sqrt(dot_product(x,x))
            if(r.le.xsp(i)) then
               insphere=i
               exit
            endif
         enddo
!
!  do the calculations
!
         if(insphere.eq.0) then
!
!  outside a sphere: field = scattered
!
            do i=1,nsphere
               x=xg(:)-rpos(:,i)
               ri0=(1.d0,0.d0)
               call vwhcalc(x,ri0,nodr(i),3,vwh(1:3,noff_nf(i)+1:noff_nf(i)+nblk_nf(i)))
            enddo
            efield(:)=matmul(vwh(:,1:neqns),amnp_nf(1:neqns))
            hfield(:)=matmul(vwh(:,1:neqns),amn3mp_nf(1:neqns))/dcmplx(0.d0,1.d0)
         else
!
!  inside a sphere: field = internal
!
            i=insphere
            if(abs(ri(1,i)-ri(2,i)).eq.0) then
               x=xg(:)-rpos(:,i)
               call vwhcalc(x,ri(1,i),nodr(i),1,vwh)
               efield(:)=matmul(vwh(:,1:nblk_nf(i)), &
                   cmnp_nf(noff_nf(i)+1:noff_nf(i)+nblk_nf(i)))
               hfield(:)=matmul(vwh(:,1:nblk_nf(i)), &
                   cmn3mp_nf(noff_nf(i)+1:noff_nf(i)+nblk_nf(i)))*ri(1,i)/dcmplx(0.d0,1.d0)
            else
               nblki=nodr(i)*(nodr(i)+2)
               allocate(vwhleft(3,2,nblki),vwhright(3,2,nblki))
               x=xg(:)-rpos(:,i)
               call vwhcalc(x,ri(1,i),nodr(i),1,vwhleft)
               call vwhcalc(x,ri(2,i),nodr(i),1,vwhright)
               efield=0.d0
               hfield=0.d0
               do n=1,nblki
                  cn1=cmnp_nf(noff_nf(i)+2*n-1)
                  cn2=cmnp_nf(noff_nf(i)+2*n)
                  efield(:)=efield(:)+(vwhleft(:,1,n)+vwhleft(:,2,n))*cn1 &
                         +(vwhright(:,1,n)-vwhright(:,2,n))*cn2
                  hfield(:)=hfield(:)+((vwhleft(:,2,n)+vwhleft(:,1,n))*cn1*ri(1,i) &
                         +(vwhright(:,2,n)-vwhright(:,1,n))*cn2*ri(2,i))/dcmplx(0.d0,1.d0)
               enddo
               deallocate(vwhleft,vwhright)
            endif
         endif
         end subroutine nearfieldspherepart
!
!  nearfieldincidentpart calculates the incident field at point xg using a regular
!  vswh expansion
!
!
!  last revised: 15 January 2011
!
         subroutine nearfieldincidentpart(xg,nodrpw,efield,hfield)
         use specialfuncs
         use numconstants
         implicit none
         integer :: nblkpw,nodrpw
         real(8) :: xg(3),r,epspw
         complex(8) :: vwhpw(3,nodrpw*(nodrpw+2)*2),vwhpwaxial(3,4*nodrpw), &
                       efield(3),hfield(3)
!
!  oblique incidence: use the full expansion
!
         if(axialinc.eq.0) then
            call vwhcalc(xg,(1.d0,0.d0),nodrpw,1,vwhpw)
            nblkpw=nodrpw*(nodrpw+2)*2
            efield(:)=matmul(vwhpw(:,1:nblkpw),pmnp_nf(1:nblkpw))
            hfield(:)=matmul(vwhpw(:,1:nblkpw),pmn3mp_nf(1:nblkpw))/dcmplx(0.d0,1.d0)
         else
!
!  axial incidence: use the shortcut
!
            call vwhaxialcalc(xg,(1.d0,0.d0),nodrpw,1,vwhpwaxial)
            nblkpw=4*nodrpw
            efield(:)=matmul(vwhpwaxial(:,1:nblkpw),pmnp_nf(1:nblkpw))
            hfield(:)=matmul(vwhpwaxial(:,1:nblkpw),pmn3mp_nf(1:nblkpw))/dcmplx(0.d0,1.d0)
         endif
         end subroutine nearfieldincidentpart
!
!  nearfieldincidentcoef generates the reshaped array of incident field coefficients
!
!
!  last revised: 15 January 2011
!
         subroutine nearfieldincidentcoef(nodrpw,alpha,beta,gamma,cbeam,epspw)
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: m,n,p,nn1,mn,mnp,nodrpw
         real (8) :: alpha,beta,cbeam,gamma,epspw,cgamma,sgamma
         complex(8), allocatable :: pmnp0(:,:,:,:)
         allocate(pmnp0(0:nodrpw+1,nodrpw,2,2))
         if(allocated(pmnp_nf)) deallocate(pmnp_nf,pmn3mp_nf)
         if(beta.ne.0.d0) then
            axialinc=0
            ndimpw=2*nodrpw*(nodrpw+2)
         else
            axialinc=1
            ndimpw=4*nodrpw
         endif
         allocate(pmnp_nf(ndimpw),pmn3mp_nf(ndimpw))
         if(cbeam.eq.0.d0) then
            call planewavecoef(alpha,beta,nodrpw,pmnp0)
         else
            call gaussianbeamcoef(alpha,beta,cbeam,nodrpw,pmnp0)
         endif
         cgamma=cos(gamma)
         sgamma=sin(gamma)
         if(axialinc.eq.0) then
            do n=1,nodrpw
               nn1=n*(n+1)
               do p=1,2
                  do m=-n,-1
                     mn=nn1+m
                     mnp=2*(mn-1)+p
                     pmnp_nf(mnp)=pmnp0(n+1,-m,p,1)*cgamma+pmnp0(n+1,-m,p,2)*sgamma
                     pmn3mp_nf(mnp)=pmnp0(n+1,-m,3-p,1)*cgamma+pmnp0(n+1,-m,3-p,2)*sgamma
                  enddo
                  do m=0,n
                     mn=nn1+m
                     mnp=2*(mn-1)+p
                     pmnp_nf(mnp)=pmnp0(m,n,p,1)*cgamma+pmnp0(m,n,p,2)*sgamma
                     pmn3mp_nf(mnp)=pmnp0(m,n,3-p,1)*cgamma+pmnp0(m,n,3-p,2)*sgamma
                  enddo
               enddo
            enddo
         else
            do n=1,nodrpw
               do p=1,2
                  mnp=4*(n-1)+p
                  pmnp_nf(mnp)=pmnp0(n+1,1,p,1)*cgamma+pmnp0(n+1,1,p,2)*sgamma
                  pmn3mp_nf(mnp)=pmnp0(n+1,1,3-p,1)*cgamma+pmnp0(n+1,1,3-p,2)*sgamma
                  pmnp_nf(mnp+2)=pmnp0(1,n,p,1)*cgamma+pmnp0(1,n,p,2)*sgamma
                  pmn3mp_nf(mnp+2)=pmnp0(1,n,3-p,1)*cgamma+pmnp0(1,n,3-p,2)*sgamma
               enddo
            enddo
         endif
         deallocate(pmnp0)
         end subroutine nearfieldincidentcoef
!
!  nearfieldpointcalc: if newcalc = 1, generates the reshaped incident, scattered, and
!                      internal field coefficients, and returns with newcalc=0
!                      if newcalc = 0, generates the field at point xg
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
         subroutine nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    gamma,epspw,xg,newcalc,efield,hfield)
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: nsphere,neqns,nodr(nsphere),i,j,k,m,n,p,nn1,mn,nodrpw,newcalc, &
                    insphere
         real (8) :: alpha,beta,cbeam,xsp(nsphere),rpos(3,nsphere),xg(3),xgp(3), &
                     gamma,epspw,rplot,cgamma,sgamma
         complex(8) :: amnp(neqns,2),ri(2,nsphere),efield(3),hfield(3),einc(3),hinc(3),ri0, &
                       ct1,ct2
         complex(8), allocatable :: pmnp0(:,:,:,:),cnmie(:,:,:),amnpt(:,:),cmnpt(:,:)
!
!  initialization operations: newcalc=1
!
         if(newcalc.eq.1) then
            if(allocated(amnp_nf)) deallocate(amnp_nf,cmnp_nf,amn3mp_nf,cmn3mp_nf, &
                     noff_nf,nblk_nf)
            allocate(amnp_nf(neqns),cmnp_nf(neqns),amn3mp_nf(neqns),cmn3mp_nf(neqns), &
                     noff_nf(nsphere),nblk_nf(nsphere))
            noff_nf(1)=0
            do i=1,nsphere
               nblk_nf(i)=2*nodr(i)*(nodr(i)+2)
               if(i.lt.nsphere) noff_nf(i+1)=noff_nf(i)+nblk_nf(i)
            enddo
            cgamma=cos(gamma)
            sgamma=sin(gamma)
            do i=1,nsphere
               ri0=2.d0/(1.d0/ri(1,i)+1.d0/ri(2,i))
               allocate(pmnp0(0:nodr(i)+1,nodr(i),2,2),cnmie(2,2,nodr(i)),&
                        amnpt(2,nodr(i)*(nodr(i)+2)),cmnpt(2,nodr(i)*(nodr(i)+2)))
               call getmiedata(which_sphere=i,sphere_int_mie_coefficients=cnmie)
               do p=1,2
                  pmnp0(0:nodr(i)+1,1:nodr(i),1:2,p) &
                      =reshape(amnp(noff_nf(i)+1:noff_nf(i)+nblk_nf(i),p),(/nodr(i)+2,nodr(i),2/))
               enddo
               if(abs(ri(1,i)-ri(2,i)).gt.1.d-10) then
                  do n=1,nodr(i)
                     nn1=n*(n+1)
                     do p=1,2
                        do m=-n,-1
                           mn=nn1+m
                           amnpt(p,mn)=pmnp0(n+1,-m,p,1)*cgamma+pmnp0(n+1,-m,p,2)*sgamma
                        enddo
                        do m=0,n
                           mn=nn1+m
                           amnpt(p,mn)=pmnp0(m,n,p,1)*cgamma+pmnp0(m,n,p,2)*sgamma
                        enddo
                     enddo
                     do m=-n,n
                        mn=nn1+m
                        ct1=amnpt(1,mn)*cnmie(1,1,n)+amnpt(2,mn)*cnmie(1,2,n)
                        ct2=amnpt(1,mn)*cnmie(2,1,n)+amnpt(2,mn)*cnmie(2,2,n)
                        cmnpt(1,mn)=ct1
                        cmnpt(2,mn)=-(0.d0,1.d0)/ri0*ct2
                     enddo
                  enddo
               else
                  do n=1,nodr(i)
                     nn1=n*(n+1)
                     do p=1,2
                        do m=-n,-1
                           mn=nn1+m
                           amnpt(p,mn)=pmnp0(n+1,-m,p,1)*cgamma+pmnp0(n+1,-m,p,2)*sgamma
                           cmnpt(p,mn)=amnpt(p,mn)*cnmie(p,p,n)
                        enddo
                        do m=0,n
                           mn=nn1+m
                           amnpt(p,mn)=pmnp0(m,n,p,1)*cgamma+pmnp0(m,n,p,2)*sgamma
                           cmnpt(p,mn)=amnpt(p,mn)*cnmie(p,p,n)
                        enddo
                     enddo
                  enddo
               endif
               amnp_nf(noff_nf(i)+1:noff_nf(i)+nblk_nf(i))= &
                      reshape(amnpt(1:2,1:nodr(i)*(nodr(i)+2)),(/nblk_nf(i)/))
               cmnp_nf(noff_nf(i)+1:noff_nf(i)+nblk_nf(i))= &
                      reshape(cmnpt(1:2,1:nodr(i)*(nodr(i)+2)),(/nblk_nf(i)/))
               amn3mp_nf(noff_nf(i)+1:noff_nf(i)+nblk_nf(i))= &
                      reshape(amnpt(2:1:-1,1:nodr(i)*(nodr(i)+2)),(/nblk_nf(i)/))
               cmn3mp_nf(noff_nf(i)+1:noff_nf(i)+nblk_nf(i))= &
                      reshape(cmnpt(2:1:-1,1:nodr(i)*(nodr(i)+2)),(/nblk_nf(i)/))
               deallocate(pmnp0,cnmie,amnpt,cmnpt)
            enddo
            rplot=sqrt(dot_product(xg,xg))
            rplotmax=rplot
            call planewavetruncationorder(rplot,epspw,nodrpw)
            nodrpwmax=nodrpw
            call nearfieldincidentcoef(nodrpw,alpha,beta,gamma,cbeam,epspw)
            newcalc=0
            return
         endif
!
!  point calculation operations: newcalc=0
!  first determine the required order of the incident field expansion, and recalculate
!  field coefficients, if necessary
!
         rplot=sqrt(dot_product(xg,xg))
         rplotmax=max(rplot,rplotmax)
         call planewavetruncationorder(rplot,epspw,nodrpw)
         if(nodrpw.gt.nodrpwmax) then
            nodrpwmax=nodrpw
            call nearfieldincidentcoef(nodrpw,alpha,beta,gamma,cbeam,epspw)
         endif
         efield=0.d0
         hfield=0.d0
!
!  calculate the sphere contribution to the field
!
         call nearfieldspherepart(xg,nsphere,xsp,rpos,ri,&
                    nodr,neqns,insphere,efield,hfield)
!
!  if the point is external to the spheres, calculate the incident field
!
         if(insphere.eq.0) then
            call nearfieldincidentpart(xg,nodrpw,einc,hinc)
            efield=efield+einc
            hfield=hfield+hinc
         endif
         end subroutine nearfieldpointcalc
!
!  nearfieldgridcalc is an MPI--enabled subroutine for calculating field points on a
!  rectangular grid.   Writes the data to nfoutunit.
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  changed so that input positions are defined relative to sphere position file origin, and
!  not the gb focal point.
!
         subroutine nearfieldgridcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    nfplane,nfplanepos0,nfplanevert0,gbfocus,deltax,gamma,nfoutunit,epspw, &
                    nfoutdata,runprintunit)
         use mpidefs
         use mpidata
         use intrinsics
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: nsphere,neqns,nodr(nsphere),nfplane,runprintunit,npoints1,npoints2, &
                    npoints,i,j,k,np23,gcoord(3),rank,numprocs,nrowperproc,nrowrem, &
                    npoints1by5,plotincfield,nsp,nfoutunit,nfoutdata,newcalc
         real (8) :: alpha,beta,cbeam,xsp(nsphere),rpos(3,nsphere),nfplanepos,&
                     nfplanevert(2,2),frowperproc,rowsum,xg(3),xgp(3),deltax,gamma,epspw, &
                     time1,time2,xplot(3,nsphere),xi0,ri0,esquare,xgpmax(3),rplot, &
                     gbfocus(3),nfplanepos0,nfplanevert0(2,2)
         complex(8) :: amnp(neqns,2),ri(2,nsphere),efield(3),hfield(3)
         integer, allocatable :: efindex(:),efnum(:)
         complex(8), allocatable :: efieldrow(:,:),efieldrowt(:,:),hfieldrow(:,:),hfieldrowt(:,:)
         call ms_mpi(mpi_command='size',mpi_size=numprocs)
         call ms_mpi(mpi_command='rank',mpi_rank=rank)
!
!  determine the plane
!
         if(nfplane.eq.1) then
            gcoord=(/2,3,1/)
         elseif(nfplane.eq.2) then
            gcoord=(/3,1,2/)
         else
            gcoord=(/1,2,3/)
         endif
!
!  shift the coordinates to gb focal origin
!
         nfplanevert(1,1)=nfplanevert0(1,1)-gbfocus(gcoord(1))
         nfplanevert(1,2)=nfplanevert0(1,2)-gbfocus(gcoord(1))
         nfplanevert(2,1)=nfplanevert0(2,1)-gbfocus(gcoord(2))
         nfplanevert(2,2)=nfplanevert0(2,2)-gbfocus(gcoord(2))
         nfplanepos=nfplanepos0-gbfocus(gcoord(3))
         xg(gcoord(3))=nfplanepos
!
!  determine the number of points
!
         npoints1=nint((nfplanevert(1,2)-nfplanevert(1,1))/deltax)+1
         npoints2=nint((nfplanevert(2,2)-nfplanevert(2,1))/deltax)+1
         npoints=npoints1*npoints2
!
!  find the maximum point-to-target origin distance and initialize the field calculation
!
         xgp(3)=nfplanepos
         rplotmax=0.d0
         xgpmax=0.d0
         do i=1,npoints1
            xgp(1)=nfplanevert(1,1)+deltax*dble(i-1)
            do j=1,npoints2
               xgp(2)=nfplanevert(2,1)+deltax*dble(j-1)
               rplot=sqrt(dot_product(xgp,xgp))
               if(rplot.gt.rplotmax) then
                  rplotmax=rplot
                  xgpmax=xgp
               endif
            enddo
         enddo
         newcalc=1
         call nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    gamma,epspw,xgpmax,newcalc,efield,hfield)
!
!  determine the intersecting spheres
!
         nsp=0
         do i=1,nsphere
            xi0=abs(rpos(gcoord(3),i)-xg(gcoord(3)))
            if(xi0.le.xsp(i)) then
               nsp=nsp+1
               xplot(1,nsp)=rpos(gcoord(1),i)+gbfocus(gcoord(1))
               xplot(2,nsp)=rpos(gcoord(2),i)+gbfocus(gcoord(2))
               ri0=xsp(i)*xsp(i)-xi0*xi0
               if(ri0.ne.0.) ri0=sqrt(ri0)
               xplot(3,nsp)=ri0
            endif
         enddo
!
!  report to runprintunit
!
         if(rank.eq.0) then
            write(runprintunit,'('' near field calculations'')')
            write(runprintunit,'('' plane, position:'',i5,f9.3)') nfplane,nfplanepos0
            write(runprintunit,'('' rectangular plot vertices:'')')
            write(runprintunit,'('' min:'',3f9.3)') nfplanevert0(1:2,1)
            write(runprintunit,'('' max:'',3f9.3)') nfplanevert0(1:2,2)
            write(runprintunit,'('' number of plotting points, step size:'',i8,f8.3)') npoints, deltax
            write(runprintunit,'('' max plane wave order:'',i5)') nodrpwmax
         endif
!
!  determine the distribution of work among the processors
!
         allocate(efindex(0:numprocs-1),efnum(0:numprocs-1), &
                  efieldrow(3,npoints2),efieldrowt(3,npoints2), &
                  hfieldrow(3,npoints2),hfieldrowt(3,npoints2))
         np23=3*npoints2
         frowperproc=dble(npoints2)/dble(numprocs)
         rowsum=0.
         do i=0,numprocs-1
            efindex(i)=floor(rowsum)
            rowsum=rowsum+frowperproc
         enddo
         do i=0,numprocs-2
            efnum(i)=efindex(i+1)-efindex(i)
         enddo
         efnum(numprocs-1)=npoints2-efindex(numprocs-1)
         npoints1by5=int(npoints1/5.+.5)
!
!  do the calculations and write the results to the file
!
         if(rank.eq.0) then
            write(nfoutunit,*) npoints1,npoints2
            write(nfoutunit,*) nsp
            do i=1,nsp
               write(nfoutunit,'(3e13.5)') xplot(1,i),xplot(2,i),xplot(3,i)
            enddo
            time1=mytime()
         endif
         xg(gcoord(3))=nfplanepos
         newcalc=0
         do i=1,npoints1
            xg(gcoord(1))=nfplanevert(1,1)+deltax*dble(i-1)
            xgp(gcoord(1))=nfplanevert0(1,1)+deltax*dble(i-1)
            efieldrowt=0.d0
            efieldrow=0.d0
            hfieldrowt=0.d0
            hfieldrow=0.d0
            do j=efindex(rank)+1,efindex(rank)+efnum(rank)
               xg(gcoord(2))=nfplanevert(2,1)+deltax*dble(j-1)
               call nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    gamma,epspw,xg,newcalc,efieldrowt(:,j),hfieldrowt(:,j))
            enddo
            call ms_mpi(mpi_command='barrier')
            call ms_mpi(mpi_command='reduce',mpi_send_buf_dc=efieldrowt,mpi_recv_buf_dc=efieldrow,&
                 mpi_number=np23,mpi_rank=0,mpi_operation=ms_mpi_sum)
            if(nfoutdata.ge.2) then
               call ms_mpi(mpi_command='reduce',mpi_send_buf_dc=hfieldrowt,mpi_recv_buf_dc=hfieldrow,&
                    mpi_number=np23,mpi_rank=0,mpi_operation=ms_mpi_sum)
            endif
            if(rank.eq.0) then
               do j=1,npoints2
                  xgp(gcoord(2))=nfplanevert0(2,1)+deltax*dble(j-1)
                  if(nfoutdata.eq.0) then
                     esquare=dot_product(efieldrow(:,j),efieldrow(:,j))
                     write(nfoutunit,'(2f9.4,e13.5)') xgp(gcoord(1)),xgp(gcoord(2)),esquare
                  elseif(nfoutdata.eq.1) then
                     write(nfoutunit,'(2f9.4,6e13.5)') xgp(gcoord(1)),xgp(gcoord(2)),efieldrow(:,j)
                  else
                     write(nfoutunit,'(2f9.4,12e13.5)') xgp(gcoord(1)),xgp(gcoord(2)),efieldrow(:,j), &
                                                        hfieldrow(:,j)
                  endif
               enddo
               if(mod(i,npoints1by5).eq.0) then
                  k=i/npoints1by5
                  time2=(mytime()-time1)*dble(5-k)/dble(k)
                  call timewrite(runprintunit,' estimated time remaining:',time2)
               endif
            endif
         enddo
         deallocate(efindex,efnum,efieldrow,efieldrowt,hfieldrow,hfieldrowt)
         end subroutine nearfieldgridcalc
!
!  nearfieldaverage is an MPI--enabled subroutine for calculating average field
!  values along a line.
!
         subroutine nearfieldaverage(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    nfplane,nfplaneposstart0,nfplaneposend0,numberplanes,nfplanevert0,gbfocus, &
                    deltax,gamma,epspw,runprintunit,efieldavez,hfieldavez,svecavez)
         use mpidefs
         use mpidata
         use intrinsics
         use specialfuncs
         use spheredata
         use miecoefdata
         use numconstants
         implicit none
         integer :: nsphere,neqns,nodr(nsphere),nfplane,runprintunit,npoints1,npoints2, &
                    npoints,i,j,k,np23,gcoord(3),rank,numprocs,nrowperproc,nrowrem, &
                    npoints1by5,plotincfield,nsp,nfoutunit,nfoutdata,newcalc,nsend
         integer :: numberplanes
         real (8) :: alpha,beta,cbeam,xsp(nsphere),rpos(3,nsphere),nfplanepos,&
                     nfplanevert(2,2),frowperproc,rowsum,xg(3),xgp(3),deltax,gamma,epspw, &
                     time1,time2,xplot(3,nsphere),xi0,ri0,esquare,xgpmax(3),rplot, &
                     gbfocus(3),nfplanepos0,nfplanevert0(2,2),nfplaneposstart,nfplaneposend, &
                     deltanfplane,nfplaneposstart0,nfplaneposend0,svec(3),svecave(3)
         real (8) :: svecavez(3,numberplanes)
         real(8), allocatable :: xypoint(:,:)
         complex(8) :: amnp(neqns,2),ri(2,nsphere),efield(3),hfield(3),hfieldave(3),efieldave(3)
         complex(8) :: efieldavez(3,numberplanes),hfieldavez(3,numberplanes)
         integer, allocatable :: efindex(:),efnum(:)
         call ms_mpi(mpi_command='size',mpi_size=numprocs)
         call ms_mpi(mpi_command='rank',mpi_rank=rank)
!
!  determine the plane
!
         if(nfplane.eq.1) then
            gcoord=(/2,3,1/)
         elseif(nfplane.eq.2) then
            gcoord=(/3,1,2/)
         else
            gcoord=(/1,2,3/)
         endif
!
!  shift the coordinates to gb focal origin
!
         nfplanevert(1,1)=nfplanevert0(1,1)-gbfocus(gcoord(1))
         nfplanevert(1,2)=nfplanevert0(1,2)-gbfocus(gcoord(1))
         nfplanevert(2,1)=nfplanevert0(2,1)-gbfocus(gcoord(2))
         nfplanevert(2,2)=nfplanevert0(2,2)-gbfocus(gcoord(2))
         nfplaneposstart=nfplaneposstart0-gbfocus(gcoord(3))
         nfplaneposend=nfplaneposend0-gbfocus(gcoord(3))
         xg(gcoord(3))=nfplanepos
!
!  determine the number of points
!
         npoints1=nint((nfplanevert(1,2)-nfplanevert(1,1))/deltax)+1
         npoints2=nint((nfplanevert(2,2)-nfplanevert(2,1))/deltax)+1
         npoints=npoints1*npoints2
         allocate(efindex(0:numprocs-1),efnum(0:numprocs-1),xypoint(2,npoints))
         frowperproc=dble(npoints)/dble(numprocs)
         rowsum=0.
         do i=0,numprocs-1
            efindex(i)=floor(rowsum)
            rowsum=rowsum+frowperproc
         enddo
         do i=0,numprocs-2
            efnum(i)=efindex(i+1)-efindex(i)
         enddo
         efnum(numprocs-1)=npoints-efindex(numprocs-1)

         xgp(gcoord(3))=max(abs(nfplaneposstart),abs(nfplaneposend))
         rplotmax=0.d0
         xgpmax=0.d0
         k=0
         do i=1,npoints1
            xgp(gcoord(1))=nfplanevert(1,1)+deltax*dble(i-1)
            do j=1,npoints2
               xgp(gcoord(2))=nfplanevert(2,1)+deltax*dble(j-1)
               k=k+1
               xypoint(1,k)=xgp(gcoord(1))
               xypoint(2,k)=xgp(gcoord(2))
               rplot=sqrt(dot_product(xgp,xgp))
               if(rplot.gt.rplotmax) then
                  rplotmax=rplot
                  xgpmax=xgp
               endif
            enddo
         enddo
         newcalc=1
         call nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    gamma,epspw,xgpmax,newcalc,efield,hfield)
         newcalc=0

         do k=1,numberplanes
            if(numberplanes.eq.1) then
               nfplanepos=nfplaneposstart
            else
               nfplanepos=nfplaneposstart+(nfplaneposend-nfplaneposstart)*(k-1)/dble(numberplanes-1)
            endif
!
!  find the maximum point-to-target origin distance and initialize the field calculation
!
            xg(3)=nfplanepos
            efieldave=0.
            hfieldave=0.
            svecave=0.
            do i=efindex(rank)+1,efindex(rank)+efnum(rank)
               xg(gcoord(1))=xypoint(1,i)
               xg(gcoord(2))=xypoint(2,i)
               efield=0.d0
               hfield=0.d0
               call nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    gamma,epspw,xg,newcalc,efield,hfield)
               efieldave=efieldave+efield
               hfieldave=hfieldave+hfield
               hfield=conjg(hfield)/2.d0
               svec(1)=-efield(3)*hfield(2)+efield(2)*hfield(3)
               svec(2)=efield(3)*hfield(1)-efield(1)*hfield(3)
               svec(3)=-efield(2)*hfield(1)+efield(1)*hfield(2)
               svecave=svecave+svec
            enddo
            call ms_mpi(mpi_command='barrier')
            efield=0.d0
            hfield=0.d0
            svec=0.d0
            call ms_mpi(mpi_command='reduce',mpi_send_buf_dc=efieldave,mpi_recv_buf_dc=efield,&
                    mpi_number=3,mpi_rank=0,mpi_operation=ms_mpi_sum)
            call ms_mpi(mpi_command='reduce',mpi_send_buf_dc=hfieldave,mpi_recv_buf_dc=hfield,&
                    mpi_number=3,mpi_rank=0,mpi_operation=ms_mpi_sum)
            call ms_mpi(mpi_command='reduce',mpi_send_buf_dp=svecave,mpi_recv_buf_dp=svec,&
                    mpi_number=3,mpi_rank=0,mpi_operation=ms_mpi_sum)
            call ms_mpi(mpi_command='barrier')
            if(rank.eq.0) then
               efield=efield/dble(npoints)
               hfield=hfield/dble(npoints)
               svec=svec/dble(npoints)
               efieldavez(:,k)=efield
               hfieldavez(:,k)=hfield
               svecavez(:,k)=svec
               i=gcoord(3)
               write(runprintunit,'('' plane:'',i5,f9.3,2e12.4)') k,nfplanepos, &
                   svec(i)
               call flush(runprintunit)
            endif
         enddo
         nsend=3*numberplanes
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=efieldavez, &
              mpi_number=nsend,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=hfieldavez, &
              mpi_number=nsend,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=svecavez, &
              mpi_number=nsend,mpi_rank=0)
         call ms_mpi(mpi_command='barrier')

         deallocate(efindex,efnum,xypoint)
         end subroutine nearfieldaverage


      end module nearfield
!
!  module solver: subroutines for solving interaction equations for fixed orientation
!  and T matrix problems
!
!
!  last revised: 15 January 2011
!
      module solver
      implicit none

      contains
!
!  tmatrixsoln: calculation of T matrix via solution of interaction equations for
!  a generalized plane wave expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: call for sphereqeff changed
!
         subroutine tmatrixsoln(neqns,nsphere,nodr,nodrt,xsp,rpos,epssoln,epscon,niter,&
                    calctmatrix,tmatrixfile,fftranpresent,niterstep,qext,qabs,qsca,istat)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         use miecoefdata
         use spheredata
         use translation
         use scatprops
         implicit none
         integer :: iter,niter,neqns,nsphere,nodr(nsphere),ntran(nsphere),nodrt, &
                    nodrmax,i,ierr,istat,m,n,p,k,l,q,noff,nblk,ma,na,ka,la,mn,istart,iunit, &
                    rank,iexit(1),calctmatrix,lt,kt,qt,nt,mt,it,nodrtt,lstart(1),numsolns,isoln, &
                    isolnstart,igroup,ngroup,rgrank,lold,grank,nsend,nodrta(1), &
                    fftranpresent,niterstep
         integer, allocatable :: lindex(:),kindex(:)
         real(8) :: eps,err,qext(nsphere),qabs(nsphere),qsca(nsphere),xsp(nsphere),xv, &
                    rpos(3,nsphere),xij(3),qabsklq(nsphere),qscaklq(nsphere),qextklq(nsphere), &
                    f2,qexttot,qabstot,qscatot,qextold(1),qscaold(1),errqe,errqs, &
                    timetran,timesolve,time1,time2,epssoln,epscon,dtemp(4),timeorder,&
                    at1,at2,at3,at4
         real(8) :: qextl(nsphere),qabsl(nsphere),qscal(nsphere)
         real(8), allocatable :: qextgroup(:,:),qabsgroup(:,:),qscagroup(:,:)
         complex(8) :: amnp(neqns),pmnp(neqns),pmnpan(neqns)
         complex(8), allocatable :: pmnp0(:,:,:),ac(:,:,:,:),pmnpt(:,:,:),amnp0(:,:,:), &
                    amnp0group(:,:,:)
         character*30 :: tmatrixfile
         character*4 :: timeunit
         data istart,iexit/1,0/
         rank=base_rank
         rgrank=root_group_rank
         grank=group_rank
         ngroup=number_groups
         call getrunparameters(run_print_unit=iunit)
         xv=(sum(xsp**3.d0))**(1.d0/3.d0)
         nodrmax=maxval(nodr)
         qext=0.d0
         qabs=0.d0
         qsca=0.d0
         qextold=0.d0
         qscaold=0.d0
!
!  perform T matrix file operations as needed
!
         if(rank.eq.0) then
            if(calctmatrix.eq.1) then
               open(3,file=tmatrixfile)
               write(3,'(i4)') nodrt
               lstart(1)=1
            else
               open(3,file=tmatrixfile)
               write(iunit,'('' finding end of record to file '',a)') tmatrixfile
               read(3,*) nodrtt
               do l=1,nodrt
                  do k=-l,l
                     do q=1,2
                        read(3,'(3i5)',end=20,err=20) lt,kt,qt
                        do n=1,l
                           do m=-n,n
                              read(3,'(2i5,4e17.9)',end=20,err=20) nt,mt,at1,at2,at3,at4
                           enddo
                        enddo
                     enddo
                  enddo
                  do i=1,nsphere
                     read(3,'(i5,3e17.9)',end=20,err=20) it,qextl(i),qabsl(i),qscal(i)
                  enddo
                  qext=qext+qextl
                  qabs=qabs+qabsl
                  qsca=qsca+qscal
               enddo
20             close(3)
               open(3,file=tmatrixfile)
               qextold(1)=0.d0
               qabstot=0.d0
               do i=1,nsphere
                  qextold=qextold+qext(i)*xsp(i)*xsp(i)/xv/xv
                  qabstot=qabstot+qabs(i)*xsp(i)*xsp(i)/xv/xv
               enddo
               qscaold(1)=qextold(1)-qabstot
               lstart(1)=lt
               read(3,*) nodrtt
               write(iunit,'('' calculations begin with order '',i5)') lstart(1)
               do l=1,lstart(1)-1
                  do k=-l,l
                     do q=1,2
                        read(3,'(3i5)') lt,kt,qt
                        do n=1,l
                           do m=-n,n
                              read(3,'(2i5,4e17.9)') nt,mt,at1,at2,at3,at4
                           enddo
                        enddo
                     enddo
                  enddo
                  do i=1,nsphere
                     read(3,'(i5,3e17.9)') it,at1,at2,at3
                  enddo
               enddo
            endif
         endif
         call ms_mpi(mpi_command='bcast',mpi_send_buf_i=lstart,mpi_number=1,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qextold,mpi_number=1,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qscaold,mpi_number=1,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qext,mpi_number=nsphere,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qabs,mpi_number=nsphere,mpi_rank=0)
         call ms_mpi(mpi_command='bcast',mpi_send_buf_dp=qsca,mpi_number=nsphere,mpi_rank=0)
         call ms_mpi(mpi_command='barrier')
         allocate(amnp0group(2,nodrt*(nodrt+2),0:ngroup-1),qextgroup(nsphere,0:ngroup-1), &
                  qabsgroup(nsphere,0:ngroup-1),qscagroup(nsphere,0:ngroup-1))
         numsolns=2*nodrt*(nodrt+2)
         allocate(lindex(numsolns),kindex(numsolns))
!
!  find the starting point
!
         i=0
         do l=1,nodrt
            do k=-l,l
               do q=1,2
                  i=i+1
                  lindex(i)=l
                  kindex(i)=k
               enddo
            enddo
         enddo
         do i=1,numsolns
            if(lindex(i).eq.lstart(1).and.kindex(i).eq.-lstart(1)) exit
         enddo
         isolnstart=i
         qextl=0.d0
         qabsl=0.d0
         qscal=0.d0
         lold=0
!
!  begin the loop over RHS of the interaction equations.   The solutions are distributed
!  among ngroup groups of processors
!
         do isoln=isolnstart,numsolns,ngroup
            if(rank.eq.0) timeorder=mytime()
            do igroup=0,ngroup-1
               l=lindex(isoln+igroup)
               k=kindex(isoln+igroup)
               q=mod(isoln+igroup-1,2)+1
!
!  calculate the RHS
!
               if(l.eq.-k.and.q.eq.1) then
                  if(allocated(ac)) deallocate(ac,amnp0)
                  allocate(ac(2,nodrmax*(nodrmax+2),-l:l,nsphere),amnp0(0:l+1,l,2))
                  do i=1,nsphere
                     xij=rpos(:,i)
                     call gentrancoef(1,xij,(1.d0,0.d0),1,nodrmax,l,l,0,0,ac(1,1,-l,i))
                  enddo
               endif
               if(igroup.eq.rgrank) then
                  if(k.le.-1) then
                     ka=l+1
                     la=-k
                  else
                     ka=k
                     la=l
                  endif
                  noff=0
                  do i=1,nsphere
                     nblk=nodr(i)*(nodr(i)+2)*2
                     allocate(pmnp0(0:nodr(i)+1,nodr(i),2))
                     do p=1,2
                        do n=1,nodr(i)
                           do m=-n,n
                              mn=n*(n+1)+m
                              if(m.le.-1) then
                                 ma=n+1
                                 na=-m
                              else
                                 ma=m
                                 na=n
                              endif
                              pmnp0(ma,na,p)=ac(abs(p-q)+1,mn,k,i)
                           enddo
                        enddo
                     enddo
                     pmnp(noff+1:noff+nblk)=reshape(pmnp0,(/nblk/))
                     deallocate(pmnp0)
                     noff=noff+nblk
                  enddo
!
!  multiply RHS by mie coefficients
!
                  call miecoeffmult(1,nsphere,neqns,pmnp,pmnpan)
                  amnp=pmnpan
!
!  call the solver
!
                  if(fftranpresent.eq.1) then
                     call cbicgff(neqns,nsphere,niter,epssoln,pmnpan,amnp,0, &
                          niterstep,iter,err)
                  else
                     call cbicg(neqns,nsphere,niter,epssoln,pmnpan,amnp,0,iter,err)
                  endif
                  if(iter.gt.niter.or.err.gt.epssoln) istat=1
                  call sphereqeff(nsphere,neqns,nodr,nodrmax,xsp,amnp,amnp, &
                         pmnp,pmnp,qextklq,qabsklq,qscaklq)
                  qextgroup(1:nsphere,igroup)=qextklq(1:nsphere)
                  qabsgroup(1:nsphere,igroup)=qabsklq(1:nsphere)
                  qscagroup(1:nsphere,igroup)=qscaklq(1:nsphere)
!
!  compute the target-based expansion
!
                  ntran=l
                  amnp0=0.d0
                  call amncommonorigin(neqns,nsphere,nodr,ntran,l,rpos, &
                                       amnp,amnp0)
                  do n=1,l
                     do m=-n,n
                        if(m.le.-1) then
                           ma=n+1
                           na=-m
                        else
                           ma=m
                           na=n
                        endif
                        mn=n*(n+1)+m
                        do p=1,2
                           amnp0group(p,mn,igroup)=amnp0(ma,na,p)
                        enddo
                     enddo
                  enddo
               endif
            enddo
!
!  send the solutions to the rank 0 processor
!
            call ms_mpi(mpi_command='barrier')
            if(grank.eq.0) then
               if(rank.ne.0) then
                  l=lindex(isoln+rgrank)
                  nblk=l*(l+2)
                  nsend=2*nblk
                  call ms_mpi(mpi_command='send',mpi_send_buf_dc=amnp0group(1,1,rgrank),&
                       mpi_number=nsend,mpi_rank=0,mpi_comm=root_group_comm)
                  call ms_mpi(mpi_command='send',mpi_send_buf_dp=qextgroup(1,rgrank),&
                       mpi_number=nsphere,mpi_rank=0,mpi_comm=root_group_comm)
                  call ms_mpi(mpi_command='send',mpi_send_buf_dp=qabsgroup(1,rgrank),&
                       mpi_number=nsphere,mpi_rank=0,mpi_comm=root_group_comm)
                  call ms_mpi(mpi_command='send',mpi_send_buf_dp=qscagroup(1,rgrank),&
                       mpi_number=nsphere,mpi_rank=0,mpi_comm=root_group_comm)
               else
                  do igroup=1,ngroup-1
                     l=lindex(isoln+igroup)
                     nblk=l*(l+2)
                     nsend=2*nblk
                     call ms_mpi(mpi_command='recv',mpi_recv_buf_dc=amnp0group(1,1,igroup),&
                          mpi_number=nsend,mpi_rank=igroup,mpi_comm=root_group_comm)
                     call ms_mpi(mpi_command='recv',mpi_recv_buf_dp=qextgroup(1,igroup),&
                          mpi_number=nsphere,mpi_rank=igroup,mpi_comm=root_group_comm)
                     call ms_mpi(mpi_command='recv',mpi_recv_buf_dp=qabsgroup(1,igroup),&
                          mpi_number=nsphere,mpi_rank=igroup,mpi_comm=root_group_comm)
                     call ms_mpi(mpi_command='recv',mpi_recv_buf_dp=qscagroup(1,igroup),&
                          mpi_number=nsphere,mpi_rank=igroup,mpi_comm=root_group_comm)
                  enddo
               endif
            endif
            call ms_mpi(mpi_command='barrier')
!
!  write results, check for convergence
!
            if(rank.eq.0) then
               do igroup=0,ngroup-1
                  l=lindex(isoln+igroup)
                  k=kindex(isoln+igroup)
                  q=mod(isoln+igroup-1,2)+1
                  qextl=qextl+qextgroup(1:nsphere,igroup)
                  qabsl=qabsl+qabsgroup(1:nsphere,igroup)
                  qscal=qscal+qscagroup(1:nsphere,igroup)
                  qext=qext+qextgroup(1:nsphere,igroup)
                  qabs=qabs+qabsgroup(1:nsphere,igroup)
                  qsca=qsca+qscagroup(1:nsphere,igroup)
                  write(3,'(3i5)') l,k,q
                  do n=1,l
                     do m=-n,n
                        mn=n*(n+1)+m
                        write(3,'(2i5,4e17.9)') n,m,amnp0group(1,mn,igroup), &
                                                amnp0group(2,mn,igroup)
                     enddo
                  enddo
                  if(istart.eq.1.and.igroup.eq.0) then
                     time1=mytime()-timeorder
                     call timewrite(iunit,' time per group solution:',time1)
                     time2=time1*dble(numsolns-isolnstart)/dble(ngroup)
                     call timewrite(iunit,' estimated t matrix calcuation time:',time2)
                     write(iunit,'(''  n   # its  qext         qabs'',&
                         &''         qsca      error     est. time rem.'')')
                     call flush(iunit)
                     istart=0
                  endif
                  if(igroup.eq.0) then
                     timeorder=mytime()-timeorder
                     time2=timeorder*dble(numsolns-isoln)/dble(ngroup)
                     if(time2.gt.3600.d0) then
                        time2=time2/3600.d0
                        timeunit=' hrs'
                     elseif(time2.gt.60.d0) then
                        time2=time2/60.d0
                        timeunit=' min'
                     else
                        timeunit=' sec'
                     endif
                  endif
                  iexit(1)=0
                  if(k.eq.l.and.q.eq.2) then
                     qexttot=0.d0
                     qabstot=0.d0
                     do i=1,nsphere
                        qexttot=qexttot+qext(i)*xsp(i)*xsp(i)/xv/xv
                        qabstot=qabstot+qabs(i)*xsp(i)*xsp(i)/xv/xv
                        write(3,'(i5,3e17.9)') i,qextl(i),qabsl(i),qscal(i)
                     enddo
                     qextl=0.d0
                     qabsl=0.d0
                     qscal=0.d0
                     qscatot=qexttot-qabstot
                     errqe=qexttot-qextold(1)
                     errqs=qscatot-qscaold(1)
                     err=max(errqe,errqs)
                     write(iunit,'(i4,i5,4e13.5,f8.2,a4)') l,iter,qexttot,qabstot, &
                        qscatot,err,time2,timeunit
                     call flush(iunit)
                     qextold(1)=qexttot
                     qscaold(1)=qscatot
                     if(err.le.epscon) iexit(1)=1
                  endif
                  if(iexit(1).eq.1) then
                     nodrt=l
                     exit
                  endif
               enddo
            endif
            call ms_mpi(mpi_command='bcast',mpi_send_buf_i=iexit,mpi_number=1,mpi_rank=0)
            call ms_mpi(mpi_command='barrier')
            if(iexit(1).eq.1) then
!
!  solution has converged
!
               deallocate(amnp0group,qextgroup,qabsgroup,qscagroup,lindex,kindex,ac,amnp0)
               nodrta(1)=nodrt
               call ms_mpi(mpi_command='bcast',mpi_send_buf_i=nodrta,mpi_number=1,mpi_rank=0)
               nodrt=nodrta(1)
               if(rank.eq.0) then
                  write(iunit,'('' T matrix converged, order:'',i5)') nodrt
                  close(3)
                  open(3,file=tmatrixfile,form='formatted',access='direct',recl=4)
                  write(3,'(i4)',rec=1) nodrt
                  close(3)
               endif
               return
            endif
         enddo
         deallocate(amnp0group,qextgroup,qabsgroup,qscagroup,lindex,kindex,ac,amnp0)
         if(rank.eq.0) then
            write(*,'('' T matrix did not converge to set epsilon'')')
            close(3)
         endif
         end subroutine tmatrixsoln
!
!  solution of interaction equations for a fixed orientation
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: modification of efficiency calculation, to calculate
!           polarized components
!  30 March 2011: took out gbfocus argument: this is not needed since positions are defined
!  relative to the gb focus.
!  20 April 2011: used 2-group MPI formulation
!
         subroutine fixedorsoln(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,&
                    eps,epstran,niter,amnp,qext,qabs,qsca,maxerr,maxiter,iterwrite, &
                    fftranpresent,niterstep,istat)
         use mpidefs
         use mpidata
         use intrinsics
         use numconstants
         use specialfuncs
         use miecoefdata
         use translation
         use scatprops
         implicit none
         integer :: iter,niter,neqns,nodrmax,k,nsphere,i,ierr,istat,rank,maxiter,iterwrite
         integer :: nodr(nsphere),m1,n1,p,rgrank,grank,ngroup,sendrank,numprocs, &
                    fftranpresent,niterstep
         real(8) :: alpha,beta,eps,err,qext(nsphere,3),maxerr,&
                    qabs(nsphere,3),qsca(nsphere,3),cbeam,gbfocus(3),epstran
         real(8) :: xsp(nsphere), rpos(3,nsphere),maxerra(1)
         complex(8) :: amnp(neqns,2)
         complex(8), allocatable :: pmnp(:,:),pmnpan(:)
         rank=base_rank
         rgrank=root_group_rank
         grank=group_rank
         ngroup=number_groups
         numprocs=number_proc
         sendrank=numprocs/2
         nodrmax=maxval(nodr)
         allocate(pmnp(neqns,2))
         gbfocus=0.d0
         if(cbeam.eq.0.d0) then
            call sphereplanewavecoef(nsphere,neqns,nodr,nodrmax,alpha,beta,rpos,pmnp)
         else
            call spheregaussianbeamcoef(nsphere,neqns,nodr,alpha,beta,cbeam, &
                    rpos,gbfocus,epstran,pmnp)
         endif
         istat=0
         maxiter=0
         maxerr=0.
!
!  calculate the two solutions
!
         allocate(pmnpan(neqns))
         if(ngroup.eq.1) then
            do k=1,2
               call miecoeffmult(1,nsphere,neqns,pmnp(1,k),pmnpan)
               amnp(1:neqns,k)=pmnpan(1:neqns)
               if(fftranpresent.eq.1) then
                  call cbicgff(neqns,nsphere,niter,eps,pmnpan,amnp(1,k),iterwrite, &
                               niterstep,iter,err)
               else
                  call cbicg(neqns,nsphere,niter,eps,pmnpan,amnp(1,k),iterwrite, &
                               iter,err)
               endif
               maxiter=max(iter,maxiter)
               maxerr=max(err,maxerr)
               if(iter.gt.niter.or.err.gt.eps) istat=1
               call ms_mpi(mpi_command='barrier')
            enddo
         else
            k=rgrank+1
            call miecoeffmult(1,nsphere,neqns,pmnp(1,k),pmnpan)
            amnp(1:neqns,k)=pmnpan(1:neqns)
            if(fftranpresent.eq.1) then
               call cbicgff(neqns,nsphere,niter,eps,pmnpan,amnp(1,k),iterwrite, &
                            niterstep,iter,err)
            else
               call cbicg(neqns,nsphere,niter,eps,pmnpan,amnp(1,k),iterwrite, &
                            iter,err)
            endif
            maxiter=max(iter,maxiter)
            maxerr=max(err,maxerr)
            call ms_mpi(mpi_command='barrier')
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=amnp(1,2),&
                  mpi_number=neqns,mpi_rank=sendrank)
         endif
         deallocate(pmnpan)
!
!  efficiency factor calculations
!
         call sphereqeff(nsphere,neqns,nodr,nodrmax,xsp,amnp(1,1),amnp(1,1), &
                      pmnp(1,1),pmnp(1,1),qext(1,1),qabs(1,1),qsca(1,1))
         call sphereqeff(nsphere,neqns,nodr,nodrmax,xsp,amnp(1,2),amnp(1,2), &
                      pmnp(1,2),pmnp(1,2),qext(1,2),qabs(1,2),qsca(1,2))
         call sphereqeff(nsphere,neqns,nodr,nodrmax,xsp,amnp(1,1),amnp(1,2), &
                      pmnp(1,1),pmnp(1,2),qext(1,3),qabs(1,3),qsca(1,3))
         call ms_mpi(mpi_command='barrier')
         deallocate(pmnp)
         end subroutine fixedorsoln
!
! hybrid bcgm, using far field translation
! november 2011
!
         subroutine cbicgff(neqns,nsphere,niter,eps,pnp,anp,iterwrite,niterstep,iter,err)
         use mpidefs
         use mpidata
         use intrinsics
         use spheredata
         use miecoefdata
         use numconstants
         use specialfuncs
         use translation
         implicit none
         integer :: neqns,niter,iter,nsphere,writetime,nodr(nsphere),noff(nsphere+1),&
                    nblk(nsphere),rank,ierr,iexit,i,j,m,n,p,iunit,iterwrite,ip1,ip2, &
                    np1,np2,nsend,numprocs,grank,istore,itermax,istep,niterstep
         real(8) :: eps,err,erra(1),enorm,time1,time2,epsstep,errstep
         complex(8)  :: pnp(neqns),anp(neqns),gnp(neqns),gnpold(neqns),pgnp(neqns), &
                        cr(neqns)
         data writetime/0/
         rank=base_rank
         grank=group_rank
         numprocs=proc_per_group
         call getmiedata(sphere_order=nodr,sphere_block=nblk,sphere_block_offset=noff)
         if(rank.eq.0) then
            call getrunparameters(run_print_unit=iunit)
         endif
         err=0.d0
         iter=0
         enorm=dot_product(pnp,pnp)
         gnpold=0.d0
         if(enorm.eq.0.d0) return
         gnp=0.d0
         ip1=mpi_sphere_index(grank)+1
         ip2=mpi_sphere_index(grank)+mpi_sphere_number(grank)
         do i=ip1,ip2
            do j=1,nsphere
               if(i.ne.j) then
                  call fftranslationerror(anp(noff(j)+1:noff(j)+nblk(j)), &
                   gnp(noff(i)+1:noff(i)+nblk(i)),j,i,nodr(j),nodr(i))
               endif
            enddo
         enddo
         do i=0,numprocs-1
            ip1=mpi_sphere_index(i)+1
            ip2=mpi_sphere_index(i)+mpi_sphere_number(i)
            np1=noff(ip1)+1
            np2=noff(ip2)+nblk(ip2)
            nsend=np2-np1+1
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=gnp(np1:np2),mpi_number=nsend, &
                 mpi_rank=i,mpi_comm=group_comm)
         enddo
         call miecoeffmult(1,nsphere,neqns,gnp,gnp)

         iter=0
         epsstep=eps
         istep=0
         do
            istep=istep+1
            gnpold=gnp
            pgnp=pnp+gnp
            call cbicg(neqns,nsphere,niterstep,epsstep,pgnp,anp,0,itermax,errstep)
            iter=iter+min(itermax,niterstep)
            gnp=0.d0
            ip1=mpi_sphere_index(grank)+1
            ip2=mpi_sphere_index(grank)+mpi_sphere_number(grank)
            do i=ip1,ip2
               do j=1,nsphere
                  if(i.ne.j) then
                     call fftranslationerror(anp(noff(j)+1:noff(j)+nblk(j)), &
                      gnp(noff(i)+1:noff(i)+nblk(i)),j,i,nodr(j),nodr(i))
                  endif
               enddo
            enddo
            do i=0,numprocs-1
               ip1=mpi_sphere_index(i)+1
               ip2=mpi_sphere_index(i)+mpi_sphere_number(i)
               np1=noff(ip1)+1
               np2=noff(ip2)+nblk(ip2)
               nsend=np2-np1+1
               call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=gnp(np1:np2),mpi_number=nsend, &
                    mpi_rank=i,mpi_comm=group_comm)
            enddo
            call miecoeffmult(1,nsphere,neqns,gnp,gnp)
            err=dot_product(gnp-gnpold,gnp-gnpold)/enorm
            if(rank.eq.0.and.iterwrite.eq.1) then
               write(iunit,'('' step,iteration,bcgm err,correc err:'',2i5,2e13.5)') &
                           istep,iter,errstep,err
               call flush(iunit)
            endif
            epsstep=eps
            err=max(err,errstep)
            if((err.lt.eps).or.iter.gt.niter) exit
         enddo
         end subroutine cbicgff

!
! iteration solver
! generalized complex biconjugate gradient method
! original code: Piotr Flatau, although not much remains.
! specialized to the multiple sphere problem
!
!
!  last revised: 15 January 2011
!  october 2011: translation calls modified
!
         subroutine cbicg(neqns,nsphere,niter,eps,pnp,anp,iterwrite,iter,err)
         use mpidefs
         use mpidata
         use intrinsics
         use spheredata
         use miecoefdata
         use numconstants
         use specialfuncs
         use translation
         implicit none
         integer :: neqns,niter,iter,nsphere,writetime,nodr(nsphere),noff(nsphere+1),&
                    nblk(nsphere),rank,ierr,iexit,i,j,m,n,p,iunit,iterwrite,ip1,ip2, &
                    np1,np2,nsend,numprocs,grank,istore
         real(8) :: eps,err,erra(1),enorm,time1,time2
         complex(8)  :: pnp(neqns),anp(neqns)
         complex(8) :: cak(1),csk,cbk,csk2(1)
         complex(8) :: cr(neqns),cp(neqns),cw(neqns),cq(neqns),cap(neqns),caw(neqns), &
                       crt(neqns),capt(neqns),cawt(neqns),ccw(neqns)
         data writetime/0/
         rank=base_rank
         grank=group_rank
         numprocs=proc_per_group
         call getmiedata(sphere_order=nodr,sphere_block=nblk,sphere_block_offset=noff)
         ip1=mpi_sphere_index(grank)+1
         ip2=mpi_sphere_index(grank)+mpi_sphere_number(grank)
         np1=noff(ip1)+1
         np2=noff(ip2)+nblk(ip2)
         nsend=np2-np1+1
         crt=0.d0
         iexit=0
         if(rank.eq.0) then
            call getrunparameters(run_print_unit=iunit)
         endif
         err=0.d0
         iter=0
         enorm=dot_product(pnp,pnp)
         cr=0.d0
         if(enorm.eq.0.d0) return
         do i=ip1,ip2
            do j=1,nsphere
               if(i.ne.j) then
                  call rottranjtoi(anp(noff(j)+1:noff(j)+nblk(j)), &
                       cr(noff(i)+1:noff(i)+nblk(i)),j,i,nodr(j),nodr(i),1,1)
               endif
            enddo
         enddo
         do i=0,numprocs-1
            ip1=mpi_sphere_index(i)+1
            ip2=mpi_sphere_index(i)+mpi_sphere_number(i)
            np1=noff(ip1)+1
            np2=noff(ip2)+nblk(ip2)
            nsend=np2-np1+1
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=cr(np1:np2),mpi_number=nsend, &
                 mpi_rank=i,mpi_comm=group_comm)
         enddo
         call miecoeffmult(1,nsphere,neqns,cr,cr)
         cr=pnp-anp+cr
         cq=conjg(cr)
         cw=cq
         cp=cr
         csk=dot_product(conjg(cr),cr)
         if(cdabs(csk).eq.0.d0) return
!
!  here starts the main iteration loop
!
         do iter=1,niter
            ip1=mpi_sphere_index(grank)+1
            ip2=mpi_sphere_index(grank)+mpi_sphere_number(grank)
            np1=noff(ip1)+1
            np2=noff(ip2)+nblk(ip2)
            nsend=np2-np1+1
            cak(1)=(0.d0,0.d0)
            cawt=(0.d0,0.d0)
            capt=(0.d0,0.d0)
            if(rank.eq.0) then
               if(writetime.eq.0) time1=mytime()
            endif
            ccw=conjg(cw)
            cap=0.d0
            caw=0.d0
            call miecoeffmult(1,nsphere,neqns,ccw,ccw)
            do i=ip1,ip2
               do j=1,nsphere
                  if(i.ne.j) then
                     call rottrantwojtoi(cp(noff(j)+1:noff(j)+nblk(j)), &
                       ccw(noff(j)+1:noff(j)+nblk(j)), &
                       cap(noff(i)+1:noff(i)+nblk(i)), &
                       caw(noff(i)+1:noff(i)+nblk(i)), &
                       j,i,nodr(j),nodr(i))
!                     call rottranjtoi(cp(noff(j)+1:noff(j)+nblk(j)), &
!                       cap(noff(i)+1:noff(i)+nblk(i)),j,i,nodr(j),nodr(i),1,1)
!                     call rottranjtoi(ccw(noff(j)+1:noff(j)+nblk(j)), &
!                       caw(noff(i)+1:noff(i)+nblk(i)),j,i,nodr(j),nodr(i),-1,-1)
                  endif
               enddo
            enddo
            call miecoeffmult(ip1,ip2,neqns,cap,cap)
            cap(np1:np2)=cp(np1:np2)-cap(np1:np2)
            caw(np1:np2)=cw(np1:np2)-conjg(caw(np1:np2))
            cak(1)=dot_product(cw(np1:np2),cap(np1:np2))
            call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dc=cak,mpi_number=1, &
                    mpi_operation=ms_mpi_sum,mpi_comm=group_comm)
            cak(1)=csk/cak(1)
            anp(np1:np2)=anp(np1:np2)+cak(1)*cp(np1:np2)
            cr(np1:np2)=cr(np1:np2)-cak(1)*cap(np1:np2)
            cq(np1:np2)=cq(np1:np2)-conjg(cak(1))*caw(np1:np2)
            csk2(1)=dot_product(cq(np1:np2),cr(np1:np2))
            err=dot_product(cr(np1:np2),cr(np1:np2))
            erra(1)=err
            call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dc=csk2,mpi_number=1, &
                    mpi_operation=ms_mpi_sum,mpi_comm=group_comm)
            call ms_mpi(mpi_command='allreduce',mpi_recv_buf_dp=erra,mpi_number=1, &
                    mpi_operation=ms_mpi_sum,mpi_comm=group_comm)
            err=erra(1)
            err=err/enorm
            if(err.lt. eps) exit
            cbk=csk2(1)/csk
            cp(np1:np2)=cr(np1:np2)+cbk*cp(np1:np2)
            cw(np1:np2)=cq(np1:np2)+conjg(cbk)*cw(np1:np2)
            csk=csk2(1)
            do i=0,numprocs-1
               ip1=mpi_sphere_index(i)+1
               ip2=mpi_sphere_index(i)+mpi_sphere_number(i)
               np1=noff(ip1)+1
               np2=noff(ip2)+nblk(ip2)
               nsend=np2-np1+1
               call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=cp(np1:np2),mpi_number=nsend, &
                    mpi_rank=i,mpi_comm=group_comm)
               call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=cw(np1:np2),mpi_number=nsend, &
                    mpi_rank=i,mpi_comm=group_comm)
            enddo
            if(rank.eq.0.and.iter.eq.1.and.writetime.eq.0) then
               time2=mytime()-time1
               call timewrite(iunit,' time per iteration:',time2)
               writetime=1
            endif
            if(rank.eq.0.and.iterwrite.eq.1) then
                write(iunit,'('' iter, err:'',i5,e13.5)') iter,err
                call flush(iunit)
            endif
         enddo
!
!  arrive here with a converged solution
!
         do i=0,numprocs-1
            ip1=mpi_sphere_index(i)+1
            ip2=mpi_sphere_index(i)+mpi_sphere_number(i)
            np1=noff(ip1)+1
            np2=noff(ip2)+nblk(ip2)
            nsend=np2-np1+1
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=anp(np1:np2),mpi_number=nsend, &
                    mpi_rank=i,mpi_comm=group_comm)
         enddo
         end subroutine cbicg

      end module solver
