!
!  generic mstm calling program.
!
!  this code is intended for user modification.   It does not employ input files, command line arguments, etc.
!
!  all parameters must be hardwired.    read the code for the basic ideas.
!
      program main
      use mpidefs
      use mpidata
      use intrinsics
      use spheredata
      use numconstants
      use specialfuncs
      use miecoefdata
      use translation
      use solver
      use scatprops
      use nearfield
      implicit none
      integer :: nsphere,neqns,nodrmax,nodrt,i,k,niter,istat,numtheta, &
                 nblkt,nodrg,m,n,p,l,q,mn,kl,m1,n1,l1,k1,q1,w,klm,mnm,ikm, &
                 fixedorrandom,numargs,calctmatrix,maxiter,nodrta(1),calcnf, &
                 calcamn,ip1,ip2,ma,na,nsend,nfplane,nfoutunit,nfoutdata, &
                 maxmbperproc,trackiterations,nonactive,normalizesm, &
                 storetranmat,calcsm,niterstep,fftranpresent
      integer, allocatable :: nodr(:),ntran(:),sphereblk(:),sphereoff(:)
      real (8) :: alphadeg,betadeg,alpha,beta,epsmie,epstran,epssoln, &
                  qexttot,qabstot,xv,scalefac,qscatot,asymparm, &
                  rireal,riimag,phideg,theta1d,theta2d,thetad,costheta,phi, &
                  sm(4,4),time1,time2,fc1,fc2,fc3,fc4,epstcon,qabslm,absrat, &
                  cbeam,gbfocus(3),maxerr,nfplanepos,nfplanevert(2,2), &
                  deltax,gammadeg,epspw,gamma,qexttotpar,qexttotper, &
                  qabstotpar,qabstotper,qscatotpar,qscatotper,cphi,sphi,s11, &
                  nfdistance
      real(8), allocatable :: xsp(:), rpos(:,:),qext(:,:),qabs(:,:), &
               qsca(:,:),smc(:,:,:),smt(:,:,:)
      complex(8) :: sa(4),chiralfactor
      complex(8), allocatable :: amnp(:,:),amnp0(:,:,:,:),ri(:,:), &
                  gmn(:),amnp1(:,:,:),amnp2(:,:,:)
      character*30 :: inputfile,spherefile,parmfile,outfile,tmatrixfile,&
                      amnfile,nfoutfile,runfile
      complex(8), allocatable :: pmnp0(:,:,:,:)
      integer :: ierr,rank,printinputdata,runprintunit,numprocs
!
! this main program was set up to perform a loop of T matrix calculations,
! all involving scaled sphere positions from the file 'ran200fvp5.pos'
!
!
! extra variable declarations go here
!
      integer :: case,numcases
      real(8) :: xsp0,rpos0(3,200),rotang,volfrac,volfrac0,volfrac1,posscale
      complex(8) :: ri0
!
! define the sphere properties and run parameters
!
      nsphere=200
      xsp0=3.d0
      ri0=(1.31d0,0.0d0)
      chiralfactor=0.d0
      volfrac0=0.1d0
      volfrac1=0.5d0

      allocate(xsp(nsphere),rpos(3,nsphere),nodr(nsphere),ntran(nsphere), &
               ri(2,nsphere),sphereblk(nsphere),sphereoff(nsphere+1))
      spherefile='ran200fvp5.pos'
      open(1,file=spherefile)
      do i=1,nsphere
         read(1,*) xsp(i),rpos0(:,i)
      enddo
      close(1)

      outfile='test.dat'
      runfile=' '
      xsp=xsp0
      xv=xsp0*dble(nsphere)**.3333d0
      ri=ri0
      epsmie=1.d-4
      epstran=1.d-6
      epssoln=1.d-10
      niter=2000
      fixedorrandom=1
      theta1d=0.
      theta2d=180.
      numtheta=181
      maxmbperproc=1500.
      storetranmat=1
      nfdistance=10000.
      normalizesm=0
      calcamn=1
      amnfile='amn-temp.dat'
      phideg=0.d0
      alphadeg=0.d0
      betadeg=0.d0
      trackiterations=0
      calcnf=0
      cbeam=0.d0
      nfplane=1
      nfplanepos=0.d0
      deltax=0.2d0
      nfplanevert=reshape((/-40.d0,-40.d0,40.d0,40.d0/),(/2,2/))
      gammadeg=0.d0
      nfoutfile='nftest2b.dat'
      nfoutunit=2
      nfoutdata=1
      epspw=0.01d0
      calctmatrix=1
      tmatrixfile='tmatrix-temp.dat'
      epstcon=1.d-6
!
! this erases the output file: the code runs a loop where
! output values are appended to the file
!
       open(1,file=outfile)
       close(1,status='delete')
       if(runfile.ne.' ') then
         runprintunit=4
         open(runprintunit,file=runfile)
       else
         runprintunit=6
       endif
       call setrunparameters(run_print_unit=runprintunit)
!
! initialize mpi
!
      call ms_mpi(mpi_command='init')
!
! this is the main variable loop.    In this example the volume fraction of the spheres
! is changed with each case.
!
      numcases=20
      do case=1,numcases-1
      volfrac=volfrac0+(volfrac1-volfrac0)*dble(case-1)/dble(numcases-1)
      posscale=(volfrac1/volfrac)**.33333d0
      rpos=rpos0*posscale*xsp0
      write(runprintunit,'('' case, fv:'',i5,f8.2)') case, volfrac
!
! the rest of the code is basically the same as the mstm
! main program, without the getrunparameter calls.
!
      if(numtheta.gt.0) then
         if(allocated(smt)) deallocate(smt)
         allocate(smt(4,4,numtheta))
      endif
!
!  determine if optical activity is present
!
      nonactive=1
      do i=1,nsphere
         if(cdabs(ri(1,i)-ri(2,i)).gt.1.d-10) then
            nonactive=0
            exit
         endif
      enddo
!
!  calculation of sphere mie coefficients, order limits
!
      call miecoefcalc(nsphere,xsp,ri,epsmie)
      call getmiedata(sphere_order=nodr,max_order=nodrmax,number_equations=neqns, &
           sphere_block=sphereblk,sphere_block_offset=sphereoff)
!
!  determine the size of the parallel run and set it up
!
      call ms_mpi(mpi_command='size',mpi_size=numprocs)
      call ms_mpi(mpi_command='rank',mpi_rank=rank)
      call ms_mpi(mpi_command='barrier')
      call mpisetup(nsphere,nodr,rpos,fixedorrandom,maxmbperproc,storetranmat, &
           nfdistance,fftranpresent,runprintunit)
      call ms_mpi(mpi_command='barrier')
!
! this was an option for moving the GB focal point
!
      gbfocus=(/0.d0,dble(case-1),0.d0/)*.5d0
      gbfocus=0.d0
!
! translation matrix calculation
!
      call mpirottranmtrxsetup(nsphere,nodr,rpos,(1.d0,0.d0),storetranmat, &
           nfdistance,runprintunit)
      call ms_mpi(mpi_command='barrier')
!
!  determine orders required to expand scattered fields about target origin
!
      call tranorders(nsphere,nodr,rpos,epstran,ntran,nodrt)
!
!  report the size of the run
!
      if(rank.eq.0) then
         write(runprintunit,'('' maximum sphere order:'',i5)') nodrmax
         write(runprintunit,'('' estimated T matrix order:'',i5)') nodrt
         write(runprintunit,'('' number of equations:'',i9)') neqns
         call flush(runprintunit)
      endif
!
! the main calculations
!
      if(fixedorrandom.eq.1) then
!
!  random orientation option
!
         if(allocated(qext)) deallocate(qext,qabs,qsca)
         allocate(qext(nsphere,1), qabs(nsphere,1), qsca(nsphere,1))
         if(calctmatrix.ge.1) then
!
!  this option calculates the T matrix either from the beginning or where left off
!
            if(rank.eq.0) time1=mytime()
            call tmatrixsoln(neqns,nsphere,nodr,nodrt,xsp,rpos,epssoln,epstcon,niter,&
                 calctmatrix,tmatrixfile,fftranpresent,niterstep,qext,qabs,qsca,istat)
            if(rank.eq.0) then
               time2=mytime()-time1
               call timewrite(runprintunit,' execution time:',time2)
            endif
            call rottranmtrxclear()
         else
!
!  and this has the T matrix already calculated and stored in the file.
!
!  read the order of the T matrix and broadcast to the processors.
!
            if(rank.eq.0) then
               open(3,file=tmatrixfile)
               read(3,*) nodrt
               close(3)
               write(runprintunit,'('' t matrix order:'',i5)') nodrt
               call flush(runprintunit)
            endif
            nodrta(1)=nodrt
            call ms_mpi(mpi_command='bcast',mpi_send_buf_i=nodrta,mpi_number=1,mpi_rank=0)
            nodrt=nodrta(1)
            call ms_mpi(mpi_command='barrier')
         endif
!
!  the T matrix is available; calculate the random orientation scattering matrix
!
         nblkt=nodrt*(nodrt+2)
         nodrg=nodrt*2
         if(allocated(smc)) deallocate(smc)
         allocate(smc(4,4,0:nodrg))
         call ranorientscatmatrix(xv,nsphere,nodrt,nodrg,cbeam,tmatrixfile,smc,qext, &
              qabs,qsca)
         if(rank.eq.0) then
            qexttot=sum(qext(:,1)*xsp*xsp)/xv/xv
            qabstot=sum(qabs(:,1)*xsp*xsp)/xv/xv
            qscatot=qexttot-qabstot
            asymparm=dble(smc(1,1,1)/smc(1,1,0))/3.d0
            call ranorienscatmatrixcalc(numtheta,theta1d,theta2d,1,smc,nodrg,smt)
         endif
      else
!
!  fixed orientation option
!
         alpha=alphadeg*pi/180.d0
         beta=betadeg*pi/180.d0
         phi=phideg*pi/180.d0
         if(allocated(amnp)) deallocate(amnp)
         allocate(amnp(neqns,2))
         if(allocated(qext)) deallocate(qext,qabs,qsca)
         allocate(qext(nsphere,3), qabs(nsphere,3), qsca(nsphere,3))
         if(calcamn.eq.1) then
!
!  this option calculates the scattering coefficients
!
            if(rank.eq.0) time1=mytime()
            call fixedorsoln(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,epssoln,&
                epstran,niter,amnp,qext,qabs,qsca,maxerr,maxiter,trackiterations, &
                fftranpresent,niterstep,istat)
!
!  write the scattering coefficients to the file
!
            if(rank.eq.0) then
               time2=mytime()-time1
               write(runprintunit,'('' max iterations, soln error:'',i6,e13.5)') &
                          maxiter,maxerr
               call timewrite(runprintunit,' execution time:',time2)
               open(3,file=amnfile)
               do i=1,nsphere
                  write(3,'(6e13.5)') qext(i,:),qabs(i,:),qsca(i,:)
                  allocate(amnp1(0:nodr(i)+1,nodr(i),2),amnp2(0:nodr(i)+1,nodr(i),2))
                  ip1=sphereoff(i)+1
                  ip2=sphereoff(i)+sphereblk(i)
                  amnp1=reshape(amnp(ip1:ip2,1),(/nodr(i)+2,nodr(i),2/))
                  amnp2=reshape(amnp(ip1:ip2,2),(/nodr(i)+2,nodr(i),2/))
                  do n=1,nodr(i)
                     do m=-n,n
                        if(m.le.-1) then
                           ma=n+1
                           na=-m
                        else
                           ma=m
                           na=n
                        endif
                        write(3,'(4e17.9)') amnp1(ma,na,1),amnp2(ma,na,1)
                        write(3,'(4e17.9)') amnp1(ma,na,2),amnp2(ma,na,2)
                     enddo
                  enddo
                  deallocate(amnp1,amnp2)
               enddo
               close(3)
            endif
         else
!
!  this option reads the scattering coefficients from the file
!
            if(rank.eq.0) then
               open(3,file=amnfile)
               do i=1,nsphere
                  read(3,'(6e13.5)') qext(i,:),qabs(i,:),qsca(i,:)
                  allocate(amnp1(0:nodr(i)+1,nodr(i),2),amnp2(0:nodr(i)+1,nodr(i),2))
                  do n=1,nodr(i)
                     do m=-n,n
                        if(m.le.-1) then
                           ma=n+1
                           na=-m
                        else
                           ma=m
                           na=n
                        endif
                        read(3,'(4e17.9)') amnp1(ma,na,1),amnp2(ma,na,1)
                        read(3,'(4e17.9)') amnp1(ma,na,2),amnp2(ma,na,2)
                     enddo
                  enddo
                  ip1=sphereoff(i)+1
                  ip2=sphereoff(i)+sphereblk(i)
                  amnp(ip1:ip2,1)=reshape(amnp1(0:nodr(i)+1,1:nodr(i),1:2),(/sphereblk(i)/))
                  amnp(ip1:ip2,2)=reshape(amnp2(0:nodr(i)+1,1:nodr(i),1:2),(/sphereblk(i)/))
                  deallocate(amnp1,amnp2)
               enddo
               close(3)
            endif
!
!  broadcast the scattering coefficients to the other processors
!
            nsend=neqns*2
            call ms_mpi(mpi_command='bcast',mpi_send_buf_dc=amnp,mpi_number=nsend,mpi_rank=0)
         endif
!
!  calculate the efficiency factors
!
         cphi=cos(phi)
         sphi=sin(phi)
         qexttotpar=sum((qext(:,1)*cphi*cphi+2.d0*qext(:,3)*cphi*sphi+qext(:,2)*sphi*sphi) &
                    *xsp*xsp)/xv/xv
         qexttotper=sum((qext(:,1)*sphi*sphi-2.d0*qext(:,3)*cphi*sphi+qext(:,2)*cphi*cphi) &
                    *xsp*xsp)/xv/xv
         qabstotpar=sum((qabs(:,1)*cphi*cphi+2.d0*qabs(:,3)*cphi*sphi+qabs(:,2)*sphi*sphi) &
                    *xsp*xsp)/xv/xv
         qabstotper=sum((qabs(:,1)*sphi*sphi-2.d0*qabs(:,3)*cphi*sphi+qabs(:,2)*cphi*cphi) &
                    *xsp*xsp)/xv/xv
         qscatotpar=qexttotpar-qabstotpar
         qscatotper=qexttotper-qabstotper
         qexttot=(qexttotpar+qexttotper)*.5d0
         qabstot=(qabstotpar+qabstotper)*.5d0
         qscatot=(qscatotpar+qscatotper)*.5d0
         qext(:,1)=(qext(:,1)+qext(:,2))*.5d0
         qabs(:,1)=(qabs(:,1)+qabs(:,2))*.5d0
         qsca(:,1)=(qsca(:,1)+qsca(:,2))*.5d0
         call rottranmtrxclear()
!
!  calculate the target-based expansion and rotate to the incident field frame
!
         allocate(amnp0(0:nodrt+1,nodrt,2,2),pmnp0(0:nodrt+1,nodrt,2,2))
         do k=1,2
            call amncommonorigin(neqns,nsphere,nodr,ntran,nodrt,rpos, &
                 amnp(1:neqns,k),amnp0(0:,1:,1:,k))
            call rotvec(alpha,beta,0.d0,nodrt,nodrt,amnp0(0:,1:,1:,k),1)
         enddo
!
!  calculate the asymmetry parameter and the scattering matrix
!
         allocate(gmn(0:2))
         call s11expansion(amnp0,nodrt,0,1,gmn)
         asymparm=dble(gmn(1)/gmn(0))/3.d0
         do i=1,numtheta
            thetad=theta1d+(theta2d-theta1d)*(i-1)/max(1.d0,dble(numtheta-1))
            costheta=cos(thetad*pi/180.d0)
            call scatteringmatrix(amnp0,nodrt,xv,costheta,phi,sa,smt(:,:,i))
         enddo
         deallocate(amnp0,pmnp0,gmn)
      endif
!
!  output file operations
!
      if(rank.eq.0) then
         open(1,file=outfile,position='append')
         if(nonactive.eq.0) then
            write(1,'('' sphere S.P., pos. (x,y,z), ref. index (L,R), Qext, Qsca, Qabs, Qabs/Qabs,LM'')')
         else
            write(1,'('' sphere S.P., pos. (x,y,z), ref. index, Qext, Qsca, Qabs, Qabs/Qabs,LM'')')
         endif
         do i=1,nsphere
            call getmiedata(which_sphere=i,sphere_qabs=qabslm)
            if(dimag(ri(1,i)).eq.0.d0.and.dimag(ri(2,i)).eq.0.d0) then
               absrat=1.d0
            else
               absrat=qabs(i,1)/qabslm
            endif
            if(nonactive.eq.0) then
               write(1,'(i5,4f10.4,4f10.6,3e13.5,f8.4)') i, xsp(i),rpos(:,i)+gbfocus, ri(:,i), &
                  qext(i,1),qsca(i,1),qabs(i,1),absrat
            else
               write(1,'(i5,4f10.4,2f10.6,3e13.5,f8.4)') i, xsp(i),rpos(:,i)+gbfocus, ri(1,i), &
                  qext(i,1),qsca(i,1),qabs(i,1),absrat
            endif
         enddo
         if(fixedorrandom.eq.1) then
            write(1,'('' total ext, abs, scat efficiencies, w.r.t. xv, and asym. parm'')')
            write(1,'(6e13.5)') qexttot,qabstot,qscatot,asymparm
         else
            write(1,'('' unpolarized total ext, abs, scat efficiencies, w.r.t. xv, and asym. parm'')')
            write(1,'(6e13.5)') qexttot,qabstot,qscatot,asymparm
            write(1,'('' parallel total ext, abs, scat efficiencies'')')
            write(1,'(6e13.5)') qexttotpar,qabstotpar,qscatotpar
            write(1,'('' perpendicular total ext, abs, scat efficiencies'')')
            write(1,'(6e13.5)') qexttotper,qabstotper,qscatotper
         endif

         write(1,'('' scattering matrix elements'')')
         if(normalizesm.eq.0) then
            write(1,'('' theta    s11         s22         s33'',&
                   &''         s44'',&
                   &''         s21         s32         s43         s31'',&
                   &''         s42         s41'')')
            do i=1,numtheta
               thetad=theta1d+(theta2d-theta1d)*(i-1)/max(1.d0,dble(numtheta-1))
               write(1,'(f8.2,10e12.4)') thetad,smt(1,1,i),smt(2,2,i),smt(3,3,i), &
                     smt(4,4,i),smt(1,2,i),smt(2,3,i),smt(3,4,i),smt(1,3,i), &
                     smt(2,4,i),smt(1,4,i)
            enddo
         else
            write(1,'('' theta    s11         s22/s11     s33'',&
                   &''/s11     s44'',&
                   &''/s11     s21/s11     s32/s11     s43/s11     s31'',&
                   &''/s11     s42/s11     s41/s11'')')
            do i=1,numtheta
               thetad=theta1d+(theta2d-theta1d)*(i-1)/max(1.d0,dble(numtheta-1))
               s11=smt(1,1,i)
               write(1,'(f8.2,10e12.4)') thetad,smt(1,1,i),smt(2,2,i)/s11,smt(3,3,i)/s11, &
                     smt(4,4,i)/s11,smt(1,2,i)/s11,smt(2,3,i)/s11,smt(3,4,i)/s11,smt(1,3,i)/s11, &
                     smt(2,4,i)/s11,smt(1,4,i)/s11
            enddo
         endif
         if(fixedorrandom.eq.1) then
            write(1,'('' scattering matrix expansion coefficients'')')
            write(1,'(''    w  a11         a22         a33         '',&
             &''a23         a32         a44         a12         '',&
             &''a34         a13         a24         a14'')')
            do w=0,nodrg
               write(1,'(i5,11e12.4)') w,smc(1,1,w),smc(2,2,w),&
                 smc(3,3,w),smc(2,3,w),smc(3,2,w),smc(4,4,w),&
                 smc(1,2,w),smc(3,4,w),smc(1,3,w),smc(2,4,w),&
                 smc(1,4,w)
            enddo
         endif
         close(1)
      endif
!
!  near field calculation options
!
      if(fixedorrandom.eq.0.and.calcnf.eq.1) then
!
! this was a modification to the main: the near field file
! is opened for appending
!
         if(rank.eq.0) then
            open(nfoutunit,file=nfoutfile,position='append')
         endif
         gamma=gammadeg*pi/180.d0
         call nearfieldgridcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri,amnp, &
                    nfplane,nfplanepos,nfplanevert,gbfocus,deltax,gamma,nfoutunit,epspw, &
                    nfoutdata,runprintunit)
         if(rank.eq.0) then
            close(nfoutunit)
         endif
      endif
!
!  all done!
!
!
! and this ends the main variable loop.
!
      enddo
!
! time to go home
!
      call ms_mpi(mpi_command='finalize')
      end
