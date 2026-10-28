c
c same as differ 2 but no implicit declarations
      implicit none

c
c
c functions
      real*8 dot33
c
      real*8 dtr, rtd, conv
C arrays :
      real*8 v(8),dv(8),vw(8),grad(8),sig(0:8),kin(1:5,0:8)
      real*8 kgrad(5,8)
c 4-vectors
      real*8 k_i(4), k_f(4), p_f(4)
c 3-vectors
      real*8 ue_i(3), ue_f(3), up_f(3), uq(3), n_e(3), n_p(3)
      real*8 up_r(3), pv_r(3), qv(3), pv_f(3)

c local variables
      real*8 hbarc, pi
      real*8 ph_b, th_b, ph_e, th_e, ph_p, th_p
      real*8 phi, phi_0, phi_b0, phi_e0, phi_p0
      real*8 theta, theta_0, the_b0, the_e0, the_p0, the_q
      real*8 theta_prq, thee, thnp
      real*8 cph_b, cph_e, cph_p, cphi
      real*8 cth_b, cth_e, cth_p, cthe_q, cthee, ctheta, cthnp
      real*8 sph_b, sph_e, sph_p, sphi
      real*8 sth_b, sth_e, sth_p, sthe_q, sthee, stheta, sthnp
      real*8 de, def, dphi_b, dphi_e, dphi_p
      real*8 dth_b, dth_e, dth_p
      real*8 e0, ef0, ef, ei, enp, ep, er
      
      real*8 omega, omega0, q, qcm, rec
      real*8 pi_mag, pipf, ppi, ppp, pr_mag
      real*8 w

      real*8 rsa, sigma
      integer i1, i2, indep, iunit
      integer i, j, k, n, ib
c real variables
      real*8 k_fact,mn,mp
      real rho_paris
c
      character infil*80,comnt*80,strfil*80,outfil*80
c 
      logical ifread,lcros,lcc1,laren,lerr,lperc
c


c common blocks for thephi
      common/electn/th_e,ph_e,cth_e,cph_e,sth_e,sph_e,w
      common/proton/th_p,ph_p,cth_p,cph_p,sth_p,sph_p
      common/beam/  th_b,ph_b,cth_b,cph_b,sth_b,sph_b,ei
c
c
c constants
      data pi/3.141592654/
      data hbarc/197.33/
      data mn/939.6/
      data mp/938.3/
      data rsa/5.e-4/
      data lcros/.false./,laren/.false./,lcc1/.false./
c
c conversion factors
c
      dtr=pi/180
      rtd=180/pi
c open files
c      call chread('input file ?',infil,i1,i2)
      open (15,file='differ.in',status='old')
c      call chread('output file ?',outfil,i1,i2)
      open (16,file='differ.out',status='unknown')
c
c calculation of cc1 cross sections not yet included
c
c            lcc1=ifread('use PWIA ? ','Yy')
      lcc1 = .true.
      laren = .false.
      if (.not.lcc1) then
         laren=.true.
         call chread ('structure function file ?',strfil,i1,i2)
         open (18,file=strfil(i1:i2),status='old')
      endif
c     
c      lperc=ifread('all variations in % ?','yY')
      lperc = .true.
c
c read the input file
c
c electron side
100   read (15,*,end=999) ib
      read (15,'(a)',end=999) comnt
      read (15,*,end=999) the_e0,phi_e0 ,ef0 
      read (15,*,end=999) dth_e,dphi_e,def

      the_e0=the_e0*dtr
      phi_e0=phi_e0*dtr
      dth_e=dth_e/1000.
      dphi_e=dphi_e/1000.

c proton side
      read (15,'(a)',end=999) comnt
      read (15,*,end=999) the_p0,phi_p0
      read (15,*,end=999) dth_p,dphi_p

      the_p0=the_p0*dtr
      phi_p0=phi_p0*dtr
      dth_p=dth_p/1000.
      dphi_p=dphi_p/1000.

c beam side
      read (15,'(a)',end=999) comnt
      read (15,*,end=999) the_b0,phi_b0,e0
      read (15,*,end=999) dth_b,dphi_b,de

      the_b0=the_b0*dtr
      phi_b0=phi_b0*dtr
      dth_b=dth_b/1000.
      dphi_b=dphi_b/1000.

c indep defines whether omega is independent of E
c indep=0 : not independent
c indep=1 : independent variation of omega and E
      read (15,'(a)',end=999) comnt
      read (15,*,end=999) indep
c end of input file

c
c fill the array's v and dv (vectors with the central values and step 
c sizes)
       v(1)=the_e0
      dv(1)=dth_e
       v(2)=phi_e0
      dv(2)=dphi_e
       v(3)=ef0
      dv(3)=def
c
       v(4)=the_p0
      dv(4)=dth_p
       v(5)=phi_p0
      dv(5)=dphi_p
c
       v(6)=the_b0
      dv(6)=dth_b
       v(7)=phi_b0
      dv(7)=dphi_b
       v(8)=e0
      dv(8)=de
c      
c calculate the different cross sections and put them into array sig
c
      do 1000 j=0,8
c
c put the current angle and energy setting into the common block
c
         if (j.eq.0) then
            do 1100 k=1,8
1100           vw(k)=v(k)
         else
            do 1110 k=1,8
               if (k.eq.j) then
      	          vw(k)=v(k)+dv(k)
               else
      		  vw(k)=v(k)
      	       endif
1110        continue
         endif
c
c special case : incident energy and omega are not indepentdent
c without an energy loss system
c
	 if (j.eq.8.and.indep.eq.0) then
              vw(3)=vw(3)+dv(j)
         endif

         ei=vw(8)
         ef = vw(3)
c
         th_b=vw(6)
         ph_b=vw(7)
         cth_b=cos(th_b)
         cph_b=cos(ph_b)
         sth_b=sin(th_b)
         sph_b=sin(ph_b)
c
         w = ei - ef
c
         th_e=vw(1)
         ph_e=vw(2)
         cth_e=cos(th_e)
         cph_e=cos(ph_e)
         sth_e=sin(th_e)
         sph_e=sin(ph_e)
c
         th_p=vw(4)
         ph_p=vw(5)
         cth_p=cos(th_p)
         cph_p=cos(ph_p)
         sth_p=sin(th_p)
         sph_p=sin(ph_p)
c Use cartesian coordinates : z-direction in beam 0-direction
c					y-direction normal to hall floor
c					x-direction normal to yz plane 
c form unit vectors
c
c incident beam
         ue_i(1)=ei*cth_b*sph_b
         ue_i(2)=-ei*sth_b
         ue_i(3)=ei*cth_b*cph_b

c scattered electron
         ue_f(1)=ef*cph_e*sth_e
         ue_f(2)=ef*sph_e
         ue_f(3)=ef*cph_e*cth_e

c outgoing proton
         up_f(1)=-cph_p*sth_p
         up_f(2)=sph_p
         up_f(3)=cph_p*cth_p
c
c momentum transfer
c
         qv(1) = ue_i(1) - ue_f(1)
         qv(2) = ue_i(2) - ue_f(2)
         qv(3) = ue_i(3) - ue_f(3)

         uq(1) = qv(1)
         uq(2) = qv(2)
         uq(3) = qv(3)

c momentum transfer
         q = sqrt(dot33(uq,uq))
c
c create unit vectors
c
         call normalize (ue_i)
         call normalize (ue_f)
         call normalize (up_f)
         call normalize (uq)

c electron scattering angle :
         cthee = dot33(ue_i,ue_f)
         if (cthee.gt.1) then
            thee = 0
         elseif (cthee.lt.-1) then
            thee = pi
         else
            thee = acos(cthee)
         endif

c theta q
         cthe_q = dot33(uq,ue_i)
         if (cthe_q.gt.1) then
            the_q = 0
         elseif (cthe_q.lt.-1) then
            the_q = pi
         else
            the_q = acos(cthe_q)
         endif

c theta
         ctheta = dot33(uq,up_f)
         if (ctheta.gt.1) then
            theta = 0
         elseif (ctheta.lt.-1) then
            theta = pi
         else
            theta = acos(ctheta)
         endif


c electron plane
         call cross3(ue_i,ue_f,n_e)

c proton plane
         call cross3(uq,up_f,n_p)

c norms    
         call normalize (n_e)
         call normalize (n_p)

c phi
         cphi = dot33 (n_e, n_p)
         if (cphi.gt.1) then
            phi = 0
         elseif (cphi.lt.-1) then
            phi = pi
         else
            phi = acos(cphi)
         endif
c
c call the subroutine
c

600   format (1x,a,1pe10.3,5x,a,1pe10.3)
c          
c
c cross section calculation :
c
c transformation to the center off mass system
c
c

         lerr=.false.
         call relkin(w,q,ctheta,conv,cthnp,enp,qcm,ppp,ppi,lerr)
c
c check for unphysical kinematics
c
         if (lerr) then
            write (6,*) ' unphysical kinematics '
            stop
         endif
         if (dabs(cthnp).gt.1.)cthnp= dsign(1.d0,cthnp)
         thnp=acos(cthnp)

c
c final proton momentum vector
c
         pv_f(1) = ppp*up_f(1)
         pv_f(2) = ppp*up_f(2)
         pv_f(3) = ppp*up_f(3)
c
c recoil momentum vector
c
         pv_r(1) = qv(1) - pv_f(1)
         pv_r(2) = qv(2) - pv_f(2)
         pv_r(3) = qv(3) - pv_f(3)
         
         up_r(1) = pv_r(1)
         up_r(2) = pv_r(2)
         up_r(3) = pv_r(3)

c
c central value for check only
c         
         if (j .eq. 0) then
            pr_mag=sqrt(dot33(pv_r,pv_r))
         endif

         call normalize (up_r)

c theta_prq
         ctheta = dot33(uq,up_r)
         if (ctheta.gt.1) then
            theta_prq = 0
         elseif (ctheta.lt.-1) then
            theta_prq = pi
         else
            theta_prq = acos(ctheta)
         endif



c 
c store selected kinematical variables into the array kin
c
         kin(1,j)=ppi
         kin(2,j)=ppp
         kin(3,j)=q
         kin(4,j)=the_q/dtr
c         kin(5,j)=enp
         kin(5,j)=theta_prq/dtr
c
c store theta_cm and phi into varialbles theta_0, phi_0
c
	   if (j.eq.0) then
		theta_0 = thnp/dtr
	      phi_0 = phi/dtr
	   endif
c
c calculate the d(ee'p)n crossection according to arenhoevel
c
         if (laren) then
            call crosec(ei,w,enp,thnp,phi,thee,sigma)
            sigma=sigma*conv
         endif
c
c calculate the d(ee'p)n cross section in PWIA using the CC1 cross sect.
c
         if (lcc1) then
c
c Use cartesian coordinates : z-direction in beam 0-direction
c					y-direction normal to hall floor
c					x-direction normal to yz plane 
c
c form the 4-vectors k_i, k_f, p_f
c
c incident beam
      	    k_i(1)=ei*cth_b*sph_b
            k_i(2)=-ei*sth_b
            k_i(3)=ei*cth_b*cph_b
            k_i(4)=ei
c scattered electron
            k_f(1)=ef*cph_e*sth_e
            k_f(2)=ef*sph_e
            k_f(3)=ef*cph_e*cth_e
            k_f(4)=ef
c outgoing proton
            p_f(1)=-ppp*cph_p*sth_p
            p_f(2)=ppp*sph_p
            p_f(3)=ppp*cph_p*cth_p
            p_f(4)=sqrt(ppp**2+mp**2)
c
            call subcc1(k_i,k_f,p_f,sigma,k_fact,pi_mag)
c
c recoil factor
c
            er=sqrt(pi_mag**2+mn**2)
            ep=sqrt(ppp**2+mp**2)
            pipf=(q**2-(ppp**2+pi_mag**2))/2.
            rec=1.-ep/er*pipf/ppp**2
c
	    sigma=sigma*k_fact*dble(rho_paris(sngl(ppi)))/rec
c
c check for a negative cross section
c
            if (sigma.lt.0) then
               write (6,*) ' unphysical cross section', sigma, k_fact,
     >       rec,  rho_paris(sngl(ppi)), ppi
               
               stop
            endif
c
	 endif
c
c fill array sig
         sig(j)=sigma

c
c calculate derivatives
c
         if (j.gt.0) then

            grad(j)=(sig(j)-sig(0))/dv(j)/sig(0)*100

            do 1200 n=1,5
      	       if (lperc) then
                  kgrad(n,j)=(kin(n,j)-kin(n,0))/dv(j)/kin(n,0)*100
       	       else
                  kgrad(n,j)=(kin(n,j)-kin(n,0))/dv(j)
               endif
1200        continue
         endif
c
1000   continue
c
c
      iunit=6
6000    write (iunit,'(/a)')
     &  '-------------------------------------------------------------'
        write (iunit,605) 'central cross section : ',sig(0)
        write (iunit,605) 'incident energy : ', e0
        	write (iunit,630) 'scattering angle :',thee/dtr,
     &                  'omega :',e0-ef0
        	write (iunit,630) 'initial momentum :',ppi,
     &                  'final momentum :',ppp
        	write (iunit,630) 'Theta CMS :',theta_0,
     &                  'Phi CMS :',phi_0
        	write (iunit,630) 'q_lab :',q,
     &                  'theta_q :',the_q/dtr
        	write (iunit,630) 'pr :',pr_mag,
     &                  'theta_prq :',kin(5,0)
        write (iunit,*) 'ei: ',ei,'w :',w,'enp :',enp,'thnp :',thnp,
     &                  'phi :',phi,'thee :',thee
      if (lperc) then
         write (iunit,'(/1x,a)') 'all derivatives are in %'
      else
         write (iunit,'(/1x,a)') 'cross section variation is in %'
      endif
      write (iunit,'(1x,a/)') 'all angle derives  are in 1/degree '
      
c
      write (iunit,620)
620   format(7x,4x,'dsig',8x,'dpi',9x,'dpf',10x,'dq',9x,
     &       'dth_q',8x,'dth_prq'/)

      write (iunit,*) 
c
      write (iunit,610) 'dth_e  ',grad(1)*dtr,(kgrad(n,1)*dtr,n=1,5)
      write (iunit,610) 'dph_e  ',grad(2)*dtr,(kgrad(n,2)*dtr,n=1,5)
      write (iunit,610) 'def    ',grad(3),(kgrad(n,3),n=1,5)
c
      write (iunit,*) '  '
c
      write (iunit,610) 'dth_p  ',grad(4)*dtr,(kgrad(n,4)*dtr,n=1,5)
      write (iunit,610) 'dph_p  ',grad(5)*dtr,(kgrad(n,5)*dtr,n=1,5)
c
      write (iunit,*) '  '
c
      write (iunit,610) 'dth_b  ',grad(6)*dtr,(kgrad(n,6)*dtr,n=1,5)
      write (iunit,610) 'dph_b  ',grad(7)*dtr,(kgrad(n,7)*dtr,n=1,5)
      write (iunit,610) 'dE     ',grad(8),(kgrad(n,8),n=1,5)
      if (iunit.eq.6) then
         iunit=16
         goto 6000
      endif
c
605   format(a,6(1pe10.3,2x))
610   format(a,6(1pe10.3,2x))
630   format(2(a,2x,1pe10.3,2x))
c
	goto 100
c
999   close (15)
      close (16)
      stop
      end
