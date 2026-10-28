c
c no impliciti statments here
c
c these are a set of routines used by the programs PHASPA, DIFFER,
c EVNTCRE and EVNTREAD.
c
      subroutine relate (e,qvec,omega,cthe_b,
     &                   the_e,cthe_e,sthe_e)
c
c
      implicit none
c input
      real*8 e, qvec, omega, cthe_b
c output
      real*8 the_e, cthe_e, sthe_e
c local
      real*8 arg, thee0, cthee0

      arg=sqrt((qvec**2-omega**2)/(4.*e*(e-omega)))
      if (dabs(arg).gt.1.) arg= dsign(1.d0,arg)
      thee0=2.*asin(arg)
      cthee0=cos(thee0)
      cthe_e=cthee0/cthe_b
      if (dabs(cthe_e).gt.1.) cthe_e= dsign(1.d0,cthe_e)
      the_e=acos(cthe_e)
      sthe_e=sin(the_e)
      return
      end
c
c
       real*8 function ang(sint,cost)

       implicit none
c
c input
       real*8 sint, cost
c local
       real*8 pi, arg1
       
       logical q1,q2,q3,q4
       data pi/3.141592654/
c
       if (dabs(cost).gt.1.) cost=dsign(1.d0,cost)
       if (dabs(sint).gt.1.) sint=dsign(1.d0,sint)
       arg1=acos(cost)
c
c find the quadrant
c
       q1=sint.ge.0.and.cost.ge.0
       q2=sint.ge.0.and.cost.le.0
       q3=sint.le.0.and.cost.le.0
       q4=sint.le.0.and.cost.ge.0
c
       if (q1.or.q2) ang=arg1
       if (q3.or.q4) ang=2.*pi-arg1 
       return
       end
c
c
      subroutine relkin(omega,qlab,cth, 
     &                 conv,cthnpc,enp,qc,ppp,ppi,lerr)
c      implicit double precision (a-h,o-z)
       implicit none
c
c input
       real*8 omega, qlab, cth
c output
       real*8 cthnpc,enp,qc,ppp,ppi, conv
       logical lerr
c local
       real*8 md,mn,mp,ml,kl2sq 
       real*8 pi, amu, dtr, rtd
       real*8 ql, e0, a0, a, b, c, discr
       real*8 beta, cthnpl, fak1, fak2, fact2, el, et
       real*8 gamma, pcm, ppcm, wc
c 
c calculates the relevant kinematic variables in the center of 
c mass system and the final and initial proton momenta in the 
c laboratory frame. 
c 
c the input varable cth = cos (theta) where theta is the angle 
c between q and the out coming proton. 
c cthnpc is the cosine of the theta-np angle in the cms frame 
c 

      data pi/3.141592654/ 
      data amu/931.48/ 
      dtr=pi/180. 
      rtd=180./pi 
      md=2.014*amu 
      mn=1.0087*amu 
      mp=1.0073*amu 
c 
c calcualate final momentum 
c 
      ql=qlab
      e0=md+omega 
      a0=ql**2-e0**2-mp**2+mn**2 
      a=4*e0**2-4*ql**2*cth**2 
      b=4*ql*a0*cth 
      c=4*e0**2*mp**2-a0**2 
      discr=b**2-4*a*c
      if (discr.lt.0) then
         lerr=.true.
         return
      endif
      ppp=(-b+sqrt(b**2-4*a*c))/2/a 
c 
c calculate theta-np lab and the initial momentum 
c 
      fak1=-0.5*(ql-2*ppp*cth)
      fak2=sqrt(ppp**2+(ql/2)**2-ppp*ql*cth)
      cthnpl=fak1/fak2
      kl2sq=ppp**2+(ql/2)**2-ppp*ql*cth
      ppi=sqrt(kl2sq+(ql/2)**2-sqrt(kl2sq)*ql*cthnpl)
c
c make a Lorentz tranformation to the CMS frame
c
      beta=ql/(md+omega)
      gamma= 1/sqrt(1-beta**2)
c
      qc=gamma*(ql-beta*omega)
      wc=gamma*(omega-beta*ql)
c
      cthnpc=md*cthnpl/sqrt(md**2+qc**2-qc**2*cthnpl**2)
c 
c calculate invariant mass
c
      et=sqrt((omega+md)**2-ql**2)
      enp=et-(mn+mp)
c
c     cross section transformation factor cms->lab
c
      pcm=sqrt(enp*(enp+2.*(mn+mp)))
      ppcm= pcm/2.
      ml=enp+mn+mp
      el=sqrt(ml**2+ql**2)
      fact2=(ppp/ppcm)**3*ml/el
      conv=fact2/(1.+ql*ml*cthnpc/(2*ppcm*el))
c
      return
      end
c
c
      subroutine crosec(einc,omega,enp,thcm,phi,thee,sigma) 
  
      implicit none
c input
      real*8 einc, omega, enp, thcm, phi, thee
c output
      real*8 sigma
c local
      real*8 pi, amu, hbarc, alpha
      real*8 qnue2, qlab2, xi, eta
      real*8 k1, k2
      real*8 md,mn,mp,ml,kl2sq 
      real*8 rho00, rho11, rho01, rho1_1
      real*8 r, c, dtr, rtd
      real*8 f00, f11, f01, f1_1

      data pi/3.141592654/,hbarc/197.33/
      data amu/931.48/ 
      dtr=pi/180. 
      rtd=180./pi 
      md=2.014*amu 
      mn=1.0087*amu 
      mp=1.0073*amu 
      alpha=1./137. 
c 
      k1=einc 
      k2=(einc-omega) 
c 
      qnue2=4*k1*k2*(sin(thee/2))**2 
      qlab2=k1**2+k2**2-2*k1*k2*cos(thee) 
c 
      xi=qnue2/qlab2 
      eta=(tan(thee/2))**2 
c 
      rho00=qnue2*xi**2/(2*eta) 
      rho11=0.5*qnue2*(1+xi/(2*eta)) 
      rho01=qnue2*xi/eta*sqrt((xi+eta)/8) 
      rho1_1=-qnue2*xi/(4*eta) 
c 
      r=(1+(k2-k1*cos(thee))/(omega+md))**(-1.) 
      c=alpha/(6*pi**2)*r/qnue2**2*k2/k1
c
c get the structure functions
c
       call strfun(enp,thcm,f00,f11,f01,f1_1)
c
c cross section
c
      sigma=c*(rho00*f00+rho11*f11+
     &         rho01*f01*cos(phi)+rho1_1*f1_1*cos(2*phi))*hbarc
c
c convert fm2 to nb
c
      sigma=sigma*(1.e+7)
c
c
      return
      end
c
      subroutine strfun(enp,thcm,fl,ft,flt,ftt)  

c returns the 5 response functions for a selected value of E_np
c and interpolates the value of thcm
c
      implicit double precision (a-h,o-z)
c
c  
      dimension f00(12,100),f11(12,100),f01(12,100),f_11(12,100)
	dimension f01p(12,100) 

      dimension energy(12),th(100)  
c  
      logical first  
      data first/.true./  
      data pi/3.141592654/  
      dtr=pi/180.  
c
c read the structure function file
c  
      if (first) then
	   read(18,*) nset  
         do  j=1,nset  
            read (18,*) energy(j)  
            do  k=1,49  
               read (18,1500) th(k),f00(j,k),f11(j,k),f01(j,k),f_11(j,k)
     &				,f01p(j,k)  
               th(k)=th(k)*dtr  
1500           format (10x,f5.0,36x,5f11.0)  
	      enddo
	   enddo
         first=.false.  
      endif  
c  
c interpolation in angles only  
c
c find the closest energy
c  
      do i=1,nset-1  
         if (enp.ge.energy(i).and.enp.lt.energy(i+1)) then
c 
c found the interval
c
	      delta = enp - energy(i)
		enp_bin = energy(i+1)-energy(i)
c
c closest energy
c
		if (delta.le.(enp_bin/2.)) then
		   j1 = i
		else
	         j1 = i+1
		endif
         endif
	enddo
c
      do k=1,49 
         if (thcm.ge.th(k)) then
            k1=k
            k2=k+1
         endif
	enddo
c 
      fth=(thcm-th(k1))/(th(k2)-th(k1))  
c 
      fl=f00(j1,k1)+(f00(j1,k2)-f00(j1,k1))*fth 
c 
      flt=f01(j1,k1)+(f01(j1,k2)-f01(j1,k1))*fth 
c 
      ft=f11(j1,k1)+(f11(j1,k2)-f11(j1,k1))*fth 
c 
      ftt=f_11(j1,k1)+(f_11(j1,k2)-f_11(j1,k1))*fth 
c 
      fltp=f01p(j1,k1)+(f01p(j1,k2)-f01p(j1,k1))*fth 
c
      return 
      end 


c
c
      subroutine thephi(theta,ctheta,phi,thee,the_q,q)
c
c calculates the angles theta phi theta-q and the momentum transfer for 
c a given proton, beam and electron direction and for a given combination 
c of the incident energy and omega
c
c the arguments are transferred via common blocks
c theta and phi are the results
c
      implicit double precision(a-h,o-z)
c
c 
c common blocks for the arguments
c 
      logical first, lcros, center, lvax

      common/electn/the_e,phi_e,cthe_e,cphi_e,sthe_e,sphi_e,omega
      common /proton/the_p,phi_p,cthe_p,cphi_p,sthe_p,sphi_p
      common /beam/  the_b,phi_b,cthe_b,cphi_b,sthe_b,sphi_b,ei
c constants
      data pi/3.141592654/
      data hbarc/197.33/
      data rsa/5.e-4/
      data thmin/1.e+30/,thmax/-1.e+30/
      data phmin/1.e+30/,phmax/-1.e+30/
      data first/.true./,lcros/.false./,center/.true./,lvax/.false./
c
c conversion factors
c
      dtr=pi/180
      rtd=180/pi
c
c begin the calculation
c
         cosa=cthe_b*cphi_b
         cosb=cthe_p*cphi_p
         cosc=cthe_e*cphi_e
c
c
         a=acos(cosa)
         b=acos(cosb)
         c=acos(cosc)
c
         sina=sin(a)
         sinb=sin(b)
         sinc=sin(c)
c
         if (dabs(the_b).lt.rsa.and.dabs(phi_b).lt.rsa) then
            sinbet=0.
            cosbet=0.
         else
            sinbet=sthe_b/sina
            cosbet=sphi_b/sina
         endif
         if (dabs(the_e).lt.rsa.and.dabs(phi_e).lt.rsa) then
            sinalp=0.
         else
            sinalp=sphi_e/sinc
         endif
         if (dabs(phi_e).lt.rsa) then
            cosalp=1.
         elseif(dabs(phi_e).gt.rsa.and.the_e.lt.rsa) then
            cosalp=0.
         else
            cosalp=(cphi_e-cosc*cthe_e)/(sinc*sthe_e)
         endif
         alp=ang(sinalp,cosalp)
         bet=ang(sinbet,cosbet)
         gam=alp+bet
c
         cthee=cosc*cosa+sinc*sina*cos(gam)
         if (dabs(cthee).gt.1.)cthee= dsign(1.d0,cthee)
         thee=acos(cthee)
         sthee=sin(thee)
c
         if (dabs(the_p).lt.rsa.and.dabs(phi_p).lt.rsa) then
            sinep=0.
         else
            sinep=sphi_p/sinb
         endif
         if (dabs(phi_p).lt.rsa) then
            cosep=1.
         elseif (dabs(phi_p).gt.rsa.and.dabs(the_p).lt.rsa )then
            cosep=0.
         else
            cosep=(cphi_p-cosb*cthe_p)/(sinb*sthe_p)
         endif
         ep=ang(sinep,cosep)
         delt=pi-(alp+ep)
         cdelt=cos(delt)
         sdelt=sin(delt)
         cose=cosb*cosc+sinb*sinc*cdelt
         if (dabs(cose).gt.1.)cose= dsign(1.d0,cose)
         ae=acos(cose)
         sine=sin(ae)
c
         sinzet=sina*sin(gam)/sin(thee)
         if (dabs(sinzet).gt.1.)sinzet= dsign(1.d0,sinzet)
         if (dabs(a).lt.rsa) then
            coszet=1.
         elseif (dabs(c).lt.rsa) then
            coszet=1.
         else
            coszet=(cosa-cosc*cthee)/(sinc*sthee)
         endif
         zet=ang(sinzet,coszet)
         sineta=sinb*sdelt/sine
         if (dabs(sineta).gt.1.)sineta= dsign(1.d0,sineta)
         if (dabs(c).lt.rsa) then
            coseta=1.
         else
            coseta=(cosb-cose*cosc)/(sine*sinc)
         endif
         eta=ang(sineta,coseta)
         chi=zet+eta
c
c calculate theta-q
c
         q=sqrt(omega**2+
     &          2*ei*(ei-omega)*(1-cthee))
         cthe_q=1./q*(ei-(ei-omega)*cthee)
         if (dabs(cthe_q).gt.1.)cthe_q= dsign(1.d0,cthe_q)       

         the_q=acos(cthe_q)
c
c calculate the final angles
c
         ctheta=cose*cos(the_q+thee)+sine*sin(the_q+thee)*cos(chi)
         if (dabs(ctheta).gt.1.)ctheta= dsign(1.d0,ctheta)
         theta=acos(ctheta)
         stheta=sin(theta)
c
         if (dabs(theta).lt.rsa) then
            cosphi=1.
            sphi=0.
         else
            cosphi=(cose-ctheta*cos(the_q+thee))/
     &             (stheta*sin(the_q+thee))
            sinphi=sin(chi)/stheta*sine
         endif
         phi=ang(sinphi,cosphi)
c
c  end of calculation
c
      return
      end
c
c
      subroutine earray(theta0,phi0,dth,dph,icol,a,b,ivar1,omega0,
     &                romega,domega,ivar2,theta,phi,omega,nang,nomega)
      implicit double precision (a-h,o-z)
c
c
C this routine sets the theta-phi array of the electron directions
C according to the choosen collimator :
C     collimators are : ICOL=1 rectangular
C                            2 spherical
C                            3 elliptical
C     theta range   : A
C     phi range     : B
C
C NANG is the highest index in the theta-phi arrays
C
C in addition the array OMEGA is also set
C
C NOMEGA is the highest index in OMEGA
C
      parameter (mnp=10200)
      dimension THETA(0:mnp),PHI(0:mnp),OMEGA(0:100)
C
      data pi/3.141592654/
c
      dtr=pi/180.
c
      theta0=theta0*dtr
      phi0=phi0*dtr
c
      if (ivar1.lt.1) then
         ia=0
         theta(ia)=theta0
         phi(ia)=phi0
         nang=0
         goto 1000
      endif
c
      nt=2*idnint(a/2./dth)
      np=2*idnint(b/2./dph)
      if (nt.gt.100) nt=100
      if (np.gt.100) np=100
      rt=float (nt)
      rp=float (np)
c
      tstart=-rt/2.*dth
      pstart=-rp/2.*dph
c
      ia=1
      theta(0)=theta0
      phi(0)=phi0
      do 100 j=0,np
         phil=pstart+float(j)*dph
         do 100 k=0,nt
            thel=tstart+float(k)*dth
            goto (10,20,20) icol
10          if ((abs(thel).le.a/2).and.(abs(phil).le.b/2)) then
               if (ia.gt.mnp) then
                   write (6,*) ' electron angle array too big '
                   stop
               endif
               theta(ia)=thel+theta0
               phi(ia)=phil+phi0
               ia=ia+1 
               goto 100
            endif
            goto 100
20          test=4*thel**2/a**2+4*phil**2/b**2
            if (test.le.1) then
               if (ia.gt.mnp) then
                   write (6,*) ' electron angle array too big '
                   stop
               endif
               theta(ia)=thel+theta0
               phi(ia)=phil+phi0
               ia=ia+1
            endif
100   continue
      nang=ia-1
1000  if (ivar2.lt.1) then
         nomega=0
         omega(nomega)=omega0
         goto 999
      endif
      no=2*idnint(romega/domega/2.)
      if (no+1.gt.100) then
          write (6,*) ' omega array too big '
          stop
      endif
      r0=float(no)
      ostart=omega0-r0/2.*domega
      do 30 k=1,no+1
         omega(k)=ostart+float(k-1)*domega
30    continue
      omega(0)=omega0
      nomega=no+1
999   return
      end
               
c
c
      subroutine parray(theta0,phi0,dth,dph,icol,a,b,ivar,
     &                  theta,phi,nang)
      implicit double precision (a-h,o-z)
c
C this routine sets the theta-phi array of the proton directions
C according to the choosen collimator :
C     collimators are : ICOL=1 rectangular
C                            2 spherical
C                            3 elliptical
C     theta range   : A
C     phi range     : B
C
C NANG is the highest index in the theta-phi arrays
      parameter (mnp=10200)
C
      dimension THETA(0:mnp),PHI(0:mnp)
C
      data pi/3.141592654/
c
      dtr=pi/180.
c
      theta0=theta0*dtr
      phi0=phi0*dtr
c
      if (ivar.lt.1) goto 999
c
      nt=2*idnint(a/2./dth)
      np=2*idnint(b/2./dph)
      if (nt.gt.100) nt=100
      if (np.gt.100) np=100
      rt=float (nt)
      rp=float (np)
c
      tstart=-rt/2.*dth
      pstart=-rp/2.*dph
c
      ia=1
      theta(0)=theta0
      phi(0)=phi0
      do 100 j=0,np
         phil=pstart+float(j)*dph
         do 100 k=0,nt
            thel=tstart+float(k)*dth
            goto (10,20,20) icol
10          if ((abs(thel).le.a/2).and.(abs(phil).le.b/2)) then
               if (ia.gt.mnp) then
                   write (6,*) ' proton angle array too big '
                   stop
               endif
               theta(ia)=thel+theta0
               phi(ia)=phil+phi0
               ia=ia+1 
               goto 100
            endif
            goto 100
20          test=4*thel**2/a**2+4*phil**2/b**2
            if (test.le.1) then
               if (ia.gt.mnp) then
                   write (6,*) ' proton angle array too big '
                   stop
               endif
               theta(ia)=thel+theta0
               phi(ia)=phil+phi0
               ia=ia+1
            endif
100   continue
      nang=ia-1
      return
999   ia=0
      theta(ia)=theta0
      phi(ia)=phi0
      nang=ia
      return
      end
               
c
c
      subroutine barray(theta0,phi0,dth,dph,icol,a,b,ivar1,e0,
     &                 re,de,ivar2,theta,phi,e,nang,ne)
      implicit double precision (a-h,o-z)
c
C this routine sets the theta-phi array of the electron beam directions
C according to the choosen beam spot form :
C     collimators are : ICOL=1 rectangular
C                            2 spherical
C                            3 elliptical
C     theta range   : A
C     phi range     : B
C
C NANG is the number of theta-phi angles pairs created
C
C in addition the array of the beam energies E is also set
C
C NE is the highest index in E
c
      parameter (mnp=10200)
      dimension THETA(0:mnp),PHI(0:mnp),E(0:50)
c
      data pi/3.141592654/
c
      dtr=pi/180.
c
      theta0=theta0*dtr
      phi0=phi0*dtr
c
      if (ivar1.lt.1) then
         ia=0
         theta(ia)=theta0
         phi(ia)=phi0
         nang=0
         goto 1000
      endif
c
      nt=2*idnint(a/2./dth)
      np=2*idnint(b/2./dph)
      if (nt.gt.100) nt=100
      if (np.gt.100) np=100
      rt=float (nt)
      rp=float (np)
c
      tstart=-rt/2.*dth
      pstart=-rp/2.*dph
c
      ia=1
      theta(0)=theta0
      phi(0)=phi0
      do 100 j=0,np
         phil=pstart+float(j)*dph
         do 100 k=0,nt
            thel=tstart+float(k)*dth
            goto (10,20,20) icol
10          if ((abs(thel).le.a/2).and.(abs(phil).le.b/2)) then
               if (ia.gt.mnp) then
                   write (6,*) ' beam angle array too big '
                   stop
               endif
               theta(ia)=thel+theta0
               phi(ia)=phil+phi0
               ia=ia+1 
               goto 100
            endif
            goto 100
20          test=4*thel**2/a**2+4*phil**2/b**2
            if (test.le.1) then
               if (ia.gt.mnp) then
                   write (6,*) ' beam angle array too big '
                   stop
               endif
               theta(ia)=thel+theta0
               phi(ia)=phil+phi0
               ia=ia+1
            endif
100   continue
      nang=ia-1
1000  if (ivar2.lt.1) then
         ne=0
         e(ne)=e0
         goto 999
      endif
      ne=2*idnint(re/de/2.)
      if (ne+1.gt.50) then
         write (6,*) ' incident energy array too big '
         stop
      endif
      if (ne.gt.200) ne=200
      rre=float(ne)
      estart=e0-rre/2.*de
      do 30 j=1,ne+1
         e(j)=estart+float(j-1)*de
30    continue
      e(0)=e0
999   return
      end
               
c
c
      subroutine cccos(theta,ctheta,nang)
      implicit double precision (a-h,o-z)
c
c
      parameter (mnp=10200)
      dimension THETA(0:mnp),CTHETA(0:mnp)
      do 10 j=0,nang
        ctheta(j)=cos(theta(j))
10    continue
      return
      end
c
c
      subroutine ccsin(theta,stheta,nang)
      implicit double precision (a-h,o-z)
c
c
      parameter (mnp=10200)
      dimension THETA(0:mnp),STHETA(0:mnp)
      do 10 j=0,nang
         stheta(j)=sin(theta(j))
10    continue
      return 
      end
c
c
c
c	subroutine subcc1
c
c	author:	R.W. Lourie
c	date:	V1.0	18-AUG-1986
c	modifications:	W. Boeglin
c	        Instead of angles , the four vectors of the incident 
c		beam, scattered electron and the outgoing proton are
c		entered in the argument list.
c		floor: x,z plane
c		       y-axis points out if the plane
c		       x-axis points to the electron side
c		       z-axis points in the beam direction
c	purpose:
C		Calculates elementary OFF-SHELL (e,p) cross-section
C		using the prescription given by DeForest in
C		Nucl. Phys.  A392 (1983)
C
c	variables:
c	input:
C		k_i : incident 4-vector
c		k_f : scattered 4-vector
c		p_f : outgoing proton 4-vector
C	OUTPUT:
C		sigma_ep: The off-shell cross section.
C		k_fact:   Kinematic factor ( = E_p*P_F)
C		pi_mag:   Magnitude of initial momentum (determined as
C			  PF - Q componentwise).
C
C	ROUTINES CALLED:  
C		form_fact
C
C
        subroutine subcc1 (k_i,k_f,p_f,  sigma_ep,k_fact,pi_mag)
	implicit real*8 (a-h,k,m,o-z)
	real*8 k_i(4),k_f(4),q(4),p_i(4),p_f(4),n_e(4),n_p(4)
	character*7 fftype
c Proton mass in MeV
	parameter (mp = 938.28)	
	parameter (pi = 3.14159)
	parameter (hbarc = 197.3286)
c fine structure constant
	parameter (alpha = 7.297e-3)		
c Proton anomalous moment
	parameter (kappa = 1.793)			
	data a,z,m_eff /1.,1.,1./
	data fftype /'mainz  '/
c
        dtr=pi/180.
        e_p_f=p_f(4)
        e_f=k_f(4)
c
C	Form components of momentum transfer 4-vector.
C	Form components of proton initial 4-momentum vector.
C
	do  10 i=1,4
	   q(i) = k_i(i) - k_f(i)
	   p_i(i) = p_f(i) - q(i)
10	continue
c
c vector normal to electron scattering plane
c
        call cross(k_i,k_f,n_e)
c
c vector normal to q p_f plane
c
        call cross(q,p_f,n_p)
c
c calculate out of plane angle phi_p_f
c
        rn_e=sqrt(dot3(n_e,n_e))
        rn_p=sqrt(dot3(n_p,n_p))
        cosphi=dot3(n_e,n_p)/(rn_e*rn_p)
        phi_p_f=acos(cosphi)
c
C	Form all dot products.
C
	qmu2  = -dot(q,q)
        pf=sqrt(p_f(4)**2-dot(p_f,p_f))
c
c calculate the electron scattering angle
c
        sin_phi=sqrt(qmu2/(4.*k_i(4)*k_f(4)))
        phi_e=2.*asin(sin_phi)
c
	pf_dot_q = 0.0
	pi_mag = 0.0
	q_mag = 0.0
	do 20 i=1,3
	   pi_mag = pi_mag + p_i(i)**2
	   q_mag = q_mag + q(i)**2
	   pf_dot_q = pf_dot_q + p_f(i)*q(i)
20	continue
	pi_mag = sqrt(pi_mag)
	q_mag = sqrt(q_mag)
c Angle between P_F and Q
	gamma = acos(pf_dot_q/(pf*q_mag))    
c
C	Set up kinematic factors etc.
C
c Initial proton energy
	e_bar = sqrt(pi_mag**2+mp**2)		
	omega_bar = e_p_f - e_bar
	q_mu_bar = sqrt(q_mag**2 - omega_bar**2)
	ebp = e_bar*e_p_f
	x = q_mu_bar**2/4./mp**2
	q_mu = (abs(qmu2))**.5
	qq = (q_mu/q_mag)**2
	tan_e2 = tan(phi_e/2.)**2
C
C	Put together various response functions
C
	call form_fact(fftype,a,z,q_mu,m_eff,ge2,
     #			     gm2,f1p,f2p,f1n,f2n)
	w_c = ((e_bar+e_p_f)**2*(f1p**2+kappa**2*x*f2p**2)-q_mag**2
     #      *(f1p+kappa*f2p)**2)/(4.*ebp)
c
	w_t = (f1p+kappa*f2p)**2*q_mu_bar**2/2./ebp
c
	w_s = (f1p**2+kappa**2*x*f2p**2)*pf**2*sin(gamma)**2/ebp
c
	w_i = -(e_bar+e_p_f)*(f1p**2+kappa**2*x*f2p**2)*pf
     #      * sin(gamma)/ebp
C
C	Calculate Off-shell cross-section.
C
	sigma_mott = 4.*((alpha*cos(phi_e/2.)*e_f)/qmu2)**2
	sigma_ep = sigma_mott*(qq**2*w_c + (qq/2.+tan_e2)*w_t
     #           + qq*sqrt(qq+tan_e2)*w_i*cos(phi_p_f)
     #           + (qq*cos(phi_p_f)**2+tan_e2)*w_s)*hbarc**2.*10000.
	k_fact = e_p_f*pf
c
        return
	end
C
C	Function DOT(A,B)
C
C	Purpose:
C
C		Takes dot product of two 4-vectors.
C		Metric:  DOT(A,A) = A_time**2 - A_space .dot. A_space
C
C
	function dot(a,b)
	implicit real*8 (a-h,k,o-z)
	real*8 a(4),b(4)
	space = 0.
	do 10 i = 1,3
	   space = space + a(i)*b(i)
10	continue
	time = a(4)*b(4)
	dot = time - space
	return
	end
c
C	SUBROUTINE CROSS(A,B,C)
C
C	Purpose:
C
C		Takes the cross product of the space part of two 
C	        4-vectors.
C		CROSS(A,B,C): C = A X B
C
C
	subroutine cross(a,b,c)
	implicit real*8 (a-h,k,o-z)
	real*8 a(4),b(4),c(4)
c
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        c(4)=0.
c
	return
	end
c
C	Function DOT3(A,B)
C
C	Purpose:
C
C		Takes 3-dot product of two 4-vectors.
C
C
	function dot3(a,b)
	implicit real*8 (a-h,k,o-z)
	real*8 a(4),b(4)
	space = 0.
	do 10 i = 1,3
	   space = space + a(i)*b(i)
10	continue
	dot3 = space
	return
	end
c
C	Function DOT33(A,B)
C
C	Purpose:
C
C		Takes 3-dot product of two 3-vectors.
C
C
	function dot33(a,b)
	implicit real*8 (a-h,k,o-z)
	real*8 a(3),b(3)
	space = 0.
	do 10 i = 1,3
	   space = space + a(i)*b(i)
10	continue
	dot33 = space
	return
	end
c
C	SUBROUTINE CROSS3(A,B,C)
C
C	Purpose:
C
C		Takes the cross product of the space part of two 
C	        3-vectors.
C		CROSS(A,B,C): C = A X B
C
C
	subroutine cross3(a,b,c)
	implicit real*8 (a-h,k,o-z)
	real*8 a(3),b(3),c(3)
c
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
c
	return
	end
c
c
c subroutine normalize : replace a vector by its unit vector
c
c call normalize (vec)
c
	subroutine normalize(vec)
	implicit real*8 (a-h,o-z)
	real*8 vec(3), norm
c
	norm = sqrt(dot33(vec,vec))
	do i=1,3
	   vec(i) = vec(i)/norm
	enddo
	return
	end

c
C
C -------------------------------------------------------------------
C     Deuteron momentum distribution.
C
C     Momenta are in GeV/c.
C     Constants are from Krautschneider's Ph.D. thesis
C -------------------------------------------------------------------
C
      double precision function h2spec(prmag)
      implicit real*8 (a-h,o-z)
      parameter (pi = 3.14159)
      if(prmag .le. 300.)then
c No rescattering.
        term = 0.         
      else
	term = 4.
      endif
c Convert to GeV/c
      pr = prmag/1000.  
      pr2 = pr**2
      t1 = pr2 + 0.002088
      t2 = pr2 + 0.0676
c
C     The constant TERM from K. Mueller Ph.D. thesis simulates
C     rescattering contribution for high
C     momenta to reach agreement with experimental data.
C     Normal Value:  TERM=4.
c
      h2spec = (1./t1-1./t2)**2 + term
      h2spec = 4.*pi*h2spec * 1.0e-12
      return
      end
c
c
c       subroutine form_fact
C
C
C	AUTHOR:	B.H.Cottman
C	DATE:	V1.0	23-SEP-83
C	MODIFICATIONS:	
C		Lots of additions by RWL: Janssens, Blatnik and Zovko,
C		and Gari/Krumpelmann formfactors.
C	PURPOSE:
C		Calculate dipole, MAINZ fit or FIT 8.2 mean formfactors.
C	VARIABLES:
C	INPUT:
C		a:	Number of nucleons
C		z:	Number of protons
C		q_mu:	4-momentum
C		m_eff:	Factor to scale nucleon mass.	
C		
C	OUTPUT:
C		ge2:	Electric formfactor squared	
C		gm2:	Magnetic formfactor squared	
C		f1p,n:	Dirac formfactor for proton and neutron
C		f2p,n:	Pauli formactor for proton and neutron
C	ROUTINES CALLED:  
C			None
C
	subroutine form_fact(fftype,a,z,q_mu,m_eff,ge2,gm2,
     #			     f1p,f2p,f1n,f2n)
c
	implicit real*8 (a-h,k,m,o-z)
	character*7 fftype
	parameter (pi = 3.1415927)
	parameter (mp = 938.2796)
	parameter (hbarc = 197.3286)
c
C  dipole fit
	gdp(qq)=1./(1.+qq/0.71e+6)**2				
c
C mainz fit
	g_e_m(t)=0.312/(1+t/6.0)+				
     &		     1.312/(1+t/15.02)+
     &               (-0.709)/(1+t/44.08)+
     &               0.085/(1+t/154.2)
	g_m_m(t)=( 0.694/(1+t/8.5)+
     &		       0.719/(1+t/15.02)+
     &                 (-0.418)/(1+t/44.08)+
     &                 0.005/(1+t/355.4) )
c
	q2 = q_mu**2
	q2_f=q2/hbarc**2
	mn = m_eff*mp
C neutron mag. moment
	mun=-1.913					
C  proton mag. moment
	mup=2.793					
	kappan=mun
	kappap=mup-1.
	mpi = 0.1395
	x=1.+q2/4./mp**2
c
C	Dipole fit
C
	if(fftype .eq. 'dipole ') then
C  dipole chosen
	 f1p=(1.+q2*mup/4./mp**2)*gdp(q2)/x		
	 f1n=(q2*mun/2./mp**2)*gdp(q2)/x
	 f2p=gdp(q2)/x
	 f2n=(2.-x)*gdp(q2)/x
	elseif(fftype .eq. 'mainz  ') then
c
C mainz fit
c
	 f1p=(g_e_m(q2_f)+(x-1.)*mup*g_m_m(q2_f))/x		
	 f1n=(q2*mun/2./mp**2)*g_m_m(q2_f)/x
	 f2p=(mup*g_m_m(q2_f)-g_e_m(q2_f))/(kappap*x)
	 f2n=(2.-x)*g_m_m(q2_f)/x
	elseif(fftype .eq. 'fit8_2 ') then
c
C fit 8.2 chosen
c
	 t=-q2/1.e+6/MPI**2					
	 f1rho = (.4775+.045/(1-t*.01948/.355)**2) / (1-t*.01948/.536)
	 f2rho = (2.6675+.481/(1-t*.01948/.268)) / (1-t*.01948/.603)
	 f1v = f1rho + (0.03144/(1-t/75.736)) - (0.08575/(1-t/308.179))
     1   + (.03164/(1-t/446.642))
	 f1s = (-.61668/(1-t/53.37)) + (1.15702/(1-t/31.468))
     1   - (.040333/(1-t/166.648))
	 f2v = f2rho - (1.3512/(1-t/75.736)) + (.03295/(1-t/308.179))
     1   + (.02151/(1-t/446.642))
	 f2s = (.12097/(1-t/53.37)) + (-.17534/(1-t/31.468))
     1   - (.00363/(1-t/166.648))
	 f1p = (f1v+f1s)
	 f1n = f1s-f1v
	 f2p = (f2v+f2s)/kappap
	 f2n = (f2s-f2v)/kappan
	elseif(fftype .eq. 'gari   ')then
c
c Gari
c
	  call gari(q2/1.0e+06,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	elseif(fftype .eq. 'blatnik')then
c
c Blatnik
c
	  call b_and_z(-q2/1.0e+06,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	elseif(fftype .eq. 'janssen')then
c
c Janssen
c
	  call janssens(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	  goto 777
	endif
c
	xx=  q2/4.0/mn/mn
C  neutron and proton...
	gen=f1n-xx*kappan*f2n			
C  elec. and mag. form..
	gep=f1p-xx*kappap*f2p			
C  factors
	gmn=f1n+kappan*f2n			
	gmp=f1p+kappap*f2p
C  nuclear elec. ff per nucleon
777	ge2=(z*gep**2+(a-z)*gen**2)		
C nuclear mag. ff per nucleon
	gm2=(z*gmp**2+(a-z)*gmn**2)		
c
	return
	end
c
	subroutine gari(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
c
c   program by gond yen, november 20, 1986
c
C   nucleon form factor parametrization by M. Gari and
C   W. Kruempelmann:
C   --- Z. Phys. A322(1985)689-693
C       Phys. Lett. B173(1986)10-14
C
C   the unit used in this code is GeV
C
C   dipole form: FD
C   FD=1.0/(Q2+0.71)**2
C
      implicit real*8 (a-h,l,k,m,o-z)
      data mn,mr,mw /0.938e0,0.776e0,0.784e0/
      data k,kv,ks,kr,kw /1.793e0,3.706e0,-0.12e0,6.62e0,0.163e0/
      data cr,cw /0.377e0,0.411e0/
      data lam1,lam2,lmqcd /0.795e0,2.27e0,0.29e0/
      mr2=mr**2
      mw2=mw**2
      mn2=mn**2
      lm1SQ=lam1**2
      lm2SQ=lam2**2
      lmqSQ=lmqcd**2
      qt2=q2*dlog((lm2SQ+q2)/lmqSQ)/dlog(lm2SQ/lmqSQ)
      f1=lm1SQ/(lm1SQ+qt2)*lm2SQ/(lm2SQ+qt2)
      f2=lm1SQ/(lm1SQ+qt2)*(lm2SQ/(lm2SQ+qt2))**2
      f1v=f1*(cr*mr2/(mr2+q2)+1.0-cr)
      f1s=f1*(cw*mw2/(mw2+q2)+1.0-cw)
      kf2v=f2*(cr*kr*mr2/(mr2+q2)+kv-cr*kr)
      kf2s=f2*(cw*kw*mw2/(mw2+q2)+ks-cw*kw)
C  Dirac and Pauli form factors: F1, F2's
      f1n=0.5*(f1s-f1v)
      f2n=0.5*(kf2s-kf2v)/((ks-kv)/2.)
      f1p=0.5*(f1s+f1v)
      f2p=0.5*(kf2s+kf2v)/k
C  Sachs form factors: GE, GM's
      gmn=f1n+(ks-kv)/2.*f2n
      gen=f1n-(ks-kv)/2.*q2/4.0/mn2*f2n
      gmp=f1p+k*f2p
      gep=f1p-k*q2/4.0/mn2*f2p
c
      return
      end

	subroutine b_and_z(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	implicit real*8 (a-h,k,m,o-z)
	data bs,bv,mu_s,mu_v,mp /-0.91,-1.10,-0.06,1.853,0.93828/
	data t_r,t_o,t_p,t_rp,t_rpp,t_op /0.585,0.614,1.039,1.3,
     #	     2.1,1.4/

c
	r_s(t) = 1./((1.-t/t_o)*(1.-t/t_p)*(1.-t/t_op))	
	r_v(t) = 1./((1.-t/t_r)*(1.-t/t_rp)*(1.-t/t_rpp))	
c
	x = q2/4./mp**2
	mu_p = mu_s + mu_v + 1.
	mu_n = mu_s - mu_v
	gep = (0.5 + x*(mu_s + 2.*mp**2*bs))*r_s(q2)
     #	    + (0.5 + x*(mu_v + 2.*mp**2*bv))*r_v(q2)
	gen = (0.5 + x*(mu_s + 2.*mp**2*bs))*r_s(q2)
     #	    - (0.5 + x*(mu_v + 2.*mp**2*bv))*r_v(q2)
	gmp = (0.5 + mu_s + bs*q2/2.)*r_s(q2)
     #      + (0.5 + mu_v + bv*q2/2.)*r_v(q2)
	gmn = (0.5 + mu_s + bs*q2/2.)*r_s(q2)
     #      - (0.5 + mu_v + bv*q2/2.)*r_v(q2)
	f1p = (gep - x*gmp)/(1.-x)
	f2p = (gmp - gep)/(mu_p-1.)/(1.-x)
	f1n = (gen - x*gmn)/(1.-x)
	f2n = (gmn - gen)/mu_n/(1.-x)
	return
	end
c
	subroutine janssens(q2,f1p,f2p,f1n,f2n,gep,gmp,gen,gmn)
	implicit real*8 (a-h,k,m,o-z)
	data mp,mu_p,mu_n,hbarc /938.28,2.793,-1.913,197.3286/
c
	ge_s(t) = (2.5/(1.+t/15.7)-1.6/(1.+t/26.7)+0.1)*0.5
	gm_s(t) = (3.33/(1.+t/15.7)-2.77/(1.+t/26.7)+0.44)*0.44
	ge_v(t) = (1.16/(1.+t/8.19)-0.16)*0.5
	gm_v(t) = (1.11/(1.+t/8.19)-0.11)*2.353
c
	tt = q2/hbarc/hbarc
	gep = ge_s(tt) + ge_v(tt)
	gen = ge_s(tt) - ge_v(tt)
	gmp = gm_s(tt) + gm_v(tt)
	gmn = gm_s(tt) - gm_v(tt)
c
	x = q2/4./mp**2
	f1p = (gep+x*gmp)/(1.+x)
	f2p = (gmp-gep)/(1.+x)/(mu_p-1.)
	f1n = (gen+x*gmn)/(1.+x)
	f2n = (gmn-gen)/(1.+x)/mu_n
	return
	end

      real function rho_paris(pm)
c 
c return paris momentum distribution, does a polynomial interpolation
c between data points
c
c requires : polint and locate
c
      implicit none
      integer MAXDAT 
      parameter (MAXDAT=100)
      
      real p(MAXDAT), rho(MAXDAT)
      real pm, rho_, drho_

      integer j,k,m

      data m/4/  ! 4 point interpolation
      
c----------------------------------------------------------------------
c intialzie data: Paris momentum distribution
c
       data p(  1)/ 9.866E+00/,rho(  1)/ 1.540E+03/
       data p(  2)/ 1.973E+01/,rho(  2)/ 1.188E+03/
       data p(  3)/ 2.960E+01/,rho(  3)/ 8.177E+02/
       data p(  4)/ 3.947E+01/,rho(  4)/ 5.299E+02/
       data p(  5)/ 4.933E+01/,rho(  5)/ 3.358E+02/
       data p(  6)/ 5.920E+01/,rho(  6)/ 2.129E+02/
       data p(  7)/ 6.907E+01/,rho(  7)/ 1.366E+02/
       data p(  8)/ 7.893E+01/,rho(  8)/ 8.918E+01/
       data p(  9)/ 8.880E+01/,rho(  9)/ 5.933E+01/
       data p( 10)/ 9.867E+01/,rho( 10)/ 4.021E+01/
       data p( 11)/ 1.085E+02/,rho( 11)/ 2.773E+01/
       data p( 12)/ 1.184E+02/,rho( 12)/ 1.944E+01/
       data p( 13)/ 1.283E+02/,rho( 13)/ 1.383E+01/
       data p( 14)/ 1.381E+02/,rho( 14)/ 9.979E+00/
       data p( 15)/ 1.480E+02/,rho( 15)/ 7.289E+00/
       data p( 16)/ 1.579E+02/,rho( 16)/ 5.385E+00/
       data p( 17)/ 1.677E+02/,rho( 17)/ 4.021E+00/
       data p( 18)/ 1.776E+02/,rho( 18)/ 3.032E+00/
       data p( 19)/ 1.875E+02/,rho( 19)/ 2.308E+00/
       data p( 20)/ 1.973E+02/,rho( 20)/ 1.773E+00/
       data p( 21)/ 2.072E+02/,rho( 21)/ 1.374E+00/
       data p( 22)/ 2.171E+02/,rho( 22)/ 1.074E+00/
       data p( 23)/ 2.269E+02/,rho( 23)/ 8.461E-01/
       data p( 24)/ 2.368E+02/,rho( 24)/ 6.723E-01/
       data p( 25)/ 2.467E+02/,rho( 25)/ 5.387E-01/
       data p( 26)/ 2.565E+02/,rho( 26)/ 4.352E-01/
       data p( 27)/ 2.664E+02/,rho( 27)/ 3.546E-01/
       data p( 28)/ 2.763E+02/,rho( 28)/ 2.914E-01/
       data p( 29)/ 2.861E+02/,rho( 29)/ 2.416E-01/
       data p( 30)/ 2.960E+02/,rho( 30)/ 2.020E-01/
       data p( 31)/ 3.059E+02/,rho( 31)/ 1.703E-01/
       data p( 32)/ 3.157E+02/,rho( 32)/ 1.448E-01/
       data p( 33)/ 3.256E+02/,rho( 33)/ 1.242E-01/
       data p( 34)/ 3.355E+02/,rho( 34)/ 1.074E-01/
       data p( 35)/ 3.453E+02/,rho( 35)/ 9.357E-02/
       data p( 36)/ 3.552E+02/,rho( 36)/ 8.213E-02/
       data p( 37)/ 3.651E+02/,rho( 37)/ 7.259E-02/
       data p( 38)/ 3.749E+02/,rho( 38)/ 6.456E-02/
       data p( 39)/ 3.848E+02/,rho( 39)/ 5.777E-02/
       data p( 40)/ 3.947E+02/,rho( 40)/ 5.196E-02/
       data p( 41)/ 4.045E+02/,rho( 41)/ 4.697E-02/
       data p( 42)/ 4.144E+02/,rho( 42)/ 4.263E-02/
       data p( 43)/ 4.243E+02/,rho( 43)/ 3.884E-02/
       data p( 44)/ 4.341E+02/,rho( 44)/ 3.551E-02/
       data p( 45)/ 4.440E+02/,rho( 45)/ 3.255E-02/
       data p( 46)/ 4.539E+02/,rho( 46)/ 2.991E-02/
       data p( 47)/ 4.637E+02/,rho( 47)/ 2.754E-02/
       data p( 48)/ 4.736E+02/,rho( 48)/ 2.540E-02/
       data p( 49)/ 4.835E+02/,rho( 49)/ 2.346E-02/
       data p( 50)/ 4.933E+02/,rho( 50)/ 2.169E-02/
       data p( 51)/ 5.032E+02/,rho( 51)/ 2.007E-02/
       data p( 52)/ 5.131E+02/,rho( 52)/ 1.858E-02/
       data p( 53)/ 5.229E+02/,rho( 53)/ 1.722E-02/
       data p( 54)/ 5.328E+02/,rho( 54)/ 1.595E-02/
       data p( 55)/ 5.427E+02/,rho( 55)/ 1.478E-02/
       data p( 56)/ 5.525E+02/,rho( 56)/ 1.370E-02/
       data p( 57)/ 5.624E+02/,rho( 57)/ 1.270E-02/
       data p( 58)/ 5.723E+02/,rho( 58)/ 1.176E-02/
       data p( 59)/ 5.821E+02/,rho( 59)/ 1.089E-02/
       data p( 60)/ 5.920E+02/,rho( 60)/ 1.009E-02/
       data p( 61)/ 6.019E+02/,rho( 61)/ 9.333E-03/
       data p( 62)/ 6.117E+02/,rho( 62)/ 8.632E-03/
       data p( 63)/ 6.216E+02/,rho( 63)/ 7.979E-03/
       data p( 64)/ 6.315E+02/,rho( 64)/ 7.370E-03/
       data p( 65)/ 6.413E+02/,rho( 65)/ 6.804E-03/
       data p( 66)/ 6.512E+02/,rho( 66)/ 6.276E-03/
       data p( 67)/ 6.611E+02/,rho( 67)/ 5.785E-03/
       data p( 68)/ 6.709E+02/,rho( 68)/ 5.327E-03/
       data p( 69)/ 6.808E+02/,rho( 69)/ 4.902E-03/
       data p( 70)/ 6.907E+02/,rho( 70)/ 4.507E-03/
       data p( 71)/ 7.005E+02/,rho( 71)/ 4.140E-03/
       data p( 72)/ 7.104E+02/,rho( 72)/ 3.799E-03/
       data p( 73)/ 7.203E+02/,rho( 73)/ 3.483E-03/
       data p( 74)/ 7.301E+02/,rho( 74)/ 3.190E-03/
       data p( 75)/ 7.400E+02/,rho( 75)/ 2.918E-03/
       data p( 76)/ 7.499E+02/,rho( 76)/ 2.667E-03/
       data p( 77)/ 7.597E+02/,rho( 77)/ 2.434E-03/
       data p( 78)/ 7.696E+02/,rho( 78)/ 2.220E-03/
       data p( 79)/ 7.795E+02/,rho( 79)/ 2.021E-03/
       data p( 80)/ 7.893E+02/,rho( 80)/ 1.839E-03/
       data p( 81)/ 7.992E+02/,rho( 81)/ 1.670E-03/
       data p( 82)/ 8.091E+02/,rho( 82)/ 1.516E-03/
       data p( 83)/ 8.189E+02/,rho( 83)/ 1.373E-03/
       data p( 84)/ 8.288E+02/,rho( 84)/ 1.242E-03/
       data p( 85)/ 8.387E+02/,rho( 85)/ 1.123E-03/
       data p( 86)/ 8.485E+02/,rho( 86)/ 1.013E-03/
       data p( 87)/ 8.584E+02/,rho( 87)/ 9.123E-04/
       data p( 88)/ 8.683E+02/,rho( 88)/ 8.205E-04/
       data p( 89)/ 8.781E+02/,rho( 89)/ 7.367E-04/
       data p( 90)/ 8.880E+02/,rho( 90)/ 6.603E-04/
       data p( 91)/ 8.979E+02/,rho( 91)/ 5.908E-04/
       data p( 92)/ 9.077E+02/,rho( 92)/ 5.277E-04/
       data p( 93)/ 9.176E+02/,rho( 93)/ 4.704E-04/
       data p( 94)/ 9.275E+02/,rho( 94)/ 4.184E-04/
       data p( 95)/ 9.373E+02/,rho( 95)/ 3.715E-04/
       data p( 96)/ 9.472E+02/,rho( 96)/ 3.291E-04/
       data p( 97)/ 9.571E+02/,rho( 97)/ 2.908E-04/
       data p( 98)/ 9.669E+02/,rho( 98)/ 2.564E-04/
       data p( 99)/ 9.768E+02/,rho( 99)/ 2.255E-04/
       data p(100)/ 9.867E+02/,rho(100)/ 1.978E-04/
c----------------------------------------------------------------------
      
c from NR Vol2 p.113      

      call locate (p, MAXDAT, pm, j)
      k = min( max(j-(m-1)/2,1), MAXDAT+1-m)

c interpolate

      call polint(p(k), rho(k), m, pm, rho_, drho_)

      rho_paris = rho_/1.0d9     ! convert from GeV^-3 to MeV^-3

      return
      end




      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=100) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)STOP ' PROBLEM in POLINT'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      SUBROUTINE LOCATE(XX,N,X,J)
      DIMENSION XX(N)

      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END

