c FORTRAN Code: Monte Carlo Simulation for fraction of muons showing an RID event
c Author: Robert J. Nemiroff
c Computes fraction of HAWC events that show RID.
       program fractionRID
c  
c      integer*4, parameter :: n=10
       integer*4 i, ntot, nRID, n, nmax
       real*8 c, pi, seed, pin, eps, Cangle
       real*8 R, H, Rd
       real*8 xd(4), yd(4), zd(4), th(4)
       real*8 xg, yg, zg, rg
       real*8 xs, ys, zs, rs
       real*8 xm, ym, zm, rm
       real*8 xe, ye, ze, re
       real*8 x_ge, y_ge, z_ge
       real*8 x_de(4), y_de(4), z_de(4), r_de(4)
       real*8 top, bot
       character*3 before, now, RID
c                             Constants
       nmax=1000000000
       c=3.0d8
       pi=dacos(-1.0d0)
       seed=2.8d0
       eps=0.1d0
       Cangle=dacos(1.0d0/1.33d0)*(180.0d0/pi)
       write (*,*) ' Input seed, please: '
       read (*,*) seed
c                             R=radius of tank (meters)
c                             H=height of tank (m)
c Source: https://www.sciencedirect.com/
c science/article/pii/S0273117713001695?via%3Dihub
       R=7.30d0/2.0d0
       H=4.5d0
c                             Ground dist. to detectors (check)
       Rd=7.3d0/4.0d0
       Rd=1.850d0
c                             Location of four ground detectors
       xd(1)=0.0d0
       yd(1)=0.0d0
       zd(1)=0.0d0
       xd(2)=-Rd
       yd(2)=0.0d0
       zd(2)=0.0d0
       xd(3)=Rd*cos(pi/3.0d0)
       yd(3)=Rd*sin(pi/3.0d0)
       zd(3)=0.0d0
       xd(4)=Rd*cos(pi/3.0d0)
       yd(4)=-Rd*sin(pi/3.00)
       zd(4)=0.0d0
c
       ntot=0
       nRID=0
c ********************
c                             Grand loop over muons
c ********************
       do 700 n=1, nmax
c                             Location muon hits ground
       xg=3.0d0*R*(2.0d0*pin(seed) - 1.0d0)
       yg=3.0d0*R*(2.0d0*pin(seed) - 1.0d0)
       zg=0.0d0
       rg=sqrt(xg**2 + yg**2)
       if (rg.ge.5.0d0*R) goto 700
c                             Direction muon arrived on sky
 200   xs=(2.0d0*pin(seed) - 1.0d0)
       ys=(2.0d0*pin(seed) - 1.0d0)
       zs=pin(seed)
       phi=atan(ys,xs)
c                             Angular direction of muon
       phi=2.0d0*pi*pin(seed)
       theta=acos(zs)
c                             Horizon dimming
       if (pin(seed).gt.zs**2) goto 200
c
c ********************** Find where muon enters tank ***
c
       xm=xg
       ym=yg
       zm=zg
       rm=sqrt(xm**2 + ym**2)
c                             Is muon in the tank?
       now='out'
       if (rm.lt.R.and.zm.lt.H) now='in'
c                             Step the photon back along track
 300   xm=xm + xs*eps
       ym=ym + ys*eps
       zm=zm + zs*eps
       rm=sqrt(xm**2 + ym**2)
c      write (*,310) xm/R, ym/R, zm/H, rm/R
c310   format (2x, 4(f7.4, 3x))
c                             If muon missed tank try again
c      if (rm.gt.3.0d0*R) write (*,*) ' Muon missed tank: wide'
       if (rm.gt.3.0d0*R) goto 700
c      if (zm.gt.H) write (*,*) ' Muon missed tank: high'
       if (zm.gt.H) goto 700
       before=now
c
       now='out'
       if (rm.lt.R.and.zm.lt.H) now='in'
c      if (before.eq.'out'.and.now.eq.'in') write (*,*) "In"
c                             If muon left tank then analyze
c      if (before.eq.'in'.and.now.eq.'out') write (*,*) "Out"
       if (before.eq.'in'.and.now.eq.'out') goto 500
       goto 300

c                             The muon has just entered the tank
c                             Three points: g, e, d
c                             For ground, entry, detector
c                             Find the vertex angles
 500   xe=xm
       ye=ym
       ze=zm
       re=rm
c      write (*,*) 
c      write (*,*) ' xm / R = ', xm/R
c      write (*,*) ' ym / R = ', ym/R
c      write (*,*) ' zm / H = ', zm/H
c      write (*,*) ' rm / R = ', rm/R
c      write (*,*) ' rg / R = ', rg/R
c
       x_ge=xg - xe
       y_ge=yg - ye
       z_ge=zg - ze
       r_ge=sqrt(x_ge**2 + y_ge**2 + z_ge**2)
c      write (*,*) '   phi = ', phi*(180.0d0/pi)
c      write (*,*) ' theta = ', theta*(180.0d0/pi)
c      write (*,*) 
c
       RID='no'
       ntot=ntot + 1
       do 600 i=1, 4
       x_de(i)=xd(i) - xe
       y_de(i)=yd(i) - ye
       z_de(i)=zd(i) - ze
       r_de(i)=sqrt(x_de(i)**2 + y_de(i)**2 + z_de(i)**2)
c 
       top=x_ge*x_de(i) + y_ge*y_de(i) + z_ge*z_de(i)
       bot=r_ge * r_de(i)
       th(i)=acos(top/bot)*(180.0d0/pi)
c      write (*,*) ' th(i) = ', th(i)
       if (th(i).lt.Cangle) RID='yes'
 600   continue
       if (RID.eq.'yes') nRID=nRID + 1
 700   continue
c
       write (*,*) ' nRID, ntot = ', nRID, ntot
       write (*,*) ' RID fraction = ', float(nRID)/float(ntot)
       
c                             Stop the madness.
 999   stop
       end
c
c*******************************************************
c
c  Subroutine that uses pi**n algorithm.
       real*8 function pin(pseed)
       implicit real*8 (a-h,o-z)
c
       pi=dacos(-1.0d0)
       pseed=pseed*pi
       if (pseed.gt.10.0d0) pseed=pseed/10.0d0
       pin=pseed*10000.0d0 - int(pseed*10000.0d0)
c
       return
       end
