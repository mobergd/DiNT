c
c   Dint – version 2.0  is licensed under the Apache License, Version 2.0 (the "License");
c   you may not use Dint – version 2.0 except in compliance with the License.
c   You may obtain a copy of the License at
c       http://www.apache.org/licenses/LICENSE-2.0
c   The license is also given in the LICENSE file.
c   Unless required by applicable law or agreed to in writing, software
c   distributed under the License is distributed on an "AS IS" BASIS,
c   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
c   See the License for the specific language governing permissions and limitations under the License.
c
c -------------------------------------------------------------------------------------------
c  Dint : Direct Nonadiabatic Trajectories A code for non-Born–Oppenheimer molecular dynamics 
c  
c  version 2.0                                    
c
c  A. W. Jasper                  
c  Argonne National Laboratory     
c
c  Rui Ming Zhang                 
c  Tsinghua University
c               
c  and                  
c    
c  D. G. Truhlar                 
c  University of Minnesota
c
c  copyright  2020
c  Donald G. Truhalar and Regents of the University of Minnesota 
c----------------------------------------------------------------------------------------------

      subroutine getrelbz2(rcom,eint,nbounce,tdelay,td0,tdf)

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'

      integer i,j,k,arr,index,nprod,nbounce
      double precision rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf),tmp,tmpi,
     & rr(3),jp(9),jq(9),jm(9),xv,xj,rin,rout,evib,erotx,rinner,
     & ediatot,kerel,kedia,xjx,xjy,xjz,eint,mmtot,rcom,xm1,xm2,
     & td0,tdf,rm,rmm,rxx,bb,dot,d1,d2

c     temp
      integer wellindex1,wellindex2
      double precision xwell,ptmp

c     for unit 31
      integer ifrag,nfrag(2),fragind(mnat,2),itmp,itmp1,itmp2,nh
      double precision rr1,rr2,ppfrag(3,mnat),xxfrag(3,mnat),
     &  mmfrag(mnat),bigjag(3,2),bigjtotag(2),comxag(3,2),compag(3,2),
     &  relmu,kerelfrag(2),jfragcom(2),erotag(3,2),erottotag(2),
     &  jfragrot(3,2),jfragrottot,rf,comx(3),comp(3),lorb(3),lorbtot,
     &  eorb,erel,tmp3j(3),tmp3e(3),mmfragtot(2),tempfrag,kefrag(2),
     &  evibag(2),tdelay

      double precision tmpprint(50)
      common/tmp/tmpprint

      save rmm,rm

c HARDCODED
      nfrag(1) = 4
      nfrag(2) = 4
c      nh=0
c      do i=1,nfrag(1)
c      if (symbol(i).eq."C".or.symbol(i).eq."c") nh=nh+1
c      enddo
c      rinner=3.1952d0*nh**0.2453d0
c      rinner=rinner/autoang
      rinner=0.d0

C FOR DISSOCIATION REACTIONS   ! for bz + bz AJ
      xm1=0.d0
      xm2=0.d0
      do i=1,nfrag(1)
        xm1=xm1+mm(i)
      enddo
      do i=nfrag(1)+1,nfrag(1)+nfrag(2)
        xm2=xm2+mm(i)
      enddo
c      mmtot = mmag(1)+mmag(2)
c      mmfragtot(1) = mmag(1)
c      mmfragtot(2) = mmag(2)
      mmtot = xm1+xm2
      mmfragtot(1) = xm1
      mmfragtot(2) = xm2
      do i=1,3
      comx(i)=0.d0
      comp(i)=0.d0
      comxag(i,1)=0.d0
      comxag(i,2)=0.d0
      compag(i,1)=0.d0
      compag(i,2)=0.d0
      enddo

      ifrag=1
      do i=1,nfrag(1)
        do k=1,3
          comx(k)=comx(k) + xx(k,i)*mm(i)                  ! center of mass
          comp(k)=comp(k) + pp(k,i)                        ! overall momentum
          comxag(k,ifrag)=comxag(k,ifrag) + xx(k,i)*mm(i)  ! center of mass for fragment
          compag(k,ifrag)=compag(k,ifrag) + pp(k,i)        ! overall momentum for fragment
        enddo
      enddo
      do k=1,3
        comx(k)=comx(k)/mmtot                                 ! overall center of mass
      enddo

      ifrag=2
      do i=nfrag(1)+1,nfrag(1)+nfrag(2)
        do k=1,3
          comx(k)=comx(k) + xx(k,i)*mm(i)                  ! center of mass
          comp(k)=comp(k) + pp(k,i)                        ! overall momentum
          comxag(k,ifrag)=comxag(k,ifrag) + xx(k,i)*mm(i)  ! center of mass for fragment
          compag(k,ifrag)=compag(k,ifrag) + pp(k,i)        ! overall momentum for fragment
        enddo
      enddo

      do k=1,3
        comx(k)=comx(k)/mmtot                                 ! overall center of mass
      enddo

      do j=1,2
      do k=1,3
        comxag(k,j)=comxag(k,j)/mmfragtot(j)-comx(k)          ! center of mass for fragment
        compag(k,j)=compag(k,j)-comp(k)/mmtot*mmfragtot(j)    ! overall momentum for fragment
      enddo
      enddo

c     relative reduced mass
      relmu = mmfragtot(1)*mmfragtot(2)/(mmfragtot(1)+mmfragtot(2))
      rf = 0.
      do i=1,3
        rf = rf + (comxag(i,1)-comxag(i,2))**2
      enddo
      rf = dsqrt(rf)  ! final fragment separation
      lorb(1) = (comxag(2,1)-comxag(2,2))*compag(3,1)
     &         -(comxag(3,1)-comxag(3,2))*compag(2,1)
      lorb(2) = (comxag(3,1)-comxag(3,2))*compag(1,1)
     &         -(comxag(1,1)-comxag(1,2))*compag(3,1)
      lorb(3) = (comxag(1,1)-comxag(1,2))*compag(2,1)
     &         -(comxag(2,1)-comxag(2,2))*compag(1,1)
      lorbtot = dsqrt(lorb(1)**2+lorb(2)**2+lorb(3)**2)
      eorb = lorbtot**2/(2*relmu*rf**2)
c      write(6,*)'XX 1',(comxag(i,1),i=1,3)
c      write(6,*)'XX 2',(comxag(i,2),i=1,3)
c      write(6,*)'XX 1',(compag(i,1),i=1,3)
c      write(6,*)'XX 2',(compag(i,2),i=1,3)
      eint = (compag(1,1)**2+compag(2,1)**2+compag(3,1)**2)/(2.d0*relmu)
c      erel = eint - eorb + tmpprint(11)
      erel = eint - eorb

      rcom=rf

c angle
      dot=0.d0
      d1=0.d0
      d2=0.d0
      do i=1,3
      dot=dot+(compag(i,1)-compag(i,2))*(comxag(i,1)-comxag(i,2))
      d1=d1+(compag(i,1)-compag(i,2))**2
      d2=d2+(comxag(i,1)-comxag(i,2))**2
      enddo
      dot=dot/dsqrt(d1*d2)  ! cos theta
      dot=dsqrt(1.d0-dot**2)  ! sin theta
      bb=dot*rcom
c      write(6,*)bb*autoang,rcom*autoang,bqci*autoang

c bounces
      if (istep.eq.0) nbounce=0
      if (istep.gt.2) then
      if (rm.lt.rf.and.rm.lt.rmm) nbounce=nbounce+1
      endif
c      write(6,*)istep,rmm*autoang,rm*autoang,rf*autoang,nbounce
      rmm=rm
      rm=rf

c t_delay
c      rxx = dsqrt(max(rcom**2-bqci**2,0.d0))
      rxx = dsqrt(max(rcom**2-bb**2,0.d0))
      if (istep.eq.0) td0 = (rxx-rinner)/dsqrt(erel*2.d0/relmu)
      tdf = (rxx-rinner)/dsqrt(erel*2.d0/relmu)
      tdelay = time - td0 - tdf

c      write(6,*)'XX',relmu,rinner,rcom,eorb,erel,tdelay,td0,tdf

c      if (mod(istep,nprint).eq.0.and.lwrite(95)) 
c     &  write(95,*)rcom*autoang,tdelay*autofs,nbounce
      return
      end
