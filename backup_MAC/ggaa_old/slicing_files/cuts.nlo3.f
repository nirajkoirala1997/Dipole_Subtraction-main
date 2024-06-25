c
c      subroutine cuts3(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,
c     &     r34,r35,r45,Et5,QT,ipass)      
c      implicit double precision (a-h,o-z)
c
c      common/machine/mid
c      common/cone/ET_iso,r0,rgg
c      common/nviso/niso
c      
c	mid = 1
c
c      ipass=0
c
c      IF(mid.eq.0)THEN         
cc         r0=0.4d0
cc         rgg=0.3d0
cc         ET_iso=1.0d0
c         
c         if( (abs(y1).ge.0.9d0 .or. abs(y2).ge.0.9d0) ) return          
c         if( (pt1.lt.14d0 .or. pt2.lt.13d0) 
c     &        .and. (pt1.lt.13d0 .or. pt2.lt.14d0) )return
c         
c         if(r34.lt.rgg)return         
c         if(r35.gt.r0 .and. r45.gt.r0)then
c            ipass=1
c            return
c         endif                  
c
c         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
c         IF(r45.le.r0 .and. Et5.gt.ET_iso) return                 
c         
c
c        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
c        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return
c               
c         
c         ipass=1               
c         return
c
c      ELSEIF(mid.eq.1)THEN
cc         r0=0.4d0
cc         rgg=0.4d0
cc         ET_iso=15d0
c
c         if( (abs(y1).ge.2.50d0 .or. abs(y2).ge.2.50d0) ) return
c         if( (pt1.lt.40d0 .or. pt2.lt.25d0) 
c     &        .and. (pt1.lt.25d0 .or. pt2.lt.40d0) )return
c         
c         if(r34.lt.rgg)return         
c         if(r35.gt.r0 .and. r45.gt.r0)then
c            ipass=1
c            return
c         endif            
c
c         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
c         IF(r45.le.r0 .and. Et5.gt.ET_iso) return                 
c         
c
c        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
c        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return
c        
c        ipass=1               
c         return
c      ENDIF
c      end
      

      subroutine cuts3(p1,p2,p3,p4,p5,ipass)
      implicit double precision (a-h,o-z)

      common/machine/mid
      common/cone/ET_iso,r0,rgg
      common/nviso/niso
      common/amf/am3,am4,am5

      ipass=0

c     yy1: rapidity of photons 3 
      y1 = rapid(p3)

c     yy2: rapiditiy of photon 4
      y2 = rapid(p4)

c     ppt1 ppt2
      pt1 = eperp(p3)
      pt2 = eperp(p4)

c     cone
      r34 = rjet(p3,p4)
      r35 = rjet(p3,p5)
      r45 = rjet(p4,p5)

c     transverse energy
      Et5 = eperp(p5)  ! for a massless particle Et=pt

      IF(mid.eq.0)THEN

         if( (abs(y1).ge.0.9d0 .or. abs(y2).ge.0.9d0) ) return
         if( (pt1.lt.14d0 .or. pt2.lt.13d0)
     &        .and. (pt1.lt.13d0 .or. pt2.lt.14d0) )return

         if(r34.lt.rgg)return
         if(r35.gt.r0 .and. r45.gt.r0)then
            ipass=1
            return
         endif

         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
         IF(r45.le.r0 .and. Et5.gt.ET_iso) return


        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return


         ipass=1
         return

      ELSEIF(mid.eq.1)THEN
c         r0=0.4d0
c         rgg=0.4d0
c         ET_iso=15d0

         if( (abs(y1).ge.2.50d0 .or. abs(y2).ge.2.50d0) ) return
         if( (pt1.lt.40d0 .or. pt2.lt.25d0)
     &        .and. (pt1.lt.25d0 .or. pt2.lt.40d0) )return

         if(r34.lt.rgg)return
         if(r35.gt.r0 .and. r45.gt.r0)then
            ipass=1
            return
         endif

         IF(r35.le.r0 .and. Et5.gt.ET_iso) return
         IF(r45.le.r0 .and. Et5.gt.ET_iso) return


        IF(r35.le.r0 .and. Et5.gt.hh(r35)) return
        IF(r45.le.r0 .and. Et5.gt.hh(r45)) return

        ipass=1
         return
      ENDIF
      end
