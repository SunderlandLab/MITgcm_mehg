C $Header: /u/gcmpack/MITgcm/pkg/hg/hg_bioaccum.F,v 1.10 2010/03/16 00:14:47 jmc Exp $
C $Name: checkpoint64g $

C To calculate the chemistry of Hg species in the ocean, yxzhang, 5/30/13

#include "GCHEM_OPTIONS.h"
#include "HG_OPTIONS.h"

CStartOfInterFace
      SUBROUTINE HG_BIOACCUM(  Ptr,
     &                         bi,bj,imin,imax,jmin,jmax,
     &                         myIter,myTime,myThid)

C     /==========================================================\
C     | SUBROUTINE HG_FOODWEB                                    |
C     | o Calculate the changes to Hg species because            |
C     |   of foodweb dynamics                                    |
C     |   based on Asif Quesi 2013 paper                         |
C     |==========================================================|
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "HG_SIZE.h"
#include "HG.h"

C     == Routine arguments ==
      INTEGER myIter
      _RL myTime
      INTEGER myThid
      _RL  Ptr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      INTEGER bi, bj, imin, imax, jmin, jmax
#ifdef FOODW_MODE
C     == Local variables ==            
      INTEGER i,j,k
      INTEGER np,nz
      _RL total_mmhg
      _RL vcfo(npmax)  !volume concentration factor, m3/m3 
      _RL fphyo(npmax) !fraction of MMhg on phytoplankton
      _RL total_fphyo, phytopo(npmax), phygrowo(npmax), doco
      _RL cwo, cwnewo    !MMHg concentration
      _RL phymmhgo(npmax), growo(npmax)
      INTEGER active(npmax)     !active species
      _RL tiny
      PARAMETER(tiny=1. _d -32)
      
C     == Phytoplankton variables
      _RL dm(npmax)   !diameter of phytoplankton, um
      _RL svr         !surface area to volume ratio, um-1, 3/r or 6/d
      _RL fcarbon(npmax)   !carbon content of wet phytoplankton
      DATA dm /1.2 _d 1,1.0 _d 1,1.8 _d 0,6.0 _d -1,1.0 _d 1,5.0 _d 0/
C      DATA fcarbon /5.0 _d -2/  ! 90% water, Hammerschmidt et al. 2006
C                                ! 50% carbon over dry weight
      DATA fcarbon /0.1236, 0.1357,0.3272,0.5749,0.1357,0.1937/ 
                                 !Moal et al. 1987; y=0.396x^0.487
C     == Partitioning adjusting time scale
      _RL t_bioaccum(npmax)
      DATA t_bioaccum / 4.0 _d 0,3.3 _d 0,6.0 _d -1,2.0 _d -1,3.3 _d 0,1.7 _d 0 /  ! dm/3 hour, dm: um Kim 2014

#ifdef ALLOW_DIAGNOSTICS      
      CHARACTER*8   diagName
      _RL  mmhgbc_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,npmax)  ! MMHG bioaccumulation from sea water
#endif

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_HG

C     == Partition MMHg onto phytoplankton
      DO j=jMin,jMax
       DO i=iMin,iMax
        DO k=1,Nr
         IF (HFacC(i,j,k,bi,bj).gt.0. _d 0) THEN
           ! Localize variables
           DO np=1,npmax
               phymmhgo(np) = max(Ptr(i,j,k,bi,bj,iPhy+np-1),0. _d 0) !mol/m3
               phytopo(np) = max(phytoP(i,j,k,bi,bj,np),0. _d 0) !mmol C/m3
               phygrowo(np) = max(phygrow(i,j,k,bi,bj,np),0. _d 0) !mmol C/m3/s
               IF ( phytopo(np) .gt. tiny ) THEN
                   active(np) = 1  ! active means update
                   growo(np) = phygrowo(np) / phytopo(np)  !s-1
               ELSE
                   active(np) = 0  ! non-active means keeps c not change
               ENDIF
           ENDDO
c
           doco = max(doc(i,j,k,bi,bj),0. _d 0)  ! mmol C/m3 -> umol C/L (uM)
           DO np=1,npmax
            IF ( active(np) .eq. 1 ) THEN
#ifdef BC_Schartup2017
             !vcfo(np) = 2.76 _d 6 / dm(np) * dexp( -0.01 * doco )
             vcfo(np) = 2.832 _d 6 / dm(np) * dexp( -0.008 * doco )
#else
             vcfo(np) = 1.8959 _d 6 * 6.0 _d 0 / dm(np) ! - 1.5823 _d 6  ! unitless, Lee 2015, 0 intercept
#ifdef BC_Mason1996
     &       * ( 1. _d 0 - dexp( -1. _d 0 / ( growo(np) 
     &       * t_bioaccum(np) ) ) )  ! biodilution, 1/grow as mean age
#endif
#ifdef BC_Hammerschmidt2013
     &       / 5.0 _d 1 ! corrected by LIS & NW Atlantic ocean measurements, logBCFdiatom=4.3
#endif
     &       * dexp( - 5. _d -3 * doco )
#endif 
            ENDIF
           ENDDO
           cwo = max(Ptr(i,j,k,bi,bj,immhg),0 _d 0) ! MMHg seawater concentration, mol/m3
c
           ! calculate partitioning
           total_mmhg = cwo          
           total_fphyo = 1. _d 0
           DO np=1,npmax
             IF ( active(np) .eq. 1 ) THEN
               fphyo(np) = vcfo(np) * phytopo(np)
     &                     * 1.2 _d -8 / fcarbon(np)  !1e-3 * 12 *1e-3 /fcarbon / (rou = 1000)     
                                                      ! mmol C / m3 -> m3 cell/ m3 water
               IF ( phymmhgo(np) .ge. ( cwo * fphyo(np) ) ) THEN
                  active(np) = 0    ! phytoplankton not losing their MMHg
               ELSE
                  total_mmhg = total_mmhg + phymmhgo(np)
                  total_fphyo = total_fphyo + fphyo(np)
               ENDIF
             ENDIF
           ENDDO

           ! update ptracers
           IF ( total_mmhg .gt. tiny ) THEN
c             dynamic update approach
c             cwnewo = cwo + 
c     &        (total_mmhg / total_fphyo - cwo)
c     &        * ( 1. _d 0 - dexp( - deltaTclock / 
c     &          t_bioaccum(np) ) ) 
c             Ptr(i,j,k,bi,bj,immhg) = cwnewo
c             cwnewo = total_mmhg
             DO np=1,npmax
              IF ( active(np) .eq. 1 ) THEN
c               dynamic update approach
c               Ptr(i,j,k,bi,bj,iPhy+np-1) =
c     &          ( total_mmhg - cwnewo ) * 
c     &          fphyo(np) / max( tiny, total_fphyo - 1. _d 0 )
c               Ptr(i,j,k,bi,bj,iPhy+np-1) =
c     &          phymmhgo(np) +
c     &          ( total_mmhg * fphyo(np) / total_fphyo -
c     &            phymmhgo(np) ) *
c     &          ( 1. _d 0 - dexp( - deltaTclock / 
c     &           t_bioaccum(np) ) )
c               cwnewo = cwnewo -
c     &           Ptr(i,j,k,bi,bj,iPhy+np-1)
c               equilibrium approach
                Ptr(i,j,k,bi,bj,iPhy+np-1) =
     &            total_mmhg / total_fphyo *
     &            fphyo(np)           
              ENDIF
             ENDDO
             Ptr(i,j,k,bi,bj,immhg) = total_mmhg / total_fphyo
           ENDIF

#ifdef ALLOW_DIAGNOSTICS
           DO np = 1,npmax
            IF ( active(np) .eq. 1 ) THEN
             mmhgbc_diag(i,j,k,np) =
     &           ( Ptr(i,j,k,bi,bj,iPhy+np-1) - phymmhgo(np) )
     &           / deltaTclock
            ELSE
             mmhgbc_diag(i,j,k,np) = 0. _d 0
            ENDIF
           ENDDO
#endif
C
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      
c write diagnostics
#ifdef ALLOW_DIAGNOSTICS
        DO np=1,npmax
         WRITE(diagName,'(A6,I2.2)') 'MMHGBC',np
         CALL DIAGNOSTICS_FILL( mmhgbc_diag(1-Olx,1-Oly,1,np), diagname,
     &                          0,Nr,2,bi,bj,myThid )
        ENDDO
#endif

#endif
#endif
#endif
       RETURN
       END
