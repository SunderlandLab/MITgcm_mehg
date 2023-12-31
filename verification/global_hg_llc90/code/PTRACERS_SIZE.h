C $Header: /u/gcmpack/MITgcm/verification/tutorial_cfc_offline/code/PTRACERS_SIZE.h,v 1.1 2005/10/01 02:28:30 edhill Exp $
C $Name: checkpoint64g $

#ifdef ALLOW_PTRACERS

CBOP
C    !ROUTINE: PTRACERS_SIZE.h
C    !INTERFACE:
C #include PTRACERS_SIZE.h
 
C    !DESCRIPTION:
C Contains passive tracer array size (number of tracers).

C PTRACERS_num defines how many passive tracers are allocated/exist.
C  and is set here (default 1)
C
C     Number of tracers
      INTEGER PTRACERS_num
      PARAMETER(PTRACERS_num = 14 )

CEOP

#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER    maxpass
      PARAMETER( maxpass     = PTRACERS_num + 2 )
#endif

#endif /* ALLOW_PTRACERS */
