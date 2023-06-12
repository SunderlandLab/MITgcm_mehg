C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C                                                                  
C Parameter Header File                                            
C                                                                  
C Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  
C       (http://www.cs.vt.edu/~asandu/Software/KPP)                
C KPP is distributed under GPL, the general public licence         
C       (http://www.gnu.org/copyleft/gpl.html)                     
C (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           
C (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            
C     With important contributions from:                           
C        M. Damian, Villanova University, USA                      
C        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
C                                                                  
C File                 : chem_hg_Parameters.h                      
C Time                 : Thu Mar 27 11:36:36 2014                  
C Working directory    : /n/home05/yxzhang/kpp/kpp-2.2.3/mehg      
C Equation file        : chem_hg.kpp                               
C Output root filename : chem_hg                                   
C                                                                  
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




C NSPEC - Number of chemical species                               
      INTEGER NSPEC
      PARAMETER ( NSPEC = 14 )
C NVAR - Number of Variable species                                
      INTEGER NVAR
      PARAMETER ( NVAR = 14 )
C NVARACT - Number of Active species                               
      INTEGER NVARACT
      PARAMETER ( NVARACT = 4 )
C NFIX - Number of Fixed species                                   
      INTEGER NFIX
      PARAMETER ( NFIX = 1 )
C NREACT - Number of reactions                                     
      INTEGER NREACT
      PARAMETER ( NREACT = 10 )
C NVARST - Starting of variables in conc. vect.                    
      INTEGER NVARST
      PARAMETER ( NVARST = 1 )
C NFIXST - Starting of fixed in conc. vect.                        
      INTEGER NFIXST
      PARAMETER ( NFIXST = 15 )
C NONZERO - Number of nonzero entries in Jacobian                  
      INTEGER NONZERO
      PARAMETER ( NONZERO = 31 )
C LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
      INTEGER LU_NONZERO
      PARAMETER ( LU_NONZERO = 31 )
C CNVAR - (NVAR+1) Number of elements in compressed row format     
      INTEGER CNVAR
      PARAMETER ( CNVAR = 15 )
C CNEQN - (NREACT+1) Number stoicm elements in compressed col format
      INTEGER CNEQN
      PARAMETER ( CNEQN = 11 )
C NHESS - Length of Sparse Hessian                                 
      INTEGER NHESS
      PARAMETER ( NHESS = 1 )
C NLOOKAT - Number of species to look at                           
      INTEGER NLOOKAT
      PARAMETER ( NLOOKAT = 0 )
C NMONITOR - Number of species to monitor                          
      INTEGER NMONITOR
      PARAMETER ( NMONITOR = 0 )
C NMASS - Number of atoms to check mass balance                    
      INTEGER NMASS
      PARAMETER ( NMASS = 1 )

C Index declaration for variable species in C and VAR              
C   VAR(ind_spc) = C(ind_spc)                                      

      INTEGER ind_RR1
      PARAMETER ( ind_RR1 = 1 )
      INTEGER ind_RR2
      PARAMETER ( ind_RR2 = 2 )
      INTEGER ind_RR3
      PARAMETER ( ind_RR3 = 3 )
      INTEGER ind_RR4
      PARAMETER ( ind_RR4 = 4 )
      INTEGER ind_RR5
      PARAMETER ( ind_RR5 = 5 )
      INTEGER ind_RR6
      PARAMETER ( ind_RR6 = 6 )
      INTEGER ind_RR7
      PARAMETER ( ind_RR7 = 7 )
      INTEGER ind_RR8
      PARAMETER ( ind_RR8 = 8 )
      INTEGER ind_RR9
      PARAMETER ( ind_RR9 = 9 )
      INTEGER ind_RR10
      PARAMETER ( ind_RR10 = 10 )
      INTEGER ind_DMHg
      PARAMETER ( ind_DMHg = 11 )
      INTEGER ind_MMHg
      PARAMETER ( ind_MMHg = 12 )
      INTEGER ind_Hg0_aq
      PARAMETER ( ind_Hg0_aq = 13 )
      INTEGER ind_Hg2_aq
      PARAMETER ( ind_Hg2_aq = 14 )

C Index declaration for fixed species in C                         
C   C(ind_spc)                                                     


C Index declaration for dummy species                              


C Index declaration for fixed species in FIX                       
C    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)                 


C NJVRP - Length of sparse Jacobian JVRP                           
      INTEGER NJVRP
      PARAMETER ( NJVRP = 10 )

C NSTOICM - Length of Sparse Stoichiometric Matrix                 
      INTEGER NSTOICM
      PARAMETER ( NSTOICM = 30 )
