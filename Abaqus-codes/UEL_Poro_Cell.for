************************************************************************
! Plane strain user element for large deformation, fluid pressure, and astrocyte
************************************************************************
! Author: Fatemeh Jafarabadi, Apr 2023, last edited June 2026
! 
! This subroutine is for plane strain elements
! 4-node element with 4pt Gauss integration
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
***********************************************************************
! User element statement in the input file: 
! 
! 2D elements
!  *User Element,Nodes=4,Type=U1,Iproperties=4,Properties=20,Coordinates=2,
!  Variables=<nVars>,Unsymm, 1,2,11,12
!***********************************************************************
      module global
      !  globalSdv(X,Y,Z)
      !  X - element pointer
      !  Y - integration point pointer
      !  Z - SDV pointer
      !
      !  numElem: Total number of elements in the real mesh
      !  ElemOffset: Offset between element num on real mesh and dummy mesh
      IMPLICIT NONE
      
      integer numElem,elemOffset,err

      parameter(numElem=18781) 
      parameter(ElemOffset=100000) 

      real*8, allocatable :: globalSdv(:,:,:)

      end module global
***********************************************************************
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used for visualization on dummy mesh 

      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      UVAR(1) = globalSdv(NOEL-elemOffset,NPT,1)   ! swelling volume change 
      UVAR(2) = globalSdv(NOEL-elemOffset,NPT,2)   ! total volume change 
      UVAR(3) = globalSdv(NOEL-elemOffset,NPT,3)   ! elastic volume change 
      UVAR(4) = globalSdv(NOEL-elemOffset,NPT,4)   ! x-comp of Darcy's velocity
      UVAR(5) = globalSdv(NOEL-elemOffset,NPT,5)   ! y-comp of Darcy's velocity
      UVAR(6) = globalSdv(NOEL-elemOffset,NPT,6)   ! magnitude of Darcy's velocity
      UVAR(7) = globalSdv(NOEL-elemOffset,NPT,7)   ! x-comp of astrocyte flux      
      UVAR(8) = globalSdv(NOEL-elemOffset,NPT,8)   ! y-comp of astrocyte flux      
      UVAR(9) = globalSdv(NOEL-elemOffset,NPT,9)   ! astrocyte source      
      UVAR(10) = globalSdv(NOEL-elemOffset,NPT,10) ! astrocyte sink              
      UVAR(11) = globalSdv(NOEL-elemOffset,NPT,11) ! pressure source      
      UVAR(12) = globalSdv(NOEL-elemOffset,NPT,12) ! total cytokine concentration

      return
      end subroutine UVARM
****************************************************************************
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,COORDS,MCRD,NNODE,UALL,DUALL,VEL,ACCN,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      ! This subroutine calls the main element subroutine (UPE4)

      use global

      IMPLICIT NONE

      ! Variables defined in UEL, passed back to Abaqus
      real(8)  RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),ENERGY(8)

      ! Variables passed to UEL
      integer NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC
      integer JELEM,NDLOAD,JDLTYP(MDLOAD,*),NPREDF,LFLAGS(*),MLVARX
      integer NJPROP,JPROPS(NJPROP),MDLOAD

      real(8) COORDS(MCRD,NNODE),UALL(NDOFEL),DUALL(MLVARX,*)
      real(8) VEL(NDOFEL),ACCN(NDOFEL)
      real(8) TIME(2),PARAMS(*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*)
      real(8) PREDEF(2,NPREDF,NNODE),PROPS(*),DTIME,PNEWDT,PERIOD

      integer nDim,nInt
      !  
      !----------------------------------------------------------------
      ! Initial checks
      !----------------------------------------------------------------
      ! Do nothing if a ``dummy'' step
      !
      if(DTIME.eq.0.0) return
      !----------------------------------------------------------------
      ! Call the element subroutine 
      !----------------------------------------------------------------
      nInt = JPROPS(3)  ! number of integration points
      nDim = JPROPS(4)  ! dimension
      !
      call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,COORDS,MCRD,NNODE,UALL,DUALL,VEL,ACCN,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,nDim,nInt)
      !
      ! Done with this element, RHS and AMATRX passed to Abaqus
      !----------------------------------------------------------------
      return
      end subroutine UEL
************************************************************************ 
!                         2D Element Subroutine
************************************************************************
      subroutine UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,COORDS,MCRD,NNODE,UALL,DUALL,VEL,ACCN,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,nDim,nInt)

      ! This subroutine performs the main element level calculations
      ! In this subroutine, "_t" denotes a quantity evaluated at 
      ! the start of the increment, and "_tau" denotes the same 
      ! quantity evaluated at the end of the increment. 

      use global

      IMPLICIT NONE

      ! Variables defined in UEL, passed back to Abaqus
      real(8)  RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),ENERGY(8)

      ! Variables passed to UEL
      integer NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC
      integer JELEM,NDLOAD,JDLTYP(MDLOAD,*),NPREDF,LFLAGS(*),MLVARX
      integer NJPROP,JPROPS(NJPROP),MDLOAD

      real(8) COORDS(MCRD,NNODE),UALL(NDOFEL),DUALL(MLVARX,*)
      real(8) VEL(NDOFEL),ACCN(NDOFEL)
      real(8) TIME(2),PARAMS(*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*)
      real(8) PREDEF(2,NPREDF,NNODE),PROPS(*),DTIME,PNEWDT,PERIOD

      integer i,j,k,nIntPt,nDim,intPt
      integer nInt,jj,stat,ngSdv
      integer nlSdv,kk

      real*8 u(nNode,nDim),du(nNode,NDOFEL)
      real*8 U_t(nDim,1),U_tau(nDim,1),dudt(nNode,nDim)
      real*8 uOld(nNode,NDOFEL),div_v
      real*8 coordsC(mcrd,nNode),coordint(nDim,1)

      real*8 POld(nNode),dP(nNode),PNew(nNode),velocity_D(nDim,1)
      real*8 P_t,P_tau,dPdX(nDim,1),dPdt,dP_dt(nNode),alpha_f
      real*8 source_P,K_h,flux_D,biot,K_f,Sf0
      
      real*8 aNew(nNode),da(nNode),aOld(nNode)
      real*8 a_tau,a_t,dadX(nDim,1),dadt,da_dt(nNode)
      real*8 qaflux(nDim,1),source_a,sink_a
      real*8 delta_a,Diff_a,dqada_adv(nDim,1)

      real*8 dshxi(nNode,nDim),xi(nInt,nDim),dshC(nNode,nDim)
      real*8 sh0(nNode),detMapJ0,dsh0(nNode,nDim),dshC0(nNode,nDim)
      real*8 Del_N(nDim,nNode),detMapJ0C,Nvec(1,nNode),Iden(3,3)
      real*8 sh(nNode),detMapJ,dsh(nNode,nDim),detMapJC
      
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),Je_tau
      real*8 F_tau(3,3),F_t(3,3),detF_tau,detF,detF_t,detFs,dJdt,J_tau
      real*8 c_total

      real*8 Ru(nDim*nNode,1),Rp(nNode,1),Ra(nNode,1)
      real*8 Kuu(nDim*nNode,nDim*nNode),Kup(nDim*nNode,nNode),Kua(nDim*nNode,nNode)
      real*8 Kpu(nNode,nDim*nNode),Kpp(nNode,nNode),Kpa(nNode,nNode)
      real*8 Kau(nNode,nDim*nNode),Kap(nNode,nNode),Kaa(nNode,nNode)

      real*8 T_tau_s(3,3)
      real*8 Smat(3,1),Bmat(3,nDim*nNode),Gmat(4,nDim*nNode)
      real*8 Qmat(4,4),G0mat(4,nDim*nNode)
      real*8 Amat(4,4),AmatUP(3,1),AmatUA(3,1),AmatAU(2,4)
      real*8 D_uu(3,3,3,3),D_up(3,3),D_ua(3,3),D_au(3,3,3)
      real*8 zero,one,half
      parameter(zero=0.d0,one=1.d0,half=1.d0/2.d0)

      ! Element parameters
      nlSdv = JPROPS(1) ! number of local sdvs per integ point
      ngSdv = JPROPS(2) ! number of global sdvs per integ point
      
      ! Allocate memory for the globalSdvs
      if (.not. allocated(globalSdv)) then
         err = 0
         allocate(globalSdv(numElem,nInt,ngSdv), stat=err)
         if (err .ne. 0) then
            write(*,*) '-------------------------------------------------'
            write(*,*) 'error allocating globalSdv'
            write(*,*) '-------------------------------------------------'
            write(*,*) '   err=',err
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- ElementType = UPE4'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- numElem =',numElem
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- nDim =',nDim
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- ngSdv =',ngSdv
         write(*,*) '-------------------------------------------------'
      endif

      call ONEM(Iden)
      
      !----------------------------------------------------------------
      ! Obtain parameters
      !----------------------------------------------------------------
      biot          = PROPS(3)
      K_h           = PROPS(4) 
      K_f           = PROPS(5)
      Sf0           = PROPS(6)
      delta_a       = PROPS(8)
      Diff_a        = PROPS(9)  
      alpha_f       = PROPS(12) 

      !----------------------------------------------------------------
      ! Initialize the residual and tangent matrices to zero
      !----------------------------------------------------------------
      Ru  = zero; Rp  = zero; Ra  = zero
      Kuu = zero; Kup = zero; Kua = zero
      Kpu = zero; Kpp = zero; Kpa = zero
      Kau = zero; Kap = zero; Kaa = zero
      
      !----------------------------------------------------------------
      ! Obtain nodal displacements, pressure, and astrocyte density 
      !----------------------------------------------------------------
      u = zero; du = zero; uOld = zero; coordsC = zero
      PNew = zero; dP = zero; POld = zero; dP_dt = zero
      aNew = zero; da = zero; aOld = zero; da_dt = zero

      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = UALL(k)
            du(i,j) = DUALL(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
            ! Obtain current nodal coordinates
            coordsC(j,i) = COORDS(j,i) + u(i,j)
         enddo
         ! 
         k = k + 1
         PNew(i) = UALL(k)
         dP(i)   = DUALL(k,1)
         POld(i) = PNew(i) - dP(i)
         dP_dt(i) = VEL(k)
         ! 
         k = k + 1
         aNew(i) = UALL(k)
         da(i) = DUALL(k,1)
         aOld(i) = aNew(i) - da(i) 
         da_dt(i) = VEL(k)
      enddo

      !----------------------------------------------------------------
      ! Shape functions
      !----------------------------------------------------------------
      ! At the element centroid
      call calcShape2DLinear(one,zero,one,sh0,dshxi)

      ! Map shape functions from local to global reference coordinate system
      call mapShape2Da(nNode,dshxi,COORDS,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif

      ! Map shape functions from local to global current coordinate system
      call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif
      !----------------------------------------------------------------
      ! F-bar method
      !----------------------------------------------------------------
      ! Calculate the deformation gradient at the element centroid
      Fc_tau = Iden; Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! 2D plane-strain implementation detF
      detFc_t = Fc_t(1,1)*Fc_t(2,2) - Fc_t(1,2)*Fc_t(2,1)
      detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
      !
      !----------------------------------------------------------------
      ! Calculations at the integration points
      !----------------------------------------------------------------
      ! Integration points
      call xint2D4pt(xi,w,nIntPt) 
      !
      ! Loop over integration points
      jj = 0 ! jj is used for tracking the state variables
      do intPt=1,nIntPt
         !
         ! Shape functions and their local gradients
         call calcShape2DLinear(nIntPt,xi,intPt,sh,dshxi)
         !
         ! Map shape functions from local to global reference coordinate system
         call mapShape2Da(nNode,dshxi,COORDS,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         ! Map shape functions from local to global current coordinate system
         call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         ! Pressure & Astrocyte and their derivatives at this intPt
         P_tau = zero; P_t = zero; dPdX = zero; dPdt = zero
         a_tau = zero; a_t = zero; dadX = zero; dadt = zero

         do k=1,nNode
            ! Pressure calculations
            P_tau = P_tau + PNew(k)*sh(k)
            P_t   = P_t + POld(k)*sh(k)
            dPdt  = dPdt + dP_dt(k)*sh(k) 
            ! Astrocyte density  calculations
            a_tau = a_tau + aNew(k)*sh(k)
            a_t   = a_t + aOld(k)*sh(k)
            dadt  = dadt + da_dt(k)*sh(k)             
            do i=1,nDim
               ! Pressure gradient
               dPdX(i,1) = dPdX(i,1) + PNew(k)*dshC(k,i)
               ! Astrocytes gradient
               dadX(i,1) = dadX(i,1) + aNew(k)*dshC(k,i)
            enddo
         enddo

         ! Displacement and its time derivatives at this integration point
         U_tau = zero; U_t = zero; dudt = zero; div_v = zero
         do j=1,nDim
            do i=1,nNode
               U_tau(j,1) = U_tau(j,1) + sh(i)*u(i,j)
               U_t(j,1)   = U_t(j,1) + sh(i)*uOld(i,j)
               dudt(i,j) = (u(i,j)-uOld(i,j))/DTIME
               div_v = div_v + dshC(i,j) * dudt(i,j)
            enddo      
         enddo
         !
         ! Deformation gradient at this integration point
         F_tau = Iden; F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! 2D plane-strain implementation
            detF_t = F_t(1,1)*F_t(2,2) - F_t(1,2)*F_t(2,1)
            detF_tau = F_tau(1,1)*F_tau(2,2) - F_tau(1,2)*F_tau(2,1)
            do i=1,nDim
               do j=1,nDim
                  F_tau(i,j) =((detFc_tau/detF_tau)**half)*F_tau(i,j)
                  F_t(i,j) = ((detFc_t/detF_t)**half)*F_t(i,j)
               enddo
            enddo
         endif
         !
         ! Calculate the deformation determinant after applying the F-bar method
         call MDET(F_tau,detF)
         !----------------------------------------------------------------
         Nvec = zero
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo 
         !
         coordint  = matmul(COORDS,Nvec)  
         !----------------------------------------------------------------
         ! Calculate cytokine concentration 
         !----------------------------------------------------------------
         call cytokineFunction(PROPS,NPROPS,TIME,nDim,coordint,c_total)
         
         !----------------------------------------------------------------
         ! Perform the constitutive time integration at this integration point 
         !----------------------------------------------------------------
         call timeIntegral(PROPS,NPROPS,nDim,DTIME,coordint,
     +              T_tau_s,F_tau,F_t,P_tau,dPdX,detFs,dJdt,
     +              J_tau,Je_tau,velocity_D,source_P,flux_D,
     +              a_tau,a_t,dadX,source_a,sink_a,qaflux,
     +              dqada_adv,c_total,D_uu,D_up,D_ua,D_au)      

         !----------------------------------------------------------------
         ! Save the state variables at this integration point
         !----------------------------------------------------------------
         jj = jj + nlSdv

         globalSdv(JELEM,intPt,1) = detFs   
         globalSdv(JELEM,intPt,2) = J_tau
         globalSdv(JELEM,intPt,3) = Je_tau 
         globalSdv(JELEM,intPt,4) = velocity_D(1,1) 
         globalSdv(JELEM,intPt,5) = velocity_D(2,1)       
         globalSdv(JELEM,intPt,6) = flux_D    
         globalSdv(JELEM,intPt,7) = qaflux(1,1)     
         globalSdv(JELEM,intPt,8) = qaflux(2,1)          
         globalSdv(JELEM,intPt,9) = source_a
         globalSdv(JELEM,intPt,10) = sink_a
         globalSdv(JELEM,intPt,11) = source_P
         globalSdv(JELEM,intPt,12) = c_total

         !----------------------------------------------------------------
         ! Required matrices for residuals and tangents
         !----------------------------------------------------------------
         Smat(1,1) = T_tau_s(1,1)
         Smat(2,1) = T_tau_s(2,2)
         Smat(3,1) = T_tau_s(1,2)
         !
         Bmat = zero; G0mat = zero; Gmat = zero
         do kk=1,nNode
            ! Bmat
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,2+nDim*(kk-1)) = dshC(kk,1)
            ! G0mat
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(4,2+nDim*(kk-1)) = dshC0(kk,2)
            ! Gmat
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dshC(kk,2)
         enddo
         !
         Del_N = transpose(dshC)
         !
         AmatUP = zero
         AmatUP(1,1) = D_up(1,1)
         AmatUP(2,1) = D_up(2,2)
         AmatUP(3,1) = D_up(1,2)
         !
         AmatUA = zero
         AmatUA(1,1) = D_ua(1,1)
         AmatUA(2,1) = D_ua(2,2)
         AmatUA(3,1) = D_ua(1,2)
         !
         AmatAU = zero
         AmatAU(1,1) = D_au(1,1,1)
         AmatAU(1,2) = D_au(1,2,1)
         AmatAU(1,3) = D_au(1,1,2)
         AmatAU(1,4) = D_au(1,2,2)
         AmatAU(2,1) = D_au(2,1,1)
         AmatAU(2,2) = D_au(2,2,1)
         AmatAU(2,3) = D_au(2,1,2)
         AmatAU(2,4) = D_au(2,2,2)
         !
         Amat = zero
         Amat(1,1) = D_uu(1,1,1,1)
         Amat(1,2) = D_uu(1,1,2,1)
         Amat(1,3) = D_uu(1,1,1,2)
         Amat(1,4) = D_uu(1,1,2,2)
         Amat(2,1) = D_uu(2,1,1,1)
         Amat(2,2) = D_uu(2,1,2,1)
         Amat(2,3) = D_uu(2,1,1,2)
         Amat(2,4) = D_uu(2,1,2,2)
         Amat(3,1) = D_uu(1,2,1,1)
         Amat(3,2) = D_uu(1,2,2,1)
         Amat(3,3) = D_uu(1,2,1,2)
         Amat(3,4) = D_uu(1,2,2,2)
         Amat(4,1) = D_uu(2,2,1,1)
         Amat(4,2) = D_uu(2,2,2,1)
         Amat(4,3) = D_uu(2,2,1,2)
         Amat(4,4) = D_uu(2,2,2,2)
         !
         Qmat = zero
         Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau_s(1,1)
         Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau_s(1,2)
         Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau_s(1,2)
         Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau_s(2,2)
         Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau_s(1,1)
         Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau_s(1,2)
         Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau_s(1,2)
         Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau_s(2,2)
         !
         !----------------------------------------------------------------
         ! Residuals 
         !----------------------------------------------------------------
         ! Displacement residual vector
         Ru = Ru + detMapJC*w(intPt)*
     +        (
     +       -matmul(transpose(Bmat),Smat)
     +       +biot*P_tau*matmul(transpose(Bmat),AmatUP)
     +        )

         ! Pressure residual vector
         Rp = Rp + detMapJC*w(intPt)*
     +        (
     +       -matmul(dshC,velocity_D)
     +       +biot*transpose(Nvec)*div_v
     +       +(one/K_f)*transpose(Nvec)*dPdt
     +       -transpose(Nvec)*source_P
     +        )  

         ! Astrocyte density residual vector
         Ra = Ra + detMapJC*w(intPt)*
     +        (
     +       +transpose(Nvec)*(dadt + a_tau*dJdt/detF)
     +       +matmul(dshC,qaflux)
     +       -transpose(Nvec)*source_a
     +       +transpose(Nvec)*sink_a
     +        ) 
     
         !----------------------------------------------------------------
         ! Tangents 
         !----------------------------------------------------------------
         ! Displacement tangent matrix
         Kuu = Kuu + detMapJC*w(intPt)*
     +        (
     +        +matmul(matmul(transpose(Gmat),Amat),Gmat)
     +        +matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +        )

         ! Displacement - Pressure tangent matrix
         Kup = Kup + detMapJC*w(intPt)*
     +        (
     +        -biot*matmul(matmul(transpose(Bmat),AmatUP),Nvec)
     +        )
         ! Displacement - Astrocytes tangent matrix
         Kua = Kua + detMapJC*w(intPt)*
     +        (
     +        +matmul(matmul(transpose(Bmat),AmatUA),Nvec)
     +        )
         ! Pressure tangent matrix
         Kpp = Kpp + detMapJC*w(intPt)*
     +        (
     +        -K_h*matmul(dshC,Del_N)
     +        -one/(DTIME*K_f)*matmul(transpose(Nvec),Nvec)
     +        )
         ! Pressure  - Astrocytes tangent matrix
         Kpa = Kpa + detMapJC*w(intPt)*
     +        (
     +        +(Sf0*c_total*alpha_f)*matmul(transpose(Nvec),Nvec)
     +        )  
         ! Astrocytes  - Displacement tangent matrix 
         Kau = Kau + detMapJC*w(intPt)*
     +        (
     +        +matmul(matmul(dshC,AmatAU),Gmat)
     +        )
         ! Astrocytes tangent matrix
         Kaa = Kaa + detMapJC*w(intPt)*
     +        (
     +        -(one/DTIME + dJdt/detF)*matmul(transpose(Nvec),Nvec)
     +        -Diff_a*matmul(dshC,transpose(dshC))
     +        -matmul(matmul(dshC,dqada_adv),Nvec)
     +        -delta_a*matmul(transpose(Nvec),Nvec)
     +        )

      enddo

      ! Pressure  - Displacement tangent matrix
      Kpu = (one/DTIME)*transpose(Kup)

      ! Astrocytes  - Pressure tangent matrix
      Kap = zero
      ! End the loop over body integration points
      !----------------------------------------------------------------
      ! Return RHS vector and the Stiffness matrix to Abaqus
      !----------------------------------------------------------------
      call assembleElement(nDim,nNode,NDOFEL,
     +     Ru,Rp,Ra,Kuu,Kup,Kua,Kpu,Kpp,Kpa,Kau,Kap,Kaa,RHS,AMATRX)     
      ! 
      !----------------------------------------------------------------

      return
      end subroutine UPE4
************************************************************************ 
! Main subroutines 
************************************************************************ 
      subroutine timeIntegral(PROPS,NPROPS,nDim,DTIME,coordint,
     +              T_tau_s,F_tau,F_t,P_tau,dPdX,detFs,dJdt,
     +              J_tau,Je_tau,velocity_D,source_P,flux_D,
     +              a_tau,a_t,dadX,source_a,sink_a,qaflux,
     +              dqada_adv,c_total,D_uu,D_up,D_ua,D_au) 

      ! This subroutine computes everything required for the time integration
      ! Similar to UPE4, "_t" denotes a quantity at the start of the increment 
      ! and "_tau" denotes a quantity at the end of the increment. 
      
      IMPLICIT NONE

      integer i,k,m,NPROPS,stat,nDim

      real(8) DTIME
      real*8 Iden(3,3),PROPS(NPROPS),lambda,mu
      real*8 F_t(3,3),F_tau(3,3),dJdt,Finv(3,3)
      real*8 Fe_tau(3,3),Je_tau,J_tau,J_t,detF
      real*8 T_tau(3,3),T_tau_s(3,3),T_tau_f(3,3),Be_tau(3,3)
      real*8 detFs,Fs_tau(3,3),Fsinv(3,3),dFsda(3,3),fac_s

      real*8 P_tau,K_h,source_P,Sf0,alpha_f
      real*8 flux_D,dPdX(nDim,1),velocity_D(nDim,1),biot,K_f

      real*8 a_tau,a_t,dadX(nDim,1),qaflux(nDim,1),theta_s,k_s,inita
      real*8 source_a,Diff_a,Sa0,delta_a,velocity_a
      real*8 sink_a,dqadFFT(3,3,3),qaIden(3,3,3)
      real*8 dqada_adv(nDim,1),advection_a(nDim,1),diffusion_a(nDim,1)

      real*8 coordint(nDim,1),N_ref(3,1)
      real*8 xcoord,ycoord,initR,c_total,x_center,y_center
      real*8 N_ref_norm,n_cur_norm,n_cur(3,1)
      real*8 D_uu(3,3,3,3),D_up(3,3),D_ua(3,3),D_au(3,3,3)

      real*8 zero,one,two,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,third=1.d0/3.d0)

      call ONEM(Iden)

      ! Obtain material properties
      lambda        = PROPS(1)  
      mu            = PROPS(2) 
      biot          = PROPS(3)
      K_h           = PROPS(4) 
      K_f           = PROPS(5)
      Sf0           = PROPS(6)
      Sa0           = PROPS(7)
      delta_a       = PROPS(8)
      Diff_a        = PROPS(9)
      velocity_a    = PROPS(10)
      k_s           = PROPS(11)  
      alpha_f       = PROPS(12) 
      inita         = PROPS(13)
      initR         = PROPS(14) 
      x_center      = PROPS(18)
      y_center      = PROPS(19)

      !----------------------------------------------------------------
      ! Calculate migration vector for astrocyte advection term 
      !----------------------------------------------------------------
      ! Coordinates at this integration point 
      xcoord = coordint(1,1)
      ycoord = coordint(2,1)

      ! Normal vector in reference coordinates along cytokine gradient
      N_ref(1,1) = xcoord - x_center
      N_ref(2,1) = ycoord - y_center 
      N_ref(3,1) = zero

      N_ref_norm = sqrt(N_ref(1,1)**two + N_ref(2,1)**two + N_ref(3,1)**two)

      if (N_ref_norm .gt. 1.d-12) then
         N_ref = N_ref/N_ref_norm
      else
         N_ref = zero
      endif  

      ! Normal vector in spatial coordinates along cytokine gradient
      n_cur = matmul(F_tau,N_ref)
      n_cur_norm = sqrt(n_cur(1,1)**two + n_cur(2,1)**two + n_cur(3,1)**two)

      if (n_cur_norm .gt. 1.d-12) then
         n_cur = n_cur/n_cur_norm
      else
         n_cur = zero
      endif 

      !----------------------------------------------------------------
      ! Assemble astrocyte flux and its derivatives
      !----------------------------------------------------------------
      advection_a = a_tau*velocity_a*n_cur
      diffusion_a = Diff_a*dadX
      qaflux = advection_a + diffusion_a
      dqada_adv = velocity_a*n_cur
      !
      dqadFFT = zero
      do i=1,3
         do k = 1,3
            do m = 1,3          
               dqadFFT(i,k,m) = dqadFFT(i,k,m)
     +         +a_tau*velocity_a*n_cur(m,1)*(Iden(i,k)-n_cur(i,1)*n_cur(k,1))     
            enddo
         enddo
      enddo

      qaIden = zero
      do i = 1,3     
         do k = 1,3   
            do m = 1,3 
               qaIden(i,k,m) = qaIden(i,k,m) 
     +         +Iden(i,k)*qaflux(m,1)-qaflux(i,1)*Iden(k,m)
            enddo
         enddo
      enddo

      !----------------------------------------------------------------
      ! Calculate source & sink terms for astrocytes and pressure
      !----------------------------------------------------------------
      source_P = Sf0*c_total*(one + alpha_f*a_tau)
      source_a = Sa0*c_total
      sink_a = delta_a*a_tau

      !----------------------------------------------------------------
      ! Calculate swelling deformation 
      !----------------------------------------------------------------
      ! Swelling deformation gradient 
      theta_s = (one + k_s*(a_tau-inita))**third
      Fs_tau = theta_s*Iden

      ! Calculate dFsda
      fac_s = (k_s*third)*(one + k_s*(a_tau-inita))**(-two*third)
      dFsda = fac_s*Iden
      
      ! Inverse of Fs, its determinant, and its transpose
      call matInv3D(Fs_tau,Fsinv,detFs,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detFs.lt.zero'
         call xit
      endif   

      ! Inverse of F_tau, its determinant, and its transpose  
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif

      !----------------------------------------------------------------
      ! Calculate Cauchy stress and jacobian 
      !----------------------------------------------------------------
      ! Elastic deformation gradient and elastic Jacobian 
      Fe_tau = matmul(F_tau,Fsinv)
      call MDET(Fe_tau,Je_tau)
      Be_tau = matmul(Fe_tau,transpose(Fe_tau))

      T_tau_s  = (one/Je_tau)*((lambda*dlog(Je_tau)-mu)*Iden + mu*Be_tau)
      T_tau_f  = biot*P_tau*Iden              
      T_tau = T_tau_s - T_tau_f

      ! Total jacobian 
      call MDET(F_tau,J_tau)
      call MDET(F_t,J_t)

      ! Time derivative of Jacobian 
      dJdt = (J_tau - J_t)/DTIME

      !----------------------------------------------------------------
      ! Calculate Darcys Flux 
      !----------------------------------------------------------------
      flux_D = zero       
      velocity_D = -K_h*dPdX
      do i=1,nDim
        flux_D = flux_D + velocity_D(i,1)**two
      enddo
      flux_D = sqrt(flux_D)

      !----------------------------------------------------------------
      ! Calculate mechanical constitutive matrix
      !----------------------------------------------------------------
      call tangents(PROPS,NPROPS,F_tau,Finv,detF,
     +       Be_tau,Fe_tau,Je_tau,detFs,Fsinv,dFsda,
     +       qaIden,dqadFFT,D_uu,D_up,D_ua,D_au)
 
      !----------------------------------------------------------------
      return
      end subroutine timeIntegral

************************************************************************
      subroutine tangents(PROPS,NPROPS,F_tau,Finv,detF,
     +       Be_tau,Fe_tau,Je_tau,detFs,Fsinv,dFsda,
     +       qaIden,dqadFFT,D_uu,D_up,D_ua,D_au)
      !
      ! This subroutine calculates all constitutive tangents
      !
      IMPLICIT NONE

      integer i,j,k,l,m,n,NPROPS,a,b
      !
      real*8 PROPS(NPROPS),Iden(3,3)
      real*8 AR_mat1(3,3,3,3),AR_mat2(3,3,3,3),AR_mat3(3,3,3,3)
      real*8 AR_mat(3,3,3,3),A_mat(3,3,3,3)
      real*8 D_uu(3,3,3,3),D_up(3,3),D_ua(3,3),D_au(3,3,3)
      real*8 dqadFFT(3,3,3),qaIden(3,3,3)

      real*8 dTdFs(3,3,3,3),dTdFs1(3,3,3,3),dTdFs2(3,3,3,3)
      real*8 dTda(3,3),Be_tau(3,3)
      real*8 Fe_tau(3,3),Je_tau,Fsinv(3,3),dFsda(3,3),detFs

      real*8 mu,lambda
      real*8 F_tau(3,3),Finv(3,3),detF
      real*8 zero,one
      parameter(zero=0.d0,one=1.d0)
      !
      
      call ONEM(Iden)

      lambda        = PROPS(1)  
      mu            = PROPS(2)

      ! Calculate the AR_mat  
      AR_mat1 = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  AR_mat1(i,j,k,l) = AR_mat1(i,j,k,l)
     +            +detFs*(lambda*dlog(Je_tau) - mu + lambda)*Finv(j,i)*Finv(l,k)
     +            +detFs*(mu - lambda*dlog(Je_tau))*Finv(l,i)*Finv(j,k)
               enddo
            enddo
         enddo
      enddo 


      AR_mat2 = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  do a = 1,3
                     do b = 1,3          
                        AR_mat2(i,j,k,l) = AR_mat2(i,j,k,l)
     +                   +detFs*mu*F_tau(i,a)*Finv(l,k)*Fsinv(a,b)*Fsinv(j,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo 

      
      AR_mat3 = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3 
                  do b = 1,3         
                     AR_mat3(i,j,k,l) = AR_mat3(i,j,k,l)
     +                +detFs*mu*Iden(i,k)*Fsinv(l,b)*Fsinv(j,b)
                  enddo 
               enddo
            enddo
         enddo
      enddo 

      AR_mat = AR_mat1 + AR_mat2 + AR_mat3

      ! Map AR_mat to spatial one 
      A_mat = zero 

      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  do m = 1,3
                     do n = 1,3          
                       A_mat(i,j,k,l) = A_mat(i,j,k,l)
     +                 +(one/detF)*F_tau(j,m)*F_tau(l,n)*AR_mat(i,m,k,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo            
          
      ! Calculate dTdFs
      dTdFs1 = zero     
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  dTdFs1(i,j,k,l) = dTdFs1(i,j,k,l)
     +            +(one/Je_tau)
     +            *(
     +             (lambda*dlog(Je_tau)-mu-lambda)*Iden(i,j)*Fsinv(l,k)
     +            +mu*Be_tau(i,j)*Fsinv(k,l)
     +             ) 
               enddo             
            enddo
         enddo
      enddo   
      
      dTdFs2 = zero     
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3
                  do a = 1,3
                  dTdFs2(i,j,k,l) = dTdFs2(i,j,k,l)
     +            +(one/Je_tau)
     +            *(
     +            -mu*Fe_tau(i,k)*Fe_tau(j,a)*Fsinv(l,a)
     +            -mu*Fe_tau(i,a)*Fe_tau(j,k)*Fsinv(l,a)
     +            )
                  enddo
               enddo             
            enddo
         enddo
      enddo 

      dTdFs = dTdFs1 + dTdFs2
      
      ! Calculate dTda 
      dTda = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3         
                  dTda(i,j) = dTda(i,j)
     +            +dTdFs(i,j,l,k)*dFsda(k,l)
               enddo             
            enddo
         enddo
      enddo      
     
      ! spatial modulus tangent
      D_uu = A_mat
      D_up = Iden
      D_ua = dTda
      D_au = qaIden-dqadFFT

      return
      end subroutine tangents

************************************************************************
      subroutine cytokineFunction(PROPS,NPROPS,TIME,nDim,coordint,c_total)

      ! This subroutine calculates the cytokine concentration 

      IMPLICIT NONE
      integer NPROPS,nDim

      real(8) TIME(2)
      real*8 c_total,c1,c2
      real*8 epsilon,beta,tp
      real*8 PROPS(NPROPS),initR
      real*8 xcoord,ycoord,x_center,y_center
      real*8 coordint(nDim,1),coordDiff,radius
      real*8 one,two
      parameter(one=1.d0,two=2.d0)

      initR     = PROPS(14) 
      epsilon   = PROPS(15) 
      beta      = PROPS(16)
      tp        = PROPS(17)   
      x_center  = PROPS(18)
      y_center  = PROPS(19)

      ! Extract coordinates at this integration point 
      xcoord = coordint(1,1)
      ycoord = coordint(2,1)

      ! Calculate the radius based on the polar coordinates 
      radius = sqrt((xcoord - x_center)**two + (ycoord - y_center)**two)
      
      ! Calculate the spatial component of cytokines
      coordDiff = (radius-initR)/initR
      c1 = (exp(-beta*coordDiff))/(one + exp(-beta*coordDiff))

      ! Calculate the temporal component of cytokines
      c2 = epsilon**two / ((TIME(2) - tp)**two + epsilon**two)

      ! Calculate the total normalized concentration of pro-inflammatory 
      c_total = c1*c2

      return 
      end subroutine cytokineFunction
************************************************************************   
      subroutine assembleElement(nDim,nNode,NDOFEL,
     +     Ru,Rp,Ra,Kuu,Kup,Kua,Kpu,Kpp,Kpa,Kau,Kap,Kaa,RHS,AMATRX)  

      !
      ! This subroutine assembles the local elements residuals and tangents
      !
      IMPLICIT NONE

      integer i,j,A11,A12,B11,B12,nDim,nNode,NDOFEL,nDofN

      real*8 Ru(nDim*nNode,1),Rp(nNode,1),Ra(nNode,1)
      real*8 Kuu(nDim*nNode,nDim*nNode),Kup(nDim*nNode,nNode),Kua(nDim*nNode,nNode)
      real*8 Kpu(nNode,nDim*nNode),Kpp(nNode,nNode),Kpa(nNode,nNode)
      real*8 Kau(nNode,nDim*nNode),Kap(nNode,nNode),Kaa(nNode,nNode)
      real*8 RHS(NDOFEL,1),AMATRX(NDOFEL,NDOFEL)
      real*8 zero
      parameter(zero=0.d0)

      ! Total number of degrees of freedom per node
      nDofN = NDOFEL/nNode


      ! Initialize RHS and AMATRX 
      RHS(:,1) = zero
      AMATRX = zero

      
      ! Element level residual
      do i=1,nNode
         A11 = nDofN*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! Displacement
         RHS(A11,1) = Ru(A12,1)
         RHS(A11+1,1) = Ru(A12+1,1)
         !
         ! Pressure
         RHS(A11+2,1) = Rp(i,1)
         !
         ! Astrocytes
         RHS(A11+3,1) = Ra(i,1)
      enddo
      !
      ! Element level tangent matrix
      do i=1,nNode
         do j=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = nDofN*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! Displacement
            AMATRX(A11,B11) = Kuu(A12,B12)
            AMATRX(A11,B11+1) = Kuu(A12,B12+1)
            AMATRX(A11+1,B11) = Kuu(A12+1,B12)
            AMATRX(A11+1,B11+1) = Kuu(A12+1,B12+1)
            !
            ! Pressure
            AMATRX(A11+2,B11+2) = Kpp(i,j)
            !
            ! Astrocytes
            AMATRX(A11+3,B11+3) = Kaa(i,j)
            !
            ! Displacement - Pressure
            AMATRX(A11,B11+2) = Kup(A12,j)
            AMATRX(A11+1,B11+2) = Kup(A12+1,j)
            !
            ! Displacement - Astrocytes
            AMATRX(A11,B11+3) = Kua(A12,j)
            AMATRX(A11+1,B11+3) = Kua(A12+1,j)
            !
            ! Pressure - Displacement
            AMATRX(A11+2,B11) = Kpu(i,B12)
            AMATRX(A11+2,B11+1) = Kpu(i,B12+1)
            !
            ! Pressure - Astrocytes 
            AMATRX(A11+2,B11+3) = Kpa(i,j)
            !
            ! Astrocytes  - Displacement
            AMATRX(A11+3,B11) = Kau(i,B12)
            AMATRX(A11+3,B11+1) = Kau(i,B12+1)
            !
            ! Astrocytes  - Pressure
            AMATRX(A11+3,B11+2) = Kap(i,j)
            !
         enddo
      enddo

      return
      end subroutine assembleElement
!************************************************************************
!     Helper subroutines
!************************************************************************
      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine gets the integration point locations
      !
      IMPLICIT NONE
      !
      integer nIntPt
      !
      real*8 xi(4,2), w(4)
      real*8 zero,one,third
      parameter(zero=0.d0,one=1.d0,third=1.d0/3.d0)

      ! Initialize
      w = zero
      xi = zero

      ! Number of Gauss points
      nIntPt = 4

      ! Gauss weights
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one
      

      ! Gauss pt locations in master element
      xi(1,1) = -dsqrt(third)
      xi(1,2) = -dsqrt(third)
      xi(2,1) = dsqrt(third)
      xi(2,2) = -dsqrt(third)
      xi(3,1) = -dsqrt(third)
      xi(3,2) = dsqrt(third)
      xi(4,1) = dsqrt(third)
      xi(4,2) = dsqrt(third)

      return
      end subroutine xint2D4pt

!************************************************************************
      subroutine calcShape2DLinear(nIntPt,xi_int,intPt,sh,dshxi)
      !
      ! This subroutine calculates the shape functions and their derivatives 
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      ! sh(i) = shape function of node i at the intPt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      IMPLICIT NONE
      integer intPt,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 one,fourth
      parameter(one=1.d0,fourth=1.d0/4.d0)
      
      ! Location in the master element
      xi = xi_int(intPt,1)
      eta = xi_int(intPt,2)
      
      
      ! The shape functions
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)

      return
      end subroutine calcShape2DLinear
!*************************************************************************
      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! This subroutine map derivatives of shape fns from xi-eta-zeta to x-y-z 
      !
      IMPLICIT NONE
      integer i,j,k,nNode,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detMapJ
      !
      real*8 zero
      parameter(zero=0.d0)


      ! Calculate the mapping Jacobian matrix:
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)

      ! Calculate first derivatives wrt x, y, z
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))

      return
      end subroutine mapShape2Da
!************************************************************************
      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! This subroutine returns the inverse, and the determinant for 3D 
      !
      IMPLICIT NONE
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv
      real*8 zero,one
      parameter(zero=0.d0,one=1.d0)

      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. zero) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = one/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D
!************************************************************************
      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! This subroutine returns the inverse, and the determinant for 2D 
      !
      IMPLICIT NONE
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv
      real*8 zero,one
      parameter(zero=0.d0,one=1.d0)
      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. zero) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = one/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!************************************************************************
      subroutine MDET(A,det)
      !
      ! This subroutine calculates the determinant of a 3 by 3 matrix [A]
      !
      IMPLICIT NONE
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +        + A(1,2)*A(2,3)*A(3,1)
     +        + A(1,3)*A(2,1)*A(3,2)
     +        - A(3,1)*A(2,2)*A(1,3)
     +        - A(3,2)*A(2,3)*A(1,1)
     +        - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine MDET
!************************************************************************
      subroutine ONEM(A)
      !
      ! This subroutine creates a 3 by 3 identity matrix [A]
      !
      IMPLICIT NONE
      !
      integer i,j
      real*8 A(3,3)
      real*8 zero,one
      parameter(zero=0.d0,one=1.d0)

      do i=1,3
         do j=1,3
            if (i .eq. j) then
               A(i,j) = one
            else
               A(i,j) = zero
            end if
         end do
      end do


      return
      end subroutine ONEM
!****************************************************************************
