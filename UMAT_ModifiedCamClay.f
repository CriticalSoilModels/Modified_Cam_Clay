*USER SUBROUTINES
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      implicit double precision (a-h, o-z)
!     implicit real(8) (a-h,o-z)


      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

	! Userdefined parameters
	
      	Integer :: iounit = 0, i, ios
      	Save iounit
      	Character(255) :: PrjDir, dbg_file
      	
      	real(8) :: xphics, xNu, xkappa, xlambda, xe0, zeta
      	real(8) :: pp0,p0,q0,pp,p,q, p_trial, q_trial, f,p1,q1
      	real(8) :: dEpsV, dEpsD, xJ2, dEps(6), theta
      	real(8) :: xK, xG, r
      	real(8) :: df_dp,df_dq,xA,d_xlambda,ddlambda
      	real(8) :: T, dT, dT_new dsubepsp, dsubepsq
      	real(8) :: dp, dq, d5p,d5q
      	
      	real(8) :: F1,F2
      	real(8) :: dSig(6), Sig0(6), Sig(6), D(6,6)
      	
      	real(8) :: yield
      	
      	real(8) :: nu_u
      	
      	Integer :: j
          Integer, parameter :: MAXITER=100

	!From usr_lib and usr_add
	Integer :: iOpt
	real(8) :: S1, S2, S3, xN1(3), xN2(3), xN3(3)
	
      
	Integer :: sub
    ! Get parameters from Props and STATEV
      xphics = Props(1)
      xNu = Props(2)
      xkappa = Props(3)
      xlambda = Props(4)
      xe0 = Props(5)
      pp = STATEV(1)
    ! Accumulated plastic strains
      do i = 1,NTENS
        EpsP(i) = STATEV(1+i) 
      end do

      zeta = (xe0 + 1) / (xlambda - xkappa) !helping parameter

    ! Calculate constitutive stresses
      Sig=stress
      dEps=dstran

    ! Set tolerance for yield surface
    ! Reccommended tolerance error (10-6 to 10-9)
      FTOL = 1e-6   
    
    ! Do the predictor corrector scheme. The subroutine calculate the plastic and elastic parts and returns the updated stress and state variables
      MaxIter = 100000
      call implicit_predictor_corrector_integration(xkappa,XNu,&
      xe0,dEps,xphics,FTOL,MaxIter,zeta,Sig,EpsP,dEpsP,pp)    
      
    ! update state variables    
      STATEV(1) = pp
      do i = 1,NTENS
        STATEV(1+i) = EpsP(i)
      end do

    ! if (isundr  ==  1) then        Calculation of pore pressure not needed because done outside the subroutine
    !   Swp = Swp0 - BulkW*(dEpsV)
    ! end if
    
    ! update stress   
      stress=Sig

    ! Calculate effective/elastic D-matrix
      Call FormDEMCC(stress, xkappa, xNu, xe0, DDSDDE, 6, xG, xK)   ! also updates K and G 

    End SUBROUTINE UMAT
      
      
      Subroutine FormDEMCC(Sig0, xkappa,xNu, xe0, D, Id, xG, xK)
C***********************************************************************
C
C    Function:  To form the elastic material stiffness matrix for MCC model (Hooke)

C
C I  xkappa     : slope of the U/R line in e-ln(p') plane
C I  xNu        : Poisson's ratio
C I  xe0        : initial void ratio
C O   D(i,j)    : Resulting matrix
C I   Id        : (First) dimension of D
C O   xG         : Shear modulus
C O   xK         : Bulk modulus
C                                                  
C                              D1  D2  D2 o  o  o  
C     Structure of             D2  D1  D2 o  o  o  
C     elastic D matrix         D2  D2  D1 o  o  o  
C                              o   o   o  G  o  o  
C                              o   o   o  o  G  o  
C                              o   o   o  o  o  G  
C                                                  
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension D(Id,Id)
      Dimension Sig0(6)

      D = 0.0
       P = MAX(-(Sig0(1)+Sig0(2)+Sig0(3))/3., 1d0)
          xK = (xe0 + 1)/xkappa * P  !bulk modulus at start of time step, assumed constant
          r = 3. * ( 1. - 2.*xNu) / ( 2. * (1.+xNu))
          xG = r*xK
      FAC= 2*xG / (1D0 - 2*xNU)
      D1 = FAC * (1D0 - xNU)
      D2 = FAC *   xNU

      Do I=1,3
        Do J=1,3
          D(I,J)=D2
        End Do
        D(I,I)=D1
      End Do
      Do I=4,6
        D(I,I)=xG
      End Do

      End
      
      !-----------------------------------------------------------------------
! Function computing the yield function
!-----------------------------------------------------------------------
	function yield(p, j, pp, theta, xphics)
      !-----------------------------------------------------------------------
      !
      !  function: computing the yield function
      !
      !  input:	p, q, pp, g_theta
      !  output:	yield
      !
      !-----------------------------------------------------------------------

	implicit none

        real(8) :: p, j, pp, g_theta
	  real(8) :: yield
        
      g_theta = cos(theta) + ( (sin(theta) * sin(xphics)) / sqrt(3) )
      g_theta = sin(xphics) / g_theta  
      
        yield = ( j / (p*g_theta) )**2 - ((pp/p) - 1)

      end function yield


!-----------------------------------------------------------------------
! Subroutine computing the derivatives
!-----------------------------------------------------------------------
	subroutine derivatives(sig,p,j,xphics,theta,dgdp,dfdsig,dgdsig)
!-----------------------------------------------------------------------
!	input:	sig,p,j,xphics,theta
!	output:	dfdsig,dgdsig
!-----------------------------------------------------------------------
	  implicit none
	  real(8), intent(in)  :: p, j, xphics, theta
	  real(8), intent(in), dimension(6)  :: sig 
	  real(8), intent(out), dimension(6) :: dgdp,dfdsig,dgdsig
        
        ! Local Variables 
        real(8) :: g_theta,dfdp,dfdj,dfdtheta,
        real(8) :: dgdj,dgdtheta,dets
        real(8), dimension(6)  :: dpdsig,djdsig,d_dets_dsig,dthetadsig
        
       g_theta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt(3.0d0))
       g_theta = sin(xphics) / g_theta 
        
        dfdp = (1/p) * ( 1 - ( ( j / (p*g_theta) )**2) )
        dfdj = (2*j) / ((p*g_theta)**2)
        dfdtheta = (2*(j**2)) / (sqrt(3)*(p**2)*g_theta*(sin(xphics)))
        dfdtheta = dfdtheta * ( (cos(theta)*sin(xphics)) - sin(theta))

        
        dgdp = (1/p) * ( 1 - ( j / (p*g_theta) ) )
        dgdj = (2*j) / ((p*g_theta)**2)
        dgdtheta = 0 
        
        dpdsig = 0.33333333d0*[1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0]
        djdsig = (1/(2*j))*[sig(1)-p, sig(2)-p, sig(3)-p, 2*sig(4), &
                     2*sig(5), 2*sig(6)]
        d_dets_dsig = [ (((sig(2)-p)*(sig(3)-p)) - (sig(5)**2) ), &
               (((sig(1)-p)*(sig(3)-p)) - (sig(6)**2) ), &
               (((sig(1)-p)*(sig(2)-p)) - (sig(4)**2) ), &
               (( 2* sig(4) * (p-sig(3)) ) + (2*sig(5)*sig(6)) ), &
               (( 2* sig(5) * (p-sig(1)) ) + (2*sig(4)*sig(6)) ), &
               (( 2* sig(6) * (p-sig(2)) ) + (2*sig(4)*sig(5)) ) ]
        dets = ((sig(1)-p)*(sig(2)-p)*(sig(3)-p))  & 
             -  ((sig(1)-p)*(sig(5)**2)) 
            - ((sig(2)-p)*(sig(6)**2)) - ((sig(3)-p)*(sig(4)**2)) & 
                + (2*sig(4)*sig(5)*sig(6))
        dthetadsig = (((dets/j)*djdsig) - d_dets_dsig) * (sqrt(3)/2)
        dthetadsig = dthetadsig / (cos((3.0d0)*theta)*(j**3))
        
        dfdsig = (dfdp * dpdsig) + (dfdj * djdsig) + &
            (dfdtheta * dthetadsig)
        dgdsig = (dgdp * dpdsig) + (dgdj * djdsig) + &
            (dgdtheta * dthetadsig)
        
        
	end subroutine derivatives

	
!-----------------------------------------------------------------------
! Subroutine computing the material stiffness parameters
!-----------------------------------------------------------------------
	subroutine stiffnessMCC(p,xe0,xkappa,xNu,xG,xK)
!-----------------------------------------------------------------------
!	input:	p,xe0,xkappa,xNu
!	output:	xG,xK
!	local:	r (xG-xK ratio)
!-----------------------------------------------------------------------
	  implicit none
	  real(8), intent(in) :: p,xe0,xkappa,xNu
	  real(8), intent(out) :: xG,xK
	  real(8) :: r, abs_p
	  
	  abs_p = abs(p)
	  xK = (xe0 + 1)/xkappa * abs_p  !bulk modulus at start of time step, assumed constant
	  r = 3. * ( 1. - 2.*xNu) / ( 2. * (1.+xNu))
	  xG = r*xK !shear modulus
	end subroutine stiffnessMCC
	  
	
	
!-----------------------------------------------------------------------
! Function computing dlambda
!-----------------------------------------------------------------------
	function dlambda(dfdp,dfdq,dfdtheta,dgdp,dgdq,dgdtheta,xA,xK,xG,&
                      dep,deq)
!-----------------------------------------------------------------------
!	input:	dfdp,dfdq,dfdtheta,dgdp,dgdq,dgdtheta
!		xK,xG
!		dep,deq
!	output:	dlambda
!-----------------------------------------------------------------------
	  implicit none
	  real(8) :: dfdp,dfdq,dfdtheta,dgdp,dgdq,dgdtheta,xA,xK,xG
	  real(8) :: dep,deq
	  real(8) :: dlambda
	  
        dpdsig(1) = 1/3 
        dpdsig(2) = 1/3 
        dpdsig(3) = 1/3 
        dpdsig(4) = 0 
        dpdsig(5) = 0 
        dpdsig(6) = 0 
        
        dqdsig(1) = 1/  * (S3 + S1) 
        dqdsig(2) = 1/q * (S1 + S3)
        dqdsig(3) = 1/q * (S1 + S3)
        dqdsig(4) = 0 
        dqdsig(5) = 0 
        dqdsig(6) = 0
        
        
	  dlambda = (dfdp*xK*dep + dfdq*3.*xG*deq) /
     &			(dfdp**2*xK + dfdq**2*3.*xG + xA)
     	  return
     	end function dlambda
	


	subroutine error(dp, dq, p, q, d5p, d5q, dT,dT_new,T,sub)
      !-----------------------------------------------------------------------
      !	
      !  function: estimating error and calculating dT_new if necessary
      !
      !  input:   dp,dq,d5p,d5q
      !	          p,q
      !		      dT, T
      !  output:  dT_new
      !		      sub
      !
      !-----------------------------------------------------------------------

      implicit none

        ! arguments
        real(8), intent(in) :: dp, dq, d5p, d5q, p, q
	  real(8), intent(in) :: dT,T
        real(8), intent(out):: dT_new
        integer, intent(out) :: sub			! control parameter

        ! local variables
	  real(8) :: R, beta
	  real(8), parameter :: SSTOL = 1e-6	! not userdefined yet

	  ! estimating error
	  R = sqrt( (d5p - dp)**2 + (d5q - dq)**2 ) 
        R = R / sqrt( (p + dp)**2 + (q +dq)**2 )

        ! evaluating need for substepping
	  if (R  >  SSTOL) then  
	    sub = 1 !the increment will be rejected
      else
	    sub = 0 !the increment is accurate enough and will be accepted
      end if  
        
        !compute new size of substep
        beta = 0.8 * sqrt(SSTOL / R)
	  if (beta  <  0.1) beta = 0.1
	  if (beta  >  2.) beta = 2.

        dT_new = beta * dT
        ! checking the accumulated substeps
                ! checking the accumulated substeps
	  if ((sub == 0).and.(T + dT + dT_new > 1.)) then 
		dT_new = 1. - (T+dT)
      end if


      end subroutine error
!-----------------------------------------------
!	Subroutine with RKF45 stress estimate
!-----------------------------------------------
      subroutine rk(p,q,pp,&
     2		 zeta, xM,phics,theta, xK, xG,
     3	            dep, deq,
     4	      	    dp,dq,dpstar,dqstar,d_xlambda)
!---------------------------------------------------------------
!	input:	zeta, xM, xK, xG	material properties
!		p,q,pp			stress state
!		dep,deq			substep strain incr
!	output: dp,dq,d5p,d5q		stress increments RK4 and 5
!		d_xlambda
!
!	local:	xT,xN
!		dep1-6,deq1-6,dp1-6,dq1-6
!-------------------------------------------------------------------
      implicit none
      real(8),intent(in) :: p,q,pp,zeta,xM,phics,theta
      real(8), intent(in) :: dep,deq
      real(8), intent(in) :: xK,xG
      
      real(8) :: dfdp,dfdq,xA
      real(8) :: dlambda
      
      real(8) :: d_xlambda,dp,dq,dpstar,dqstar
      real(8) :: dptemp,dqtemp
      
      real(8),dimension(2,7) :: dsig
      real(8),dimension(7) :: gamma
      real(8),dimension(6,6) :: koeff
      
      integer :: i,j
      
      
	DO j = 1,6
	  dsig(1,j) = 0.
	  dsig(2,j) = 0.
	END DO
      
      !Initialising the coefficient vectors of DOPRI54
        gamma = (/5179./57600., 0., 7571./16695., 393./640.,
     2            -92097./339200., 187./2100.,1./40./) !coefficients for order 5 stress update
     	!coefficients for order 4 not needed
      
        koeff(1,:) = (/ 1./5.,       	0.,          0.,         
     2        		0.,		0.,	     0./)
        koeff(2,:) = (/ 3./40.,      9./40.,         0.,
     2			0.,		0.,	     0. /)
        koeff(3,:) = (/ 44./45., -56./15.,32./9.,
     2			0.,		0.,	     0. /)
        koeff(4,:) = (/ 19372./6561.,   -25360./2187., 64448./6561.,
     2              -212./729., 0.,0./)
        koeff(5,:) = (/ 9017./3168.,     -355./33.,   46732./5247.,
     2              49./176., -5103./18656.,0. /)
        koeff(6,:) = (/ 35./384.,     0.,   500./1113.,
     2              125./192., -2187./6784.,11./84. /)
      
      ! Computing the partial stress increments
        DO i = 1,7
          IF (i  >  1) THEN
            dptemp = sum(koeff(i-1,:)*dsig(1,1:i-1))
            dqtemp = sum(koeff(i-1,:)*dsig(2,1:i-1))
          ELSE
            dptemp = 0.
            dqtemp = 0.
          END IF
          
     1     CALL derivatives(p+dptemp,q+dqtemp,pp,zeta,xM,phics,theta,
     2                    dfdp,dfdq,dfdtheta,dgdp,dgdq,dgdtheta,xA)     
          
          if (i  ==  1) then
     1       d_xlambda = dlambda(dfdp,dfdq,dfdtheta,dgdp,dgdq,
     2         dgdtheta,xA,xK,xG,dep,deq)        
          end if
          
          dsig(1,i) = xK*dep-d_xlambda*xK*dfdp
          dsig(2,i) = 3.*xG*deq-d_xlambda*3.*xG*dfdq
        END DO
        
        dp = sum(dsig(1,:)*gamma)
        dq = sum(dsig(2,:)*gamma)
        
        !for error estimating
        dpstar = dptemp 
        dqstar = dqtemp


      end subroutine rk

      Subroutine GetModelCount(nMod)
      !
      ! Return the maximum model number (iMod) in this DLL
      !
      Integer (Kind=4) nMod

      nMod = 1 ! Maximum model number (iMod) in current DLL

      Return
      End ! GetModelCount

      Subroutine GetModelName( iMod , ModelName )
      !
      ! Return the name of the different models
      !
      Integer  iMod
      Character (Len= * ) ModelName
      Character (Len=255) tName

      tName = 'MCC_Explicit'

      LT = Len_Trim(tName)
      ModelName= Char(lt) // tName(1:Lt)
      Return
      End ! GetModelName

      Subroutine GetParamCount( iMod , nParam )
      !
      ! Return the number of parameters of the different models
      !
      nParam = 5
      
      
      Return
      End ! GetParamCount

      Subroutine GetParamName( iMod , iParam, ParamName )
      !
      ! Return the parameters name of the different models
      !
      Character (Len=255) ParamName, Units
      Call GetParamAndUnit(iMod,iParam,ParamName,Units)
      Return
      End

      Subroutine GetParamUnit( iMod , iParam, Units )
      !
      ! Return the units of the different parameters of the different models
      !
      Character (Len=255) ParamName, Units
      Call GetParamAndUnit(iMod,iParam,ParamName,Units)
      Return
      End

      Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
      !
      ! Return the parameters name and units of the different models
      !
      ! Units: use F for force unit
      !            L for length unit
      !            T for time unit
      !
      Character (Len=255) ParamName, Units, tName
      Select Case (iMod)
        Case (1)
          ! ModName = 'DP'
          Select Case (iParam)
            Case (1)
              ParamName = 'M'     ; Units     = '-'
            Case (2)
              ParamName = 'nu'      ; Units     = '-'
            Case (3)
              ParamName = 'kappa'   ; Units     = '-'
            Case (4)
              ParamName = 'lambda'  ; Units     = '-'
            Case (5)
              ParamName = 'e_0'       ; Units    = '-'
            Case Default
              ParamName = 'xxx'     ; Units     = '-'
          End Select
        Case Default
          ! model not in DLL
          ParamName = ' N/A '        ; Units     = ' N/A '
      End Select
      tName    = ParamName
      LT       = Len_Trim(tName)
      ParamName= Char(lt) // tName(1:Lt)

      tName = Units
      LT    = Len_Trim(tName)
      Units = Char(lt) // tName(1:Lt)

      Return
      End ! GetParamAndUnit



      Subroutine MZEROI(I,K)
C
C***********************************************************************
C
C     Function: To make an integre array I with Dimension K to zero
C
C***********************************************************************
C
      Dimension I(*)

      Do J=1,K
        I(J)=0
      End Do

      Return
      End

      Subroutine SETRVAL(R,K,V)
C
C***********************************************************************
C
C     Function: To fill a real array R with Dimension K with value V
C
C***********************************************************************
C
      Implicit real(8) (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J)=V
      End Do

      Return
      End

      Subroutine SETIVAL(I,K,IV)
C
C***********************************************************************
C
C     Function: To fill an integer array I with Dimension K with value IV
C
C***********************************************************************
C
      Implicit real(8) (A-H,O-Z)
      Dimension I(*)

      Do J=1,K
        I(J)=IV
      End Do

      Return
      End

      Subroutine COPYIVEC(I1,I2,K)
C
C***********************************************************************
C
C     Function: To copy an integer array I1 with Dimension K to I2
C
C***********************************************************************
C
      Implicit real(8) (A-H,O-Z)
      Dimension I1(*),I2(*)

      Do  J=1,K
        I2(J)=I1(J)
      End Do

      Return
      End
     
C***********************************************************************
      Subroutine MulVec(V,F,K)
C***********************************************************************
C
C     Function: To multiply a real vector V with dimension K by F
C
C***********************************************************************
C
      IMPLICIT real(8) (A-H,O-Z)
      DIMENSION V(*)

      Do J=1,K
        V(J)=F*V(J)
      End Do

      Return
      End     ! Subroutine Mulvec
C***********************************************************************

      
      Subroutine MatMatSq(n, xMat1, xMat2, xMatR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2 for square matrices, size n
C
C I   n     : Dimension of matrices
C I   xMat1 : Matrix (n,*)
C I   xMat2 : Matrix (n,*)
C O   xMatR : Resulting matrix (n,*)
C
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension xMat1(n,*),xMat2(n,*),xMatR(n,*)
C**********************************************************************

      Do I=1,n
        Do J=1,n
          X=0
          Do K=1,n
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMatSq

C***********************************************************************
      Subroutine WriVal ( io, C , V )
C***********************************************************************
C
C Write (Double) value to file unit io (when io>0)
C
C***********************************************************************
C
      Implicit real(8) (A-H,O-Z)
      Character C*(*)

      If (io <= 0) Return

      Write(io,*) C,V
    1 Format( A,3x, 1x,1p,e12.5)
      Return
      End
C***********************************************************************
      Subroutine WriIVl ( io, C , I )
C***********************************************************************
C
C Write (integer) value to file unit io (when io>0)
C
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Character C*(*)

      If (io <= 0) Return

      Write(io,*) C,I
    1 Format( A,3x, 1x,I6)
      Return
      End
C***********************************************************************
      Subroutine WriIVc ( io, C , iV , n )
C***********************************************************************
C
C Write (integer) vector to file unit io (when io>0)
C
C***********************************************************************
      Character C*(*)
      Dimension iV(*)

      If (io <= 0) Return

      Write(io,*) C
      Write(io,1) (iv(i),i=1,n)
    1 Format( ( 2(3x,5i4) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec ( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io <= 0) Return

      If (Len_Trim(C) <= 6) Then
        Write(io,2) C,( V(i),i=1,n)
      Else
        Write(io,*) C
        Write(io,1) ( V(i),i=1,n)
      End If
    1 Format( ( 2(1x, 3(1x,1p,e10.3) ) ) )
    2 Format( A, ( T7, 2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec5( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 5 values per line
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io <= 0) Return

      Write(io,*) C
      Write(io,1) ( V(i),i=1,n)
    1 Format( 5(1x,1p,e12.5) )
      Return
      End
C***********************************************************************
      Subroutine WriMat ( io, C , V , nd, nr, nc )
C***********************************************************************
C
C Write (Double) matrix to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Character C*(*)
      Dimension V(nd,*)

      If (io <= 0) Return

      Write(io,*) C
      Do j=1,nr
        Write(io,1) j,( V(j,i),i=1,nc)
      End Do
    1 Format(i4, (  T7,2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************



      subroutine setveclen(xn,n,xl)
C**********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension xN(*)
      x=0
      do i=1,n
        x=x+xn(i)**2
      end do
      if (x /= 0) Then
        f=xl/sqrt(x)
        do i=1,3
          xn(i)=xn(i)*f
        end do
      end if
      return
      end ! setveclen

C***********************************************************************
      Subroutine MatVec(xMat,IM,Vec,N,VecR)
C***********************************************************************
C
C     Calculate VecR = xMat*Vec
C
C I   xMat  : (Square) Matrix (IM,*)
C I   Vec   : Vector
C I   N     : Number of rows/colums
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine MatVec

C***********************************************************************
      Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
C***********************************************************************
C
C     Calculate VecR() = R1*Vec1()+R2*Vec2()
C
C I   Vec1,
C I   Vec2  : Vectors
C I   R1,R2 : Multipliers
C I   N     : Number of rows
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec
C
C***********************************************************************

      Subroutine CarSig(S1,S2,S3,xN1,xN2,xN3,SNew)
C***********************************************************************
C
C     Returns the Cartesian stresses using the principal stresses S1..S3
C     and the principal directions
C
C I   S1..S3   : Principal stresses
C I   xN1..xN3 : Principal directions (xNi for Si)
C
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*),SNew(*)
      Dimension SM(3,3),T(3,3),TT(3,3),STT(3,3)
C***********************************************************************
C
C**** Fill transformation (rotation) matrix
C
      Do I=1,3
        T(I,1) = xN1(I)
        T(I,2) = xN2(I)
        T(I,3) = xN3(I)
        TT(1,I) = T(I,1)
        TT(2,I) = T(I,2)
        TT(3,I) = T(I,3)
      End Do
!      Call MatTranspose(T,3,TT,3,3,3)

      Call MZeroR(SM,9)
      SM(1,1) = S1
      SM(2,2) = S2
      SM(3,3) = S3
C
C**** SMnew = T*SM*TT
C
      
      Call MatMat(SM ,3,  TT,3 , 3,3,3 ,STT,3)
      Call MatMat( T ,3, STT,3 , 3,3,3 ,SM ,3)
!     Call MatMatSq(3, SM,  TT, STT )   ! STT = SM*TT
!     Call MatMatSq(3,  T, STT, SM  )   ! SM  =  T*STT
C
C**** Extract cartesian stress vector from stress matrix
C
      Do I=1,3
        SNew(I) = SM(I,I)
      End Do
      SNew(4) = SM(2,1)
      SNew(5) = SM(3,2)
      SNew(6) = SM(3,1)

      Return
      End     ! Subroutine CarSig
C**********************************************************************

            subroutine PrincipalSig(IOpt, S, xN1, xN2, xN3, S1, S2, S3,
     1 P, Q,J,theta)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions
      !            from cartesian stress vector
      !
      !  IOpt            I   I     flag to calculate principal direction (IOpt = 1)
      !  IntGlo          I   I     global ID of Gauss point or particle 
      !  S               I   R()   cartesian stress
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      !  J               O   R     sqrt J2
      !  theta           O   R     lode angle
      ! 
      !-------------------------------------------------------------------
      
      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt
        real(8), intent(in) :: S(6)
        real(8), intent(out) :: xN1(3), xN2(3), xN3(3),
     &                                   S1, S2, S3, P, Q

        if (IOpt .eq. 1) then
          call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q,J,theta) ! Calculate principal direction
        else
          call Eig_3a(0,S,S1,S2,S3,P,Q) ! Do not calculate principal direction
        end if
     
      end subroutine PrincipalSig
C**********************************************************************
      subroutine Eig_3(iOpt, St, xN1, xN2, xN3, S1, S2, S3, P, Q, &
          J, theta)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions 
      !            from cartesian stress vector
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      ! 
      !  IOpt            I   I     flag for output writing (IOpt = 1) 
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      !  J               O   R     sqrt J2
      !  theta           O   R     lode angle
      ! 
      !-------------------------------------------------------------------

      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt
        real(8), intent(in) :: St(6)
        real(8), intent(out) :: xN1(3), xN2(3), xN3(3),
     &                                   S1, S2, S3, P, Q, J, theta
        
        ! local variables
        real(8) :: A(3,3), V(3,3)
        real(8) :: abs_max_s, tol
        real(8) :: tau, sign_tau, t, c, s
        real(8) :: temp1, temp2, temp3
        integer :: i, k, it, itmax, ip, iq
        integer :: iS1, iS2, iS3

        ! Put cartesian stress vector into matrix A
        A(1,1) = St(1) ! xx
        A(1,2) = St(4) ! xy = yx
        A(1,3) = St(6) ! zx = xz

        A(2,1) = St(4) ! xy = yx
        A(2,2) = St(2) ! yy
        A(2,3) = St(5) ! zy = yz

        A(3,1) = St(6) ! zx = xz
        A(3,2) = St(5) ! zy = yz
        A(3,3) = St(3) ! zz

        ! Set V to unity matrix
        V(1,1) = 1
        V(2,1) = 0
        V(3,1) = 0

        V(1,2) = 0
        V(2,2) = 1
        V(3,2) = 0

        V(1,3) = 0
        V(2,3) = 0
        V(3,3) = 1

        ! get maximum value of cartesian stress vector
        abs_max_s = 0.0
        do i = 1,6
          if (abs(St(i)) .gt. abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-16 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it .lt. itmax) .and. 
     &             (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) .gt. tol) )
        
          it = it + 1
          do k = 1,3
                
            if (k .eq. 1) then
              ip = 1
              iq = 2
            else if (k .eq.2) then
              ip = 2
              iq = 3
            else
              ip = 1
              iq = 3
            end if
            
            if (abs(A(ip,iq)) .gt. 1d-50) then
                
              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
              if (tau .ge. 0.0) then
                sign_tau = 1.0
              else
                sign_tau = -1.0
              end if
              
              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
              c = 1.0 / sqrt(1.0 + t**2)
              s = t * c
              
              temp1 = c * A(1, ip) - s * A(1, iq)
              temp2 = c * A(2, ip) - s * A(2, iq)
              temp3 = c * A(3, ip) - s * A(3, iq)
              A(1, iq) = s * A(1, ip) + c * A(1, iq)
              A(2, iq) = s * A(2, ip) + c * A(2, iq)
              A(3, iq) = s * A(3, ip) + c * A(3, iq)
              A(1, ip) = temp1
              A(2, ip) = temp2
              A(3, ip) = temp3

              temp1 = c * V(1, ip) - s * V(1, iq)
              temp2 = c * V(2, ip) - s * V(2, iq)
              temp3 = c * V(3, ip) - s * V(3, iq)
              V(1, iq) = s * V(1, ip) + c * V(1, iq)
              V(2, iq) = s * V(2, ip) + c * V(2, iq)
              V(3, iq) = s * V(3, ip) + c * V(3, iq)
              V(1, ip) = temp1
              V(2, ip) = temp2
              V(3, ip) = temp3

              temp1 = c * A(ip, 1) - s * A(iq, 1)
              temp2 = c * A(ip, 2) - s * A(iq, 2)
              temp3 = c * A(ip, 3) - s * A(iq, 3)
              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
              A(ip, 1) = temp1
              A(ip, 2) = temp2
              A(ip, 3) = temp3
            end if ! A(ip,iq)<>0
            
          end do ! k
          

          
        end do ! while
     
        ! get principal stresses from diagonal of A
        S1 = A(1, 1)
        S2 = A(2, 2)
        S3 = A(3, 3)
      
        ! derived invariants
        P = (S1 + S2 + S3) / 3.
        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )
        J = Q /sqrt(3).
        theta = atan( (1/sqrt(3)) * ( ( (2*(S2-S3)) / (S1-S3) ) - 1. ) )

        ! Sort eigenvalues S1 <= S2 <= S3
        iS1 = 1
        iS2 = 2
        iS3 = 3
        
        if (S1 .gt. S2) then
          t   = S2
          S2  = S1
          S1  = t
          it  = iS2
          iS2 = iS1
          iS1 = it
        end if
        
        if (S2 .gt. S3) then
          t   = S3
          S3  = S2
          S2  = t
          it  = iS3
          iS3 = iS2
          iS2 = it
        end if
        
        if (S1 .gt. S2) then
          t   = S2
          S2  = S1
          S1  = t
          it  = iS2
          iS2 = iS1
          iS1 = it
        end if
        
        ! get corresponding principal directions from V
        do i = 1,3
          xN1(i) = V(i, is1)
          xN2(i) = V(i, is2)
          xN3(i) = V(i, is3)
        end do
      
        ! optional output writing


      end subroutine Eig_3

C**********************************************************************
      
      subroutine Eig_3a(iOpt, St, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses from cartesian stress vector
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      ! 
      !  IOpt            I   I     flag for output writing (IOpt = 1) 
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------

      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt
        real(8), intent(in) :: St(6)
        real(8), intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        real(8) :: A(3,3)
        real(8) :: abs_max_s, tol
        real(8) :: tau, sign_tau, t, c, s
        real(8) :: temp1, temp2, temp3
        integer :: i, k, it, itmax, ip, iq

        ! Put cartesian stress vector into matrix A
        A(1,1) = St(1) ! xx
        A(1,2) = St(4) ! xy = yx
        A(1,3) = St(6) ! zx = xz

        A(2,1) = St(4) ! xy = yx
        A(2,2) = St(2) ! yy
        A(2,3) = St(5) ! zy = yz

        A(3,1) = St(6) ! zx = xz
        A(3,2) = St(5) ! zy = yz
        A(3,3) = St(3) ! zz

        ! get maximum value of cartesian stress vector
        abs_max_s = 0.0
        do i = 1,6
          if (abs(St(i)) .gt. abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-20 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it .lt. itmax) .and. 
     &             (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) .gt. tol) )
        
          it = it + 1
          do k = 1,3
            if (k .eq. 1) then
              ip = 1
              iq = 2
            else if (k .eq.2) then
              ip = 2
              iq = 3
            else
              ip = 1
              iq = 3
            end if
            
            if (abs(A(ip,iq)) .gt. 1d-50) then
                
              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
              if (tau .ge. 0.0) then
                sign_tau = 1.0
              else
                sign_tau = -1.0
              end if
              
              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
              c = 1.0 / sqrt(1.0 + t**2)
              s = t * c

              temp1 = c * A(1, ip) - s * A(1, iq)
              temp2 = c * A(2, ip) - s * A(2, iq)
              temp3 = c * A(3, ip) - s * A(3, iq)
              A(1, iq) = s * A(1, ip) + c * A(1, iq)
              A(2, iq) = s * A(2, ip) + c * A(2, iq)
              A(3, iq) = s * A(3, ip) + c * A(3, iq)
              A(1, ip) = temp1
              A(2, ip) = temp2
              A(3, ip) = temp3

              temp1 = c * A(ip, 1) - s * A(iq, 1)
              temp2 = c * A(ip, 2) - s * A(iq, 2)
              temp3 = c * A(ip, 3) - s * A(iq, 3)
              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
              A(ip, 1) = temp1
              A(ip, 2) = temp2
              A(ip, 3) = temp3
  
            end if ! A(ip,iq)<>0
            
          end do ! k
          
          ! optional output writing

          
        end do ! while
     
        ! get principal stresses from diagonal of A
        S1 = A(1, 1)
        S2 = A(2, 2)
        S3 = A(3, 3)
      
        ! derived invariants
        P = (S1 + S2 + S3) / 3.
        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )

        ! Sort eigenvalues S1 <= S2 <= S3
        if (S1 .gt. S2) then
          t   = S2
          S2  = S1
          S1  = t
        end if
        
        if (S2 .gt. S3) then
          t   = S3
          S3  = S2
          S2  = t
        end if
        
        if (S1 .gt. S2) then
          t   = S2
          S2  = S1
          S1  = t
        end if

        ! optional output writing

        
      end subroutine Eig_3a      
C************************************
     
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2
C
C I   xMat1 : Matrix (Id1,*)
C I   xMat2 : Matrix (Id2,*)
C I   nR1   : Number of rows in resulting matrix    (= No rows in xMat1)
C I   nC2   : Number of columns in resulting matrix (= No cols in xMat2)
C I   nC1   : Number of columns in matrix xMat1
C             = Number  rows    in matrix xMat2
C O   xMatR : Resulting matrix (IdR,*)
C
C***********************************************************************
      Implicit real(8) (A-H,O-Z)
      Dimension xMat1(Id1,*),xMat2(Id2,*),xMatR(IdR,*)
C**********************************************************************

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMat
C***********************************************************************

      subroutine MZEROR(R, K)
      !-----------------------------------------------
      !
      !  Function: set real array to zero values
      !
      !  R   I/O   R()    Real array
      !  K   I     I      Length of real array
      !
      !-----------------------------------------------
      
      implicit none
      
        ! arguments
        integer, intent(in) :: K
        real(8), dimension(K), intent(inout) :: R

        ! set each entry to zero
        R(1:K) = 0.0

      end subroutine MZEROR
      
C**********************************************************************
      
!-----------------------------------------------
!	Subroutine for implicit stress estimate
!-----------------------------------------------
      subroutine implicit_predictor_corrector_integration(xkappa,XNu,&
        xe0,dEps,xphics,FTOL,MaxIter,zeta,Sig,EpsP,dEpsP,pp)
!---------------------------------------------------------------
!	input
!        xkappa:     
!           xNu:  
!           xe0:
!          dEps:                 
!        xphics:                                                                            
!          FTOL:         Tolerence on the yield surface 
!       MaxIter:         Maximum iterations for the algorithm                      
!          zeta:         (v / lambda + kappa)      
!---------------------------------------------------------------
!	input/output: 
!            Sig:        Updated stress (Voight notation) 
!           EpsP:        Plastic strain
!          dEpsP:        Incremental plastic strain 
!             pp:        State variable - preconsolidation pressure 
!---------------------------------------------------------------
!	output: 
!          dEpsP:        Incremental plastic strain 
!---------------------------------------------------------------
!	local:	
!           Sigu:        Local variables for calculation - Current stress 
!          EpsPu:        Local variables for calculation - Plastic strain 
!         dEpsPu:        Local variables for calculation - Incremental plastic strain                  
!              D:        Elastic D Matrix        
!             xG:        Shear Mod
!             xK:        Bulk Mod
!          dSigu:        Stress increment (total)
!           iOpt:        Integer for eigen value options
!  xN1, xN2, xN3:        Eigen directions (not used directly)
!     S1, S2, S3:        Principal stresses
!    p,q,j,theta:        Invariants 
!              F:        Yield surface value        
!        counter:        Counter for the implicit algorithm (exits if higher) 
!           dgdp:        Derivative of plastic potential with pressure = plastic volumetric strain
!         dfdsig:        Derivative of yield potential with current stress state (n vector)                
!         dgdsig:        Derivative of plastic potential with current stress state (m vector)                
!          sqrt3:        Square root of 3  
!         gtheta:        used in 3D algorithm in yield surface instead of M                     
!        result1:        dummy variable for calculations 
!        result2:        dummy variable for calculations 
!        dlambda:        Lambda dot  
!              A:        Hardening/Softening parameter  

!-------------------------------------------------------------------
      implicit none
      ! Input variables
      real(8),intent(in) :: xkappa,xNu,xe0,xphics,FTOL,MaxIter,zeta
      real(8), intent(in), dimension(6)  :: dEps
      
      ! Input/Output variables
      real(8), intent(inout), dimension(6)  :: Sig,EpsP,dEpsP
      real(8), intent(inout)  :: pp

      ! Output variables
      real(8), intent(inout), dimension(6)  :: dEpsP

      ! Local variables 
      real(8)  :: xK,xG,xN1,xN2,xN3,S1,S2,S3
      real(8)  :: p,q,j,theta,F,dgdp,dfdsig,dgdsig,sqrt3
      real(8)  :: gtheta,result1,result2,dlambda,A
      real(8), dimension(6)  :: Sigu,EpsPu,dEpsPu,dSigu
      real(8),dimension(6,6) :: D
      integer :: iOpt,counter
      
    ! Initialization  
      DEpsP = 0.0d0
      F = 0.0d0

    !Store variables for updating
      Sigu = Sig
      EpsPu = EpsP
      
    ! Update G,K and evaluate D
      Call FormDEMCC(Sigu, xkappa, xNu, xe0, D, 6, xG, xK)
      
      call MatVec(D, 6, dEps, 6, dSigu)                                
      Sigu = Sigu + dSigu
      
      
    ! Calculate F for the updated Sig - find variants first
      iOpt = 1
      call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p, q, &
          j, theta)
      
      F = yield(p, j, pp, theta, xphics)
   
      if (abs(F) < FTOL) then
              ! Prediction of the stress and strain values are correct and the values can be updated and returned
              ! Update Sig, EpsP, dEpsP
              Sig = Sigu
              EpsP(:) = 0
              dEpsP(:) = 0
      
              ! Update state parameters values
              pp = pp
      
              ! Exit the subroutine
              return
      end if
      
    ! Max number of plastic descent iterations
      MaxIter = 100000
      counter = 0
          
      do while(abs(F) >= FTOL .and. counter <= MaxIter)
              
        ! Calc n_vec, m_vec
          call derivatives(Sigu,p,j,xphics,theta,dgdp,dfdsig,dgdsig)    
          
        ! Calculate A (make function later)
          sqrt3 = sqrt(3.0d0)
          g_theta = cos(theta) + ((sin(theta) * sin(xphics)) / sqrt3)
          g_theta = sin(xphics) / g_theta 
          A = zeta * (pp/(p**2))* (1 - (( j / (p*g_theta) )**2) )
          
        ! n * D * m = dfdsig * D * dgdsig
          Call FormDEMCC(Sigu, xkappa, xNu, xe0, D, 6, xG, xK)
          result1 = matmul(D, dgdsig)
          result2 = dot_product(dfdsig, result1)
          
        ! dlambda
          dlambda = F / (result2 + A)
          
        ! Update the stress 
          Sigu = Sigu - (dlambda*result1)
          
        ! Acc plastic strain
          EpsPu = EpsPu + (dlambda * dgdsig)

        ! Update the state parameters (pp)
          pp = pp * exp(zeta * dlambda * dgdp) 
          
        ! Calc the yield function value
          iOpt = 1
          call PrincipalSig(iOpt, Sigu, xN1, xN2, xN3, S1, S2, S3, p,&
              q,j,theta)
          F = yield(p, j, pp, theta, xphics)         
          
        ! Update the counter
          counter = counter + 1
               
      end do 
      
    ! Return the integrated parameters
      Sig = Sigu
      dEpsP = EpsPu-EpsP
      EpsP = EpsPu

      end subroutine implicit_predictor_corrector_integration