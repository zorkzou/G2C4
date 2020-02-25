!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! G2C4: an interface between CFour and Gaussian via the "External" keyword of Gaussian.
!
! Wenli Zou,  Email: qcband@gmail.com
! Institute of Modern Physics, Northwest University, Xi'an, China
!
! Feb. 23, 2020
!
! 1. CFour ver. 2.00beta and 2.1 have been tested. Old versions do no work.
! 2. The Cartesian geometry saved in CFour's ZMAT file is in Bohr, thus UNITS=1,COORDINATES=1 should always be set up.
! 3. Numerical gradient and numerical Hessian calculations by CFour are not supported. They may be set up in Gaussian (e.g. freq=
!    numer), then CFour will be requested to do energy or analytic gradient calculations.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program G2C4
  implicit real(kind=8) (a-h,o-z)
  parameter (iGIN=41, iGOU=42, iCTP=51, iCIN=52, iCOU=61, iC_GRD=62, iC_DIP=63, iC_FCM=64, iC_APT=65, iC_POL=66)
  character*200 :: ctmp,tag
  allocatable   :: IZAG(:), XYZG(:), DIPG(:), POLG(:), APTG(:), GRDG(:), FCMG(:)
  allocatable   :: IZAC(:), XYZC(:), DIPC(:), POLC(:), APTC(:), GRDC(:), FCMC(:)
  allocatable   :: RMAT(:), Scr(:)

!---------------------------------------------------------------------------------------------------------------------------------
! 1. read arguments
!---------------------------------------------------------------------------------------------------------------------------------
  call RdArg(Imode,iGIN,iGOU,iCTP,iCIN,iCOU,iC_GRD,iC_DIP,iC_FCM,iC_APT,iC_POL,ctmp)

!---------------------------------------------------------------------------------------------------------------------------------
! 2. read *.EIn
!---------------------------------------------------------------------------------------------------------------------------------
  call RdEIn_Line1(iGIN,natom,nder,icharge,multi)
  na3 = 3*natom
  allocate(IZAG(natom), XYZG(na3))
  call RdEIn_XYZ(iGIN,natom,IZAG,XYZG)

!---------------------------------------------------------------------------------------------------------------------------------
! 3. set array length
!---------------------------------------------------------------------------------------------------------------------------------
  ! a special case: nder = 0 if natom = 1
  if(natom == 1) nder = 0
  ntt = na3*(na3+1)/2
  nss = na3*na3

  if(nder == 0) then
    lengrd = 1
    lenfcm = 1
    lenapt = 1
  else if(nder == 1) then
    lengrd = na3
    lenfcm = 1
    lenapt = 1
  else if(nder == 2) then
    lengrd = na3
    lenfcm = nss
    lenapt = 3*na3
  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 4. read data from CFour output & scratch files (Imode = 1)
!---------------------------------------------------------------------------------------------------------------------------------
  if(Imode == 1) then
    allocate(IZAC(natom), XYZC(na3), DIPC(3), POLC(9), GRDC(lengrd), FCMC(lenfcm), APTC(lenapt))
    XYZC = 0.0d0
    DIPC = 0.0d0
    POLC = 0.0d0
    GRDC = 0.0d0
    FCMC = 0.0d0
    APTC = 0.0d0

    call RdC4Ene(iCOU,Energy,ctmp,tag)

    if(nder > 0 .and. natom > 1) then
      call RdC4GRD(iC_GRD,natom,IZAC,XYZC,GRDC)
      call RdC4DIP(iC_DIP,DIPC)
    end if

    if(nder > 1 .and. natom > 1) then
      call RdC4FCM(iC_FCM,natom,FCMC)
      call RdC4APT(iC_APT,natom,APTC)
      call RdC4POL(iC_POL,POLC)
    end if

  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 5. Prepare new data for Gaussian / CFour
!---------------------------------------------------------------------------------------------------------------------------------
  if(Imode == 1) then
    allocate(DIPG(3), GRDG(lengrd), FCMG(lenfcm), APTG(lenapt), POLG(9))
    DIPG = DIPC
    GRDG = GRDC
    FCMG = FCMC
    APTG = APTC
    POLG = POLC

    if(nder > 0) then

      allocate(RMAT(9), Scr(ntt))

      ! Check: symmetry-equivalent atoms may be reordered by CFour if symmetry being used, which may lead to wrong results in
      ! gradient and Hessian calculations.
      call ckorder(iCOU,natom,IZAC,XYZC,Scr,ctmp,tag)

      ! rotate the dipole moment vector and gadients to initial orientation
      ! calculate rotation matrix
      call qrotmol(0,6,natom,XYZG,XYZC,RMAT,Scr)

      call rotvec(1,-1,RMAT,DIPG,Scr)
      call rotvec(natom,-1,RMAT,GRDG,Scr)

      if(nder > 1) then
        ! rotate hessian, apt, and polar to initial orientation
        call rothess(natom,-1,RMAT,FCMG,Scr)
        call Sq2Tr(na3,FCMG,Scr)
        call rotmat(natom,-1,RMAT,APTG,Scr)
        call rotmat(1,-1,RMAT,POLG,Scr)
        call Sq2Tr(3,POLG,Scr)
      end if

      deallocate(RMAT, Scr)
    end if

  else
    allocate(IZAC(natom), XYZC(na3))
    IZAC = IZAG
    XYZC = XYZG

  end if

!---------------------------------------------------------------------------------------------------------------------------------
! 6. Generate Gaussian's *.EOu file (Imode = 1) or CFour's ZMAT file (Imode = 0)
!---------------------------------------------------------------------------------------------------------------------------------
  if(Imode == 1) then
    call WrtEOU(iGOU,nder,natom,Energy,DIPG,GRDG,FCMG,APTG,POLG)

  else
    call WrtCar(iCTP,iCIN,natom,nder,icharge,multi,IZAC,XYZC,ctmp)

  end if

end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Check whether symmetry-equivalent atoms have been reordered by CFour.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ckorder(iCOU,natom,IZA,XYZ,Scr,ctmp,tag)
  implicit real(kind=8) (a-h,o-z)
  dimension :: IZA(natom), XYZ(3,natom), Scr(3)
  character*100 :: ctmp,tag

  rewind(iCOU)
  tag=' Z-matrix   Atomic            Coordinates (in bohr)'
  do while(.true.)
    read(iCOU,"(a100)",iostat=ist) ctmp
      if(ist /= 0) call XError("No Cartesian coordinates found in the CFour output file.")
    if(index(ctmp,tag(1:51)) > 0) exit
  end do

  read(iCOU,"(/)",iostat=ist)
    if(ist /= 0) call XError("Please check Cartesian coordinates in the CFour output file.")

  do i=1,natom
    read(iCOU,"(a100)",iostat=ist) ctmp
      if(ist /= 0) call XError("Please check Cartesian coordinates in the CFour output file.")
    read(ctmp(13:),*) IZtmp, Scr
    if(IZtmp < 1) call XError("Dummy atom or ghost atom cannot be used.")
    if(IZtmp /= IZA(i) .or. distance(XYZ(1,i),Scr) > 1.d-6) call XError("Atoms have been reordered in the CFour calculation.")
  end do

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! distance between points a and b
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distance(a,b)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: a(3),b(3),c(3)

  c = a-b
  distance = sqrt( dot_product(c,c) )

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Write CFour's ZMAT file in Cartesian coordinates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WrtCar(iCTP,iCIN,natom,nder,icharge,multi,IZA,XYZ,ctmp)
  implicit real(kind=8) (a-h,o-z)
  dimension :: IZA(natom), XYZ(3,natom)
  character*200 :: ctmp

  rewind(iCIN)
  rewind(iCTP)
  do while(.true.)
    read(iCTP,"(a200)",iostat=ist) ctmp
    if(ist /= 0) then
      if(nder == 0) then
        call XError("No $SP found in CFour templet.")
      else if(nder == 1) then
        call XError("No $GRAD found in CFour templet.")
      else
        call XError("No $FREQ found in CFour templet.")
      end if
    end if

    if(ctmp(1:1) == '!') cycle

    if(nder == 0 .and. index(ctmp,"$sp") > 0) exit
    if(nder == 1 .and. index(ctmp,"$grad") > 0) exit
    if(nder == 2 .and. index(ctmp,"$freq") > 0) exit
  end do

  loop1: do while(.true.)
    read(iCTP,"(a200)",iostat=ist) ctmp
      if(ist /= 0) call XError("No end_of_calc found in CFour templet.")

    if(ctmp(1:1) == '!') cycle

    if(index(ctmp,"file_block") > 0) then
      loop2: do while(.true.)
        read(iCTP,"(a200)",iostat=ist) ctmp
          if(ist /= 0) call XError("Please check file_block in CFour templet.")
        if(index(ctmp,"end_of_file") > 0) exit loop2
        write(iCIN,"(a)") trim(ctmp)
      end do loop2
    end if

    if(index(ctmp,"calc_block") > 0) then
      ! Title
      write(iCIN,"('Input file generated by G2C4.')")

      ! Cartesian geometry
      do i=1,natom
        call ElemZA(1,ctmp,IZA(i),ctmp)
        write(iCIN,"(a3,3f20.12)") ctmp(1:3), XYZ(:,i)
      end do

      write(iCIN,*)
      loop3: do while(.true.)
        read(iCTP,"(a200)",iostat=ist) ctmp
          if(ist /= 0) call XError("Please check calc_block in CFour templet.")
        if(index(ctmp,"end_of_calc") > 0) exit loop1

        ! replace {C}
        idxc = index(ctmp,'CHARGE={C}') + index(ctmp,'CHARGE={c}')
        if(idxc > 0) write(ctmp(idxc+7:idxc+9),"(i3.3)") icharge
        ! replace {M}
        idxm = index(ctmp,'MULT={M}') + index(ctmp,'MULT={m}')
        if(idxm > 0) write(ctmp(idxm+5:idxm+7),"(i3.3)") multi
        ! replace {R}
        idxm = index(ctmp,'REFERENCE={R}') + index(ctmp,'REFERENCE={r}')
        if(idxm > 0 .and. multi == 1) write(ctmp(idxm+10:idxm+12),"('RHF')")
        if(idxm > 0 .and. multi /= 1) write(ctmp(idxm+10:idxm+12),"('UHF')")

        write(iCIN,"(a)") trim(ctmp)
      end do loop3
    end if

  end do loop1

  write(iCIN,"(///)")

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! symmetric square matrix --> lower triangular matrix
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Sq2Tr(N,S,Scr)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: Scr(*),S(N,N)

  ntt = N*(N+1)/2
  ii=0
  Do i=1,N
    Do j=1,i
      ii=ii+1
      Scr(ii)=(S(j,i)+S(i,j))*0.5d0
    end Do
  end Do

  call Copy(ntt,Scr,S)

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! copy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Copy(n,A,B)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: A(n),B(n)

  B = A

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rotate the Hessian matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rothess(natom,mode,ROT,FCM,SCR)
  implicit double precision (a-h,o-z)
  dimension FCM(3,natom,3,natom),ROT(9),SCR(3,3,2)

  do i=1,natom
    do j=1,natom
      SCR(:,:,1) = FCM(:,j,:,i)
      call rotmat(1,mode,ROT,SCR(1,1,1),SCR(1,1,2))
      FCM(:,j,:,i) = SCR(:,:,1)
    end do
  end do

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! It performs a 3x3 rotation operation of Cartesian matrices:
! mode >=0: M' = ROT^T * M * ROT
!      < 0: M' = ROT * M * ROT^T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotmat(NM,mode,ROT,AM,SCR)
  implicit double precision (a-h,o-z)
  dimension AM(9,NM),ROT(9),SCR(9)

  if(mode >= 0) then
    do i=1,NM
      call dgemm('n','n',3,3,3,1.0d0,AM(1,i),3,ROT,3,0.0d0,SCR,3)
      call dgemm('t','n',3,3,3,1.0d0,ROT,3,SCR,3,0.0d0,AM(1,i),3)
    end do
  else
    do i=1,NM
      call dgemm('n','t',3,3,3,1.0d0,AM(1,i),3,ROT,3,0.0d0,SCR,3)
      call dgemm('n','n',3,3,3,1.0d0,ROT,3,SCR,3,0.0d0,AM(1,i),3)
    end do
  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! It performs a 3x3 rotation operation of Cartesian vectors:
! mode >=0: XYZ' = XYZ * ROT
!      < 0: XYZ' = XYZ * ROT^-1 = XYZ * ROT^T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rotvec(NV,mode,ROT,XYZ,SCR)
  implicit double precision (a-h,o-z)
  dimension XYZ(3,NV),ROT(3,3),SCR(3)

  if(mode >= 0) then
    do i=1,NV
      SCR = XYZ(:,i)
      XYZ(:,i) = 0.0d0
      do j=1,3
        do k=1,3
          XYZ(j,i)=XYZ(j,i)+ROT(k,j)*SCR(k)
        end do
      end do
    end do
  else
    do i=1,NV
      SCR = XYZ(:,i)
      XYZ(:,i) = 0.0d0
      do j=1,3
        do k=1,3
          XYZ(j,i)=XYZ(j,i)+ROT(j,k)*SCR(k)
        end do
      end do
    end do
  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This subroutine determines the rotation matrix for the best fit superimposition of two molecular Coordinates.
!
! The basic method was described in S. K. Kearsley, Acta Cryst. A45, 208 (1989), and coded in the PDBSUP program by B.Rupp and
! S.Parkin at LLNL (1996).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine qrotmol(iprint,iout,natm,cartarget,carprobe,rmat,scr)
  implicit real(kind=8) (a-h,o-z)
  dimension cartarget(3,natm), carprobe(3,natm), rmat(3,3), scr(4,4)
  allocatable :: dxp(:,:), dxm(:,:), qmat(:,:)

  allocate(dxp(3,natm), dxm(3,natm), qmat(4,4))

  ! translate the molecules to their geometric centers
  call shift(natm,scr,cartarget)
  call shift(natm,scr,carprobe)
  if(iprint > 0) then
    write(iout,"(' Cartesian coordinates of target molecule:')")
    do i=1,natm
      write(iout,"(3f20.12)")cartarget(:,i)
    end do
    write(iout,"(' Cartesian coordinates of probe molecule:')")
    do i=1,natm
      write(iout,"(3f20.12)")carprobe(:,i)
    end do
  end if

  ! coordinate differences: plus and minus
  dxp = carprobe + cartarget
  dxm = carprobe - cartarget

  ! construct Kearsley's Q-matrix
  call conqmt(iprint,iout,natm,dxp,dxm,qmat)

  ! diagonalize Q: Q * L = L * A; A --> qmat(:,1), L --> scr
  call DiagS1(iprint,iout,4,qmat,scr)

  ! construct the best fit rotation matrix using the eigenvectors
  call conrot(iprint,iout,scr,rmat)

  deallocate(dxp, dxm, qmat)

end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Constructing the best fit rotation matrix. See the 2nd equation in
! S.K.Kearsley, On the orthogonal transformation used for structural comparisons, Acta Cryst. A45, 208 (1989)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine conrot(iprint,iout,q,r)
  implicit real(kind=8) (a-h,o-z)
  dimension         :: q(4), r(3,3)

  r(1,1) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
  r(2,1) = 2 * (q(2)*q(3) + q(1)*q(4))
  r(3,1) = 2 * (q(2)*q(4) - q(1)*q(3))
  r(1,2) = 2 * (q(2)*q(3) - q(1)*q(4))
  r(2,2) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
  r(3,2) = 2 * (q(3)*q(4) + q(1)*q(2))
  r(1,3) = 2 * (q(2)*q(4) + q(1)*q(3))
  r(2,3) = 2 * (q(3)*q(4) - q(1)*q(2))
  r(3,3) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)

  if(iprint > 0) then
    write(iout,"(/,' Rotation matrix for the best fit:')")
    do i = 1, 3
      write(iout,"(4x,3f16.6)") r(i,:)
    end do
  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Solve A * L = L * e, where A is symmetric, e will be saved in the first colum of A whereas L is saved in W.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DiagS1(iprint,iout,N,A,W)
  implicit real(kind=8) (a-h,o-z)
  dimension :: A(N,N), W(N,N)

  Call CJCBJ(N,A,W)

  if(iprint > 0) then
    write(iout,"(/,' Eigenvalues and Eigenvectors of Q:')")
    write(iout,"(4x,4f16.6)") A(:,1)
    write(iout,"(4x,4('      ----------'))")
    do i = 1, N
      write(iout,"(4x,4f16.6)") W(i,:)
    end do
    write(iout,"(' where the smallest eigenvalue represents s.r.s. of the best fit.')")
  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reference BLAS level3 routine (version 3.4.0) --
! Reference BLAS is a software package provided by Univ. of Tennessee,    --
! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! November 2011
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

  ! Scalar Arguments ..
  DOUBLE PRECISION ALPHA,BETA
  INTEGER K,LDA,LDB,LDC,M,N
  CHARACTER TRANSA,TRANSB

  ! Array Arguments ..
  DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)

  ! External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME

  ! Intrinsic Functions ..
  INTRINSIC MAX

  ! Local Scalars ..
  DOUBLE PRECISION TEMP
  INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
  LOGICAL NOTA,NOTB

  ! Parameters ..
  DOUBLE PRECISION ONE,ZERO
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

  ! Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not transposed and set  NROWA, NCOLA
  ! and NROWB as the number of rows and  columns of  A  and the  number of  rows  of  B  respectively.
  NOTA = LSAME(TRANSA,'N')
  NOTB = LSAME(TRANSB,'N')
  IF (NOTA) THEN
      NROWA = M
      NCOLA = K
  ELSE
      NROWA = K
      NCOLA = M
  END IF
  IF (NOTB) THEN
      NROWB = K
  ELSE
      NROWB = N
  END IF

  ! Test the input parameters.
  INFO = 0
  IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. (.NOT.LSAME(TRANSA,'T'))) THEN
      INFO = 1
  ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. (.NOT.LSAME(TRANSB,'T'))) THEN
      INFO = 2
  ELSE IF (M<0) THEN
      INFO = 3
  ELSE IF (N<0) THEN
      INFO = 4
  ELSE IF (K<0) THEN
      INFO = 5
  ELSE IF (LDA<MAX(1,NROWA)) THEN
      INFO = 8
  ELSE IF (LDB<MAX(1,NROWB)) THEN
      INFO = 10
  ELSE IF (LDC<MAX(1,M)) THEN
      INFO = 13
  END IF
  IF (INFO/=0) THEN
      WRITE(*,"(/,' Parameter number ',i2,' is illegal.')") INFO
      call XError("DGEMM has a problem.")
  END IF

  ! Quick return if possible.
  IF ((M==0) .OR. (N==0) .OR. (((ALPHA==ZERO).OR. (K==0)).AND. (BETA==ONE))) RETURN

  ! And if  alpha==zero.
  IF (ALPHA==ZERO) THEN
      IF (BETA==ZERO) THEN
          DO J = 1,N
              DO I = 1,M
                  C(I,J) = ZERO
              END DO
          END DO
      ELSE
          DO J = 1,N
              DO I = 1,M
                  C(I,J) = BETA*C(I,J)
              END DO
          END DO
      END IF
      RETURN
  END IF

  ! Start the operations.
  IF (NOTB) THEN
      IF (NOTA) THEN

          ! Form  C := alpha*A*B + beta*C.
          DO J = 1,N
              IF (BETA==ZERO) THEN
                  DO I = 1,M
                      C(I,J) = ZERO
                  END DO
              ELSE IF (BETA/=ONE) THEN
                  DO I = 1,M
                      C(I,J) = BETA*C(I,J)
                  END DO
              END IF
              DO L = 1,K
                  IF (B(L,J)/=ZERO) THEN
                      TEMP = ALPHA*B(L,J)
                      DO I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
                      END DO
                  END IF
              END DO
          END DO
      ELSE

          ! Form  C := alpha*A**T*B + beta*C
          DO J = 1,N
              DO I = 1,M
                  TEMP = ZERO
                  DO L = 1,K
                      TEMP = TEMP + A(L,I)*B(L,J)
                  END DO
                  IF (BETA==ZERO) THEN
                      C(I,J) = ALPHA*TEMP
                  ELSE
                      C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
              END DO
          END DO
      END IF
  ELSE
      IF (NOTA) THEN

          ! Form  C := alpha*A*B**T + beta*C
          DO J = 1,N
              IF (BETA==ZERO) THEN
                  DO I = 1,M
                      C(I,J) = ZERO
                  END DO
              ELSE IF (BETA/=ONE) THEN
                  DO I = 1,M
                      C(I,J) = BETA*C(I,J)
                  END DO
              END IF
              DO L = 1,K
                  IF (B(J,L)/=ZERO) THEN
                      TEMP = ALPHA*B(J,L)
                      DO I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
                      END DO
                  END IF
              END DO
          END DO
      ELSE

          ! Form  C := alpha*A**T*B**T + beta*C
          DO J = 1,N
              DO I = 1,M
                  TEMP = ZERO
                  DO L = 1,K
                      TEMP = TEMP + A(L,I)*B(J,L)
                  END DO
                  IF (BETA==ZERO) THEN
                      C(I,J) = ALPHA*TEMP
                  ELSE
                      C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
              END DO
          END DO
      END IF
  END IF

  RETURN
END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Reference BLAS level1 routine (version 3.1) --
! Reference BLAS is a software package provided by Univ. of Tennessee,    --
! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! November 2011
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LOGICAL FUNCTION LSAME(CA,CB)

  ! Scalar Arguments ..
  CHARACTER CA,CB

  ! Intrinsic Functions ..
  INTRINSIC ICHAR

  ! Local Scalars ..
  INTEGER INTA,INTB,ZCODE

  ! Test if the characters are equal
  LSAME = CA == CB
  IF (LSAME) RETURN

  ! Now test for equivalence if both characters are alphabetic.
  ZCODE = ICHAR('Z')

  ! Use 'Z' rather than 'A' so that ASCII can be detected on Prime
  ! machines, on which ICHAR returns a value with bit 8 set.
  ! ICHAR('A') on Prime machines returns 193 which is the same as
  ! ICHAR('A') on an EBCDIC machine.
  INTA = ICHAR(CA)
  INTB = ICHAR(CB)

  IF (ZCODE==90 .OR. ZCODE==122) THEN

      ! ASCII is assumed - ZCODE is the ASCII code of either lower or upper case 'Z'.
      IF (INTA>=97 .AND. INTA<=122) INTA = INTA - 32
      IF (INTB>=97 .AND. INTB<=122) INTB = INTB - 32

  ELSE IF (ZCODE==233 .OR. ZCODE==169) THEN

      ! EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or upper case 'Z'.
      IF (INTA>=129 .AND. INTA<=137 .OR. INTA>=145 .AND. INTA<=153 .OR. INTA>=162 .AND. INTA<=169) INTA = INTA + 64
      IF (INTB>=129 .AND. INTB<=137 .OR. INTB>=145 .AND. INTB<=153 .OR. INTB>=162 .AND. INTB<=169) INTB = INTB + 64

  ELSE IF (ZCODE==218 .OR. ZCODE==250) THEN

      ! ASCII is assumed, on Prime machines - ZCODE is the ASCII code plus 128 of either lower or upper case 'Z'.
      IF (INTA>=225 .AND. INTA<=250) INTA = INTA - 32
      IF (INTB>=225 .AND. INTB<=250) INTB = INTB - 32
  END IF
  LSAME = INTA == INTB

  RETURN
END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Diagonalize a real symmetric matrix.
! The subroutines were taken from S. L. Xu's FORTRAN Algorithm book, Tsinghua, 1995.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CJCBJ(N,A,V)
  implicit real(kind=8) (a-h,o-z)
  parameter(EPS=1.0d-8,Zero=0.0d0,One=1.0d0,Two=2.0d0)
  DIMENSION A(N,N), V(N,N)
  INTEGER P,Q

  DO 20 I=1,N
    V(I,I)=One
    DO 10 J=1,N
      IF (I /= J) V(I,J)=Zero
    10  CONTINUE
  20  CONTINUE
  FF=Zero

  DO 500 I=2,N
  DO 500 J=1,I-1
  500  FF=FF+A(I,J)*A(I,J)

!ooo  FF=SQRT(Two*FF)
  FF=SQRT(Two*ABS(FF))
  205  FF=FF/(One*N)
  25  DO 30 I=2,N
  DO 30 J=1,I-1
!ooo    IF (ABS(A(I,J)) >= FF) THEN
    IF ((ABS(A(I,J)) - FF) >= EPS) THEN
      P=I
      Q=J
      GOTO 600
    END IF
  30  CONTINUE

  IF (FF >= EPS) GOTO 205
  goto 1000

  600  X=-A(P,Q)
  Y=(A(Q,Q)-A(P,P))/Two
  OMEGA=X/SQRT(X*X+Y*Y)
  IF (Y < Zero) OMEGA=-OMEGA
  SN=One+SQRT(One-OMEGA*OMEGA)
  SN=OMEGA/SQRT(Two*SN)
  CN=SQRT(One-SN*SN)
  FM=A(P,P)
  A(P,P)=FM*CN*CN+A(Q,Q)*SN*SN+A(P,Q)*OMEGA
  A(Q,Q)=FM*SN*SN+A(Q,Q)*CN*CN-A(P,Q)*OMEGA
  A(P,Q)=Zero
  A(Q,P)=Zero

  DO 60 J=1,N
    IF ((J /= P).AND.(J /= Q)) THEN
      FM=A(P,J)
      A(P,J)= FM*CN+A(Q,J)*SN
      A(Q,J)=-FM*SN+A(Q,J)*CN
    END IF
  60  CONTINUE

  DO 70 I=1,N
    IF ((I /= P).AND.(I /= Q)) THEN
      FM=A(I,P)
      A(I,P)= FM*CN+A(I,Q)*SN
      A(I,Q)=-FM*SN+A(I,Q)*CN
    END IF
  70  CONTINUE

  DO 80 I=1,N
    FM=V(I,P)
    V(I,P)= FM*CN+V(I,Q)*SN
    V(I,Q)=-FM*SN+V(I,Q)*CN
  80  CONTINUE

  GOTO 25

  1000  do i=2,N
    A(i,1) = A(i,i)
  end do
  call eigsrt(N,A,V)

  RETURN
END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Sort eigenvectors and eigenvalues.
! The subroutines were taken from G. Y. He and Y. L. Gao's FORTRAN Algorithm book, Science Press, 2002.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE eigsrt(n,d,v)
 implicit real(kind=8) (a-h,o-z)
 dimension d(n),v(n,n)

 do i=1,n-1
   k=i
   p=d(i)
   do j=i+1,n
!    if(d(j)>=p) then  ! from large to small
     if(d(j)<=p) then  ! from small to large
       k=j
       p=d(j)
     endif
   end do
   if(k/=i) then
     d(k)=d(i)
     d(i)=p
     do j=1,n
       p=v(j,i)
       v(j,i)=v(j,k)
       v(j,k)=p
     end do
   endif
 end do

 return
END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! construct Kearsley's Q-matrix. See the last equation in
! S.K.Kearsley, On the orthogonal transformation used for structural comparisons, Acta Cryst. A45, 208 (1989)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine conqmt(iprint,iout,n,xp,xm,q)
  implicit real(kind=8) (a-h,o-z)
  dimension         :: xp(3,n), xm(3,n), q(4,4)

  q = 0.0d0

  do i = 1, n
    q(1,1) = q(1,1) + xm(1,i)*xm(1,i) + xm(2,i)*xm(2,i) + xm(3,i)*xm(3,i)
    q(2,2) = q(2,2) + xp(2,i)*xp(2,i) + xp(3,i)*xp(3,i) + xm(1,i)*xm(1,i)
    q(3,3) = q(3,3) + xp(1,i)*xp(1,i) + xp(3,i)*xp(3,i) + xm(2,i)*xm(2,i)
    q(4,4) = q(4,4) + xp(1,i)*xp(1,i) + xp(2,i)*xp(2,i) + xm(3,i)*xm(3,i)

    q(2,1) = q(2,1) + xp(2,i)*xm(3,i) - xm(2,i)*xp(3,i)
    q(3,1) = q(3,1) + xm(1,i)*xp(3,i) - xp(1,i)*xm(3,i)
    q(4,1) = q(4,1) + xp(1,i)*xm(2,i) - xm(1,i)*xp(2,i)
    q(3,2) = q(3,2) + xm(1,i)*xm(2,i) - xp(1,i)*xp(2,i)
    q(4,2) = q(4,2) + xm(1,i)*xm(3,i) - xp(1,i)*xp(3,i)
    q(4,3) = q(4,3) + xm(2,i)*xm(3,i) - xp(2,i)*xp(3,i)

    q(1,2) = q(2,1)
    q(1,3) = q(3,1)
    q(1,4) = q(4,1)
    q(2,3) = q(3,2)
    q(2,4) = q(4,2)
    q(3,4) = q(4,3)
  end do

  if(iprint > 0) then
    write(iout,"(/,' Q-matrix:')")
    write(iout,"(4x,4f16.6)") q
  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! shift the coordinates to the geometric center
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Shift(n,o,xyz)
  implicit real(kind=8) (a-h,o-z)
  dimension :: xyz(3,*), o(3)

  o = 0.0d0
  do i=1,n
    o = o + xyz(:,i)
  end do
  o = o /dble(n)
  do i=1,n
    xyz(:,i) = xyz(:,i) - o
  end do

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Write Gaussian's *.EOu file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine WrtEOU(iGOU,nder,natom,Energy,DIP,GRD,FCM,APT,POL)
  Implicit Real*8(A-H,O-Z)
  dimension :: DIP(*), GRD(*), FCM(*), APT(*), POL(*)

  na3 = 3*natom
  ntt = na3*(na3+1)/2

  rewind(iGOU)
  if(nder == 0) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (0.0d0, i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (0.0d0, i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,ntt)

  else if(nder == 1) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (DIP(i), i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (GRD(i), i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (0.0d0, i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (0.0d0, i=1,ntt)

  else if(nder == 2) then
    ! energy, dipole moment
    write(iGOU,"(4d20.12)") Energy, (DIP(i), i=1,3)
    ! gradients
    write(iGOU,"(3d20.12)") (GRD(i), i=1,na3)
    ! polarizability (l.t. part)
    write(iGOU,"(3d20.12)") (POL(i), i=1,6)
    ! apt
    write(iGOU,"(3d20.12)") (APT(i), i=1,3*na3)
    ! hessian matrix (l.t. part)
    write(iGOU,"(3d20.12)") (FCM(i), i=1,ntt)

  end if

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read APT data from CFour's DIPDER file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdC4APT(iC_APT,natom,APT)
  Implicit Real*8(A-H,O-Z)
  dimension :: APT(3,3,natom)
  logical :: ifopen

  inquire(unit=iC_APT,opened=ifopen)
    if(.NOT.ifopen) return

  rewind(iC_APT)

  do ixyz = 1,3
    read(iC_APT,*,iostat=ist)
      if(ist /= 0) call XError("Please check the DIPDER file.")
    do iatom = 1,natom
      read(iC_APT,*,iostat=ist) ZA, APT(ixyz,:,iatom)
        if(ist /= 0) call XError("Please check the APT data in the DIPDER file.")
    end do
  end do

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read Hessian matrix from CFour's FCMFINAL file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdC4FCM(iC_FCM,natom,FCM)
  Implicit Real*8(A-H,O-Z)
  dimension :: FCM(*)
  logical :: ifopen

  inquire(unit=iC_FCM,opened=ifopen)
    if(.NOT.ifopen) call XError("-C_FCM is required for hessian calculations.")

  rewind(iC_FCM)
  read(iC_FCM,*,iostat=ist)
    if(ist /= 0) call XError("Please check the FCMFINAL file.")

  read(iC_FCM,*,iostat=ist) FCM(1:9*natom*natom)
    if(ist /= 0) call XError("Please check the Hessian matrix in the FCMFINAL file.")

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read Cartesian coordinates and gradients from CFour's GRD file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdC4GRD(iC_GRD,natom,IZA,XYZ,GRD)
  Implicit Real*8(A-H,O-Z)
  dimension :: IZA(natom), XYZ(3,natom), GRD(3,natom)
  logical :: ifopen

  inquire(unit=iC_GRD,opened=ifopen)
    if(.NOT.ifopen) call XError("-C_GRD is required for gradient and hessian calculations.")

  rewind(iC_GRD)
  read(iC_GRD,*,iostat=ist)
    if(ist /= 0) call XError("Please check the GRD file.")

  do i=1,natom
    read(iC_GRD,*,iostat=ist) ZA, XYZ(:,i)
      if(ist /= 0) call XError("Please check Cartesian coordinates in the GRD file.")
    IZA(i) = nint(ZA)
  end do

  do i=1,natom
    read(iC_GRD,*,iostat=ist) ZA, GRD(:,i)
      if(ist /= 0) call XError("Please check Cartesian gradients in the GRD file.")
  end do

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read polarizability tensor from CFour's POLAR file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdC4POL(iC_POL,POL)
  Implicit Real*8(A-H,O-Z)
  dimension :: POL(9)
  logical :: ifopen

  inquire(unit=iC_POL,opened=ifopen)
    if(.NOT.ifopen) return

  rewind(iC_POL)
  read(iC_POL,*,iostat=ist) POL
    if(ist /= 0) call XError("No polarizability tensor found in the POLAR file.")

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read dipole moment vector from CFour's DIPOL file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdC4Dip(iC_DIP,DIP)
  Implicit Real*8(A-H,O-Z)
  dimension :: DIP(3)
  logical :: ifopen

  inquire(unit=iC_DIP,opened=ifopen)
    if(.NOT.ifopen) return

  rewind(iC_DIP)
  read(iC_DIP,*,iostat=ist) DIP
    if(ist /= 0) call XError("No dipole moment vector found in the DIPOL file.")

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read total energy from CFour's output file
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdC4Ene(iCOU,ENE,ctmp,tag)
  Implicit Real*8(A-H,O-Z)
  character*200 :: ctmp,tag

  ENE = 0.0d0

  ! read total energy
  rewind(iCOU)
  tag="  The final electronic energy is"
  do while(.true.)
    read(iCOU,"(a64)",iostat=ist) ctmp
      if(ist /= 0) call XError("No total energy found at the end of the CFour output file.")
    if(index(ctmp,tag(1:32)) > 0) then
      read(ctmp(33:64),*) ENE
      exit
    end if
  end do

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read Cartesian coordinates (in a.u.) from Gaussian's *.EIn
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdEIn_XYZ(iGIN,natom,IZA,XYZ)
  Implicit Real*8(A-H,O-Z)
  dimension :: IZA(natom), XYZ(3,natom)

  IZA = 0
  XYZ = 0.0d0

  do i = 1, natom
    read(iGIN,*,iostat=ist) IZA(i), XYZ(:,i)
    if(ist /=0) call XError("Please check the Cartesian coordinates in *.EIn.")
  end do

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read the first line of Gaussian's *.EIn
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine RdEIn_Line1(iGIN,natom,nder,icharge,multi)
  Implicit Real*8(A-H,O-Z)

  rewind(iGIN)

  natom=-1
  nder=-1
  icharge=-1
  multi=-1

  read(iGIN,*,iostat=ist) natom,nder,icharge,multi
  if(ist /=0) call XError("Please check the first line in *.EIn.")
  if(natom < 1) call XError("Natom < 1.")
  if(nder < 0 .or. nder > 2) call XError("Nder is out of range.")
  if(multi < 1) call XError("Multi cannot be smaller than 1.")

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read arguments:  -G2C (imode = 0; default) or -C2G (imode = 1)
!
! -G2C: Gaussian's *.EIn file name (via -GIN), CFour's ZMAT file name (via -CIN), and a CFour templet file name (via -CTP) should
!       be provided.
! -C2G: Gaussian's *.EIn file name (via -GIN) and *.EOu file name (via -GOU), and CFour's output file name (via -COU to provide
!       total energy) and some of the following scratch files if nder > 0 should be provided,
!       nder = 1 : -C_GRD, and -C_DIP (optional)
!       nder = 2 : -C_GRD, -C_DIP (optional), -C_FCM, -C_APT (optional), and -C_POL (optional)
!       in which
!       -C_GRD is the GRD file with Cartesian coordinates and gradients,
!       -C_DIP is the DIPOL file with dipole moment vector,
!       -C_FCM is the FCMFINAL file with Cartesian hessian matrix,
!       -C_APT is the DIPDER file with dipole moment derivatives for IR intensity, and
!       -C_POL is the POLAR file with polarizability tensor.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RdArg(imode,iGIN,iGOU,iCTP,iCIN,iCOU,iC_GRD,iC_DIP,iC_FCM,iC_APT,iC_POL,ctmp)
  implicit real(kind=8) (a-h,o-z)
  character*200 :: ctmp
  logical :: ifopen

  imode = 0

  i = 0
  do while(.true.)
    i = i + 1
    call Get1Cmd(i,.True.,ctmp)
    istr=nonspace(ctmp)
    iend=len_trim(ctmp)

    if(iend == 0) then
      exit
    else if(ctmp(istr:iend) == '-C2G') then
      imode = 1
    else if(ctmp(istr:iend) == '-G2C') then
      imode = 0
    else if(ctmp(istr:iend) == '-GIN') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -GIN.")
      else
        open(iGIN,file=ctmp(istr:iend),status='old',iostat=ist)
        if(ist /= 0) call XError("Cannot open old file for -GIN.")
      end if
    else if(ctmp(istr:iend) == '-GOU') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -GOU.")
      else
        open(iGOU,file=ctmp(istr:iend),status='replace',iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -GOU.")
      end if
    else if(ctmp(istr:iend) == '-CTP') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -CTP.")
      else
        open(iCTP,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -CTP.")
      end if
    else if(ctmp(istr:iend) == '-CIN') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -CIN.")
      else
        open(iCIN,file=ctmp(istr:iend),status='replace',iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -CIN.")
      end if
    else if(ctmp(istr:iend) == '-COU') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -COU.")
      else
        open(iCOU,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -COU.")
      end if
    else if(ctmp(istr:iend) == '-C_GRD') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -C_GRD.")
      else
        open(iC_GRD,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -C_GRD.")
      end if
    else if(ctmp(istr:iend) == '-C_DIP') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -C_DIP.")
      else
        open(iC_DIP,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -C_DIP.")
      end if
    else if(ctmp(istr:iend) == '-C_FCM') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -C_FCM.")
      else
        open(iC_FCM,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -C_FCM.")
      end if
    else if(ctmp(istr:iend) == '-C_APT') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -C_APT.")
      else
        open(iC_APT,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -C_APT.")
      end if
    else if(ctmp(istr:iend) == '-C_POL') then
      i = i + 1
      call Get1Cmd(i,.False.,ctmp)
      istr=nonspace(ctmp)
      iend=len_trim(ctmp)
      if(iend == 0) then
        call XError("No file name provided for -C_POL.")
      else
        open(iC_POL,file=ctmp(istr:iend),iostat=ist)
        if(ist /= 0) call XError("Cannot open file for -C_POL.")
      end if
    else
      call XError("Unknown argument "//ctmp(istr:iend))
    end if
  end do

  if(Imode == 1) then
    inquire(unit=iGIN,opened=ifopen)
      if(.NOT.ifopen) call XError("-GIN is not defined!")
    inquire(unit=iGOU,opened=ifopen)
      if(.NOT.ifopen) call XError("-GOU is not defined!")
    inquire(unit=iCOU,opened=ifopen)
      if(.NOT.ifopen) call XError("-COU is not defined!")
  else
    inquire(unit=iGIN,opened=ifopen)
      if(.NOT.ifopen) call XError("-GIN is not defined!")
    inquire(unit=iCTP,opened=ifopen)
      if(.NOT.ifopen) call XError("-CTP is not defined!")
    inquire(unit=iCIN,opened=ifopen)
      if(.NOT.ifopen) call XError("-CIN is not defined!")
  end if

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an argument, in upper case if ifUP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Get1Cmd(i,ifUP,ctmp)
  Implicit Real*8(A-H,O-Z)
  logical ifUP
  character*200 ctmp

  call get_command_argument(i,ctmp)
  ! call getarg(i,ctmp)
  if(ifUP) call charl2u(ctmp)

  Return
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! position of the first non-space character in a string.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonspace(string)
  implicit double precision (a-h,o-z)
  character*(*) string

  length=LEN_TRIM(string)
  if(length <= 1) then
    i=length
  else
    do i=1,length
      if(string(i:i) /= ' ') exit
    end do
  endif

  nonspace=i

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! CHA --> cha
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charu2l(cha)
  implicit real(kind=8) (a-h,o-z)
  character*(*) :: cha
  character*1  :: U2L

  do i=1,len_trim(cha)
    cha(i:i)=U2L(cha(i:i))
  end do

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! cha --> CHA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charl2u(cha)
  implicit real(kind=8) (a-h,o-z)
  character*(*) :: cha
  character*1  :: L2U

  do i=1,len_trim(cha)
    cha(i:i)=L2U(cha(i:i))
  end do

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! L --> l
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U2L(letter)
  implicit real(kind=8) (a-h,o-z)
  character*1 :: letter,U2L

  if((ichar(letter) >= 65).and.(ichar(letter) <= 90))then
    U2L=char(ichar(letter)+32)
  else
    U2L=letter
  endif

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! l --> L
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L2U(letter)
  implicit real(kind=8) (a-h,o-z)
  character*1 :: letter,L2U

  if( ichar(letter) >= 97 .and. ichar(letter) <= 122 )then
    L2U=char(ichar(letter)-32)
  else
    L2U=letter
  endif

  return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! read an error message and stop
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine XError(inf)
  implicit real(kind=8) (a-h,o-z)
  character*(*) :: inf

  write(*,"(/,' *** Error! ',a)")trim(inf)

  stop

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Mode = 0 : returns nuclear charge zchar for an element symbol "el"; iza is not used.
!     /= 0 : returns element symbol "el" for nuclear charge iza; zchar is not used.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ElemZA(Mode,el,iza,zchar)
  implicit real(kind=8) (a-h,o-z)
  parameter (maxza=120)
  character*3 :: el,atomlib(maxza)
  data (atomlib(i),i=1,maxza) / &
   'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ',   'NA ','MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ','CA ', &
   'SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',   'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ', &
   'NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ','IN ','SN ',   'SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ', &
   'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ',   'LU ','HF ','TA ','W  ','RE ','OS ','IR ','PT ','AU ','HG ', &
   'TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',   'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ', &
   'MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ',   'RG ','CN ','NH ','FL ','MC ','LV ','TS ','OG ','UUE','UBN'/
  save atomlib

  if (Mode == 0) then

    call charl2u(el)
    zchar = 0.d0
    do i=1,maxza
      if(index(el,atomlib(i)) /= 0)then
        zchar = dble(i)
        exit
      end if
    end do

  else

    el = "???"
    if(iza > 0 .and. iza <= maxza) el = adjustl(atomlib(iza))
    call charu2l(el(2:3))

  end if

  return
end

!--- END

