! MIT License
! 
! Copyright (c) 2019 YuyaM
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  MAIN PRORGRAM  phaseAverage                             !
      !                                                          !
      !                           2019/03/08                     !
      !                           WRITTEN BY YUYA MIKI           !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !<
      !! @param ME :: maximum element number
      !! @param MP :: maximum node number
      !! @param N  :: the number of node in a element
      !! @param NODE :: node number list
      !! @param X,Y,Z :: Position
      !! @param U,V,W :: Velocity
      !! @param P     :: Pressure at element
      !! @param PN    :: Pressure at node
      !! @param numPhase :: phase average interval
      !! @param NSUM  :: Number of Flow in FLOWS
      program phaseAverage
      implicit none
!
      include 'gf2.h'
      ! parameter for GFALL
      DATA IWRITE / 2 /
      DATA INAME  / 1 /

      integer(4) :: ME,MP
      integer(4),parameter :: N=8
!
      integer(4), allocatable :: NODE(:,:)
      real(4),allocatable,dimension(:)   :: X,Y,Z,U,V,W,P,PN
      real(4),allocatable,dimension(:,:) :: UA,VA,WA,PA,PNA
!
      integer(4)    :: ISTEP,NPCHK,NEPRS,NPPRS
      real(4)       :: TIMEP,INVNUM
      character(60) :: FILEMS,FILEAR,FILEFF,FILEAV,FILEAVE,FILEFS
      integer(4)    :: numPhase,ISUM,iPhase,NSUM,IFLOWS,NFLOWS
!     [work]
      integer(4)    :: NP,NE,NDUM,IP,IE
      character(4)  :: CNUM
      ! fortran input output
      integer(4),parameter :: IUT0  = 0
      integer(4),parameter :: IUT5  = 5
      integer(4),parameter :: IUT6  = 6
      integer(4),parameter :: IUTMS = 11
      integer(4),parameter :: IUTFF = 12
      integer(4),parameter :: IUTAV = 13
      integer(4) :: IERR,FLAGOK,FLAGGFSEP,FLAGDELIM
!
! read MESH DATA
!
      write(IUT6,*)
      write(IUT6,*) '** phaseAverage: Phase averaging FLOWS data **'
!
      do
          write(IUT6,*)
          write(IUT6,*) '[In]Specify filename for MESH data'
          read (IUT5,'(A60)') FILEMS
          write(IUT6,*) '[In]Specify filename for FLOWS data'
          read (IUT5,'(A60)') FILEFF
          write(IUT6,*) '[Out]Specify filename for AVE  data'
          read (IUT5,'(A60)') FILEAV
          write(IUT6,*) 'Specify interval number for phase-average'
          read (IUT5,*) numPhase
          write(IUT6,*) 'if FLOWS IS SEPARATED    :1'
          write(IUT6,*) 'if FLOWS IS NOT SEPARATED:0'
          read (IUT5,*) FLAGGFSEP
          write(IUT6,*) 'Specify number of flowsfor phase-average'
          read (IUT5,*) NFLOWS
          write(IUT6,*) 'if want to change delimiter .P****'
          read (IUT5,*) FLAGDELIM
          IF(FLAGDELIM.eq.1)THEN
            write(IUT6,*) 'USE .P**** :0'
            write(IUT6,*) 'USE .t**** :1'
            read (IUT5,*) FLAGDELIM
          ENDIF

          write(IUT6,'(A10,A60)') "FILEMS   =", FILEMS
          write(IUT6,'(A10,A60)') "FILEFF   =", FILEFF
          write(IUT6,'(A10,A60)') "FILEAV   =", FILEAV
          write(IUT6,'(A10,I6)')  "numPhase    =", numPhase
          write(IUT6,'(A10,I6)')  "FLAGGFSEP=", FLAGGFSEP
          write(IUT6,'(A10,I6)')  "NFLOWS   =", NFLOWS
          write(IUT6,*) 'These parameters are ok? Yes:1,No:0'
          read(IUT5,*) FLAGOK
          if(FLAGOK.EQ.1) exit
      enddo
!
! check data size using MESH file
!
      call CHKMSH(FILEMS,IUTMS,IUT0,IUT6,MP,ME,IERR)
      if (IERR.NE.0) STOP
!
! allocate variables
!
      write(IUT6,*)
      write(IUT6,*)
      write(IUT6,*) 'phaseAverage: ALLOCATING VARIABLES... '
      allocate(   NODE(N,ME), STAT=LERR(01))
      allocate(        X(MP), STAT=LERR(02))
      allocate(        Y(MP), STAT=LERR(03))
      allocate(        Z(MP), STAT=LERR(04))
      allocate(        U(MP), STAT=LERR(05))
      allocate(        V(MP), STAT=LERR(06))
      allocate(        W(MP), STAT=LERR(07))
      allocate(        P(ME), STAT=LERR(08))
      allocate(       PN(MP), STAT=LERR(09))
      allocate( UA(MP,numPhase), STAT=LERR(10))
      allocate( VA(MP,numPhase), STAT=LERR(11))
      allocate( WA(MP,numPhase), STAT=LERR(12))
      allocate( PA(ME,numPhase), STAT=LERR(13))
      allocate(PNA(MP,numPhase), STAT=LERR(14))
      call CHKALC(14,LERR,IUT6,IERR) 
      if(IERR.NE.0) then
          write(IUT6,*) 'phaseAverage: allocating error       '
          STOP
      endif
      write(IUT6,*) 'phaseAverage: allocating finish       '

      write(IUT6,*)
      write(IUT6,*) 'reading MESH data...'
!
! read MESH file
!
      IACT = 1
      call GFALL(IUT0,IUT6,IUTMS,FILEMS,
     *           MCOM,NCOMFL,COMFLE,
     *           MCOM,NCOMST,COMSET,
     *           IACT,IWRITE,INAME,IRESV,  
     *           ICAST,IDATA0,IALL,ISKIP,IERR,
     *           '*GRID_3D *NODE_3D !',
     *           NAME,MP,NP,X,Y,Z,
     *           NAME,ME,N,NE,NDUM,NODE,
     *           ICHECK)     
      if(IERR.NE.0) STOP
!
! read FLOWS file
!
      IF(FLAGGFSEP.EQ.0)THEN
        IACT = 3
        ! Open FLOWS
        call GFALL(IUT0,IUT6,IUTFF,FILEFF,
     *             MCOM,NCOMFL,COMFLE,
     *             MCOM,NCOMST,COMSET,
     *             IACT,IWRITE,INAME,IRESV,  
     *             ICAST,IDATA0,IALL,ISKIP,IERR,
     *             ' !',
     *             ICHECK)     
        if(IERR.NE.0) STOP

        UA = 0.0E0
        VA = 0.0E0
        WA = 0.0E0
        PA = 0.0E0
        PNA = 0.0E0
        NSUM = 0
        NFLOWS = 0
        do
          IACT = 5
          ! Read FLOWS
          call GFALL(IUT0,IUT6,IUTFF,FILEFF,
     *               MCOM,NCOMFL,COMFLE,
     *               MCOM,NCOMST,COMSET,
     *               IACT,IWRITE,INAME,IRESV,  
     *               ICAST,IDATA0,IALL,ISKIP,IERR,
     *               '*TIME_PS *STEP_PS *VELO_3D 
     *                *PRES_3E *PRES_3D !',
     *               NAME,TIMEP,
     *               NAME,ISTEP,
     *               NAME,MP,NPCHK,U,V,W,
     *               NAME,ME,NEPRS,P,
     *               NAME,MP,NPPRS,PN,
     *               ICHECK)
          if(IERR.NE.0)   STOP
          ! NFLOWSはFLOWSの数
          NFLOWS = NFLOWS + 1
          ! ISUM = 0 - numPhase-1
          ! 1,2,3,..,0(numPhase),1,2,3..,0(numPhase)
          ISUM = mod(NFLOWS,numPhase)
          ! NSUMは平均した回数
          IF(ISUM == 0) NSUM = NSUM + 1
          ! ISUM+1 = 1 - numPhase
          do IP=1,NP
            UA (IP,ISUM+1) = UA (IP,ISUM+1) + U (IP)
            VA (IP,ISUM+1) = VA (IP,ISUM+1) + V (IP)
            WA (IP,ISUM+1) = WA (IP,ISUM+1) + W (IP)
            PNA(IP,ISUM+1) = PNA(IP,ISUM+1) + PN(IP)
          enddo
          do IE=1,NE
            PA (IE,ISUM+1) = PA (IE,ISUM+1) + P(IE)
          enddo
          write(IUT6,*) "averaged:when NFLOWS=",NFLOWS,
     *                                 ",NSUM=",NSUM,
     *                                 ",ISUM=",ISUM
          if(IACT .EQ. 7) exit
        end do

        IACT = 7
        ! Close FLOWS
        call GFALL(IUT0,IUT6,IUTFF,FILEFF,
     *             MCOM,NCOMFL,COMFLE,
     *             MCOM,NCOMST,COMSET,
     *             IACT,IWRITE,INAME,IRESV,  
     *             ICAST,IDATA0,IALL,ISKIP,IERR,
     *             ' !', 
     *             ICHECK)     
        if(IERR.NE.0) STOP
      ELSEIF(FLAGGFSEP.EQ.1)THEN
        ! SEPARATED
        write(IUT6,*) "FLAGGFSEP",FLAGGFSEP
        write(IUT6,*) "NFLOWS",NFLOWS
        UA = 0.0E0
        VA = 0.0E0
        WA = 0.0E0
        PA = 0.0E0
        PNA = 0.0E0
        NSUM = 0
        do IFLOWS=1,NFLOWS
          write(CNUM,'(I4.4)') IFLOWS
          IF(FLAGDELIM.eq.0)THEN
            FILEFS=trim(FILEFF)//".P"//CNUM
          ELSEIF(FLAGDELIM.eq.1)THEN
            FILEFS=trim(FILEFF)//".t"//CNUM
          ENDIF
          ! Read FLOWS
          IACT = 1
          call GFALL(IUT0,IUT6,IUTFF,FILEFS,
     *               MCOM,NCOMFL,COMFLE,
     *               MCOM,NCOMST,COMSET,
     *               IACT,IWRITE,INAME,IRESV,  
     *               ICAST,IDATA0,IALL,ISKIP,IERR,
     *               '*TIME_PS *STEP_PS *VELO_3D 
     *                *PRES_3E *PRES_3D !',
     *               NAME,TIMEP,
     *               NAME,ISTEP,
     *               NAME,MP,NPCHK,U,V,W,
     *               NAME,ME,NEPRS,P,
     *               NAME,MP,NPPRS,PN,
     *               ICHECK)
          if(IERR.NE.0)   STOP
          ! ISUM = 0 - numPhase-1
          ! 1,2,3,..,0(numPhase),1,2,3..,0(numPhase)
          ISUM = mod(IFLOWS,numPhase)
          ! NSUMは平均した回数
          IF(ISUM == 0) NSUM = NSUM + 1
          ! ISUM+1 = 1 - numPhase
          do IP=1,NP
            UA (IP,ISUM+1) = UA (IP,ISUM+1) + U(IP)
            VA (IP,ISUM+1) = VA (IP,ISUM+1) + V(IP)
            WA (IP,ISUM+1) = WA (IP,ISUM+1) + W(IP)
            PNA(IP,ISUM+1) = PNA(IP,ISUM+1) + PN(IP)
          end do
          do IE=1,NE
            PA (IE,ISUM+1) = PA (IE,ISUM+1) + P(IE)
          end do
          write(IUT6,*) "averaged:when NFLOWS=",IFLOWS,",NSUM=",NSUM,
     *                                 ",ISUM=",ISUM
        end do
      ELSE
         write(IUT6,*) "FLAGGFSEP",FLAGGFSEP
         STOP
      ENDIF
!
! conduct average
!
      ! ゼロ割りを防ぐ
      if(NSUM < 1)then
          WRITE(IUT6,*) "ERROR IN FLAG1"
          STOP
      endif
      INVNUM = 1.0E0 / float(NSUM)
      do iPhase=1,numPhase
        do IP=1,NP
          UA (IP,iPhase) = UA (IP,iPhase) * INVNUM
          VA (IP,iPhase) = VA (IP,iPhase) * INVNUM
          WA (IP,iPhase) = WA (IP,iPhase) * INVNUM
          PNA(IP,iPhase) = PNA(IP,iPhase) * INVNUM
        end do
        do IE=1,NE
          PA (IE,iPhase) = PA (IE,iPhase) * INVNUM
        end do
      end do

!
! output AVE file
!
      do iPhase=1,numPhase
          TIMEP = 0.0E0
          ISTEP = NSUM
          do IP=1,NP
            U  (IP)= UA (IP,iPhase)
            V  (IP)= VA (IP,iPhase)
            W  (IP)= WA (IP,iPhase)
            PN (IP)= PNA(IP,iPhase)
          end do
          do IE=1,NE
            P(IE)  = PA (IE,iPhase)
          end do
 
          IACT = 2
          write(CNUM,'(I4.4)') iPhase
          FILEAVE = trim(FILEAV)//'.P'//CNUM
          call GFALL(IUT0,IUT6,IUTAV,FILEAVE,
     *               MCOM,NCOMFL,COMFLE,
     *               MCOM,NCOMST,COMSET,
     *               IACT,IWRITE,INAME,IRESV,  
     *               ICAST,IDATA0,IALL,ISKIP,IERR, 
     *               '*TIME_PS *STEP_PS *VELO_3D 
     *                *PRES_3E *PRES_3D !',
     *               NAME,TIMEP,
     *               NAME,ISTEP,
     *               NAME,MP,NP,U,V,W,
     *               NAME,ME,NE,P,
     *               NAME,MP,NP,PN,
     *               ICHECK)
      end do

      STOP
      end program phaseAverage
