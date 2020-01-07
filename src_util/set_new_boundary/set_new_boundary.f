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
      !  MAIN PRORGRAM  set_new_boundary                         !
      !                                                          !
      !                           2020/01/07                     !
      !                           WRITTEN BY YUYA MIKI           !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !<
      !! @param ME :: maximum element number
      !! @param MP :: maximum node number
      !! @param N  :: the number of node in a element
      !! @param NODE :: node number list
      !! @param X,Y,Z :: Position
      !! @param LPBOUN :: Boundary List
      !! @version 0.1
      program set_new_boundary
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
      real(4),allocatable,dimension(:)   :: X,Y,Z
!
      integer(4)    :: ISTEP,NPCHK,NEPRS,NPPRS
      character(60) :: FILEMS,FILEBN
!     [work]
      integer(4)    :: NP,NE,NDUM,IP,IE,IN
      ! fortran input output
      integer(4),parameter :: IUT0  = 0
      integer(4),parameter :: IUT5  = 5
      integer(4),parameter :: IUT6  = 6
      integer(4),parameter :: IUTMS = 11
      integer(4),parameter :: IUTBN = 12
      integer(4) :: IERR,FLAGOK,FLAGGFSEP,FLAGDELIM

      ! [boundary]
      integer(4) :: MB,NPBOUN(18)
      integer(4),allocatable :: LPBOUN(:,:)
      integer(4),allocatable,dimension(:) :: LPINT1,LPINT2,LPINT3,
     *                                       LPSET1,LPSET2,LPSET3,
     *                                       LPDEP1,LPDEP2,
     *                                       LEHSRC

      real(4),allocatable,dimension(:) :: TEMP, HEAT, HSRC, HTRS

      ! detect
      real(4) :: z_detect_value

!
! read MESH DATA
!
      write(IUT6,*)
      write(IUT6,*) '** set new boundary **'
!
      do
          write(IUT6,*)
          write(IUT6,*) '[In]Specify filename for MESH data'
          read (IUT5,'(A60)') FILEMS
          write(IUT6,*) '[In]Specify filename for BOUN data'
          read (IUT5,'(A60)') FILEBN

          write(IUT6,'(A10,A60)') "FILEMS   =", FILEMS
          write(IUT6,'(A10,A60)') "FILEBN   =", FILEBN
          write(IUT6,*) 'These parameters are ok? Yes:1,No:0'
          read(IUT5,*) FLAGOK
          if(FLAGOK.EQ.1) exit
      enddo
!
! check data size using MESH file
!
      call CHKMSH(FILEMS,IUTMS,IUT0,IUT6,MP,ME,IERR)
      if (IERR.NE.0) STOP

      ! set maximum number of boundary nodes as MP
      MB = MP;
!
! allocate variables
!
      write(IUT6,*)
      write(IUT6,*)
      write(IUT6,*) 'set_new_boundary: ALLOCATING VARIABLES... '
      allocate(   NODE(N,ME), STAT=LERR(01))
      allocate(        X(MP), STAT=LERR(02))
      allocate(        Y(MP), STAT=LERR(03))
      allocate(        Z(MP), STAT=LERR(04))
      ! Boundary-variables allocate
      allocate(LPBOUN(MB,18), STAT=LERR(05))
      allocate(LPINT1(MB),    STAT=LERR(06))
      allocate(LPINT2(MB),    STAT=LERR(07))
      allocate(LPINT3(MB),    STAT=LERR(08))
      allocate(LPSET1(MB),    STAT=LERR(09))
      allocate(LPSET2(MB),    STAT=LERR(10))
      allocate(LPSET3(MB),    STAT=LERR(11))
      allocate(LPDEP1(MB),    STAT=LERR(12))
      allocate(LPDEP2(MB),    STAT=LERR(13))
      !
      ALLOCATE(    UWALL(MB), STAT=LERR(14))
      ALLOCATE(    VWALL(MB), STAT=LERR(15))
      ALLOCATE(    WWALL(MB), STAT=LERR(16))
      !
      ALLOCATE(    UINLT(MB), STAT=LERR(17))
      ALLOCATE(    VINLT(MB), STAT=LERR(18))
      ALLOCATE(    WINLT(MB), STAT=LERR(19))
      !
      ALLOCATE(     TEMP(MB), STAT=LERR(20))
      ALLOCATE(     HEAT(MB), STAT=LERR(21))
      !
      ALLOCATE(   LEHSRC(MB), STAT=LERR(22))
      ALLOCATE(     HSRC(MB), STAT=LERR(23))
      ALLOCATE(     HTRS(MB), STAT=LERR(24))

      call CHKALC(24,LERR,IUT6,IERR) 
      if(IERR.NE.0) then
          write(IUT6,*) 'set_new_boundary: allocating error       '
          STOP
      endif
      write(IUT6,*) 'set_new_boundary: allocating finish       '

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
      ! read BOUN file
      !
      IACT = 1
      call GFALL(IUT0,IUT6,IUTBN,FILEBN,
     *           MCOM,NCOMFL,COMFLE,
     *           MCOM,NCOMST,COMSET,
     *           IACT,IWRITE,INAME,IRESV,  
     *           ICAST,IDATA0,IALL,ISKIP,IERR,
     *           '*BC_MWAL *BC_WALL *BC_SYMT *BC_INLT *BC_FREE 
     *            *BC_BODY *BC_CYCL *BC_INTR *BC_PSET
     *            *BC_RELY *BC_WV3D *BC_IV3D 
     *            *BC_TMPN *BC_TMPV *BC_HSRN *BC_HSRV  
     *            *BC_HFXN *BC_HFXV *BC_HTRN *BC_HTRV 
     *            *BC_FORC *BC_FOIN !',
     ! movin wall
     *           NAME,MB,NPBOUN( 1),LPBOUN(1, 1), 
     ! wall
     *           NAME,MB,NPBOUN( 2),LPBOUN(1, 2), 
     ! symmetry
     *           NAME,MB,NPBOUN( 3),LPBOUN(1, 3), 
     ! inlet
     *           NAME,MB,NPBOUN( 4),LPBOUN(1, 4), 
     ! free
     *           NAME,MB,NPBOUN( 5),LPBOUN(1, 5), 
     ! body
     *           NAME,MB,NPBOUN( 6),LPBOUN(1, 6), 
     *           NAME,MB,NPBOUN( 7),LPBOUN(1, 7),LPBOUN(1, 8), 
     *           NAME,MB,NPINT     ,LPINT1,LPINT2,LPINT3,
     *           NAME,MB,NPSET     ,LPSET1,LPSET2,LPSET3,
     *           NAME,MB,NPDEP     ,LPDEP1,LPDEP2,
     *           NAME,MB,NPBOUN( 1),UWALL,VWALL,WWALL,
     *           NAME,MB,NPBOUN( 4),UINLT,VINLT,WINLT,
     *           NAME,MB,NPBOUN( 9),LPBOUN(1, 9), 
     *           NAME,MB,NPBOUN( 9),TEMP,
     *           NAME,MB,NPBOUN(10),LPBOUN(1,10),
     *           NAME,MB,NPBOUN(10),HEAT,
     *           NAME,MB,NEHSRC    ,LEHSRC,
     *           NAME,MB,NEHSRC    ,HSRC,
     *           NAME,MB,NPBOUN(11),LPBOUN(1,11),
     *           NAME,MB,NPBOUN(11),HTRS,
     *           NAME,MB,NPBOUN(15),LPBOUN(1,15),LPBOUN(1,16),
     *           NAME,MB,NPBOUN(17),LPBOUN(1,17),LPBOUN(1,18),
     *           ICHECK)     
      if(IERR.NE.0) STOP

      ! set some plane to inlet boundary
      write(iut6,*) "Inlet node number is ", NPBOUN(4)

      ! loop all elements to detect specified surface(maybe z plane.)
      
      z_detect_value = 0.4 ! please set value
      do ie = 1, ne
            do in = 1,8
                  ip = node(in,ie)
                  if(ip .eq. 0) exit
                  if(z(ip) .eq. z_detect_value)then
                        ! add this node to inlet boundary
                        NPBOUN(4) = NPBOUN(4) + 1
                        inlet_back = NPBOUN(4)
                        LPBOUN(inlet_back, 4) = ip
                  endif

            end do
      end do
      write(iut6,*) "Inlet node number is ", NPBOUN(4)
      write(iut6,*) "set_new_boundary done"

      STOP
      end program set_new_boundary