module converr
 
!$$$   module documentation block
!                .      .    .                                       .
! module:    convinfo
!   prgmmr: su          org: np2                date: 2007-03-15
! abstract:  This module contains variables and routines related
!            to the assimilation of conventional observations error
!
! program history log:
!   2007-03-15  su  - original code - move reading observation error table 
!                                     from read_prepbufr to here so all the 
!                                     processor can have the new error information 
!
! Subroutines Included:
!   sub converr_read      - allocate arrays for and read in conventional error table 
!   sub converr_destroy   - destroy conventional error arrays
!
! Variable Definitions:
!   def etabl             -  the array to hold the error table
!   def ptabl             -  the array to have vertical pressure values
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block

use kinds, only:r_kind,i_kind,r_single
use constants, only: zero
use obsmod, only : oberrflg 

  integer(i_kind) ietabl,itypex,lcount,k,m
  real(r_single),allocatable,dimension(:,:,:) :: etabl
  real(r_kind),allocatable,dimension(:)  :: ptabl

contains


  subroutine converr_read(mype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    convinfo_err      read conventional information file
!
!     prgmmr:    su    org: np2                date: 2007-03-15
!
! abstract:  This routine reads the conventional error table file
!
! program history log:
!   2008-06-04  safford -- add subprogram doc block
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

     allocate(etabl(300,33,6))

     etabl=1.e9_r_kind
      
     
     ietabl=19
     open(ietabl,file='errtable',form='formatted')
     rewind ietabl
     etabl=1.e9_r_kind
     lcount=0
     loopd : do 
        read(ietabl,100,IOSTAT=iflag) itypex
        if( iflag /= 0 ) exit loopd
100     format(1x,i3)
        lcount=lcount+1
        do k=1,33
           read(ietabl,110)(etabl(itypex,k,m),m=1,6)
110        format(1x,6e12.5)
        end do
     end do   loopd

    if(lcount.le.0 .and. mype ==0) then
       write(6,*)'CONVERR:  ***WARNING*** obs error table not available to 3dvar.'
       oberrflg=.false.
    else
       if(mype ==0) write(6,*)'CONVERR:  using observation errors from user provided table'
       allocate(ptabl(34))
       ptabl=zero
       ptabl(1)=etabl(120,1,1)
       do k=2,33
          ptabl(k)=0.5*(etabl(120,k-1,1)+etabl(120,k,1))
       enddo
       ptabl(34)=etabl(120,33,1)
    endif

     close(ietabl)

     return
  end subroutine converr_read


subroutine converr_destroy
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    converr_destroy      destroy conventional information file
!     prgmmr:    su    org: np2                date: 2007-03-15
!
! abstract:  This routine destroys arrays from converr file
!
! program history log:
!   2007-03-15  su 
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$

     deallocate(etabl,ptabl)
     return
  end subroutine converr_destroy

end module converr



