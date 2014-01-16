subroutine state2ensctl(eval,mval,grad)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    state2ensctl
!   prgmmr: kleist
!
! abstract:  Contribution from state space to ensemble control vector
!
! program history log:
!   2011-11-17  kleist - initial code
!
!   input argument list:
!     eval - Ensemble state variable variable
!     grad - Control variable
!
!   output argument list:
!     grad - Control variable
!
!$$$ end documentation block

use kinds, only: r_kind,i_kind
use control_vectors, only: control_vector,cvars3d
use gsi_4dvar, only: l4dvar,l4densvar,nobs_bins,ibin_anl
use hybrid_ensemble_parameters, only: uv_hyb_ens,dual_res,ntlevs_ens
use hybrid_ensemble_isotropic, only: ensemble_forward_model_ad
use hybrid_ensemble_isotropic, only: ensemble_forward_model_ad_dual_res
use balmod, only: strong_bk_ad
use gsi_bundlemod, only: gsi_bundlecreate
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: gsi_bundlegetvar
use gsi_bundlemod, only: gsi_bundleputvar
use gsi_bundlemod, only: gsi_bundledestroy
use gsi_bundlemod, only: assignment(=)
use gsi_bundlemod, only : self_add
use constants, only: zero,max_varname_length
use mpeu_util, only: getindex
use gsi_metguess_mod, only: gsi_metguess_get
use mod_strong, only: tlnmc_option
implicit none

! Declare passed variables
type(control_vector), intent(inout) :: grad
type(gsi_bundle)    , intent(inout) :: mval
type(gsi_bundle)    , intent(in   ) :: eval(ntlevs_ens)

! Declare local variables
character(len=*),parameter::myname='state2ensctl'
character(len=max_varname_length),allocatable,dimension(:) :: clouds
integer(i_kind) :: i,j,k,ii,jj,ic,id,istatus,nclouds

integer(i_kind), parameter :: ncvars = 5
integer(i_kind) :: icps(ncvars)
type(gsi_bundle):: wbundle_c ! work bundle
character(len=3), parameter :: mycvars(ncvars) = (/  &  ! vars from CV needed here
                               'sf ', 'vp ', 'ps ', 't  ',    &
                               'q  '/)
logical :: lc_sf,lc_vp,lc_ps,lc_t,lc_rh
real(r_kind),pointer,dimension(:,:)   :: cv_ps
real(r_kind),pointer,dimension(:,:,:) :: cv_sf,cv_vp,cv_rh,cv_tv
! Declare required local state variables
integer(i_kind), parameter :: nsvars = 5
integer(i_kind) :: isps(nsvars)
character(len=4), parameter :: mysvars(nsvars) = (/  &  ! vars from ST needed here
                               'u   ', 'v   ', 'p3d ', 'q   ', 'tsen' /)
logical :: ls_u,ls_v,ls_p3d,ls_q,ls_tsen
real(r_kind),pointer,dimension(:,:)   :: rv_ps,rv_sst
real(r_kind),pointer,dimension(:,:,:) :: rv_u,rv_v,rv_p3d,rv_q,rv_tsen,rv_tv,rv_oz
real(r_kind),pointer,dimension(:,:,:) :: rv_rank3

logical :: do_getuv,do_tv_to_tsen_ad,do_normal_rh_to_q_ad,do_getprs_ad,lstrong_bk_vars

!****************************************************************************

! Inquire about chemistry
call gsi_metguess_get('clouds::3d',nclouds,istatus)
if (nclouds>0) then
    allocate(clouds(nclouds))
    call gsi_metguess_get('clouds::3d',clouds,istatus)
endif

! Since each internal vector of grad has the same structure, pointers are
! the same independent of the subwindow jj
call gsi_bundlegetpointer (grad%step(1),mycvars,icps,istatus)
lc_sf =icps(1)>0; lc_vp =icps(2)>0; lc_ps =icps(3)>0
lc_t  =icps(4)>0; lc_rh =icps(5)>0

! Since each internal vector of grad has the same structure, pointers are
! the same independent of the subwindow jj
call gsi_bundlegetpointer (eval(1),mysvars,isps,istatus)
ls_u  =isps(1)>0; ls_v   =isps(2)>0; ls_p3d=isps(3)>0
ls_q  =isps(4)>0; ls_tsen=isps(5)>0

! Define what to do depending on what's in CV and SV
lstrong_bk_vars     =lc_sf.and.lc_vp.and.lc_ps .and.lc_t
do_getuv            =lc_sf.and.lc_vp.and.ls_u  .and.ls_v
do_tv_to_tsen_ad    =lc_t .and.ls_q .and.ls_tsen
do_normal_rh_to_q_ad=lc_t .and.lc_rh.and.ls_p3d.and.ls_q
do_getprs_ad        =lc_t .and.lc_ps.and.ls_p3d

! Initialize
mval%values=zero

do jj=1,ntlevs_ens
!  Create a temporary bundle similar to grad, and copy contents of grad into it
   call gsi_bundlecreate ( wbundle_c, grad%step(1), 'stat2ensctl work', istatus )
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble creating work bundle'
      call stop2(999)
   endif
   wbundle_c%values=zero

   call gsi_bundlegetpointer (wbundle_c,'sf' ,cv_sf ,istatus)
   call gsi_bundlegetpointer (wbundle_c,'vp' ,cv_vp ,istatus)
   call gsi_bundlegetpointer (wbundle_c,'q'  ,cv_rh ,istatus)
   call gsi_bundlegetpointer (wbundle_c,'t'  ,cv_tv, istatus)
   call gsi_bundlegetpointer (wbundle_c,'ps' ,cv_ps ,istatus)

! Get sv pointers here
!  Get pointers to required state variables
   call gsi_bundlegetpointer (eval(jj),'u'   ,rv_u,   istatus)
   call gsi_bundlegetpointer (eval(jj),'v'   ,rv_v,   istatus)
   call gsi_bundlegetpointer (eval(jj),'ps'  ,rv_ps,  istatus)
   call gsi_bundlegetpointer (eval(jj),'p3d' ,rv_p3d, istatus)
   call gsi_bundlegetpointer (eval(jj),'tv'  ,rv_tv,  istatus)
   call gsi_bundlegetpointer (eval(jj),'tsen',rv_tsen,istatus)
   call gsi_bundlegetpointer (eval(jj),'q'   ,rv_q ,  istatus)
   call gsi_bundlegetpointer (eval(jj),'oz'  ,rv_oz , istatus)
   call gsi_bundlegetpointer (eval(jj),'sst' ,rv_sst, istatus)


!  Calculate sensible temperature
   if(do_tv_to_tsen_ad) call tv_to_tsen_ad(rv_tv,rv_q,rv_tsen)

!  Adjoint to convert ps to 3-d pressure
   if(do_getprs_ad) call getprs_ad(rv_ps,rv_tv,rv_p3d)

! If calling TLNMC, already have u,v (so set last argument to true)
   if(lstrong_bk_vars) then
      if ( (tlnmc_option==3) .or. &
         (jj==ibin_anl .and. tlnmc_option==2) ) then

!  Adjoint of consistency for 3d pressure and sensible temperature
!  Calculate sensible temperature
         if(do_tv_to_tsen_ad) call tv_to_tsen_ad(rv_tv,rv_q,rv_tsen)

!  Adjoint to convert ps to 3-d pressure
         if(do_getprs_ad) call getprs_ad(rv_ps,rv_tv,rv_p3d)

!  Adjoint of strong_bk
         call strong_bk_ad(rv_u,rv_v,rv_ps,rv_tv,.true.)
      end if
   end if

   call self_add(mval,eval(jj))

!  Adjoint of control to initial state
   call gsi_bundleputvar ( wbundle_c, 't' ,  rv_tv,  istatus )
   call gsi_bundleputvar ( wbundle_c, 'ps',  rv_ps,  istatus )
   call gsi_bundleputvar ( wbundle_c, 'oz',  rv_oz,  istatus )
   call gsi_bundleputvar ( wbundle_c, 'sst', rv_sst, istatus )

!  Since cloud-vars map one-to-one, take care of them together
   do ic=1,nclouds
      id=getindex(cvars3d,clouds(ic))
      if (id>0) then
          call gsi_bundlegetpointer (eval(jj),clouds(ic),rv_rank3,istatus)
          call gsi_bundleputvar     (wbundle_c, clouds(ic),rv_rank3,istatus)
      endif
   enddo

!  Convert RHS calculations for u,v to st/vp
   if (do_getuv) then
      if(uv_hyb_ens) then
         call gsi_bundleputvar ( wbundle_c, 'sf', rv_u, istatus )
         call gsi_bundleputvar ( wbundle_c, 'vp', rv_v, istatus )
      else
         call getuv(rv_u,rv_v,cv_sf,cv_vp,1)
      end if
   end if

!  Calculate sensible temperature
   if(do_tv_to_tsen_ad) call tv_to_tsen_ad(cv_tv,rv_q,rv_tsen)

!  Adjoint of convert input normalized RH to q to add contribution of moisture
!  to t, p , and normalized rh
   if(do_normal_rh_to_q_ad) call normal_rh_to_q_ad(cv_rh,cv_tv,rv_p3d,rv_q)

!  Adjoint to convert ps to 3-d pressure
   if(do_getprs_ad) call getprs_ad(cv_ps,cv_tv,rv_p3d)

   if(dual_res) then
      call ensemble_forward_model_ad_dual_res(wbundle_c,grad%aens(1,:),jj)
   else
      call ensemble_forward_model_ad(wbundle_c,grad%aens(1,:),jj)
   end if

   call gsi_bundledestroy(wbundle_c,istatus)
   if (istatus/=0) then
      write(6,*) trim(myname),': trouble destroying work bundle'
      call stop2(999)
   endif

end do

if (nclouds>0) deallocate(clouds)


return 
end subroutine state2ensctl