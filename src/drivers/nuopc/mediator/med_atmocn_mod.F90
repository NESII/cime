module med_atmocn_mod

  use shr_kind_mod          , only: r8=>shr_kind_r8, in=>shr_kind_in
  use shr_kind_mod          , only: cs=>shr_kind_cs, cl=>shr_kind_cl
  use shr_sys_mod           , only: shr_sys_abort
  use shr_flux_mod          , only: shr_flux_atmocn, shr_flux_atmocn_diurnal
  use shr_orb_mod           , only: shr_orb_params, shr_orb_cosz, shr_orb_decl
  use shr_const_mod         , only: SHR_CONST_PI
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_getFieldN
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only: shr_nuopc_methods_ChkErr
  use seq_timemgr_mod       , only: seq_timemgr_EclockGetData
  use med_constants_mod     , only: med_constants_dbug_flag
  use ESMF
  use NUOPC

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_atmocn_init
  public med_atmocn_ocnalb
  public med_atmocn_flux

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  real(r8) , pointer :: lats        (:) ! latitudes  (degrees)
  real(r8) , pointer :: lons        (:) ! longitudes (degrees)
  integer  , pointer :: mask        (:) ! ocn domain mask: 0 <=> inactive cell
  integer  , pointer :: emask       (:) ! ocn mask on exchange grid decomp
  real(r8) , pointer :: anidr       (:) ! albedo: near infrared, direct
  real(r8) , pointer :: avsdr       (:) ! albedo: visible      , direct
  real(r8) , pointer :: anidf       (:) ! albedo: near infrared, diffuse
  real(r8) , pointer :: avsdf       (:) ! albedo: visible      , diffuse
  real(r8) , pointer :: swndr       (:) ! direct near-infrared  incident solar radiation
  real(r8) , pointer :: swndf       (:) ! diffuse near-infrared incident solar radiation
  real(r8) , pointer :: swvdr       (:) ! direct visible  incident solar radiation
  real(r8) , pointer :: swvdf       (:) ! diffuse visible incident solar radiation
  real(r8) , pointer :: uocn        (:) ! ocn velocity, zonal
  real(r8) , pointer :: vocn        (:) ! ocn velocity, meridional
  real(r8) , pointer :: tocn        (:) ! ocean temperature
  real(r8) , pointer :: zbot        (:) ! atm level height
  real(r8) , pointer :: ubot        (:) ! atm velocity, zonal
  real(r8) , pointer :: vbot        (:) ! atm velocity, meridional
  real(r8) , pointer :: thbot       (:) ! atm potential T
  real(r8) , pointer :: shum        (:) ! atm specific humidity
  real(r8) , pointer :: shum_16O    (:) ! atm H2O tracer
  real(r8) , pointer :: shum_HDO    (:) ! atm HDO tracer
  real(r8) , pointer :: shum_18O    (:) ! atm H218O tracer
  real(r8) , pointer :: roce_16O    (:) ! ocn H2O ratio
  real(r8) , pointer :: roce_HDO    (:) ! ocn HDO ratio
  real(r8) , pointer :: roce_18O    (:) ! ocn H218O ratio
  real(r8) , pointer :: dens        (:) ! atm density
  real(r8) , pointer :: tbot        (:) ! atm bottom surface T
  real(r8) , pointer :: sen         (:) ! heat flux: sensible
  real(r8) , pointer :: lat         (:) ! heat flux: latent
  real(r8) , pointer :: lwup        (:) ! lwup over ocean
  real(r8) , pointer :: evap        (:) ! water flux: evaporation
  real(r8) , pointer :: evap_16O    (:) ! H2O flux: evaporation
  real(r8) , pointer :: evap_HDO    (:) ! HDO flux: evaporation
  real(r8) , pointer :: evap_18O    (:) ! H218O flux: evaporation
  real(r8) , pointer :: taux        (:) ! wind stress, zonal
  real(r8) , pointer :: tauy        (:) ! wind stress, meridional
  real(r8) , pointer :: tref        (:) ! diagnostic:  2m ref T
  real(r8) , pointer :: qref        (:) ! diagnostic:  2m ref Q
  real(r8) , pointer :: u10         (:) ! diagnostic: 10m wind speed
  real(r8) , pointer :: duu10n      (:) ! diagnostic: 10m wind speed squared
  real(r8) , pointer :: fswpen      (:) ! fraction of sw penetrating ocn surface layer
  real(r8) , pointer :: ocnsal      (:) ! ocean salinity
  real(r8) , pointer :: lwdn        (:) ! long  wave, downward
  real(r8) , pointer :: swdn        (:) ! short wave, downward
  real(r8) , pointer :: swup        (:) ! short wave, upward
  real(r8) , pointer :: rainc       (:) ! rainc
  real(r8) , pointer :: rainl       (:) ! rainl
  real(r8) , pointer :: snowc       (:) ! snowc
  real(r8) , pointer :: snowl       (:) ! snowl
  real(r8) , pointer :: tbulk       (:) ! diurnal diagnostic: ocn bulk T
  real(r8) , pointer :: tskin       (:) ! diurnal diagnostic: ocn skin T
  real(r8) , pointer :: tskin_night (:) ! diurnal diagnostic: ocn skin T
  real(r8) , pointer :: tskin_day   (:) ! diurnal diagnostic: ocn skin T
  real(r8) , pointer :: cSkin       (:) ! diurnal diagnostic: ocn cool skin
  real(r8) , pointer :: cSkin_night (:) ! diurnal diagnostic: ocn cool skin
  real(r8) , pointer :: warm        (:) ! diurnal diagnostic: ocn warming
  real(r8) , pointer :: warmMax     (:) ! diurnal diagnostic: ocn warming, max daily value
  real(r8) , pointer :: warmMaxInc  (:) ! diurnal diagnostic: ocn warming, max daily value, increment
  real(r8) , pointer :: salt        (:) ! diurnal diagnostic: ocn salting
  real(r8) , pointer :: speed       (:) ! diurnal diagnostic: ocn speed
  real(r8) , pointer :: regime      (:) ! diurnal diagnostic: ocn regime
  real(r8) , pointer :: windMax     (:) ! diurnal diagnostic: ocn wind   , max daily value
  real(r8) , pointer :: windAvg     (:) ! diurnal diagnostic: ocn wind   , daily avg
  real(r8) , pointer :: windMaxInc  (:) ! diurnal diagnostic: ocn wind   , max daily value, increment
  real(r8) , pointer :: QsolAvg     (:) ! diurnal diagnostic: ocn Qsol   , daily avg
  real(r8) , pointer :: qSolInc     (:) ! diurnal diagnostic: ocn Qsol   , daily avg, increment
  real(r8) , pointer :: windInc     (:) ! diurnal diagnostic: ocn wind   , daily avg, increment
  real(r8) , pointer :: nInc        (:) ! diurnal diagnostic: a/o flux   , increment
  real(r8) , pointer :: ustar       (:) ! saved ustar
  real(r8) , pointer :: re          (:) ! saved re
  real(r8) , pointer :: ssq         (:) ! saved sq

  ! Fields that are not obtained via GetFldPtr
  real(r8) , pointer :: uGust       (:) ! wind gust
  real(r8) , pointer :: prec        (:) ! precip
  real(r8) , pointer :: prec_gust   (:) ! atm precip for convective gustiness (kg/m^3)

  ! TODO: put flds_wiso in config and query it - for now just hard-wire it to false
  logical :: flds_wiso = .false.

  ! Conversion from degrees to radians
  integer                :: dbug_flag = med_constants_dbug_flag
  integer                :: dbrc
  character(len=1024)    :: tmpstr
  real(r8)    ,parameter :: const_pi      = SHR_CONST_PI       ! pi
  real(r8)    ,parameter :: const_deg2rad = const_pi/180.0_r8  ! deg to rads
  character(*),parameter :: u_FILE_u = &
    __FILE__
!===============================================================================
contains
!===============================================================================

  subroutine med_atmocn_init(gcomp, aoflux_grid, FBAtmOcn, FBAtm, FBOcn, FBFrac, rc)
    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables and then use the module
    ! variables in the med_atmocn_ocnalb and med_atmocn_flux routine
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: aoflux_grid ! FBAtmOcn grid - 'atm' or 'ocn'
    type(ESMF_FieldBundle)         :: FBAtmOcn    ! Atm/Ocn flux fields
    type(ESMF_FieldBundle)         :: FBAtm       ! Atm Import fields on atmocn flux grid
    type(ESMF_FieldBundle)         :: FBOcn       ! Ocn Import fields on atmocn flux grid
    type(ESMF_FieldBundle)         :: FBfrac      ! Fraction data for various components, on their grid
    integer          , intent(out) :: rc
    !
    ! Local variables
    type(ESMF_VM)            :: vm
    integer(in)              :: iam
    type(ESMF_Field)         :: lfield
    type(ESMF_Grid)          :: lgrid
    type(ESMF_Mesh)          :: lmesh
    type(ESMF_GeomType_Flag) :: geomtype
    integer                  :: n
    integer                  :: lsize
    real(r8), pointer        :: rmask(:)  ! ocn domain mask
    real(r8), pointer        :: ofrac(:)
    real(r8), pointer        :: ifrac(:)
    integer                  :: dimCount
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    real(ESMF_KIND_R8), pointer :: ownedElemCoords(:)
    character(*),parameter   :: subName =   '(med_atmocn_init) '
    !-----------------------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! The following is for debugging
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! atm/ocn fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_avsdr', fldptr1=avsdr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_avsdf', fldptr1=avsdf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_anidr', fldptr1=anidr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_anidf', fldptr1=anidf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_tref', fldptr1=tref, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_qref', fldptr1=qref, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_ustar', fldptr1=ustar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_re', fldptr1=re, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_ssq', fldptr1=ssq, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_u10', fldptr1=u10, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_duu10n', fldptr1=duu10n, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_taux', fldptr1=taux, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_tauy', fldptr1=tauy, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_lat', fldptr1=lat, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_sen', fldptr1=sen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_evap', fldptr1=evap, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(evap)
    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_evap_16O', fldptr1=evap_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_evap_18O', fldptr1=evap_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_evap_HDO', fldptr1=evap_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(evap_16O(lsize)); evap_16O(:) = 0._r8
       allocate(evap_18O(lsize)); evap_18O(:) = 0._r8
       allocate(evap_HDO(lsize)); evap_HDO(:) = 0._r8
    end if

    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_swdn', fldptr1=swdn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_swup', fldptr1=swup, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='Faox_lwup', fldptr1=lwup, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_fswpen', fldptr1=fswpen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! fields for history output only
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_warm_diurn', fldptr1=warm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_salt_diurn', fldptr1=salt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_speed_diurn', fldptr1=speed, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_regime_diurn', fldptr1=regime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_warmmax_diurn', fldptr1=warmMax, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_windmax_diurn', fldptr1=windMax, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_warmmaxinc_diurn', fldptr1=warmMaxInc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_qsolinc_diurn', fldptr1=qSolInc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_windinc_diurn', fldptr1=windInc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_ninc_diurn', fldptr1=nInc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_tbulk_diurn', fldptr1=tbulk, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_tskin_diurn', fldptr1=tskin, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_tskin_day_diurn', fldptr1=tskin_day, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_tskin_night_diurn', fldptr1=tskin_night, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_cskin_diurn', fldptr1=cskin, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_cskin_night_diurn', fldptr1=cskin_night, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_windavg_diurn', fldptr1=windAvg, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_qsolavg_diurn', fldptr1=qSolAvg, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtmOcn, fldname='So_windmaxinc_diurn', fldptr1=windmaxinc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Ocn import fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_t', fldptr1=tocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_u', fldptr1=uocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_v', fldptr1=vocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_s', fldptr1=ocnsal, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_fswpen', fldptr1=fswpen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_16O', fldptr1=roce_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_18O', fldptr1=roce_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_HDO', fldptr1=roce_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(roce_16O(lsize)); roce_16O(:) = 0._r8
       allocate(roce_18O(lsize)); roce_18O(:) = 0._r8
       allocate(roce_HDO(lsize)); roce_HDO(:) = 0._r8
    end if

    !----------------------------------
    ! Atm import fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_z', fldptr1=zbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_u', fldptr1=ubot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_v', fldptr1=vbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_tbot', fldptr1=tbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_ptem', fldptr1=thbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum', fldptr1=shum, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_16O', fldptr1=shum_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_18O', fldptr1=shum_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_HDO', fldptr1=shum_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(shum_16O(lsize)); shum_16O(:) = 0._r8
       allocate(shum_18O(lsize)); shum_18O(:) = 0._r8
       allocate(shum_HDO(lsize)); shum_HDO(:) = 0._r8
    end if

    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_dens', fldptr1=dens, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_lwdn', fldptr1=lwdn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_rainc', fldptr1=rainc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_rainl', fldptr1=rainl, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_snowc', fldptr1=snowc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_snowl', fldptr1=snowl, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_swndr', fldptr1=swndr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_swndf', fldptr1=swndf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_swvdr', fldptr1=swvdr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_swvdf', fldptr1=swvdf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! !----------------------------------
    ! ! Get lat, lon, which are time-invariant
    ! !----------------------------------

    ! The following assumes that all fields in FBAtmOcn have the same grid - so
    ! only need to query field 1
    call shr_nuopc_methods_FB_getFieldN(FBAtmOcn, fieldnum=1, field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine if first field in either on  grid or a mesh
    call ESMF_FieldGet(lfield, geomtype=geomtype, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (geomtype == ESMF_GEOMTYPE_MESH) then
       call ESMF_LogWrite(trim(subname)//" : FBATM is on a mesh ", ESMF_LOGMSG_INFO, rc=rc)
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
       call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(lons(numOwnedElements/spatialDim))
       allocate(lats(numOwnedElements/spatialDim))
       do n = 1,numOwnedElements/spatialDim
          lons(n) = ownedElemCoords(2*n-1)
          lats(n) = ownedElemCoords(2*n)
       end do
    else if (geomtype == ESMF_GEOMTYPE_GRID) then
       call ESMF_LogWrite(trim(subname)//" : FBATM is on a grid ", ESMF_LOGMSG_INFO, rc=rc)
       call ESMF_FieldGet(lfield, grid=lgrid, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_GridGet(lgrid, dimCount=dimCount, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_GridGetCoord(lgrid, coordDim=1, farrayPtr=lons, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_GridGetCoord(lgrid, coordDim=2, farrayPtr=lats, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//": ERROR FBATM must be either on a grid or a mesh", ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    end if

    !----------------------------------
    ! Fields that are not obtained via GetFldPtr
    !----------------------------------
    allocate(uGust(lsize))     ; uGust(:) =  0.0_r8
    allocate(prec(lsize))      ; prec(:) =  0.0_r8
    allocate(prec_gust(lsize)) ; prec_gust(:) =  0.0_r8

    !----------------------------------
    ! setup the compute mask.
    !----------------------------------

    ! prefer to compute just where ocean exists, so setup a mask here.
    ! this could be run with either the ocean or atm grid so need to be careful.
    ! really want the ocean mask on ocean grid or ocean mask mapped to atm grid,
    ! but do not have access to the ocean mask mapped to the atm grid.
    ! the dom mask is a good place to start, on ocean grid, it should be what we want,
    ! on the atm grid, it's just all 1's so not very useful.
    ! next look at ofrac+ifrac in fractions.  want to compute on all non-land points.
    ! using ofrac alone will exclude points that are currently all sea ice but that later
    ! could be less that 100% covered in ice.

    ! allocate grid mask fields
    allocate(emask(lsize)); emask(:) =  0.0_r8
    allocate(mask(lsize)) ; mask(:) =  0.0_r8

    ! default compute everywhere, then "turn off" gridcells
    ! mask(:) = 1

    ! allocate(rmask(lsize))      ;
    ! allocate(ifrac(lsize))
    ! allocate(ofrac(lsize))

    ! ! use domain mask first
    ! if (trim(aoflux_grid) == 'ocn') then
    !    ! In this case, FBFrac is the FBFrac_o
    !    call shr_nuopc_methods_FB_getFldPtr(FBFrac , 'ofrac' , rmask, rc=rc)
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    rmask(:) = 1._r8
    ! end if
    ! where (rmask(:) < 0.5_r8) mask(:) = 0   ! like nint

    ! ! then check ofrac + ifrac
    ! call shr_nuopc_methods_FB_getFldPtr(FBFrac , fldname='ofrac' , fldptr1=ofrac, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call shr_nuopc_methods_FB_getFldPtr(FBFrac , fldname='ifrac' , fldptr1=ifrac, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! where (ofrac(:) + ifrac(:) <= 0.0_r8) mask(:) = 0

    ! deallocate(rmask)
    ! deallocate(ifrac)
    ! deallocate(ofrac)

    ! emask(:) = mask(:)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_atmocn_init

  !===============================================================================

  subroutine med_atmocn_ocnalb(gcomp, FBFrac_o, nextsw_cday, rc)

    !-----------------------------------------------------------------------
    ! Update albedoes on ocean grid
    !-----------------------------------------------------------------------
    !
    ! Arguments
    type(ESMF_GridComp)    :: gcomp
    type(ESMF_FieldBundle) :: FBFrac_o     ! Fraction data on ocean grid
    real(r8), intent(in)   :: nextsw_cday  ! calendar day of next atm shortwave
    integer , intent(out)  :: rc
    !
    ! Local variables
    type(ESMF_VM)     :: vm
    integer(in)       :: iam
    character(CL)     :: cvalue
    real(r8), pointer :: ofrac_o(:)
    real(r8), pointer :: ofrad_o(:)
    real(r8), pointer :: ifrac_o(:)
    real(r8), pointer :: ifrad_o(:)
    integer(in)       :: lsize      ! local size
    logical	      :: flux_albav ! flux avg option
    logical	      :: update_alb ! was albedo updated
    integer(in)	      :: n,i	    ! indices
    real(r8)	      :: rlat	    ! gridcell latitude in radians
    real(r8)	      :: rlon	    ! gridcell longitude in radians
    real(r8)	      :: cosz	    ! Cosine of solar zenith angle
    real(r8)	      :: eccen	    ! Earth orbit eccentricity
    real(r8)	      :: mvelpp	    ! Earth orbit
    real(r8)	      :: lambm0	    ! Earth orbit
    real(r8)	      :: obliqr	    ! Earth orbit
    real(r8)	      :: delta	    ! Solar declination angle  in radians
    real(r8)	      :: eccf	    ! Earth orbit eccentricity factor
    real(r8)    , parameter :: albdif = 0.06_r8 ! 60 deg reference albedo, diffuse
    real(r8)    , parameter :: albdir = 0.07_r8 ! 60 deg reference albedo, direct
    character(*), parameter :: subName =   '(med_atmocn_ocnalb) '
    !-----------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flux_albav', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flux_albav

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    ! The following is for debugging
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine indices

    update_alb = .false.
    lsize = size(anidr)

    if (flux_albav) then

       do n = 1,lsize
          anidr(n) = albdir
          avsdr(n) = albdir
          anidf(n) = albdif
          avsdf(n) = albdif

          ! Albedo is now function of latitude (will be new implementation)
          !rlat = const_deg2rad * lats(n)
          !anidr = 0.069_r8 - 0.011_r8 * cos(2._r8 * rlat)
          !avsdr = anidr
          !anidf = anidr
          !avsdf = anidr
       end do
       update_alb = .true.

    else

       !--- flux_atmocn needs swdn & swup = swdn*(-albedo)
       !--- swdn & albedos are time-aligned  BEFORE albedos get updated below ---
       do n=1,lsize

          if ( anidr(n) == 1.0_r8 ) then ! dark side of earth
             swup(n) = 0.0_r8
             swdn(n) = 0.0_r8
          else
             swup(n) = swndr(n) * (-anidr(n)) + swndf(n) * (-anidf(n)) + &
                       swvdr(n) * (-avsdr(n)) + swvdf(n) * (-avsdf(n))
             swdn(n) = swndr(n) + swndf(n) + swvdr(n) + swvdf(n)
          end if
       end do

       ! Solar declination
       ! Will only do albedo calculation if nextsw_cday is not -1.

       if (nextsw_cday >= -0.5_r8) then
          call shr_orb_decl(nextsw_cday, eccen, mvelpp,lambm0, obliqr, delta, eccf)

          ! Compute albedos
          do n = 1,lsize
             rlat = const_deg2rad * lats(n)
             rlon = const_deg2rad * lons(n)
             cosz = shr_orb_cosz( nextsw_cday, rlat, rlon, delta )
             if (cosz  >  0.0_r8) then !--- sun hit --
                anidr(n) = (.026_r8/(cosz**1.7_r8 + 0.065_r8)) +   &
                           (.150_r8*(cosz         - 0.100_r8 ) *   &
                           (cosz - 0.500_r8 ) * (cosz - 1.000_r8 )  )
                avsdr(n) = anidr(n)
                anidf(n) = albdif
                avsdf(n) = albdif
             else !--- dark side of earth ---
                anidr(n) = 1.0_r8
                avsdr(n) = 1.0_r8
                anidf(n) = 1.0_r8
                avsdf(n) = 1.0_r8
             end if
          end do
          update_alb = .true.
       endif    ! nextsw_cday

    end if   ! flux_albav

    !--- update current ifrad/ofrad values if albedo was updated
    if (update_alb) then
       call shr_nuopc_methods_FB_getFldPtr(FBFrac_o, 'ifrac', fldptr1=ifrac_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(FBFrac_o, 'ifrad', fldptr1=ifrad_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_getFldPtr(FBFrac_o, fldname='ofrac', fldptr1=ofrac_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(FBFrac_o, fldname='ofrad', fldptr1=ofrad_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ifrad_o(:) = ifrac_o(:)
       ofrad_o(:) = ofrac_o(:)
    endif

  end subroutine med_atmocn_ocnalb

  !===============================================================================

  subroutine med_atmocn_flux(gcomp, clock, ocn_prognostic, dead_comps, rc)

    !-----------------------------------------------------------------------
    ! Determine atm/ocn fluxes eother on atm or on ocean grid
    ! The module arrays are set via pointers the the mediator internal states
    ! in med_atmocn_init and are used below.
    ! gcomp (the mediator gridded component) is only needed to retreive the
    ! attributes
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_Clock)      :: clock
    logical , intent(in)  :: ocn_prognostic
    logical , intent(in)  :: dead_comps
    integer , intent(out) :: rc
    !
    ! Local variables
    character(CL) :: cvalue
    integer(in)   :: tod, dt
    integer(in)   :: n,i                     ! indices
    integer(in)   :: lsize                   ! local size
    real(r8)      :: gust_fac = huge(1.0_r8) ! wind gust factor
    logical       :: cold_start              ! .true. to initialize internal fields in shr_flux diurnal
    logical       :: read_restart            ! .true. => continue run
    logical       :: flux_diurnal            ! .true. => turn on diurnal cycle in atm/ocn fluxes
    logical,save  :: first_call = .true.
    real(r8),parameter :: albdif = 0.06_r8 ! 60 deg reference albedo, diffuse
    real(r8),parameter :: albdir = 0.07_r8 ! 60 deg reference albedo, direct
    character(*),parameter :: subName =   '(med_atmocn_flux) '
    !-----------------------------------------------------------------------

    call seq_timemgr_EClockGetData( clock, curr_tod=tod, dtime=dt )

    call NUOPC_CompAttributeGet(gcomp, name='gust_fac', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) gust_fac

    call NUOPC_CompAttributeGet(gcomp, name='flux_diurnal', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flux_diurnal

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    cold_start = .false.   ! use restart data or data from last timestep
    if (first_call) then
       if (.not.read_restart) then
          cold_start = .true.
       end if
       first_call = .false.
    end if

    ! Determine local size
    lsize = size(lons)

    ! Update ocean surface fluxes
    ! Must fabricate "reasonable" data (when using dead components)

    if (dead_comps) then
       do n = 1,lsize
          mask        (n) =   1      ! ocn domain mask            ~ 0 <=> inactive cell
          tocn        (n) = 290.0_r8 ! ocn temperature            ~ Kelvin
          uocn        (n) =   0.0_r8 ! ocn velocity, zonal        ~ m/s
          vocn        (n) =   0.0_r8 ! ocn velocity, meridional   ~ m/s
          zbot        (n) =  55.0_r8 ! atm height of bottom layer ~ m
          ubot        (n) =   0.0_r8 ! atm velocity, zonal        ~ m/s
          vbot        (n) =   2.0_r8 ! atm velocity, meridional   ~ m/s
          thbot       (n) = 301.0_r8 ! atm potential temperature  ~ Kelvin
          shum        (n) = 1.e-2_r8 ! atm specific humidity      ~ kg/kg

          !wiso note: shum_* should be multiplied by Rstd_* here?
          shum_16O    (n) = 1.e-2_r8 ! H216O specific humidity ~ kg/kg
          shum_HDO    (n) = 1.e-2_r8 ! HD16O specific humidity ~ kg/kg
          shum_18O    (n) = 1.e-2_r8 ! H218O specific humidity ~ kg/kg
          roce_16O    (n) = 1.0_r8   ! H216O surface ratio     ~ mol/mol
          roce_HDO    (n) = 1.0_r8   ! HDO   surface ratio     ~ mol/mol
          roce_18O    (n) = 1.0_r8   ! H218O surface ratio     ~ mol/mol
          dens        (n) = 1.0_r8   ! atm density             ~ kg/m^3
          tbot        (n) = 300.0_r8 ! atm temperature         ~ Kelvin
          uGust       (n) = 0.0_r8
          lwdn        (n) = 0.0_r8
          prec        (n) = 0.0_r8
          prec_gust   (n) = 0.0_r8
          fswpen      (n) = 0.0_r8
          ocnsal      (n) = 0.0_r8
          warm        (n) = 0.0_r8
          salt        (n) = 0.0_r8
          speed       (n) = 0.0_r8
          regime      (n) = 0.0_r8
          warmMax     (n) = 0.0_r8
          windMax     (n) = 0.0_r8
          qSolAvg     (n) = 0.0_r8
          windAvg     (n) = 0.0_r8
          warmMaxInc  (n) = 0.0_r8
          windMaxInc  (n) = 0.0_r8
          qSolInc     (n) = 0.0_r8
          windInc     (n) = 0.0_r8
          nInc        (n) = 0.0_r8
          tbulk       (n) = 0.0_r8
          tskin       (n) = 0.0_r8
          tskin_day   (n) = 0.0_r8
          tskin_night (n) = 0.0_r8
          cskin       (n) = 0.0_r8
          cskin_night (n) = 0.0_r8
          swdn        (n) = 0.0_r8
          swup        (n) = 0.0_r8
       enddo
    else
       do n = 1,lsize
          nInc(n) = 0._r8 ! needed for minval/maxval calculation
          if (mask(n) /= 0) then
             !--- mask missing atm or ocn data if found
             if (dens(n) < 1.0e-12 .or. tocn(n) < 1.0) then
                emask(n) = 0
                !write(logunit,*) 'aoflux tcx1',n,dens(n),tocn(n)
             endif
             !!uGust(n) = 1.5_r8*sqrt(uocn(n)**2 + vocn(n)**2) ! there is no wind gust data from ocn
             uGust(n) = 0.0_r8
             prec(n) = rainc(n) + rainl(n) + snowc(n) + snowl(n)
             prec_gust(n) = rainc(n)
             ! Note: swdn and swup are set in flux_ocnalb using data from previous timestep
          end if
       enddo
    end if

    if (flux_diurnal) then
       call shr_flux_atmocn_diurnal (lsize , zbot , ubot, vbot, thbot, &
            shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
            tocn , emask, sen , lat , lwup , &
            roce_16O, roce_HDO, roce_18O,    &
            evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
            uGust, lwdn , swdn , swup, prec, &
            fswpen, ocnsal, ocn_prognostic, flux_diurnal,    &
            lats, lons , warm , salt , speed, regime,       &
            warmMax, windMax, qSolAvg, windAvg,             &
            warmMaxInc, windMaxInc, qSolInc, windInc, nInc, &
            tbulk, tskin, tskin_day, tskin_night, &
            cskin, cskin_night, tod, dt,          &
            duu10n,ustar, re  , ssq, &
            !missval should not be needed if flux calc
            !consistent with mrgx2a fraction
            !duu10n,ustar, re  , ssq, missval = 0.0_r8, &
            cold_start=cold_start)
    else
       call shr_flux_atmocn (lsize , zbot , ubot, vbot, thbot, prec_gust, gust_fac, &
            shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
            tocn , emask, sen , lat , lwup , &
            roce_16O, roce_HDO, roce_18O,    &
            evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
            duu10n,ustar, re  , ssq)
            !missval should not be needed if flux calc
            !consistent with mrgx2a fraction
            !duu10n,ustar, re  , ssq, missval = 0.0_r8 )
    endif

    do n = 1,lsize
       if (mask(n) /= 0) then
          u10(n) = sqrt(duu10n(n))
       end if
    enddo

  end subroutine med_atmocn_flux

end module med_atmocn_mod
