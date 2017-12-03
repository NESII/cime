module cesm_init_mod

  use ESMF
  use NUOPC
  use shr_kind_mod,      only: SHR_KIND_R8, SHR_KIND_CS, SHR_KIND_CL
  use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
  use shr_const_mod,     only: shr_const_cday
  use shr_file_mod,      only: shr_file_setLogLevel, shr_file_setLogUnit
  use shr_file_mod,      only: shr_file_setIO, shr_file_getUnit
  use shr_scam_mod,      only: shr_scam_checkSurface
  use shr_map_mod,       only: shr_map_setDopole
  use shr_mpi_mod,       only: shr_mpi_bcast, shr_mpi_chkerr
  use shr_mem_mod,       only: shr_mem_init, shr_mem_getusage
  use shr_cal_mod,       only: shr_cal_date2ymd
  use shr_orb_mod,       only: shr_orb_params
  use shr_frz_mod,       only: shr_frz_freezetemp_init
  use shr_reprosum_mod,  only: shr_reprosum_setopts
  use seq_comm_mct     , only: CPLID, GLOID, logunit, loglevel
  use seq_comm_mct     , only: ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
  use seq_comm_mct     , only: num_inst_atm, num_inst_lnd, num_inst_rof
  use seq_comm_mct     , only: num_inst_ocn, num_inst_ice, num_inst_glc
  use seq_comm_mct     , only: num_inst_wav, num_inst_esp
  use seq_comm_mct     , only: num_inst_xao, num_inst_frc, num_inst_phys
  use seq_comm_mct     , only: num_inst_total, num_inst_max
  use seq_comm_mct     , only: seq_comm_iamin, seq_comm_name, seq_comm_namelen, seq_comm_iamroot
  use seq_comm_mct     , only: seq_comm_init, seq_comm_setnthreads, seq_comm_getnthreads
  use seq_comm_mct     , only: seq_comm_gloroot
  use seq_comm_mct     , only: seq_comm_getinfo => seq_comm_setptrs
  use seq_timemgr_mod  , only: seq_timemgr_type
  use seq_timemgr_mod  , only: seq_timemgr_clockInit
  use seq_timemgr_mod  , only: seq_timemgr_clockPrint
  use seq_timemgr_mod  , only: seq_timemgr_EClockGetData
  use med_infodata_mod , only: med_infodata_init1, med_infodata_init2, med_infodata
  use seq_flds_mod     , only: seq_flds_set
  use t_drv_timers_mod
  use perf_mod

  ! MV: use seq_io_mod, only : seq_io_cpl_init
  ! MV: use seq_io_mod, only : seq_io_cpl_init
  ! MV: use cplcomp_exchange_mod, only: seq_mctext_decomp

  implicit none

  private

  public :: cesm_init

#include <mpif.h>

  character(len=*) , parameter    :: u_FILE_u =  __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine cesm_init(driver, &
                       Eclock_d, Eclock_a, Eclock_l, Eclock_o, &
                       Eclock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e)

    ! USES:
    use pio             , only: file_desc_t, pio_closefile, pio_file_is_open
    use shr_pio_mod     , only: shr_pio_init1, shr_pio_init2
    use shr_const_mod   , only: shr_const_tkfrz, shr_const_tktrip
    use shr_const_mod   , only: shr_const_mwwv, shr_const_mwdair
    use shr_wv_sat_mod  , only: shr_wv_sat_set_default, shr_wv_sat_init
    use shr_wv_sat_mod  , only: ShrWVSatTableSpec, shr_wv_sat_make_tables
    use shr_wv_sat_mod,   only: shr_wv_sat_get_scheme_idx, shr_wv_sat_valid_idx
    use shr_orb_mod     , only: SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL, shr_orb_params
    use shr_file_mod    , only: shr_file_getUnit, shr_file_freeUnit
    use shr_assert_mod  , only: shr_assert_in_domain
    use seq_io_read_mod , only: seq_io_read

    ! INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)    , intent(inout) :: driver
    type(ESMF_Clock)       , intent(inout) :: EClock_d
    type(ESMF_Clock)       , intent(inout) :: EClock_a
    type(ESMF_Clock)       , intent(inout) :: EClock_l
    type(ESMF_Clock)       , intent(inout) :: EClock_o
    type(ESMF_Clock)       , intent(inout) :: EClock_i
    type(ESMF_Clock)       , intent(inout) :: EClock_g
    type(ESMF_Clock)       , intent(inout) :: EClock_r
    type(ESMF_Clock)       , intent(inout) :: EClock_w
    type(ESMF_Clock)       , intent(inout) :: EClock_e

    ! LOCAL
    ! threading control
    integer                         :: Global_Comm
    integer                         :: mpicom_GLOID          ! MPI global communicator
    integer                         :: mpicom_CPLID          ! MPI cpl communicator
    integer                         :: mpicom_OCNID          ! MPI ocn communicator for ensemble member 1
    integer                         :: iam_GLOID             ! pe number in global id
    logical                         :: iamin_CPLID           ! pe associated with CPLID
    logical                         :: iamroot_GLOID         ! GLOID masterproc
    logical                         :: iamroot_CPLID         ! CPLID masterproc
    integer                         :: nthreads_GLOID        ! OMP global number of threads
    integer                         :: nthreads_CPLID        ! OMP cpl number of threads
    integer                         :: nthreads_ATMID        ! OMP atm number of threads
    integer                         :: nthreads_LNDID        ! OMP lnd number of threads
    integer                         :: nthreads_ICEID        ! OMP ice number of threads
    integer                         :: nthreads_OCNID        ! OMP ocn number of threads
    integer                         :: nthreads_GLCID        ! OMP glc number of threads
    integer                         :: nthreads_ROFID        ! OMP glc number of threads
    integer                         :: nthreads_WAVID        ! OMP wav number of threads
    integer                         :: nthreads_ESPID        ! OMP esp number of threads
    integer                         :: pethreads_GLOID       ! OMP number of threads per task
    logical                         :: drv_threading         ! driver threading control
    character(SHR_KIND_CL)          :: cpl_seq_option        ! coupler sequencing option
    logical                         :: reprosum_use_ddpdd    ! setup reprosum, use ddpdd
    real(SHR_KIND_R8)               :: reprosum_diffmax      ! setup reprosum, set rel_diff_max
    logical                         :: reprosum_recompute    ! setup reprosum, recompute if tolerance exceeded
    logical                         :: output_perf = .false. ! require timing data output for this pe
    integer                         :: ymd                   ! Current date (YYYYMMDD)
    integer                         :: year                  ! Current date (YYYY)
    integer                         :: month                 ! Current date (MM)
    integer                         :: day                   ! Current date (DD)
    integer                         :: tod                   ! Current time of day (seconds)
    character(SHR_KIND_CL)          :: orb_mode              ! orbital mode
    character(SHR_KIND_CS)          :: tfreeze_option        ! Freezing point calculation
    integer                         :: orb_iyear             ! orbital year
    integer                         :: orb_iyear_align       ! associated with model year
    integer                         :: orb_cyear             ! orbital year for current orbital computation
    integer                         :: orb_nyear             ! orbital year associated with currrent model year
    real(SHR_KIND_R8)               :: orb_eccen             ! orbital eccentricity
    real(SHR_KIND_R8)               :: orb_obliq             ! obliquity in degrees
    real(SHR_KIND_R8)               :: orb_mvelp             ! moving vernal equinox long
    real(SHR_KIND_R8)               :: orb_obliqr            ! Earths obliquity in rad
    real(SHR_KIND_R8)               :: orb_lambm0            ! Mean long of perihelion at vernal equinox (radians)
    real(SHR_KIND_R8)               :: orb_mvelpp            ! moving vernal equinox long
    real(SHR_KIND_R8)               :: wall_time_limit       ! wall time limit in hours
    character(SHR_KIND_CS)          :: force_stop_at         ! force stop at next (month, day, etc)
    character(SHR_KIND_CS)          :: cime_model            ! currently acme or cesm
    character(SHR_KIND_CL)          :: atm_gnam              ! atm grid
    character(SHR_KIND_CL)          :: lnd_gnam              ! lnd grid
    character(SHR_KIND_CL)          :: ocn_gnam              ! ocn grid
    character(SHR_KIND_CL)          :: ice_gnam              ! ice grid
    character(SHR_KIND_CL)          :: rof_gnam              ! rof grid
    character(SHR_KIND_CL)          :: glc_gnam              ! glc grid
    character(SHR_KIND_CL)          :: wav_gnam              ! wav grid
    character(SHR_KIND_CL)          :: samegrid_ao           ! samegrid atm and ocean
    character(SHR_KIND_CL)          :: samegrid_al           ! samegrid atm and land
    character(SHR_KIND_CL)          :: samegrid_lr           ! samegrid land and rof
    character(SHR_KIND_CL)          :: samegrid_oi           ! samegrid ocean and ice
    character(SHR_KIND_CL)          :: samegrid_ro           ! samegrid runoff and ocean
    character(SHR_KIND_CL)          :: samegrid_aw           ! samegrid atm and wave
    character(SHR_KIND_CL)          :: samegrid_ow           ! samegrid ocean and wave
    character(SHR_KIND_CL)          :: samegrid_lg           ! samegrid glc and land
    character(SHR_KIND_CL)          :: samegrid_og           ! samegrid glc and ocean
    character(SHR_KIND_CL)          :: samegrid_ig           ! samegrid glc and ice
    character(SHR_KIND_CL)          :: samegrid_alo          ! samegrid atm, lnd, ocean
    logical                         :: shr_map_dopole        ! pole corrections in shr_map_mod
    logical                         :: single_column         ! scm mode logical
    real(SHR_KIND_R8)               :: scmlon                ! single column lon
    real(SHR_KIND_R8)               :: scmlat                ! single column lat
    logical                         :: atm_aero              ! atm provides aerosol data
    type(file_desc_t)               :: pioid
    integer                         :: maxthreads
    character(SHR_KIND_CS)          :: wv_sat_scheme
    real(SHR_KIND_R8)               :: wv_sat_transition_start
    logical                         :: wv_sat_use_tables
    real(SHR_KIND_R8)               :: wv_sat_table_spacing
    character(SHR_KIND_CL)          :: errstring
    character(SHR_KIND_CL)          :: cvalue
    type(ShrWVSatTableSpec)         :: liquid_spec
    type(ShrWVSatTableSpec)         :: ice_spec
    type(ShrWVSatTableSpec)         :: mixed_spec
    integer                         :: comp_id(num_inst_total)
    integer                         :: comp_comm(num_inst_total)
    integer                         :: comp_comm_iam(num_inst_total)
    logical                         :: comp_iamin(num_inst_total)
    character(len=seq_comm_namelen) :: comp_name(num_inst_total)
    logical                         :: flag
    integer                         :: i, it, n
    character(SHR_KIND_CL)          :: start_type            ! Type of startup
    logical                         :: read_restart          ! read the restart file, based on start_type
    character(SHR_KIND_CL)          :: restart_file          ! Full archive path to restart file
    character(SHR_KIND_CL)          :: restart_pfile         ! Restart pointer file
    character(SHR_KIND_CL)          :: rest_case_name        ! Short case identification
    integer                         :: unitn                 ! Namelist unit number to read
    logical                         :: exists                ! true if file exists
    integer                         :: ierr                  ! MPI error return
    integer                         :: rc                    ! return code
    logical                         :: atm_present           ! .true.  => atm is present
    logical                         :: lnd_present           ! .true.  => land is present
    logical                         :: ice_present           ! .true.  => ice is present
    logical                         :: ocn_present           ! .true.  => ocn is present
    logical                         :: glc_present           ! .true.  => glc is present
    logical                         :: glclnd_present        ! .true.  => glc is computing land coupling
    logical                         :: glcocn_present        ! .true.  => glc is computing ocean runoff
    logical                         :: glcice_present        ! .true.  => glc is computing icebergs
    logical                         :: rofice_present        ! .true.  => rof is computing icebergs
    logical                         :: rof_present           ! .true.  => rof is present
    logical                         :: flood_present         ! .true.  => rof is computing flood
    logical                         :: wav_present           ! .true.  => wav is present
    logical                         :: esp_present           ! .true.  => esp is present
    character(len=*) , parameter    :: NLFileName = "drv_in" ! input namelist filename
    integer          , parameter    :: ens1=1                ! use first instance of ensemble only
    integer          , parameter    :: fix1=1                ! temporary hard-coding to first ensemble, needs to be fixed
    real(SHR_KIND_R8), parameter    :: epsilo = shr_const_mwwv/shr_const_mwdair
    character(len=*) , parameter    :: sp_str = 'str_undefined'
    character(len=*) , parameter    :: start_type_start = "startup"
    character(len=*) , parameter    :: start_type_cont  = "continue"
    character(len=*) , parameter    :: start_type_brnch = "branch"
    character(len=*) , parameter    :: orb_fixed_year       = 'fixed_year'
    character(len=*) , parameter    :: orb_variable_year    = 'variable_year'
    character(len=*) , parameter    :: orb_fixed_parameters = 'fixed_parameters'
    character(len=*) , parameter    :: subname = '(cesm_init)'

    !----------------------------------------------------------
    !| Initialize MCT and MPI communicators and IO
    !----------------------------------------------------------

    call mpi_initialized(flag,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_initialized')
    if (.not. flag) then
       call mpi_init(ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_init')
    endif

    Global_Comm=MPI_COMM_WORLD
    comp_comm = MPI_COMM_NULL

    call shr_pio_init1(num_inst_total,NLFileName, Global_Comm)
    !
    ! If pio_async_interface is true Global_Comm is MPI_COMM_NULL on the servernodes
    ! and server nodes do not return from shr_pio_init2
    !
    !   if (Global_Comm /= MPI_COMM_NULL) then

    call seq_comm_init(Global_Comm, NLFileName)

    !--- set task based threading counts ---
    call seq_comm_getinfo(GLOID,pethreads=pethreads_GLOID,iam=iam_GLOID)
    call seq_comm_setnthreads(pethreads_GLOID)

    !--- get some general data ---
    it=1
    call seq_comm_getinfo(GLOID,mpicom=mpicom_GLOID,&
         iamroot=iamroot_GLOID,nthreads=nthreads_GLOID)

    call seq_comm_getinfo(CPLID,mpicom=mpicom_CPLID,&
         iamroot=iamroot_CPLID,nthreads=nthreads_CPLID,&
         iam=comp_comm_iam(it))

    comp_id(it)    = CPLID
    comp_comm(it)  = mpicom_CPLID
    iamin_CPLID    = seq_comm_iamin(CPLID)
    comp_iamin(it) = seq_comm_iamin(comp_id(it))
    comp_name(it)  = seq_comm_name(comp_id(it))
    do n = 1,num_inst_atm
       it=it+1
       comp_id(it)    = ATMID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ATMID(n), mpicom=comp_comm(it), nthreads=nthreads_ATMID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_lnd
       it=it+1
       comp_id(it)    = LNDID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(LNDID(n), mpicom=comp_comm(it), nthreads=nthreads_LNDID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_ocn
       it=it+1
       comp_id(it)    = OCNID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(OCNID(n), mpicom=comp_comm(it), nthreads=nthreads_OCNID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_ice
       it=it+1
       comp_id(it)    = ICEID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ICEID(n), mpicom=comp_comm(it), nthreads=nthreads_ICEID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_glc
       it=it+1
       comp_id(it)    = GLCID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(GLCID(n), mpicom=comp_comm(it), nthreads=nthreads_GLCID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_rof
       it=it+1
       comp_id(it)    = ROFID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ROFID(n), mpicom=comp_comm(it), nthreads=nthreads_ROFID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_wav
       it=it+1
       comp_id(it)    = WAVID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(WAVID(n), mpicom=comp_comm(it), nthreads=nthreads_WAVID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_esp
       it=it+1
       comp_id(it)    = ESPID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ESPID(n), mpicom=comp_comm(it), nthreads=nthreads_ESPID, iam=comp_comm_iam(it))
    enddo
    ! ESP components do not use the coupler (they are 'external')

    !----------------------------------------------------------
    !| Set logging parameters both for shr code and locally
    !----------------------------------------------------------

    if (iamroot_CPLID) then
       inquire(file='cpl_modelio.nml',exist=exists)
       if (exists) then
          logunit = shr_file_getUnit()
          call shr_file_setIO('cpl_modelio.nml',logunit)
          call shr_file_setLogUnit(logunit)
          loglevel = 1
          call shr_file_setLogLevel(loglevel)
       endif
    else
       loglevel = 0
       call shr_file_setLogLevel(loglevel)
    endif

    !----------------------------------------------------------
    ! Log info about the environment settings
    !----------------------------------------------------------

    !  When using io servers (pio_async_interface=.true.) the server tasks do not return from
    !  shr_pio_init2

    call shr_pio_init2(comp_id,comp_name,comp_iamin,comp_comm,comp_comm_iam)

    !----------------------------------------------------------
    ! Print Model heading and copyright message
    !----------------------------------------------------------

    ! MV: if (iamroot_CPLID) call seq_cesm_printlogheader()

    !----------------------------------------------------------
    !| Timer initialization (has to be after mpi init)
    !----------------------------------------------------------

    maxthreads = max(nthreads_GLOID,nthreads_CPLID,nthreads_ATMID, &
         nthreads_LNDID,nthreads_ICEID,nthreads_OCNID,nthreads_GLCID, &
         nthreads_ROFID, nthreads_WAVID, nthreads_ESPID, pethreads_GLOID )

    call t_initf(NLFileName, LogPrint=.true., mpicom=mpicom_GLOID, &
         MasterTask=iamroot_GLOID,MaxThreads=maxthreads)

    if (iamin_CPLID) then
       ! TODO: where should this be called
       ! MV: call seq_io_cpl_init()
    endif

    call t_startf('CPL:INIT')
    call t_adj_detailf(+1)

    call t_startf('CPL:cesm_init')

    !----------------------------------------------------------
    ! Memory test
    !----------------------------------------------------------

    call shr_mem_init(prt=iamroot_CPLID)

    !----------------------------------------------------------
    ! Initialize infodata
    !----------------------------------------------------------

    call med_infodata_init1(med_infodata, GLOID)

    !----------------------------------------------------------
    ! Add atm_aero to driver attributes
    !----------------------------------------------------------

    call NUOPC_CompAttributeAdd(driver, attrList=(/'atm_aero'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    ! set initial value to .false.
    call NUOPC_CompAttributeSet(driver, name='atm_aero', value='.false.', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    !----------------------------------------------------------
    ! Deterine same grid attributes
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="atm_gnam", value=atm_gnam, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeGet(driver, name="lnd_gnam", value=lnd_gnam, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeGet(driver, name="rof_gnam", value=rof_gnam, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeGet(driver, name="ice_gnam", value=ice_gnam, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeGet(driver, name="ocn_gnam", value=ocn_gnam, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeGet(driver, name="wav_gnam", value=wav_gnam, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    samegrid_ao  = '.true.'
    samegrid_al  = '.true.'
    samegrid_lr  = '.true.'
    samegrid_oi  = '.true.'
    samegrid_ro  = '.true.'
    samegrid_aw  = '.true.'
    samegrid_ow  = '.true.'
    samegrid_lg  = '.true.'
    samegrid_og  = '.true.'
    samegrid_ig  = '.true.'
    samegrid_alo = '.true.'

    ! set samegrid to true for single column
    call NUOPC_CompAttributeGet(driver, name="single_column", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) single_column
    if (.not. single_column) then
       if (trim(atm_gnam) /= trim(ocn_gnam)) samegrid_ao = '.false.'
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = '.false.'
       if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = '.false.'
       if (trim(rof_gnam) /= trim(ocn_gnam)) samegrid_ro = '.false.'
       if (trim(ocn_gnam) /= trim(ice_gnam)) samegrid_oi = '.false.'
       if (trim(atm_gnam) /= trim(wav_gnam)) samegrid_aw = '.false.'
       if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = '.false.'
       if (trim(lnd_gnam) /= trim(glc_gnam)) samegrid_lg = '.false.'
       if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_og = '.false.'
       if (trim(ice_gnam) /= trim(glc_gnam)) samegrid_ig = '.false.'
       if (samegrid_al == '.true.' .and. samegrid_ao == '.true.') then
          samegrid_alo = '.true.'
       else
          samegrid_alo = '.false.'
       end if
    endif

    call NUOPC_CompAttributeAdd(driver, attrList=(/'samegrid_ao', 'samegrid_al', 'samegrid_lr', &
         'samegrid_oi', 'samegrid_ro', 'samegrid_aw', 'samegrid_ow', 'samegrid_lg', &
         'samegrid_og', 'samegrid_ig', 'samegrid_alo'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_ao", value=samegrid_ao, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_al", value=samegrid_al, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_lr", value=samegrid_lr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_oi", value=samegrid_oi, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_ro", value=samegrid_ro, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_aw", value=samegrid_aw, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_ow", value=samegrid_ow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_lg", value=samegrid_lg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_og", value=samegrid_og, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_ig", value=samegrid_ig, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call NUOPC_CompAttributeSet(driver, name="samegrid_alo", value=samegrid_alo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    !----------------------------------------------------------
    ! Check consistency of driver attributes
    !----------------------------------------------------------

    call driver_attributes_check(driver)

    !----------------------------------------------------------
    ! Initialize coupled fields
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="cime_model", value=cime_model, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    if ( trim(cime_model) /= 'cesm') then
       call shr_sys_abort( subname//': cime_model must be set to cesm, aborting')
    end if

    call seq_flds_set(nlfilename, GLOID, cime_model)

    !----------------------------------------------------------
    ! Initialize dopole flag (as a module variable in shr_map_mod)
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="shr_map_dopole", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) shr_map_dopole

    call shr_map_setDopole(shr_map_dopole)

    !----------------------------------------------------------
    ! Initialize options for reproducible sums
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="reprosum_use_ddpdd", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) reprosum_use_ddpdd

    call NUOPC_CompAttributeGet(driver, name="reprosum_diffmax", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) reprosum_diffmax

    call NUOPC_CompAttributeGet(driver, name="reprosum_recompute", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) reprosum_recompute

    call shr_reprosum_setopts(repro_sum_use_ddpdd_in=reprosum_use_ddpdd, &
         repro_sum_rel_diff_max_in=reprosum_diffmax, repro_sum_recompute_in=reprosum_recompute)

    !----------------------------------------------------------
    ! Check cpl_seq_option
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="cpl_seq_option", value=cpl_seq_option, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    if (trim(cpl_seq_option) /= 'CESM1_ORIG' .and. &
        trim(cpl_seq_option) /= 'CESM1_ORIG_TIGHT' .and. &
        trim(cpl_seq_option) /= 'CESM1_MOD' .and. &
        trim(cpl_seq_option) /= 'CESM1_MOD_TIGHT' .and. &
        trim(cpl_seq_option) /= 'RASM_OPTION1' .and. &
        trim(cpl_seq_option) /= 'RASM_OPTION2' ) then
        call shr_sys_abort(subname//' invalid cpl_seq_option = '//trim(cpl_seq_option))
    endif

    !----------------------------------------------------------
    ! Test Threading Setup in driver happens to be valid on all pes for all IDs
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="drv_threading", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) drv_threading

    if (drv_threading) then
       if (iamroot_GLOID) write(logunit,*) ' '
       if (iamroot_GLOID) write(logunit,'(2A)    ') subname,' Test Threading in driver'
       call seq_comm_setnthreads(nthreads_GLOID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLOID = ',nthreads_GLOID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_CPLID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_CPLID = ',nthreads_CPLID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ATMID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ATMID = ',nthreads_ATMID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_LNDID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_LNDID = ',nthreads_LNDID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_OCNID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_OCNID = ',nthreads_OCNID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ICEID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ICEID = ',nthreads_ICEID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_GLCID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLCID = ',nthreads_GLCID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ROFID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ROFID = ',nthreads_ROFID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_WAVID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_WAVID = ',nthreads_WAVID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ESPID)
       if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ESPID = ',nthreads_ESPID,seq_comm_getnthreads()
       if (iamroot_GLOID) write(logunit,*) ' '
       call seq_comm_setnthreads(nthreads_GLOID)
    endif

    !-----------------------------------------------------
    ! Determine if restart is read
    !-----------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name='start_type', value=start_type, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    ! Check valid values of start type
    call NUOPC_CompAttributeGet(driver, name="start_type", value=start_type, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       call shr_sys_abort(subname//': start_type invalid = '//trim(start_type))
    end if

    read_restart = .false.
    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       read_restart = .true.
    endif

    ! Add rest_case_name and read_restart to driver attributes
    call NUOPC_CompAttributeAdd(driver, attrList=(/'rest_case_name','read_restart'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    rest_case_name = ' '
    call NUOPC_CompAttributeSet(driver, name='rest_case_name', value=rest_case_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    write(cvalue,*) read_restart
    call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(cvalue), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    !-----------------------------------------------------
    ! Read Restart (seq_io_read must be called on all pes)
    !-----------------------------------------------------

    ! Error check on restart_pfile
    call NUOPC_CompAttributeGet(driver, name="restart_pfile", value=restart_pfile, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    if ( len_trim(restart_pfile) == 0 ) then
       call shr_sys_abort( subname//': restart_pfile must be set' )
    end if

    if (read_restart) then
       !--- read rpointer if restart_file is set to sp_str ---
       if (seq_comm_iamroot(GLOID)) then

          call NUOPC_CompAttributeGet(driver, name='restart_file', value=restart_file, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

          if (trim(restart_file) == trim(sp_str)) then
             ! Read pointer file

             call NUOPC_CompAttributeGet(driver, name='restart_pfile', value=restart_file, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

             unitn = shr_file_getUnit()
             if (loglevel > 0) write(logunit,"(3A)") subname," read rpointer file ", trim(restart_pfile)
             open(unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
             if (ierr < 0) then
                call shr_sys_abort( subname//':: rpointer file open returns an'// ' error condition' )
             end if
             read(unitn,'(a)', iostat=ierr) restart_file
             if (ierr < 0) then
                call shr_sys_abort( subname//':: rpointer file read returns an'// ' error condition' )
             end if
             close(unitn)
             call shr_file_freeUnit( unitn )
             write(logunit,"(3A)") subname,' restart file from rpointer= ', trim(restart_file)
          endif
       endif
       call shr_mpi_bcast(restart_file,mpicom_GLOID)

       call NUOPC_CompAttributeSet(driver, name='restart_pfile', value=restart_file, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    endif

    !----------------------------------------------------------
    ! Initialize time manager
    !----------------------------------------------------------

    call seq_timemgr_clockInit(driver, pioid, mpicom_gloid, &
         EClock_d, EClock_a, EClock_l, EClock_o, &
         EClock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e)

    !----------------------------------------------------------
    ! Initialize infodata items which need the clocks
    !----------------------------------------------------------

    call med_infodata_init2(med_infodata)

    !----------------------------------------------------------
    ! Initialize freezing point calculation for all components
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="tfreeze_option", value=tfreeze_option, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    call shr_frz_freezetemp_init(tfreeze_option)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="orb_mode", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(driver, name="orb_iyear", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(driver, name="orb_iyear_align", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(driver, name="orb_obliq", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(driver, name="orb_eccen", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(driver, name="orb_mvelp", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) orb_mvelp

    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
          call shr_sys_abort(subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode))
       endif

    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
          call shr_sys_abort(subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode))
       endif

    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
          write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
          write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
          call shr_sys_abort(subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode))
       endif
    else
       call shr_sys_abort(subname//' ERROR: invalid orb_mode '//trim(orb_mode))
    endif

    ! Determine orbital params

    if (trim(orb_mode) == trim(orb_variable_year)) then
       call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd)
       call shr_cal_date2ymd(ymd,year,month,day)
       orb_cyear = orb_iyear + (year - orb_iyear_align)
       call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                           orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
    else
       call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, &
                           orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
    end if

    if (orb_eccen  == SHR_ORB_UNDEF_REAL .or. &
        orb_obliqr == SHR_ORB_UNDEF_REAL .or. &
        orb_mvelpp == SHR_ORB_UNDEF_REAL .or. &
        orb_lambm0 == SHR_ORB_UNDEF_REAL) then
       call shr_sys_abort(subname//': orb params incorrect')
    endif

    ! Add updated orbital params to driver attributes

    call NUOPC_CompAttributeAdd(driver, attrList=(/'orb_obliqr', 'orb_lambm0', 'orb_mvelpp'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    write(cvalue,*) orb_eccen
    call NUOPC_CompAttributeSet(driver, name="orb_eccen", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    write(cvalue,*) orb_obliqr
    call NUOPC_CompAttributeSet(driver, name="orb_obliqr", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    write(cvalue,*) orb_lambm0
    call NUOPC_CompAttributeSet(driver, name="orb_lambm0", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    write(cvalue,*) orb_mvelpp
    call NUOPC_CompAttributeSet(driver, name="orb_mvelpp", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    !----------------------------------------------------------
    ! Initialize water vapor info
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="wv_sat_scheme", value=wv_sat_scheme, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    if (.not. shr_wv_sat_valid_idx(shr_wv_sat_get_scheme_idx(trim(wv_sat_scheme)))) then
       call shr_sys_abort(subname//': "'//trim(wv_sat_scheme)//'" is not a recognized saturation vapor pressure scheme name')
    end if
    if (.not. shr_wv_sat_set_default(wv_sat_scheme)) then
       call shr_sys_abort('Invalid wv_sat_scheme.')
    end if

    call NUOPC_CompAttributeGet(driver, name="wv_sat_transition_start", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) wv_sat_transition_start

    call shr_assert_in_domain(wv_sat_transition_start, &
         ge=0._SHR_KIND_R8, le=40._SHR_KIND_R8, &
         varname="wv_sat_transition_start", msg="Invalid transition temperature range.")

    call NUOPC_CompAttributeGet(driver, name="wv_sat_use_tables", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) wv_sat_use_tables

    call NUOPC_CompAttributeGet(driver, name="wv_sat_table_spacing", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) wv_sat_table_spacing

    ! A transition range averaging method in CAM is only valid for:
    ! -40 deg C <= T <= 0 deg C
    ! shr_wv_sat_mod itself checks for values with the wrong sign, but we
    ! have to check that the range is no more than 40 deg C here. Even
    ! though this is a CAM-specific restriction, it's not really likely
    ! that any other parameterization will be dealing with mixed-phase
    ! water below 40 deg C anyway.

    call shr_wv_sat_init(shr_const_tkfrz, shr_const_tktrip, wv_sat_transition_start, epsilo, errstring)
    if (errstring /= "") then
       call shr_sys_abort('shr_wv_sat_init: '//trim(errstring))
    end if

    ! The below produces internal lookup tables in the range 175-374K for
    ! liquid water, and 125-274K for ice, with a resolution set by the
    ! option wv_sat_table_spacing.
    ! In theory these ranges could be specified in the namelist, but in
    ! practice users will want to change them *very* rarely if ever, which
    ! is why only the spacing is in the namelist.

    if (wv_sat_use_tables) then
       liquid_spec = ShrWVSatTableSpec(ceiling(200._SHR_KIND_R8/wv_sat_table_spacing), 175._SHR_KIND_R8, wv_sat_table_spacing)
       ice_spec    = ShrWVSatTableSpec(ceiling(150._SHR_KIND_R8/wv_sat_table_spacing), 125._SHR_KIND_R8, wv_sat_table_spacing)
       mixed_spec  = ShrWVSatTableSpec(ceiling(250._SHR_KIND_R8/wv_sat_table_spacing), 125._SHR_KIND_R8, wv_sat_table_spacing)
       call shr_wv_sat_make_tables(liquid_spec, ice_spec, mixed_spec)
    end if

    !----------------------------------------------------------
    ! Set single_column flags
    ! If in single column mode, overwrite flags according to focndomain file
    ! in ocn_in namelist. SCAM can reset the "present" flags for lnd,
    ! ocn, ice, rof, and flood.
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="single_column", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    read(cvalue,*) single_column

    ! NOTE: cam stand-alone aqua-planet model will no longer be supported here - only the data model aqua-planet
    ! will be supported
    if (single_column) then

       call NUOPC_CompAttributeGet(driver, name="scmlon", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
       read(cvalue,*) scmlon

       call NUOPC_CompAttributeGet(driver, name="scmlat", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
       read(cvalue,*) scmlat

       call seq_comm_getinfo(OCNID(ens1), mpicom=mpicom_OCNID)

       ! TODO: Single column mode needs to be re-implemented - previously all of the xxx_present flags were set
       ! in med_infodata calls, reset here and the put back into med_infodata - this is no longer the case
       call shr_scam_checkSurface(scmlon, scmlat, &
            OCNID(ens1), mpicom_OCNID,            &
            lnd_present=lnd_present,              &
            ocn_present=ocn_present,              &
            ice_present=ice_present,              &
            rof_present=rof_present,              &
            flood_present=flood_present,          &
            rofice_present=rofice_present)

    endif

    !----------------------------------------------------------
    ! Finalize initialization
    !----------------------------------------------------------

    if(PIO_FILE_IS_OPEN(pioid)) then
       call pio_closefile(pioid)
    endif

    call t_stopf('CPL:cesm_init')

    call t_adj_detailf(-1)
    call t_stopf('CPL:INIT')

  end subroutine cesm_init

  !===============================================================================
  subroutine driver_attributes_check( driver )

    ! !DESCRIPTION: Check that input driver config values have reasonable values

    ! !USES:

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)    , intent(INOUT) :: driver

    !----- local -----
    integer                :: lastchar       ! Last character index
    character(SHR_KIND_CL) :: start_type     ! Type of startup
    character(SHR_KIND_CL) :: rest_case_name ! Short case identification
    character(SHR_KIND_CL) :: case_name      ! Short case identification
    character(SHR_KIND_CS) :: aoflux_grid    ! grid for atm ocn flux calc
    character(SHR_KIND_CL) :: vect_map       ! vector mapping option, none, cart3d, cart3d_diag, cart3d_uvw, cart3d_uvw_diag
    character(SHR_KIND_CS) :: logFilePostFix ! postfix for output log files
    character(SHR_KIND_CL) :: outPathRoot    ! root for output log files
    character(SHR_KIND_CL) :: cvalue         ! temporary
    integer                :: rc             ! return code
    character(len=*), parameter :: u_FILE_u =  __FILE__
    character(len=*), parameter :: subname = '(driver_attributes_check) '
    !-------------------------------------------------------------------------------

    ! --- Case name ------
    call NUOPC_CompAttributeGet(driver, name="case_name", value=case_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    lastchar = len(case_name)
    if ( len_trim(case_name) == 0) then
       call shr_sys_abort( subname//': variable case_name must be set, aborting')
    end if
    if (case_name(lastchar:lastchar) /= ' ') then
       write(logunit,"(A,I4,A)")'ERROR: case_name must not exceed ', len(case_name)-1,' characters'
       call shr_sys_abort( subname//': variable case_name must be set, aborting')
    end if

    ! --- LogFile ending name -----
    call NUOPC_CompAttributeGet(driver, name="logFilePostFix", value=logFilePostFix, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort

    if ( len_trim(logFilePostFix) == 0 ) then
       call shr_sys_abort( subname//': logFilePostFix  must be set to something not blank' )
    end if

    ! --- Output path root directory -----
    call NUOPC_CompAttributeGet(driver, name="outPathRoot", value=outPathRoot, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort

    if ( len_trim(outPathRoot) == 0 ) then
       call shr_sys_abort( subname//': outPathRoot  must be set' )
    end if
    if ( index(outPathRoot, "/", back=.true.) /= len_trim(outPathRoot) ) then
       call shr_sys_abort( subname//': outPathRoot must end with a slash' )
    end if

    ! --- Grid for atm/ocean flux computations ----
    call NUOPC_CompAttributeGet(driver, name="aoflux_grid", value=aoflux_grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    if ((trim(aoflux_grid) /= 'ocn') .and. &
        (trim(aoflux_grid) /= 'atm') .and. &
        (trim(aoflux_grid) /= 'exch')) then
       write(logunit,'(2a)') 'ERROR aoflux_grid not supported = ',trim(aoflux_grid)
       call shr_sys_abort(subname//': aoflux_grid invalid = '//trim(aoflux_grid))
    endif

    ! --- Vector mapping options ----
    call NUOPC_CompAttributeGet(driver, name="vect_map", value=vect_map, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    if ((trim(vect_map) /= 'none')        .and. &
        (trim(vect_map) /= 'cart3d')      .and. &
        (trim(vect_map) /= 'cart3d_diag') .and. &
        (trim(vect_map) /= 'cart3d_uvw')  .and. &
        (trim(vect_map) /= 'cart3d_uvw_diag')) then
       write(logunit,'(2a)') 'ERROR vect_map not supported = ',trim(vect_map)
       call shr_sys_abort(subname//': vect_map invalid = '//trim(vect_map))
    endif

    ! --- Case name and restart case name ------
    ! call NUOPC_CompAttributeGet(driver, name="rest_case_name", value=rest_case_name, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()

    ! if ((trim(start_type) == start_type_cont ) .and. (trim(case_name)  /= trim(rest_case_name))) then
    !    write(logunit,'(10a)') subname,' case_name =',trim(case_name),':',' rest_case_name =',trim(rest_case_name),':'
    !    call shr_sys_abort(subname//': invalid continue restart case name = '//trim(rest_case_name))
    ! endif

  end subroutine driver_attributes_check

end module cesm_init_mod
