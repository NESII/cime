module med_infodata_mod

  ! !DESCRIPTION: A module to get, put, and store some standard scalar data

  ! !USES:
  use ESMF
  use NUOPC
  use mpi
  use shr_kind_mod, only: SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_IN, SHR_KIND_R8, SHR_KIND_I8
  use shr_sys_mod,  only: shr_sys_flush, shr_sys_abort
  use seq_comm_mct, only: logunit, loglevel
  use seq_comm_mct, only: num_inst_atm, num_inst_lnd, num_inst_rof
  use seq_comm_mct, only: num_inst_ocn, num_inst_ice, num_inst_glc
  use seq_comm_mct, only: num_inst_wav
  use seq_flds_mod, only: seq_flds_scalar_num, seq_flds_scalar_name
  use seq_flds_mod, only: seq_flds_scalar_index_nx,  seq_flds_scalar_index_ny
  use seq_flds_mod, only: seq_flds_scalar_index_flood_present
  use seq_flds_mod, only: seq_flds_scalar_index_rofice_present
  use seq_flds_mod, only: seq_flds_scalar_index_precip_fact
  use seq_flds_mod, only: seq_flds_scalar_index_nextsw_cday
  use seq_flds_mod, only: seq_flds_scalar_index_atm_aero
  use seq_flds_mod, only: seq_flds_scalar_index_dead_comps
  use shr_nuopc_methods_mod, only: shr_nuopc_methods_chkErr

  implicit none

  private  ! default private

  ! !PUBLIC TYPES:

  public :: med_infodata_type

  ! !PUBLIC MEMBER FUNCTIONS

  public :: med_infodata_init1           ! Initialize before clocks are initialized
  public :: med_infodata_init2           ! Init after clocks are initialized
  public :: med_infodata_GetData         ! Get values from infodata object
  public :: med_infodata_PutData         ! Change values in infodata object
  public :: med_infodata_CopyStateToInfodata
  public :: med_infodata_CopyInfodataToState

  ! !PUBLIC DATA MEMBERS:
  public :: med_infodata                  ! instance of infodata datatype

  ! Type to hold pause/resume signaling information
  type seq_pause_resume_type
     private
     character(SHR_KIND_CL) :: atm_resume(num_inst_atm) = ' ' ! atm resume file
     character(SHR_KIND_CL) :: lnd_resume(num_inst_lnd) = ' ' ! lnd resume file
     character(SHR_KIND_CL) :: ice_resume(num_inst_ice) = ' ' ! ice resume file
     character(SHR_KIND_CL) :: ocn_resume(num_inst_ocn) = ' ' ! ocn resume file
     character(SHR_KIND_CL) :: glc_resume(num_inst_glc) = ' ' ! glc resume file
     character(SHR_KIND_CL) :: rof_resume(num_inst_rof) = ' ' ! rof resume file
     character(SHR_KIND_CL) :: wav_resume(num_inst_wav) = ' ' ! wav resume file
     character(SHR_KIND_CL) :: cpl_resume = ' '               ! cpl resume file
  end type seq_pause_resume_type

  ! InputInfo derived type
  type med_infodata_type
     private     ! This type is opaque

     !--- set via components and held fixed after initialization ---
     logical                 :: dead_comps              ! do we have dead models
     logical                 :: rofice_present          ! does rof have iceberg coupling on
     logical                 :: rof_prognostic          ! does rof component need input data
     logical                 :: flood_present           ! does rof have flooding on
     logical                 :: ocnrof_prognostic       ! does component need rof data
     logical                 :: iceberg_prognostic      ! does the ice model support icebergs
     logical                 :: glclnd_present          ! does glc have land coupling fields on
     logical                 :: glcocn_present          ! does glc have ocean runoff on
     logical                 :: glcice_present          ! does glc have iceberg coupling on
     logical                 :: glc_coupled_fluxes      ! does glc send fluxes to other components
                                                        ! (only relevant if glc_present is .true.)
     !--- set via components and may be time varying ---
     logical                 :: atm_aero    = .false.          ! atmosphere aerosols
     real(SHR_KIND_R8)       :: nextsw_cday = -1.0_SHR_KIND_R8 ! calendar of next atm shortwave
     real(SHR_KIND_R8)       :: precip_fact =  1.0_SHR_KIND_R8 ! precip factor
     type(seq_pause_resume_type), pointer :: pause_resume => NULL()

     !--- set by driver and may be time varying
     logical                 :: glc_valid_input = .true.  ! is valid accumulated data being sent to prognostic glc
  end type med_infodata_type

  type (med_infodata_type), target :: med_infodata ! single instance for cpl and all comps

  ! used/reused in module
  integer                :: dbrc
  character(len=1024)    :: msgString
  character(*),parameter :: u_FILE_u = &
    __FILE__
  !===============================================================================

CONTAINS

  !===============================================================================
  subroutine med_infodata_init1(infodata, ID)

    ! !DESCRIPTION:
    ! Initialize pause_resume

    ! !INPUT/OUTPUT PARAMETERS:
    type(med_infodata_type), intent(INOUT) :: infodata  ! infodata object
    integer(SHR_KIND_IN),    intent(IN)    :: ID        ! seq_comm ID

    !----- local -----
    character(len=*),    parameter :: subname = '(med_infodata_Init1) '
    !-------------------------------------------------------------------------------

    if (associated(infodata%pause_resume)) then
       deallocate(infodata%pause_resume)
    end if
    nullify(infodata%pause_resume)

  end subroutine med_infodata_init1

  !===============================================================================
  subroutine med_infodata_init2(infodata)

    ! !DESCRIPTION: re-initialize pause-resume depending on the time manager setup

    ! !USES:
    use seq_timemgr_mod, only : seq_timemgr_pause_active

    ! !INPUT/OUTPUT PARAMETERS:
    type(med_infodata_type), intent(INOUT) :: infodata  ! infodata object
    !----------------------------------------------------------

    !| If pause/resume is active, initialize the resume data
    if (seq_timemgr_pause_active() .and. (.not. associated(infodata%pause_resume))) then
       allocate(infodata%pause_resume)
    end if

  end subroutine med_infodata_init2

!================================================================================

  subroutine med_infodata_CopyStateToInfodata(State, infodata, type, mpicom, rc)
    ! ----------------------------------------------
    ! Copy scalar data from State to local data on root then broadcast data
    ! to all PETs in component.
    ! ----------------------------------------------
    type(ESMF_State),  intent(in)     :: State
    type(med_infodata_type),intent(inout)  :: infodata
    character(len=*),  intent(in)     :: type
    integer,           intent(in)     :: mpicom
    integer,           intent(inout)  :: rc

    ! local variables
    integer                         :: mytask, ierr, len
    character(MPI_MAX_ERROR_STRING) :: lstring
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    real(ESMF_KIND_R8), pointer     :: farrayptr(:,:)
    real(ESMF_KIND_R8)              :: data(seq_flds_scalar_num)
    logical                         :: dead_comps
    character(len=*), parameter     :: subname='(med_infodata_CopyStateToInfodata)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, mytask, rc)
    call ESMF_StateGet(State, itemName=trim(seq_flds_scalar_name), itemType=itemType, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_NOTFOUND) then
       call ESMF_LogWrite(trim(subname)//": "//trim(seq_flds_scalar_name)//" not found", ESMF_LOGMSG_INFO, &
            line=__LINE__, file=u_FILE_u, rc=dbrc)
    else
      call ESMF_StateGet(State, itemName=trim(seq_flds_scalar_name), field=field, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (size(data) < seq_flds_scalar_num .or. size(farrayptr) < seq_flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR on data size", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
        data(1:seq_flds_scalar_num) = farrayptr(1,1:seq_flds_scalar_num)
      endif

      call MPI_BCAST(data, seq_flds_scalar_num, MPI_REAL8, 0, mpicom, rc)
      if (rc /= MPI_SUCCESS) then
        call MPI_ERROR_STRING(rc,lstring,len,ierr)
        call ESMF_LogWrite(trim(subname)//": ERROR "//trim(lstring), ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      if (type == 'atm2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        if (infodata%dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) infodata%dead_comps = .true.
        infodata%atm_aero = (nint(data(seq_flds_scalar_index_atm_aero))/=0)
        infodata%nextsw_cday = data(seq_flds_scalar_index_nextsw_cday)

      elseif (type == 'ocn2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        if (infodata%dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) infodata%dead_comps = .true.
        infodata%precip_fact=data(seq_flds_scalar_index_precip_fact)

      elseif (type == 'ice2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        if (infodata%dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) infodata%dead_comps = .true.

      elseif (type == 'lnd2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        if (infodata%dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) infodata%dead_comps = .true.

      elseif (type == 'rof2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        if (infodata%dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) infodata%dead_comps = .true.
        infodata%flood_present=(nint(data(seq_flds_scalar_index_flood_present)) /= 0)
        infodata%rofice_present=(nint(data(seq_flds_scalar_index_rofice_present)) /= 0)

      elseif (type == 'wav2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        if (infodata%dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) infodata%dead_comps = .true.

      elseif (type == 'glc2cpli') then
        write(msgString,'(2i8,2l4)') nint(data(seq_flds_scalar_index_nx)),nint(data(seq_flds_scalar_index_ny))
        call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
        call med_infodata_getData(infodata,dead_comps=dead_comps)
        if (dead_comps .or. nint(data(seq_flds_scalar_index_dead_comps))/=0) dead_comps = .true.
        call med_infodata_putData(infodata, dead_comps=dead_comps)

      elseif (type == 'atm2cpl') then
         infodata%nextsw_cday=data(seq_flds_scalar_index_nextsw_cday)

      elseif (type == 'ocn2cpl') then
         infodata%precip_fact=data(seq_flds_scalar_index_precip_fact)

      elseif (type == 'ice2cpl') then
        ! nothing

      elseif (type == 'lnd2cpl') then
        ! nothing

      elseif (type == 'rof2cpl') then
        ! nothing

      elseif (type == 'wav2cpl') then
        ! nothing

      elseif (type == 'glc2cpl') then
        ! nothing

      else
        call ESMF_LogWrite(trim(subname)//": ERROR in type = "//trim(type), ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

    endif

  end subroutine med_infodata_CopyStateToInfodata

!================================================================================

  subroutine med_infodata_CopyInfodataToState(infodata, State, type, mpicom, rc)
    ! ----------------------------------------------
    ! Copy local scalar data into State, root only,
    ! but called on all PETs in component
    ! ----------------------------------------------
    type(med_infodata_type),intent(in):: infodata
    type(ESMF_State),  intent(inout)  :: State
    character(len=*),  intent(in)     :: type
    integer,           intent(in)     :: mpicom
    integer,           intent(inout)  :: rc

    ! local variables
    integer                     :: mytask
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: ItemType
    real(ESMF_KIND_R8), pointer :: farrayptr(:,:)
    logical                     :: dead_comps
    real(ESMF_KIND_R8)          :: nextsw_cday, precip_fact
    character(len=*), parameter :: subname='(med_infodata_CopyInfodataToState)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, mytask, rc)

    call ESMF_StateGet(State, itemName=trim(seq_flds_scalar_name), itemType=itemType, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

      call ESMF_LogWrite(trim(subname)//": "//trim(seq_flds_scalar_name)//" not found", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)

    else

      call ESMF_StateGet(State, itemName=trim(seq_flds_scalar_name), field=field, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (size(farrayptr) < seq_flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR on data size", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif

        nextsw_cday = infodata%nextsw_cday
        precip_fact = infodata%precip_fact
        dead_comps = infodata%dead_comps

        farrayptr(1,seq_flds_scalar_index_nextsw_cday) = nextsw_cday
        farrayptr(1,seq_flds_scalar_index_precip_fact) = precip_fact
        if (dead_comps) then
           farrayptr(1,seq_flds_scalar_index_dead_comps) = 1._ESMF_KIND_R8
        else
           farrayptr(1,seq_flds_scalar_index_dead_comps) = 0._ESMF_KIND_R8
        end if

      endif
    endif

  end subroutine med_infodata_CopyInfodataToState

  !===============================================================================
  SUBROUTINE med_infodata_GetData( infodata, &
       flux_epbal              , &
       ocnrof_prognostic       , &
       glc_coupled_fluxes      , &
       flood_present           , &
       rofice_present          , &
       glclnd_present          , &
       glcocn_present          , &
       glcice_present          , &
       glc_valid_input         , &
       iceberg_prognostic      , &
       dead_comps              , &
       flux_epbalfact          , &
       nextsw_cday             , &
       precip_fact             , &
       atm_aero                , &
       atm_resume              , &
       lnd_resume              , &
       ocn_resume              , &
       ice_resume              , &
       glc_resume              , &
       rof_resume              , &
       wav_resume              , &
       cpl_resume              )

    implicit none

    ! !DESCRIPTION:!    Get values out of the infodata object.

    ! !INPUT/OUTPUT PARAMETERS:

    type(med_infodata_type),          intent(IN)  :: infodata                ! Input CCSM structure
    character(SHR_KIND_CL), optional, intent(IN)  :: flux_epbal     ! selects E,P,R adjustment technique
    logical,                optional, intent(OUT) :: dead_comps              ! do we have dead models
    logical,                optional, intent(OUT) :: rofice_present
    logical,                optional, intent(OUT) :: flood_present
    logical,                optional, intent(OUT) :: ocnrof_prognostic
    logical,                optional, intent(OUT) :: iceberg_prognostic
    logical,                optional, intent(OUT) :: glclnd_present
    logical,                optional, intent(OUT) :: glcocn_present
    logical,                optional, intent(OUT) :: glcice_present
    logical,                optional, intent(OUT) :: glc_coupled_fluxes
    real(SHR_KIND_R8),      optional, intent(OUT) :: nextsw_cday             ! calendar of next atm shortwave
    real(SHR_KIND_R8),      optional, intent(OUT) :: precip_fact             ! precip factor
    real(SHR_KIND_R8),      optional, intent(OUT) :: flux_epbalfact          ! adjusted precip factor
    logical,                optional, intent(OUT) :: atm_aero                ! atmosphere aerosols
    logical,                optional, intent(OUT) :: glc_valid_input
    character(SHR_KIND_CL), optional, intent(OUT) :: atm_resume(:) ! atm read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: lnd_resume(:) ! lnd read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: ice_resume(:) ! ice read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: ocn_resume(:) ! ocn read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: glc_resume(:) ! glc read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: rof_resume(:) ! rof read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: wav_resume(:) ! wav read resume state
    character(SHR_KIND_CL), optional, intent(OUT) :: cpl_resume ! cpl read resume state

    !----- local -----
    character(len=*), parameter :: subname = '(med_infodata_GetData) '

    !-------------------------------------------------------------------------------

    if ( present(dead_comps)          ) dead_comps         = infodata%dead_comps
    if ( present(rofice_present)      ) rofice_present     = infodata%rofice_present
    if ( present(flood_present)       ) flood_present      = infodata%flood_present
    if ( present(ocnrof_prognostic)   ) ocnrof_prognostic  = infodata%ocnrof_prognostic
    if ( present(iceberg_prognostic)  ) iceberg_prognostic = infodata%iceberg_prognostic
    if ( present(glclnd_present)      ) glclnd_present     = infodata%glclnd_present
    if ( present(glcocn_present)      ) glcocn_present     = infodata%glcocn_present
    if ( present(glcice_present)      ) glcice_present     = infodata%glcice_present
    if ( present(glc_coupled_fluxes)  ) glc_coupled_fluxes = infodata%glc_coupled_fluxes
    if ( present(glc_valid_input)     ) glc_valid_input    = infodata%glc_valid_input
    if ( present(atm_aero)            ) atm_aero           = infodata%atm_aero
    if ( present(nextsw_cday)         ) nextsw_cday        = infodata%nextsw_cday
    if ( present(precip_fact)         ) precip_fact        = infodata%precip_fact

    if ( present(flux_epbalfact) ) then
       if (.not. present(flux_epbal)) then
          call shr_sys_abort(subname // "Must provide flux_epbal as an input argument to determine infodata%precip_fact")
       end if
       flux_epbalfact = 1.0_SHR_KIND_R8
       if (trim(flux_epbal) == 'ocn') then
          flux_epbalfact = infodata%precip_fact
       end if
       if (flux_epbalfact <= 0.0_SHR_KIND_R8) then
          if (loglevel > 0) write(logunit,'(2a,e16.6)') &
               trim(subname),' WARNING: factor from ocn = ',flux_epbalfact
          if (loglevel > 0) write(logunit,'(2a)') &
               trim(subname),' WARNING: resetting flux_epbalfact to 1.0'
          flux_epbalfact = 1.0_SHR_KIND_R8
       end if
    endif

    if ( present(atm_resume) ) then
       if (associated(infodata%pause_resume)) then
          atm_resume(:)  = infodata%pause_resume%atm_resume(:)
       else
          atm_resume(:) = ' '
       end if
    end if
    if ( present(lnd_resume) ) then
       if (associated(infodata%pause_resume)) then
          lnd_resume(:)  = infodata%pause_resume%lnd_resume(:)
       else
          lnd_resume(:) = ' '
       end if
    end if
    if ( present(ice_resume) ) then
       if (associated(infodata%pause_resume)) then
          ice_resume(:)  = infodata%pause_resume%ice_resume(:)
       else
          ice_resume(:) = ' '
       end if
    end if
    if ( present(ocn_resume) ) then
       if (associated(infodata%pause_resume)) then
          ocn_resume(:)  = infodata%pause_resume%ocn_resume(:)
       else
          ocn_resume(:) = ' '
       end if
    end if
    if ( present(glc_resume) ) then
       if (associated(infodata%pause_resume)) then
          glc_resume(:)  = infodata%pause_resume%glc_resume(:)
       else
          glc_resume(:) = ' '
       end if
    end if
    if ( present(rof_resume) ) then
       if (associated(infodata%pause_resume)) then
          rof_resume(:)  = infodata%pause_resume%rof_resume(:)
       else
          rof_resume(:) = ' '
       end if
    end if
    if ( present(wav_resume) ) then
       if (associated(infodata%pause_resume)) then
          wav_resume(:)  = infodata%pause_resume%wav_resume(:)
       else
          wav_resume(:) = ' '
       end if
    end if
    if ( present(cpl_resume) ) then
       if (associated(infodata%pause_resume)) then
          cpl_resume     = infodata%pause_resume%cpl_resume
       else
          cpl_resume = ' '
       end if
    end if

  END SUBROUTINE med_infodata_GetData

  !===============================================================================
  SUBROUTINE med_infodata_PutData( infodata , &
       ocnrof_prognostic       , &
       glc_coupled_fluxes      , &
       flood_present           , &
       rofice_present          , &
       glclnd_present          , &
       glcocn_present          , &
       glcice_present          , &
       glc_valid_input         , &
       iceberg_prognostic      , &
       dead_comps              , &
       nextsw_cday             , &
       precip_fact             , &
       atm_aero                , &
       atm_resume              , &
       lnd_resume              , &
       ocn_resume              , &
       ice_resume              , &
       glc_resume              , &
       rof_resume              , &
       wav_resume              , &
       cpl_resume              )

    implicit none

    ! !DESCRIPTION:  Put values into the infodata object.

    ! !INPUT/OUTPUT PARAMETERS:
    type(med_infodata_type),          intent(INOUT) :: infodata                ! Input CCSM structure
    logical,                optional, intent(IN) :: dead_comps     ! do we have dead models
    logical,                optional, intent(IN) :: rofice_present
    logical,                optional, intent(IN) :: flood_present
    logical,                optional, intent(IN) :: ocnrof_prognostic
    logical,                optional, intent(IN) :: iceberg_prognostic
    logical,                optional, intent(IN) :: glclnd_present
    logical,                optional, intent(IN) :: glcocn_present
    logical,                optional, intent(IN) :: glcice_present
    logical,                optional, intent(IN) :: glc_coupled_fluxes
    logical,                optional, intent(IN) :: glc_valid_input
    logical,                optional, intent(IN) :: atm_aero       ! atm aerosols
    real(SHR_KIND_R8),      optional, intent(IN) :: nextsw_cday    ! calendar of next atm shortwave
    real(SHR_KIND_R8),      optional, intent(IN) :: precip_fact    ! precip factor
    character(SHR_KIND_CL), optional, intent(IN) :: atm_resume(:)  ! atm resume
    character(SHR_KIND_CL), optional, intent(IN) :: lnd_resume(:)  ! lnd resume
    character(SHR_KIND_CL), optional, intent(IN) :: ice_resume(:)  ! ice resume
    character(SHR_KIND_CL), optional, intent(IN) :: ocn_resume(:)  ! ocn resume
    character(SHR_KIND_CL), optional, intent(IN) :: glc_resume(:)  ! glc resume
    character(SHR_KIND_CL), optional, intent(IN) :: rof_resume(:)  ! rof resume
    character(SHR_KIND_CL), optional, intent(IN) :: wav_resume(:)  ! wav resume
    character(SHR_KIND_CL), optional, intent(IN) :: cpl_resume     ! cpl resume

    !----- local -----
    character(len=*), parameter :: subname = '(med_infodata_PutData) '
    !-------------------------------------------------------------------------------

    if ( present(dead_comps)            ) infodata%dead_comps         = dead_comps
    if ( present(rofice_present)        ) infodata%rofice_present     = rofice_present
    if ( present(flood_present)         ) infodata%flood_present      = flood_present
    if ( present(ocnrof_prognostic)     ) infodata%ocnrof_prognostic  = ocnrof_prognostic
    if ( present(iceberg_prognostic)    ) infodata%iceberg_prognostic = iceberg_prognostic
    if ( present(glclnd_present)        ) infodata%glclnd_present     = glclnd_present
    if ( present(glcocn_present)        ) infodata%glcocn_present     = glcocn_present
    if ( present(glcice_present)        ) infodata%glcice_present     = glcice_present
    if ( present(glc_coupled_fluxes)    ) infodata%glc_coupled_fluxes = glc_coupled_fluxes
    if ( present(glc_valid_input)       ) infodata%glc_valid_input    = glc_valid_input
    if ( present(atm_aero)              ) infodata%atm_aero           = atm_aero
    if ( present(nextsw_cday)           ) infodata%nextsw_cday        = nextsw_cday
    if ( present(precip_fact)           ) infodata%precip_fact        = precip_fact

    if ( present(atm_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%atm_resume(:) = atm_resume(:)
       else if (ANY(len_trim(atm_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%atm_resume(:) = atm_resume(:)
       end if
    end if
    if ( present(lnd_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%lnd_resume(:) = lnd_resume(:)
       else if (ANY(len_trim(lnd_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%lnd_resume(:) = lnd_resume(:)
       end if
    end if
    if ( present(ice_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%ice_resume(:) = ice_resume(:)
       else if (ANY(len_trim(ice_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%ice_resume(:) = ice_resume(:)
       end if
    end if
    if ( present(ocn_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%ocn_resume(:) = ocn_resume(:)
       else if (ANY(len_trim(ocn_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%ocn_resume(:) = ocn_resume(:)
       end if
    end if
    if ( present(glc_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%glc_resume(:) = glc_resume(:)
       else if (ANY(len_trim(glc_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%glc_resume(:) = glc_resume(:)
       end if
    end if
    if ( present(rof_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%rof_resume(:) = rof_resume(:)
       else if (ANY(len_trim(rof_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%rof_resume(:) = rof_resume(:)
       end if
    end if
    if ( present(wav_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%wav_resume(:) = wav_resume(:)
       else if (ANY(len_trim(wav_resume(:)) > 0)) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%wav_resume(:) = wav_resume(:)
       end if
    end if
    if ( present(cpl_resume) ) then
       if (associated(infodata%pause_resume)) then
          infodata%pause_resume%cpl_resume = cpl_resume
       else if (len_trim(cpl_resume) > 0) then
          allocate(infodata%pause_resume)
          infodata%pause_resume%cpl_resume = cpl_resume
       end if
    end if

  END SUBROUTINE med_infodata_PutData

end module med_infodata_mod
