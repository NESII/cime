module ESM

  !-----------------------------------------------------------------------------
  ! Code that specializes generic ESM Component code.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Driver, &
    driver_routine_SS             => SetServices, &
    driver_label_SetModelServices => label_SetModelServices, &
    driver_label_SetRunSequence   => label_SetRunSequence
  use MED                   , only : med_SS => SetServices
  use NUOPC_Connector       , only : cpl_SS => SetServices

  use shr_sys_mod           , only : shr_sys_abort
  use shr_kind_mod          , only : SHR_KIND_R8, SHR_KIND_CS, SHR_KIND_CL
  use shr_kind_mod          , only : CL => SHR_KIND_CL
  use shr_scam_mod          , only : shr_scam_checkSurface
  use shr_mpi_mod           , only : shr_mpi_bcast, shr_mpi_chkerr
  use shr_mem_mod           , only : shr_mem_init, shr_mem_getusage
  use shr_cal_mod           , only : shr_cal_date2ymd
  use shr_orb_mod           , only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  use shr_frz_mod           , only : shr_frz_freezetemp_init
  use shr_reprosum_mod      , only : shr_reprosum_setopts
  use shr_pio_mod           , only : shr_pio_init1, shr_pio_init2
  use shr_const_mod         , only : shr_const_tkfrz, shr_const_tktrip
  use shr_const_mod         , only : shr_const_mwwv, shr_const_mwdair
  use shr_wv_sat_mod        , only : shr_wv_sat_set_default, shr_wv_sat_init
  use shr_wv_sat_mod        , only : ShrWVSatTableSpec, shr_wv_sat_make_tables
  use shr_wv_sat_mod        , only : shr_wv_sat_get_scheme_idx, shr_wv_sat_valid_idx
  use shr_file_mod          , only : shr_file_getUnit, shr_file_freeUnit
  use shr_file_mod          , only : shr_file_setloglevel, shr_file_setlogunit, shr_file_setio
  use shr_assert_mod        , only : shr_assert_in_domain

  use seq_comm_mct          , only : CPLID, GLOID, logunit, loglevel
  use seq_comm_mct          , only : ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
  use seq_comm_mct          , only : num_inst_atm, num_inst_lnd, num_inst_rof
  use seq_comm_mct          , only : num_inst_ocn, num_inst_ice, num_inst_glc
  use seq_comm_mct          , only : num_inst_wav, num_inst_esp, num_inst_total
  use seq_comm_mct          , only : seq_comm_init, seq_comm_setnthreads, seq_comm_getnthreads
  use seq_comm_mct          , only : seq_comm_getinfo => seq_comm_setptrs, seq_comm_petlist
  use seq_comm_mct          , only : seq_comm_iamin, seq_comm_name, seq_comm_namelen, seq_comm_iamroot
  use seq_timemgr_mod       , only : seq_timemgr_clockInit, seq_timemgr_EClockGetData
  use seq_flds_mod          , only : seq_flds_set

  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_setDict_fromseqflds
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr

  use med_infodata_mod      , only : med_infodata_init1, med_infodata_init2, med_infodata
  use pio                   , only : file_desc_t, pio_closefile, pio_file_is_open
  use t_drv_timers_mod
  use perf_mod

#ifdef ESMFUSE_cam
  use  cam_comp_nuopc, only:   cam_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_pop2
  use pop2_comp_nuopc, only:   pop2_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_cice
  use cice_comp_nuopc, only:   cice_SS => SetServices
#endif
#ifdef ESMFUSE_clm
  use  clm_comp_nuopc, only:   clm_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_rtm
  use  rtm_comp_nuopc, only:   rtm_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_mosart
  use mosart_comp_nuopc, only: mosart_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_mom6
  use mom6_comp_nuopc, only:   mom6_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_ww3
  use  ww3_comp_nuopc, only:   ww3_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_cism
  use cism_comp_nuopc, only:   cism_SS => SetServices
#endif

#ifdef ESMFUSE_datm
  use datm_comp_nuopc, only: datm_SS => SetServices
#endif
#ifdef ESMFUSE_docn
  use docn_comp_nuopc, only: docn_SS => SetServices
#endif
#ifdef ESMFUSE_dice
  use dice_comp_nuopc, only: dice_SS => SetServices
#endif
#ifdef ESMFUSE_dlnd
  use dlnd_comp_nuopc, only: dlnd_SS => SetServices
#endif
#ifdef ESMFUSE_drof
  use drof_comp_nuopc, only: drof_SS => SetServices
#endif
#ifdef ESMFUSE_dwav
  use dwav_comp_nuopc, only: dwav_SS => SetServices
#endif
#ifdef ESMFUSE_dglc
  use dglc_comp_nuopc, only: dglc_SS => SetServices
#endif
#ifdef ESMFUSE_desp
  use desp_comp_nuopc, only: desp_SS => SetServices
#endif

#ifdef ESMFUSE_xatm
  use xatm_comp_nuopc, only: xatm_SS => SetServices
#endif
#ifdef ESMFUSE_xocn
  use xocn_comp_nuopc, only: xocn_SS => SetServices
#endif
#ifdef ESMFUSE_xice
  use xice_comp_nuopc, only: xice_SS => SetServices
#endif
#ifdef ESMFUSE_xlnd
  use xlnd_comp_nuopc, only: xlnd_SS => SetServices
#endif
#ifdef ESMFUSE_xrof
  use xrof_comp_nuopc, only: xrof_SS => SetServices
#endif
#ifdef ESMFUSE_xwav
  use xwav_comp_nuopc, only: xwav_SS => SetServices
#endif
#ifdef ESMFUSE_xglc
  use xglc_comp_nuopc, only: xglc_SS => SetServices
#endif
#ifdef ESMFUSE_NOTYET_xesp
  use xesp_comp_nuopc, only: xesp_SS => SetServices
#endif

  implicit none

  include 'mpif.h'
  integer, parameter      :: dbug_flag = 10
  character(len=512)      :: msgstr
  integer                 :: dbrc
  logical                 :: mastertask
  character(*), parameter :: runseq_filename = "cesm.runconfig"
  character(*), parameter :: u_FILE_u = __FILE__

  type(ESMF_Clock), target :: EClock_d
  type(ESMF_Clock), target :: EClock_a
  type(ESMF_Clock), target :: EClock_l
  type(ESMF_Clock), target :: EClock_o
  type(ESMF_Clock), target :: EClock_i
  type(ESMF_Clock), target :: EClock_g
  type(ESMF_Clock), target :: EClock_r
  type(ESMF_Clock), target :: EClock_w
  type(ESMF_Clock), target :: EClock_e

  private

  public SetServices

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    type(ESMF_Config)           :: config
    type(NUOPC_FreeFormat)      :: attrFF
    character(len=*), parameter :: subname = "(esm.F90:SetServices)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! NUOPC_Driver registers the generic methods
    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
      specRoutine=SetModelServices, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetRunSequence, &
      specRoutine=SetRunSequence, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! register an internal initialization method
    call NUOPC_CompSetInternalEntryPoint(driver, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p2"/), userRoutine=ModifyCplLists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create, open and set the config
    config = ESMF_ConfigCreate(rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ConfigLoadFile(config, runseq_filename, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_GridCompSet(driver, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    integer                     :: localrc
    type(ESMF_Time)             :: startTime
    type(ESMF_Time)             :: stopTime
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_Clock)            :: internalClock
    type(ESMF_GridComp)         :: child
#if (1 == 0)
    integer                     :: petCount, i
    integer, allocatable        :: petList(:)
#endif
    integer                     :: n, stat
    integer                     :: GlobalComm
    integer, pointer            :: petList(:)
    type(ESMF_Config)           :: config
    character(len=32), allocatable  :: compLabels(:)
    character(len=10)           :: value
    character(len=20)           :: model, prefix
    character(len=8)            :: atm_present, lnd_present, ocn_present
    character(len=8)            :: ice_present, rof_present, wav_present
    character(len=8)            :: glc_present, med_present
    integer                     :: componentCount
    type(NUOPC_FreeFormat)      :: attrFF
    type(ESMF_VM)               :: vm
    integer                     :: localPet
    character(len=*), parameter :: subname = "(esm.F90:SetModelServices)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !-------------------------------------------
    ! Read components from config file
    !-------------------------------------------

    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (localPet == 0) then
       mastertask=.true.
    else
       mastertask = .false.
    end if

    ! determine the generic component labels
    componentCount = ESMF_ConfigGetLen(config,label="CESM_component_list:", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (componentCount == 0) then
      write (msgstr, *) "No models were specified in CESM_component_list "
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, &
        file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif
    allocate(compLabels(componentCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
          line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return
    call ESMF_ConfigGetAttribute(config, valueList=compLabels, label="CESM_component_list:", &
         count=componentCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       call NUOPC_FreeFormatPrint(attrFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_info_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       call NUOPC_FreeFormatPrint(attrFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_maps_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       call NUOPC_FreeFormatPrint(attrFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_auxhist_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       call NUOPC_FreeFormatPrint(attrFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_clock_attributes::", relaxedflag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mastertask) then
       call NUOPC_FreeFormatPrint(attrFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Initialize mct and pets and cesm stuff
    !-------------------------------------------

    call esm_SetAttributes_and_InitClocks(driver, &
         Eclock_d, Eclock_a, Eclock_l, Eclock_o, &
         Eclock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e)

    call shr_nuopc_fldList_setDict_fromseqflds(rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------
    ! check component list for active components and set present flags
    !--------

    med_present = "false"
    atm_present = "false"
    lnd_present = "false"
    ocn_present = "false"
    ice_present = "false"
    rof_present = "false"
    wav_present = "false"
    glc_present = "false"

    do n=1, componentCount

      !--- construct component prefix
      prefix=trim(compLabels(n))

      !--- read in model instance name
      call ESMF_ConfigGetAttribute(config, model, label=trim(prefix)//"_model:", default="none", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--- check that there was a model instance specified
      if (trim(model) == "none") then
        ! Error condition: no model was specified
        write (msgstr, *) "No model was specified for component: ",trim(prefix)
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, &
          file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif

      if (trim(prefix) == "MED") then
        med_present = "true"
      elseif (trim(prefix) == "ATM") then
        atm_present = "true"
      elseif (trim(prefix) == "LND") then
        lnd_present = "true"
      elseif (trim(prefix) == "OCN") then
        ocn_present = "true"
      elseif (trim(prefix) == "ICE") then
        ice_present = "true"
      elseif (trim(prefix) == "ROF") then
        rof_present = "true"
      elseif (trim(prefix) == "WAV") then
        wav_present = "true"
      elseif (trim(prefix) == "GLC") then
        glc_present = "true"
      endif

    enddo

    call ESMF_LogWrite(trim(subname)//":atm_present="//trim(atm_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":lnd_present="//trim(lnd_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":ocn_present="//trim(ocn_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":ice_present="//trim(ice_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":rof_present="//trim(rof_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":wav_present="//trim(wav_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":glc_present="//trim(glc_present), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_LogWrite(trim(subname)//":med_present="//trim(med_present), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeAdd(driver, attrList=(/'atm_present','lnd_present','ocn_present', &
      'ice_present','rof_present','wav_present','glc_present','med_present'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call NUOPC_CompAttributeSet(driver, name="atm_present", value=atm_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="lnd_present", value=lnd_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="ocn_present", value=ocn_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="ice_present", value=ice_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="rof_present", value=rof_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="wav_present", value=wav_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="glc_present", value=glc_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(driver, name="med_present", value=med_present, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------
    ! determine information for each component and add to the driver
    !--------

    do n=1, componentCount

      !--- construct component prefix
      prefix=trim(compLabels(n))

      !--- read in model instance name
      call ESMF_ConfigGetAttribute(config, model, label=trim(prefix)//"_model:", default="none", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--- check that there was a model instance specified
      if (trim(model) == "none") then
        ! Error condition: no model was specified
        write (msgstr, *) "No model was specified for component: ",trim(prefix)
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, &
          file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif

#if (1 == 0)
      ! read in petList bounds
      call ESMF_ConfigGetAttribute(config, petListBounds, &
        label=trim(prefix)//"_petlist_bounds:", default=-1, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      ! handle the default situation
      if (petListBounds(1)==-1 .or. petListBounds(2)==-1) then
        petListBounds(1) = 0
        petListBounds(2) = petCount - 1
      endif
      ! set petList for this component
      allocate(petList(petListBounds(2)-petListBounds(1)+1))
      do j=petListBounds(1), petListBounds(2)
        petList(j-petListBounds(1)+1) = j ! PETs are 0 based
      enddo
#endif

      !--------
      ! ATM
      !--------

      if (trim(prefix) == "ATM") then

        if (atm_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' atm_present inconsistent = '//trim(prefix)//':'//trim(atm_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(ATMID(1),petList)
        if (trim(model) == "datm") then
#ifdef ESMFUSE_datm
          call NUOPC_DriverAddComp(driver, "ATM", datm_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xatm") then
#ifdef ESMFUSE_xatm
          call NUOPC_DriverAddComp(driver, "ATM", xatm_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "cam") then
#ifdef ESMFUSE_cam
          call NUOPC_DriverAddComp(driver, "ATM",  cam_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call esm_AddAttributes(child, driver, ATMID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read ATM attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! OCN
      !--------

      elseif (trim(prefix) == "OCN") then

        if (ocn_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' ocn_present inconsistent = '//trim(prefix)//':'//trim(ocn_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(OCNID(1),petList)
        if (trim(model) == "docn") then
#ifdef ESMFUSE_docn
          call NUOPC_DriverAddComp(driver, "OCN", docn_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xocn") then
#ifdef ESMFUSE_xocn
          call NUOPC_DriverAddComp(driver, "OCN", xocn_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "pop2") then
#ifdef ESMFUSE_NOTYET_pop2
          call NUOPC_DriverAddComp(driver, "OCN", pop2_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "mom6") then
#ifdef ESMFUSE_NOTYET_mom6
          call NUOPC_DriverAddComp(driver, "OCN", mom6_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call esm_AddAttributes(child, driver, OCNID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read OCN attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! ICE
      !--------

      elseif (trim(prefix) == "ICE") then

        if (ice_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' ice_present inconsistent = '//trim(prefix)//':'//trim(ice_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(ICEID(1),petList)
        if (trim(model) == "dice") then
#ifdef ESMFUSE_dice
          call NUOPC_DriverAddComp(driver, "ICE", dice_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "cice") then
#ifdef ESMFUSE_NOTYET_cice
          call NUOPC_DriverAddComp(driver, "ICE", cice_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xice") then
#ifdef ESMFUSE_xice
          call NUOPC_DriverAddComp(driver, "ICE", xice_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, ICEID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read ICE attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! LND
      !--------

      elseif (trim(prefix) == "LND") then

        if (lnd_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' lnd_present inconsistent = '//trim(prefix)//':'//trim(lnd_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(LNDID(1),petList)
        if (trim(model) == "dlnd") then
#ifdef ESMFUSE_dlnd
          call NUOPC_DriverAddComp(driver, "LND", dlnd_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "clm") then
#ifdef ESMFUSE_clm
          call NUOPC_DriverAddComp(driver, "LND", clm_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xlnd") then
#ifdef ESMFUSE_xlnd
          call NUOPC_DriverAddComp(driver, "LND", xlnd_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, LNDID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read LND attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! WAV
      !--------

      elseif (trim(prefix) == "WAV") then

        call seq_comm_petlist(WAVID(1),petList)
        if (trim(model) == "dwav") then
#ifdef ESMFUSE_dwav
          call NUOPC_DriverAddComp(driver, "WAV", dwav_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "ww") then
#ifdef ESMFUSE_NOTYET_ww3
          call NUOPC_DriverAddComp(driver, "WAV", ww3_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xwav") then
#ifdef ESMFUSE_xwav
          call NUOPC_DriverAddComp(driver, "WAV", xwav_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, WAVID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read WAV attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! GLC
      !--------

      elseif (trim(prefix) == "GLC") then

        call seq_comm_petlist(GLCID(1),petList)
        if (trim(model) == "cism") then
#ifdef ESMFUSE_NOTYET_cism
          call NUOPC_DriverAddComp(driver, "GLC", cism_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xglc") then
#ifdef ESMFUSE_xglc
          call NUOPC_DriverAddComp(driver, "GLC", xglc_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, GLCID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read GLC attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! ROF
      !--------

      elseif (trim(prefix) == "ROF") then

        if (rof_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' rof_present inconsistent = '//trim(prefix)//':'//trim(rof_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(ROFID(1),petList)
        if (trim(model) == "drof") then
#ifdef ESMFUSE_drof
          call NUOPC_DriverAddComp(driver, "ROF", drof_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "rtm") then
#ifdef ESMFUSE_NOTYET_rtm
          call NUOPC_DriverAddComp(driver, "ROF", rtm_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "mosart") then
#ifdef ESMFUSE_NOTYET_mosart
          call NUOPC_DriverAddComp(driver, "ROF", mosart_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xrof") then
#ifdef ESMFUSE_xrof
          call NUOPC_DriverAddComp(driver, "ROF", xrof_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, ROFID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read ROF attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! WAV
      !--------

      elseif (trim(prefix) == "WAV") then

        if (wav_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' wav_present inconsistent = '//trim(prefix)//':'//trim(wav_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(WAVID(1),petList)
        if (trim(model) == "dwav") then
#ifdef ESMFUSE_dwav
          call NUOPC_DriverAddComp(driver, "WAV", dwav_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "ww3") then
#ifdef ESMFUSE_NOTYET_ww3
          call NUOPC_DriverAddComp(driver, "WAV", ww3_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xwav") then
#ifdef ESMFUSE_xwav
          call NUOPC_DriverAddComp(driver, "WAV", xwav_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, WAVID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read WAV attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! GLC
      !--------

      elseif (trim(prefix) == "GLC") then

        if (glc_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' glc_present inconsistent = '//trim(prefix)//':'//trim(glc_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(GLCID(1),petList)
        if (trim(model) == "dglc") then
#ifdef ESMFUSE_dglc
          call NUOPC_DriverAddComp(driver, "GLC", dglc_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "cism") then
#ifdef ESMFUSE_NOTYET_cism
          call NUOPC_DriverAddComp(driver, "GLC", cism_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        elseif (trim(model) == "xglc") then
#ifdef ESMFUSE_xglc
          call NUOPC_DriverAddComp(driver, "GLC", xglc_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' model unavailable = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
#endif
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, GLCID(1), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read GLC attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", &
          relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! MED
      !--------

      elseif (trim(prefix) == "MED") then

        if (med_present /= "true") then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' med_present inconsistent = '//trim(prefix)//':'//trim(med_present), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call seq_comm_petlist(CPLID,petList)
        if (trim(model) == "cesm") then
          call NUOPC_DriverAddComp(driver, "MED", med_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg=subname//' invalid model = '//trim(prefix)//':'//trim(model), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call esm_AddAttributes(child, driver, CPLID, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read MED attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label=trim(prefix)//"_Attributes::", relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! for now read DRIVER_info_attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_info_attributes::", relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! for now read DRIVER_maps_attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_maps_attributes::", relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! for now read DRIVER_auxhist_attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_auxhist_attributes::", relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! for now read DRIVER_auxhist_attributes from config file into FreeFormat
        attrFF = NUOPC_FreeFormatCreate(config, label="DRIVER_clock_attributes::", relaxedflag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(child, attrFF, addFlag=.true., rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      else

        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
          msg=subname//' invalid component = '//trim(prefix), &
          line=__LINE__, file=u_FILE_u, rcToReturn=rc)
        return  ! bail out

      endif

    enddo

    !--------
    ! clean-up
    !--------
    deallocate(compLabels)

    !--------
    ! Set baseline clock
    !--------

    call ESMF_GridCompSet(driver, clock=Eclock_o, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_Clock_TimePrint(Eclock_o,subname//'EClock_base',rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetModelServices

  !-----------------------------------------------------------------------------

  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    integer                       :: localrc
    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Clock)              :: internalClock
    type(ESMF_Config)             :: config
    type(NUOPC_FreeFormat)        :: runSeqFF
    character(len=*), parameter :: subname = "(esm.F90:SetRunSequence)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !--------
    ! Run Sequence and Connectors
    !--------

    ! read free format run sequence from config
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mastertask) then
       call NUOPC_FreeFormatPrint(runSeqFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, autoAddConnectors=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_FreeFormatDestroy(runSeqFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------
    ! Update Clocks
    !--------

    call NUOPC_DriverSetRunSequence(driver, slot=1, clock=EClock_o, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_DriverSetRunSequence(driver, slot=2, clock=EClock_a, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_Clock_TimePrint(Eclock_o,subname//'EClock_o',rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_Clock_TimePrint(Eclock_a,subname//'EClock_a',rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostic
    call NUOPC_DriverPrint(driver, orderflag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetRunSequence

  !-----------------------------------------------------------------------------

  recursive subroutine ModifyCplLists(driver, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_CplComp), pointer     :: connectorList(:)
    integer                         :: i, j, cplListSize
    character(len=160), allocatable :: cplList(:)
    character(len=160)              :: tempString
    character(len=*), parameter :: subname = "(esm.F90:ModifyCplLists)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgstr,*) "Found ", size(connectorList), " Connectors."// &
      " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msgstr), ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do i=1, size(connectorList)
      ! query the cplList for connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
        itemCount=cplListSize, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (cplListSize>0) then
        allocate(cplList(cplListSize))
        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
          valueList=cplList, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        ! go through all of the entries in the cplList and add options
        do j=1, cplListSize
           !tempString = trim(cplList(j))//":REMAPMETHOD=bilinear"//&
           !":SrcTermProcessing=1:DUMPWEIGHTS=true:TermOrder=SrcSeq"
          tempString = trim(cplList(j))//":remapmethod=redist"
          cplList(j) = trim(tempString)
        enddo
        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), &
          name="CplList", valueList=cplList, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        deallocate(cplList)
      endif
    enddo

    deallocate(connectorList)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine ModifyCplLists

  !-----------------------------------------------------------------------------

  subroutine esm_SetAttributes_and_InitClocks(driver, &
       Eclock_d, Eclock_a, Eclock_l, Eclock_o, &
       Eclock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e)

    ! USES:

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
    character(len=*) , parameter    :: start_type_start     = "startup"
    character(len=*) , parameter    :: start_type_cont      = "continue"
    character(len=*) , parameter    :: start_type_brnch     = "branch"
    character(len=*) , parameter    :: orb_fixed_year       = 'fixed_year'
    character(len=*) , parameter    :: orb_variable_year    = 'variable_year'
    character(len=*) , parameter    :: orb_fixed_parameters = 'fixed_parameters'
    character(len=*) , parameter    :: subname = '(esm_SetAttributes_and_InitClocks)'

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

    call NUOPC_CompAttributeAdd(driver, &
         attrList=(/'samegrid_ao', 'samegrid_al', 'samegrid_lr', &
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

    call esm_CheckAttributes(driver)

    !----------------------------------------------------------
    ! Initialize cime_model
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="cime_model", value=cime_model, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    if ( trim(cime_model) /= 'cesm') then
       call shr_sys_abort( subname//': cime_model must be set to cesm, aborting')
    end if

    !----------------------------------------------------------
    ! Initialize coupled fields
    !----------------------------------------------------------

    call seq_flds_set(nlfilename, GLOID, cime_model)

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
       ! in med_infodata calls, reset here and the copied back into med_infodata - this is no longer the case
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

  end subroutine esm_SetAttributes_and_InitClocks

  !===============================================================================

  subroutine esm_CheckAttributes( driver )

    ! !DESCRIPTION: Check that input driver config values have reasonable values

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

  end subroutine esm_CheckAttributes

  !===============================================================================

  subroutine esm_AddAttributes(gcomp, driver, MCTID, rc)

    ! Add specific set of attributes to gcomp from driver attributes

    ! input/output parameters
    type(ESMF_GridComp),intent(inout) :: gcomp
    type(ESMF_GridComp),intent(in)    :: driver
    integer            ,intent(in)    :: MCTID
    integer            ,intent(inout) :: rc

    ! locals
    character(len=CL)  :: cvalue
    integer            :: n
    integer, parameter :: nattrlist = 28
    character(len=*), parameter :: attrList(nattrlist) = &
      (/ "case_name"    ,"single_column","scmlat"        ,"scmlon"               , &
         "read_restart" ,"start_type"   ,"tfreeze_option","model_version"        , &
         "orb_eccen"    ,"orb_obliqr"   ,"orb_lambm0"    ,"orb_mvelpp"           , &
         "info_debug"   ,"atm_aero"     ,"aqua_planet"   ,"brnch_retain_casename", &
         "perpetual"    ,"perpetual_ymd","hostname"      ,"username"             , &
         "atm_present"  ,"lnd_present"  ,"ocn_present"   ,"ice_present"          , &
         "rof_present"  ,"wav_present"  ,"glc_present"   ,"med_present"          /)
    character(len=*), parameter :: subname = "(esm.F90:esm_AddAttributes)"
    !-------------------------------------------

    rc = ESMF_Success

    ! First add MCTID to gcomp attributes
    write(cvalue,*) MCTID
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'MCTID'/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name='MCTID', value=trim(cvalue), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    ! Now add all the other attributes in AttrList (which have already been added to driver attributes)
    call NUOPC_CompAttributeAdd(gcomp, attrList=attrList, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,nattrlist
       call NUOPC_CompAttributeGet(driver, name=trim(attrList(n)), value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name=trim(attrList(n)), value=trim(cvalue), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine esm_AddAttributes

end module ESM
