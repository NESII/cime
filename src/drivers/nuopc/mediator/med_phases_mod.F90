module med_phases_mod

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  ! This mediator operates on two timescales and keeps two internal Clocks to
  ! do so.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_clean
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_Regrid
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_accum
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_average
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_copy
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use med_internalstate_mod , only : InternalState
  use med_internalstate_mod , only : ncomps, compname
  use med_internalstate_mod , only : compmed, compatm, complnd, compocn
  use med_internalstate_mod , only : compice, comprof, compwav, compglc
  use med_internalstate_mod , only : fldsFr, fldsAtmOcn
  use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mapfcopy
  use med_constants_mod     , only : med_constants_dbug_flag
  use med_constants_mod     , only : med_constants_statewrite_flag
  use med_constants_mod     , only : med_constants_spval_init
  use med_constants_mod     , only : med_constants_spval
  use med_constants_mod     , only : med_constants_czero
  use med_constants_mod     , only : med_constants_ispval_mask
  use med_merge_mod         , only : med_merge_auto
  use med_atmocn_mod        , only : med_atmocn_init, med_atmocn_ocnalb, med_atmocn_flux

  implicit none

  private

  integer           , parameter :: dbug_flag       = med_constants_dbug_flag
  logical           , parameter :: statewrite_flag = med_constants_statewrite_flag
  real(ESMF_KIND_R8), parameter :: spval_init      = med_constants_spval_init
  real(ESMF_KIND_R8), parameter :: spval           = med_constants_spval
  real(ESMF_KIND_R8), parameter :: czero           = med_constants_czero
  integer           , parameter :: ispval_mask     = med_constants_ispval_mask
  character(len=*)  , parameter :: ice_fraction_name = 'Si_ifrac'
  integer                :: dbrc
  type(ESMF_FieldBundle) :: FBtmp1,FBtmp2
  character(*),parameter :: u_FILE_u = &
    __FILE__

  public med_phases_prep_atm
  public med_phases_prep_ocn
  public med_phases_prep_ice
  public med_phases_prep_lnd
  public med_phases_prep_rof
  public med_phases_prep_wav
  public med_phases_prep_glc
  public med_phases_prep_map_atm2ocn
  public med_phases_accum_fast
  public med_phases_atmocn_init
  public med_phases_atmocn_ocnalb
  public med_phases_atmocn_flux

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_atm(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! This Mediator phase runs before ATM and ICE are being called and
    ! prepares the ATM and ICE import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    real(ESMF_KIND_R8), pointer :: ifrac_i(:)                 ! ice fraction on ice grid
    real(ESMF_KIND_R8), pointer :: ifrac_af(:), ifrac_afr(:)  ! ice fraction on atm grid consf map
    real(ESMF_KIND_R8), pointer :: ifrac_ad(:), ifrac_adr(:)  ! ice fraction on atm grid consd map
    real(ESMF_KIND_R8), pointer :: ifrac_ab(:), ifrac_abr(:)  ! ice fraction on atm grid bilnr map
    real(ESMF_KIND_R8), pointer :: ifrac_ap(:), ifrac_apr(:)  ! ice fraction on atm grid patch map
    real(ESMF_KIND_R8), pointer :: ocnwgt(:),icewgt(:),customwgt(:)
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_atm)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    call ESMF_TimeGet(time,timestring=timestr)
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    !tcraig, turn this off for ice2atm, use ice frac weighted mapping below
    do n1 = 1,ncomps
      n2 = compatm
      if (n1/=compice .and. is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
           call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), &
                trim(subname)//' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    ! map FBAtmOcn_o to FBAtmOcn_a
    if (is_local%wrap%med_coupling_active(compocn,compatm)) then
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBAtmOcn_a, value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_Regrid(fldsAtmOcn, is_local%wrap%FBAtmOcn_o, is_local%wrap%FBAtmOcn_a, &
         consfmap=is_local%wrap%RH(compocn,compatm,mapconsf), &
         consdmap=is_local%wrap%RH(compocn,compatm,mapconsd), &
         bilnrmap=is_local%wrap%RH(compocn,compatm,mapbilnr), &
         patchmap=is_local%wrap%RH(compocn,compatm,mappatch), &
         string='o2aatmocn', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 1) then
        call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBAtmOcn_a            , trim(subname)//' FBAtmOcn_a ', rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    endif

    !---------------------------------------
    !--- map ice to atm with frac weighting
    !---------------------------------------

    if (is_local%wrap%med_coupling_active(compice,compatm)) then
      call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compice,compatm), value=czero, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBImp(compice,compice), trim(ice_fraction_name), rc=rc)) then

        !--- tcraig, need to weight the ice2atm regrid by the ice fraction
        !--- need to compute weight by the frac mapped with the correct mapping
        !--- first compute the ice fraction on the atm grid for all active mappings

        call med_phases_map_frac(is_local%wrap%FBImp(compice,compice), trim(ice_fraction_name), &
                                 is_local%wrap%FBImp(compice,compatm), trim(ice_fraction_name), &
                                 is_local%wrap%RH(compice,compatm,mapconsf), ifrac_af, ifrac_afr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call med_phases_map_frac(is_local%wrap%FBImp(compice,compice), trim(ice_fraction_name), &
                                 is_local%wrap%FBImp(compice,compatm), trim(ice_fraction_name), &
                                 is_local%wrap%RH(compice,compatm,mapconsd), ifrac_ad, ifrac_adr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call med_phases_map_frac(is_local%wrap%FBImp(compice,compice), trim(ice_fraction_name), &
                                 is_local%wrap%FBImp(compice,compatm), trim(ice_fraction_name), &
                                 is_local%wrap%RH(compice,compatm,mapbilnr), ifrac_ab, ifrac_abr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call med_phases_map_frac(is_local%wrap%FBImp(compice,compice), trim(ice_fraction_name), &
                                 is_local%wrap%FBImp(compice,compatm), trim(ice_fraction_name), &
                                 is_local%wrap%RH(compice,compatm,mappatch), ifrac_ap, ifrac_apr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        !--- multiply FBImp(compice,compice) by ice_fraction

        call shr_nuopc_methods_FB_init(FBtmp1, fbgeom=is_local%wrap%FBImp(compice,compice), &
          fbflds=is_local%wrap%FBImp(compice,compice), name='FBtmp1_ice_x_frac', rc=rc)
        if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
        call shr_nuopc_methods_FB_reset(FBtmp1, value=czero, rc=rc)
        if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

        call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compice), trim(ice_fraction_name), fldptr1=dataPtr1, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        do n = 1,fldsFr(compice)%num
          if (fldsFr(compice)%shortname(n) /= trim(ice_fraction_name) .and. &
              shr_nuopc_methods_FB_FldChk(is_local%wrap%FBImp(compice,compice), fldsFr(compice)%shortname(n), rc=rc) .and. &
              shr_nuopc_methods_FB_FldChk(FBtmp1                              , fldsFr(compice)%shortname(n), rc=rc)) then
            call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compice) , fldsFr(compice)%shortname(n), dataPtr3, rc=rc)
            call shr_nuopc_methods_FB_GetFldPtr(FBtmp1                               , fldsFr(compice)%shortname(n), dataPtr4, rc=rc)
            ! avoid non array fields like the scalars
            if (lbound(dataptr1,1) == lbound(dataptr3,1) .and. ubound(dataptr1,1) == ubound(dataptr3,1)) then
              do i=lbound(dataptr3,1),ubound(dataptr3,1)
                dataPtr4(i) = dataPtr3(i) * dataPtr1(i)
              enddo
            endif
          endif
        enddo

        !--- regrid FBtmp1 = FBImp(compice,compice)*frac, fields with fraction multiplied

        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compice,compatm), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(compice), FBtmp1, is_local%wrap%FBImp(compice,compatm), &
           consfmap=is_local%wrap%RH(compice,compatm,mapconsf), &
           consdmap=is_local%wrap%RH(compice,compatm,mapconsd), &
           bilnrmap=is_local%wrap%RH(compice,compatm,mapbilnr), &
           patchmap=is_local%wrap%RH(compice,compatm,mappatch), &
           string=trim(compname(compice))//'2'//trim(compname(compatm)), rc=rc)
          !string='i2a', rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_clean(FBtmp1, rc=rc)
        if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

        !--- divide FBImp(compice,compatm) by ice_fraction, interpolated ice fraction

        do n = 1,fldsFr(compice)%num
          if (fldsFr(compice)%shortname(n) /= trim(ice_fraction_name) .and. &
              shr_nuopc_methods_FB_FldChk(is_local%wrap%FBImp(compice,compatm), fldsFr(compice)%shortname(n), rc=rc)) then
            call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compatm), fldsFr(compice)%shortname(n), dataPtr3, rc=rc)
            ! avoid non array fields like the scalars
            if (lbound(dataptr1,1) == lbound(dataptr3,1) .and. ubound(dataptr1,1) == ubound(dataptr3,1)) then
              if (fldsFr(compice)%mapping(n) == "conservefrac") then
                do i=lbound(dataptr3,1),ubound(dataptr3,1)
                  dataPtr3(i) = dataPtr3(i) * ifrac_afr(i)
                enddo
              elseif (fldsFr(compice)%mapping(n) == "conservedst") then
                do i=lbound(dataptr3,1),ubound(dataptr3,1)
                  dataPtr3(i) = dataPtr3(i) * ifrac_adr(i)
                enddo
              elseif (fldsFr(compice)%mapping(n) == 'bilinear') then
                do i=lbound(dataptr3,1),ubound(dataptr3,1)
                  dataPtr3(i) = dataPtr3(i) * ifrac_abr(i)
                enddo
              elseif (fldsFr(compice)%mapping(n) == 'patch') then
                do i=lbound(dataptr3,1),ubound(dataptr3,1)
                  dataPtr3(i) = dataPtr3(i) * ifrac_apr(i)
                enddo
              else
                call ESMF_LogWrite(trim(subname)//": mapping name error "//trim(fldsFr(compice)%mapping(n)), ESMF_LOGMSG_INFO, rc=rc)
                rc=ESMF_FAILURE
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
              endif
            endif
          endif
        enddo

        !--- make sure ifrac_a in the mapped bundle is correct
        !--- this is handled by ice_fraction_name check in the mult/divide phases to avoid ice_fraction weighting
!        call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compatm), trim(ice_fraction_name), dataPtr3, rc=rc)
!        do i=lbound(dataptr3,1),ubound(dataptr3,1)
!          dataPtr3(i) = ifrac_af(i)
!        enddo

        if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)) &
          deallocate(ifrac_af, ifrac_afr)
        if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapconsd), rc=rc)) &
          deallocate(ifrac_ad, ifrac_adr)
        if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapbilnr), rc=rc)) &
          deallocate(ifrac_ab, ifrac_abr)
        if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mappatch), rc=rc)) &
          deallocate(ifrac_ap, ifrac_apr)

      else
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compice,compatm), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(compice), is_local%wrap%FBImp(compice,compice), is_local%wrap%FBImp(compice,compatm), &
           consfmap=is_local%wrap%RH(compice,compatm,mapconsf), &
           consdmap=is_local%wrap%RH(compice,compatm,mapconsd), &
           bilnrmap=is_local%wrap%RH(compice,compatm,mapbilnr), &
           patchmap=is_local%wrap%RH(compice,compatm,mappatch), &
           string=trim(compname(compice))//'2'//trim(compname(compatm)), rc=rc)
          !string='i2a', rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      if (dbug_flag > 1) then
        call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(compice,compatm), trim(subname)//' FBImp(compice,compatm) ', rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    endif

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compatm), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(compatm), &
                    is_local%wrap%FBImp(compocn,compatm)   , is_local%wrap%FBfrac(compatm), 'ofrac', &
                    is_local%wrap%FBAtmOcn_a               , is_local%wrap%FBfrac(compatm), 'ofrac', &
                    is_local%wrap%FBImp(compice,compatm)   , is_local%wrap%FBfrac(compatm), 'ifrac', &
                    is_local%wrap%FBImp(complnd,compatm)   , is_local%wrap%FBfrac(compatm), 'lfrac', &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

#if (1 == 0)
    !---  ocn and ice fraction for merges

    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compatm), 'Si_ifrac', icewgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ocnwgt(lbound(icewgt,1):ubound(icewgt,1),lbound(icewgt,2):ubound(icewgt,2)))
    do j=lbound(icewgt,2),ubound(icewgt,2)
    do i=lbound(icewgt,1),ubound(icewgt,1)
      ocnwgt = 1.0_ESMF_KIND_R8 - icewgt
    enddo
    enddo

    !--- merges

    call shr_nuopc_methods_FB_FieldMerge(is_local%wrap%FBExp(compatm)  ,'surface_temperature'    ,         &
                                  is_local%wrap%FBImp(compocn,compatm) ,'sea_surface_temperature', ocnwgt, &
                                  is_local%wrap%FBImp(compice,compatm) ,'sea_ice_temperature'    , icewgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(ocnwgt)
#endif

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_atm

  !-----------------------------------------------------------------------------
  subroutine med_phases_prep_ice(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the ICE import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_ice)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    call ESMF_TimeGet(time,timestring=timestr)
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    do n1 = 1,ncomps
      n2 = compice
      if (is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(&
           fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), trim(subname)//' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compice), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(compice), &
                    FB1=is_local%wrap%FBImp(compocn,compice), FB1w=is_local%wrap%FBfrac(compice), fldw1='ofrac', &
                    FB2=is_local%wrap%FBImp(compatm,compice), FB2w=is_local%wrap%FBfrac(compice), fldw2='afrac', &
                    FB3=is_local%wrap%FBImp(compglc,compice), &
                    FB4=is_local%wrap%FBImp(comprof,compice), &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_ice

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_lnd(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the LND import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_lnd)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    call ESMF_TimeGet(time,timestring=timestr)
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    do n1 = 1,ncomps
      n2 = complnd
      if (is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), trim(subname)//' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(complnd), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(complnd), &
                    FB1=is_local%wrap%FBImp(compatm,complnd), FB1w=is_local%wrap%FBfrac(complnd), fldw1='afrac', &
                    FB2=is_local%wrap%FBImp(compglc,complnd), &
                    FB3=is_local%wrap%FBImp(comprof,complnd), &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_lnd

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_rof(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the ROF import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_rof)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    call ESMF_TimeGet(time,timestring=timestr)
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    do n1 = 1,ncomps
      n2 = comprof
      if (is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), trim(subname)//' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(comprof), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(comprof), &
                    is_local%wrap%FBImp(complnd,comprof), is_local%wrap%FBfrac(comprof), 'lfrac', &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_rof

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_wav(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the WAV import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_wav)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    call ESMF_TimeGet(time,timestring=timestr)
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    do n1 = 1,ncomps
      n2 = compwav
      if (is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), trim(subname)//' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compwav), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(compwav), &
                    FB1=is_local%wrap%FBImp(compocn,compwav), &
                    FB2=is_local%wrap%FBImp(compice,compwav), &
                    FB3=is_local%wrap%FBImp(compatm,compwav), &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_wav

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_glc(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Prepares the GLC import Fields.

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_glc)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    call ESMF_TimeGet(time,timestring=timestr)
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    do n1 = 1,ncomps
      n2 = compglc
      if (is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), trim(subname)//' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    !---------------------------------------
    !--- auto merges
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compglc), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(compglc), &
                    is_local%wrap%FBImp(complnd,compglc), is_local%wrap%FBfrac(compglc), 'lfrac', &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations
    !---------------------------------------

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_glc

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_StateItem_Flag)   :: itemType
    type(InternalState)         :: is_local
    integer                     :: i,j,n,n1,n2
    character(ESMF_MAXSTR)      :: fieldname1(10),fieldname2(10),fieldname3(10)
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:)
    real(ESMF_KIND_R8), pointer :: atmwgt(:),icewgt(:),customwgt(:)
    logical                     :: checkOK, checkOK1, checkOK2
    logical,save                :: first_call = .true.
    character(len=*),parameter  :: subname='(med_phases_prep_ocn)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- average ocn accumulator
    !---------------------------------------

    if (dbug_flag > 7) then
      call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), trim(subname)//' FBaccO_B4avg ', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call shr_nuopc_methods_FB_average(is_local%wrap%FBExpAccum(compocn), is_local%wrap%FBExpAccumCnt(compocn), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), trim(subname)//' FBaccO_avg ', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- copy to FBExp(compocn)
    !---------------------------------------

    call shr_nuopc_methods_FB_copy(is_local%wrap%FBExp(compocn), is_local%wrap%FBExpAccum(compocn), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- zero accumulator
    !---------------------------------------

    is_local%wrap%FBExpAccumCnt(compocn) = 0
    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExpAccum(compocn), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), trim(subname)//' FBaccO_AFzero ', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- update local scalar data
    !---------------------------------------

    !is_local%wrap%scalar_data(1) =

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_prep_ocn

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  subroutine med_phases_accum_fast(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is_local
    integer                     :: i,j,n,n1,n2
    logical,save                :: first_call = .true.
    character(len=*),parameter :: subname='(med_phases_accum_fast)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- mapping
    !---------------------------------------

    do n1 = 1,ncomps
      n2 = compocn
      if (is_local%wrap%med_coupling_active(n1,n2)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_Regrid(&
           fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
           consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
           consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
           bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
           patchmap=is_local%wrap%RH(n1,n2,mappatch), &
           string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
           call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(n1,n2), trim(subname) &
                //' FBImp('//trim(compname(n1))//','//trim(compname(n2))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    enddo

    !---------------------------------------
    !--- flux calculation (tcx maybe should be separate phase?)
    !---------------------------------------

    !call flux_calc(is_local%wrap%FBAtmOcn_o, rc=rc)

    !---------------------------------------
    !--- auto merges to ocn
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBExp(compocn), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call med_merge_auto(is_local%wrap%FBExp(compocn), &
                    FB1=is_local%wrap%FBImp(compatm,compocn)   , FB1w=is_local%wrap%FBfrac(compocn), fldw1='afrac', &
                    FB2=is_local%wrap%FBAtmOcn_o               , FB2w=is_local%wrap%FBfrac(compocn), fldw2='afrac', &
                    FB3=is_local%wrap%FBImp(compice,compocn)   , FB3w=is_local%wrap%FBfrac(compocn), fldw3='ifrac', &
                    FB4=is_local%wrap%FBImp(comprof,compocn)   , &
                    FB5=is_local%wrap%FBImp(compglc,compocn)   , &
                    document=first_call, string=subname, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- custom calculations to ocn
    !---------------------------------------

#if (1 == 0)

    ! atm and ice fraction
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), trim(ice_fraction_name), icewgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(atmwgt(lbound(icewgt,1):ubound(icewgt,1),lbound(icewgt,2):ubound(icewgt,2)))
    allocate(customwgt(lbound(icewgt,1):ubound(icewgt,1),lbound(icewgt,2):ubound(icewgt,2)))
    do j=lbound(icewgt,2),ubound(icewgt,2)
    do i=lbound(icewgt,1),ubound(icewgt,1)
      atmwgt = 1.0_ESMF_KIND_R8 - icewgt
    enddo
    enddo

    !-------------
    ! mean_evap_rate = mean_laten_heat_flux * (1-ice_fraction)/const_lhvap
    !-------------

!    customwgt = atmwgt / const_lhvap
!    call shr_nuopc_methods_FB_FieldMerge(is_local%wrap%FBExp(compocn),'mean_evap_rate' , &
!                                  is_local%wrap%FBImp(compatm,compocn), 'mean_laten_heat_flux' ,customwgt, rc=rc)
!    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------
    ! field_for_ocn = field_from_atm * (1-ice_fraction)
    !-------------

    call shr_nuopc_methods_FB_FieldMerge(is_local%wrap%FBExp(compocn),'mean_fprec_rate' , &
                                  is_local%wrap%FBImp(compatm,compocn), 'mean_fprec_rate' ,atmwgt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------
    ! End merges
    !-------------

    deallocate(atmwgt,customwgt)

    if (dbug_flag > 1) then
      call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compocn), trim(subname)//' FB4ocn_AFmrg ', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#endif

    !---------------------------------------
    !--- ocean accumulator
    !---------------------------------------

    call shr_nuopc_methods_FB_accum(is_local%wrap%FBExpAccum(compocn), is_local%wrap%FBExp(compocn), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    is_local%wrap%FBExpAccumCnt(compocn) = is_local%wrap%FBExpAccumCnt(compocn) + 1

    if (dbug_flag > 1) then
      call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), trim(subname)//' FBaccOcn_AFaccum ', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_accum_fast

  !-----------------------------------------------------------------------------

  subroutine med_phases_map_frac(FBin, fldin, FBout, fldout, RH, iarr, iarr_r, rc)

    type(ESMF_FieldBundle) :: FBin            ! input FB
    character(len=*)       :: fldin           ! input fieldname
    type(ESMF_FieldBundle) :: FBout           ! output FB
    character(len=*)       :: fldout          ! output fieldname
    type(ESMF_RouteHandle) :: RH              ! RH
    real(ESMF_KIND_R8), pointer :: iarr(:)    ! mapped field
    real(ESMF_KIND_R8), pointer :: iarr_r(:)  ! reciprical of mapped field
    integer :: rc                             ! return code
    !---------------------------------------

    ! local
    integer :: i
    real(ESMF_KIND_R8), pointer :: dataPtr2(:)
    character(len=*),parameter :: subname='(med_phases_map_frac)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    if (ESMF_RouteHandleIsCreated(RH, rc=rc)) then
      call shr_nuopc_methods_FB_FieldRegrid(FBin, fldin, FBout, fldout, RH, rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--- copy out the ifrac on ice grid and ifrac on atm grid
      call shr_nuopc_methods_FB_GetFldPtr(FBout, fldout, fldptr1=dataPtr2, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      allocate(iarr(lbound(dataptr2,1):ubound(dataptr2,1)))
      allocate(iarr_r(lbound(dataptr2,1):ubound(dataptr2,1)))
      do i=lbound(dataptr2,1),ubound(dataptr2,1)
        !--- compute ice fraction on atm grid and reciprocal
        if (dataPtr2(i) == 0.0_ESMF_KIND_R8) then
          iarr(i)   = 0.0_ESMF_KIND_R8
          iarr_r(i) = 0.0_ESMF_KIND_R8
        else
          iarr(i)   = dataPtr2(i)
          iarr_r(i) = 1.0_ESMF_KIND_R8/dataPtr2(i)
        endif
      enddo
    endif

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_map_frac

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_map_atm2ocn(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: n1, n2
    type(InternalState) :: is_local
    character(len=*),parameter :: subname='(med_phases_prep_atmocn_o)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    n1 = compatm
    n2 = compocn
    call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(n1,n2), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_Regrid(&
         fldsFr(n1), is_local%wrap%FBImp(n1,n1), is_local%wrap%FBImp(n1,n2), &
         consfmap=is_local%wrap%RH(n1,n2,mapconsf), &
         consdmap=is_local%wrap%RH(n1,n2,mapconsd), &
         bilnrmap=is_local%wrap%RH(n1,n2,mapbilnr), &
         patchmap=is_local%wrap%RH(n1,n2,mappatch), &
         string=trim(compname(n1))//'2'//trim(compname(n2)), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: add this to the run sequence
  end subroutine med_phases_prep_map_atm2ocn

  !-----------------------------------------------------------------------------

  subroutine med_phases_atmocn_init(gcomp, rc)
    use shr_kind_mod  , only : CL=>shr_kind_cl

    ! arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(CL)       :: cvalue
    character(CL)       :: aoflux_grid
    type(InternalState) :: is_local
    character(len=*),parameter :: subname='(med_phases_atmocn_init)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aoflux_grid

    if (trim(aoflux_grid) == 'ocn') then

       call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compatm,compocn), value=czero, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_Regrid( &
            fldsFr(compatm), is_local%wrap%FBImp(compatm,compatm), is_local%wrap%FBImp(compatm,compocn), &
            consfmap=is_local%wrap%RH(compatm,compocn,mapconsf), &
            consdmap=is_local%wrap%RH(compatm,compocn,mapconsd), &
            bilnrmap=is_local%wrap%RH(compatm,compocn,mapbilnr), &
            patchmap=is_local%wrap%RH(compatm,compocn,mappatch), &
            string=trim(compname(compatm))//'2'//trim(compname(compocn)), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(compatm,compocn), &
               trim(subname) //' FBImp('//trim(compname(compatm))//','//trim(compname(compocn))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       call med_atmocn_init(gcomp, aoflux_grid,         &
            FBAtmOcn=is_local%wrap%FBAtmOcn_o,          &
            FBAtm=is_local%wrap%FBImp(compatm,compocn), &
            FBOcn=is_local%wrap%FBImp(compocn,compocn), &
            FBfrac=is_local%wrap%FBfrac(compocn), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    elseif (trim(aoflux_grid) == 'atm') then

       call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(compocn,compatm), value=czero, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_Regrid( &
            fldsFr(compocn), is_local%wrap%FBImp(compocn,compocn), is_local%wrap%FBImp(compocn,compatm), &
            consfmap=is_local%wrap%RH(compocn,compatm,mapconsf), &
            consdmap=is_local%wrap%RH(compocn,compatm,mapconsd), &
            bilnrmap=is_local%wrap%RH(compocn,compatm,mapbilnr), &
            patchmap=is_local%wrap%RH(compocn,compatm,mappatch), &
            string=trim(compname(compocn))//'2'//trim(compname(compatm)), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(compocn,compatm), &
               trim(subname) //' FBImp('//trim(compname(compocn))//','//trim(compname(compatm))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       call med_atmocn_init(gcomp, aoflux_grid,         &
            FBAtmOcn=is_local%wrap%FBAtmOcn_a,          &
            FBAtm=is_local%wrap%FBImp(compatm,compatm), &
            FBOcn=is_local%wrap%FBImp(compocn,compatm), &
            FBfrac=is_local%wrap%FBfrac(compatm), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    else

       call ESMF_LogWrite(trim(subname)//' aoflux_grid = '//trim(aoflux_grid)//' not available', ESMF_LOGMSG_INFO, rc=dbrc)
       return

    endif

  end subroutine med_phases_atmocn_init

  !-----------------------------------------------------------------------------

  subroutine med_phases_atmocn_ocnalb(gcomp, rc)
    use seq_flds_mod  , only : seq_flds_scalar_index_nextsw_cday

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    real(ESMF_KIND_R8)         :: nextsw_cday  ! calendar day of next atm shortwave
    character(len=*),parameter :: subname='(med_phases_atmocn_ocnalb)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine nextsw_cday
    call shr_nuopc_methods_State_GetScalar(is_local%wrap%NstateImp(compatm), &
        scalar_id=seq_flds_scalar_index_nextsw_cday, value=nextsw_cday, mpicom=is_local%wrap%mpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Calculate ocean albedoes
    call med_atmocn_ocnalb(gcomp, FBfrac_o=is_local%wrap%FBfrac(compocn), nextsw_cday=nextsw_cday, rc=rc)

  end subroutine med_phases_atmocn_ocnalb

  !-----------------------------------------------------------------------------

  subroutine med_phases_atmocn_flux(gcomp, rc)
    use seq_flds_mod  , only : seq_flds_scalar_index_dead_comps

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: clock
    real(ESMF_KIND_R8)  :: rdead_comps
    logical             :: dead_comps
    logical             :: ocn_prognostic
    character(len=*),parameter :: subname='(med_phases_atmocn_flux)'
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the clock from the mediator Component
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine dead_comps
    call shr_nuopc_methods_State_GetScalar(is_local%wrap%NStateImp(compatm), &
         scalar_id=seq_flds_scalar_index_dead_comps, value=rdead_comps, mpicom=is_local%wrap%mpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    dead_comps = (rdead_comps /= 0._ESMF_KIND_r8)
    write(6,*)'DEBUG: atm rdead, dead is ',rdead_comps, dead_comps

    call shr_nuopc_methods_State_GetScalar(is_local%wrap%NStateImp(compocn), &
         scalar_id=seq_flds_scalar_index_dead_comps, value=rdead_comps, mpicom=is_local%wrap%mpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    dead_comps = (rdead_comps /= 0._ESMF_KIND_r8)
    write(6,*)'DEBUG: ocn rdead, dead is ',rdead_comps, dead_comps

    ! Determine ocn_prognostic
    ocn_prognostic = NUOPC_IsConnected(is_local%wrap%NStateImp(compocn))

    ! Calculate atm/ocn fluxes
    call med_atmocn_flux(gcomp, clock, ocn_prognostic=ocn_prognostic, dead_comps=dead_comps, rc=rc)

  end subroutine med_phases_atmocn_flux

end module med_phases_mod
