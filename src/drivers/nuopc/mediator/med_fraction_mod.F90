module med_fraction_mod

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  ! This mediator operates on two timescales and keeps two internal Clocks to
  ! do so.
  !-----------------------------------------------------------------------------

  use ESMF
  use shr_nuopc_methods_mod
  use med_internalstate_mod
  use med_constants_mod

  implicit none

  private

  integer            :: dbrc
  integer           , parameter :: dbug_flag   = med_constants_dbug_flag
  logical           , parameter :: statewrite_flag = med_constants_statewrite_flag
  real(ESMF_KIND_R8), parameter :: spval_init  = med_constants_spval_init
  real(ESMF_KIND_R8), parameter :: spval       = med_constants_spval
  real(ESMF_KIND_R8), parameter :: czero       = med_constants_czero
  character(*),parameter :: u_FILE_u = &
    __FILE__

  public med_fraction_setupflds
  public med_fraction_init
  public med_fraction_set

  integer, parameter, public :: nfracs = 5
  character(len=5),public :: fraclist(nfracs,ncomps)

  character(len=5),parameter,dimension(5) :: fraclist_a = (/'afrac','ifrac','ofrac','lfrac','lfrin'/)
  character(len=5),parameter,dimension(5) :: fraclist_o = (/'afrac','ifrac','ofrac','ifrad','ofrad'/)
  character(len=5),parameter,dimension(3) :: fraclist_i = (/'afrac','ifrac','ofrac'/)
  character(len=5),parameter,dimension(3) :: fraclist_l = (/'afrac','lfrac','lfrin'/)
  character(len=5),parameter,dimension(2) :: fraclist_g = (/'gfrac','lfrac'/)
  character(len=5),parameter,dimension(2) :: fraclist_r = (/'lfrac','rfrac'/)
  character(len=5),parameter,dimension(1) :: fraclist_w = (/'wfrac'/)

  !--- standard ---
  real(ESMF_KIND_R8),parameter :: eps_fracsum = 1.0e-02   ! allowed error in sum of fracs
  real(ESMF_KIND_R8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  real(ESMF_KIND_R8),parameter :: eps_fraclim = 1.0e-03   ! truncation limit in fractions_a(lfrac)
  logical ,parameter :: atm_frac_correct = .false. ! turn on frac correction on atm grid
  !--- standard plus atm fraction consistency ---
  !  real(ESMF_KIND_R8),parameter :: eps_fracsum = 1.0e-12   ! allowed error in sum of fracs
  !  real(ESMF_KIND_R8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  !  real(ESMF_KIND_R8),parameter :: eps_fraclim = 1.0e-03   ! truncation limit in fractions_a(lfrac)
  !  logical ,parameter :: atm_frac_correct = .true. ! turn on frac correction on atm grid
  !--- unconstrained and area conserving? ---
  !  real(ESMF_KIND_R8),parameter :: eps_fracsum = 1.0e-12   ! allowed error in sum of fracs
  !  real(ESMF_KIND_R8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  !  real(ESMF_KIND_R8),parameter :: eps_fraclim = 1.0e-20   ! truncation limit in fractions_a(lfrac)
  !  logical ,parameter :: atm_frac_correct = .true. ! turn on frac correction on atm grid

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine med_fraction_setupflds(rc)
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname='(med_fraction_setupflds)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    fraclist(:,:) = ' '
    fraclist(1:size(fraclist_a),compatm) = fraclist_a
    fraclist(1:size(fraclist_o),compocn) = fraclist_o
    fraclist(1:size(fraclist_i),compice) = fraclist_i
    fraclist(1:size(fraclist_l),complnd) = fraclist_l
    fraclist(1:size(fraclist_r),comprof) = fraclist_r
    fraclist(1:size(fraclist_w),compwav) = fraclist_w
    fraclist(1:size(fraclist_g),compglc) = fraclist_g

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_fraction_setupflds

  !-----------------------------------------------------------------------------

  subroutine med_fraction_init(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! This subroutine initializes the fractions

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr(:)
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n,n1
    character(len=*),parameter  :: subname='(med_fraction_init)'

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

    !---------------------------------------
    !--- zero out FBfracs
    !---------------------------------------

    do n1 = 1,ncomps
      if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n1),rc=rc)) then
        call shr_nuopc_methods_FB_reset(is_local%wrap%FBfrac(n1), value=czero, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo

    !---------------------------------------
    !--- Initialize fractions on atm grid/decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'afrac', dataPtr, rc=rc)
      ! If 'afrac' exists, then set it to 1.
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        dataPtr = 1.0_ESMF_KIND_R8
      endif
      if (is_local%wrap%comp_present(compocn)) then
         ! map atm 'afrac' to ocn 'afrac' conservatively
         call shr_nuopc_methods_FB_FieldRegrid(&
              is_local%wrap%FBfrac(compatm), 'afrac', &
              is_local%wrap%FBfrac(compocn), 'afrac', &
              is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      if (is_local%wrap%comp_present(compice)) then
         ! map atm 'afrac' to ice 'afrac' conservatively
         call shr_nuopc_methods_FB_FieldRegrid(&
              is_local%wrap%FBfrac(compatm), 'afrac', &
              is_local%wrap%FBfrac(compice), 'afrac', &
              is_local%wrap%RH(compatm,compice,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    endif

    !---------------------------------------
    !--- Initialize fractions on glc grid decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(compglc)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compglc), 'gfrac', dataPtr1, rc=rc)
      ! If 'gfrac' and 'frac' exists, then copy 'frac' to 'gfrac'
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), 'frac' , dataPtr2, rc=rc)
        if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
          dataPtr1 = dataPtr2
        endif
      endif
    endif

    !---------------------------------------
    !--- Initialize fractions on land grid decomp, just an initial "guess", updated later
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrin', dataPtr1, rc=rc)
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
         ! If 'lfrin' and 'Sl_lfrin' exist, then copy 'Sl_lfrin' to 'lfrin'
        ! TODO: make sure that clm sets Sl_lfrin
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(complnd,complnd) , 'Sl_lfrin' , dataPtr2, rc=rc)
        if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
          dataPtr1 = dataPtr2
        endif
        if (is_local%wrap%comp_present(compatm)) then
           ! map lnd 'lfrin' to atm 'lfrin' conservatively
           call shr_nuopc_methods_FB_FieldRegrid(&
                is_local%wrap%FBfrac(complnd), 'lfrin', &
                is_local%wrap%FBfrac(compatm), 'lfrin', &
                is_local%wrap%RH(complnd,compatm,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

           ! map atm 'afrac' to lnd 'afrac' conservatively
          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBfrac(compatm), 'afrac', &
               is_local%wrap%FBfrac(complnd), 'afrac', &
               is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    endif

    !---------------------------------------
    !--- Initialize fractions on rof grid/decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(comprof)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(comprof), 'rfrac', dataPtr1, rc=rc)
      ! If 'rfrac' exists then set it to 1.
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
         dataPtr1 = 1.0_ESMF_KIND_R8
         ! call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(comprof,comprof) , 'frac' , dataPtr2, rc=rc)
         ! if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
         !   dataPtr1 = dataPtr2
         ! endif
      endif
    endif

    !---------------------------------------
    !--- Initialize fractions on wav grid decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(compwav)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compwav), 'wfrac', dataPtr, rc=rc)
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        dataPtr = 1.0_ESMF_KIND_R8
      endif
    endif

    !---------------------------------------
    !--- Initialize fractions on ice grid/decomp (initialize ice fraction to zero)
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then
      ! determine if 'ofrac' is in Fbfrac(compice)
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', dataPtr1, rc=rc)
      ! if 'ofrac' and 'Si_imask' exist then copy 'Si_imask' to 'ofrac'
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , dataPtr2, rc=rc)
        if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
          ! If not then assume that 'Si_imask' is, then have 'Si_imask' point to the same place as 'ofrac'
          dataPtr1 = dataPtr2
        endif

        if (is_local%wrap%comp_present(compatm)) then
          ! map ice 'ofrac' to atm 'ofrac' conservatively
           call shr_nuopc_methods_FB_FieldRegrid(&
                is_local%wrap%FBfrac(compice), 'ofrac', &
                is_local%wrap%FBfrac(compatm), 'ofrac', &
                is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! map atm 'afrac' to ice 'afrac' conservatively
          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBfrac(compatm), 'afrac', &
               is_local%wrap%FBfrac(compice), 'afrac', &
               is_local%wrap%RH(compatm,compice,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    endif

    !---------------------------------------
    !--- Initialize fractions on ocean grid/decomp
    !--- These are initialized the same as for ice
    !---------------------------------------

    if (is_local%wrap%comp_present(compocn)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', dataPtr1, rc=rc)
      ! if 'ofrac' and 'Si_imask' exist then copy 'So_omask' to 'ofrac'
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compocn,compocn) , 'So_omask' , dataPtr2, rc=rc)
        if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
          dataPtr1 = dataPtr2
        endif
        if (is_local%wrap%comp_present(compatm)) then
          ! map ocn 'ofrac' to atm 'ofrac' conservatively
           call shr_nuopc_methods_FB_FieldRegrid(&
                is_local%wrap%FBfrac(compocn), 'ofrac', &
                is_local%wrap%FBfrac(compatm), 'ofrac', &
                is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! map atm 'afrac' to ocn 'afrac' conservatively
          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBfrac(compatm), 'afrac', &
               is_local%wrap%FBfrac(compocn), 'afrac', &
               is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
    endif

    !---------------------------------------
    !--- Set ofrac and lfrac on atm grid.  These should actually be mapo2a of
    !--- ofrac and lfrac but we can't map lfrac from o2a due to masked mapping
    !---  weights.  So we have to settle for a residual calculation that is
    !---  truncated to zero to try to preserve "all ocean" cells.
    !---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', dataPtr1, rc=rc)
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', dataPtr2, rc=rc)
        ! If 'ofrac' and 'lfrac' exist then set 'lfrac' to 1. - 'ofrac'
        if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
          if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compice)) then
            dataPtr2 = 1.0_ESMF_KIND_R8 - dataPtr1
            if (atm_frac_correct) then
              where (dataPtr2 < eps_fraclim) dataPtr1 = 1.0_ESMF_KIND_R8
            endif
            where (dataPtr2 < eps_fraclim) dataPtr2 = 0.0_ESMF_KIND_R8
          elseif (is_local%wrap%comp_present(complnd)) then
            dataPtr1 = 1.0_ESMF_KIND_R8 - dataPtr2
            if (atm_frac_correct) then
              where (dataPtr1 < eps_fraclim) dataPtr2 = 1.0_ESMF_KIND_R8
            endif
            where (dataPtr1 < eps_fraclim) dataPtr1 = 0.0_ESMF_KIND_R8
          endif
        endif
      endif
    endif

    !---------------------------------------
    !--- finally, set fractions_l(lfrac) from fractions_a(lfrac)
    !--- and fractions_r(lfrac) from fractions_l(lfrac)
    !--- and fractions_g(lfrac) from fractions_l(lfrac)
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then
      if (is_local%wrap%comp_present(compatm)) then
         call shr_nuopc_methods_FB_FieldRegrid(&
              is_local%wrap%FBfrac(compatm), 'lfrac', &
              is_local%wrap%FBfrac(complnd), 'lfrac', &
              is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      else
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrin', dataPtr1, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrac', dataPtr2, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        dataPtr2 = dataPtr1
      endif

      if (is_local%wrap%comp_present(comprof)) then
         call shr_nuopc_methods_FB_FieldRegrid(&
              is_local%wrap%FBfrac(complnd), 'lfrac', &
              is_local%wrap%FBfrac(comprof), 'lfrac', &
              is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      if (is_local%wrap%comp_present(compglc)) then
         ! tcraig l2g_consf does not exist yet
         ! call shr_nuopc_methods_FB_FieldRegrid(&
         !      is_local%wrap%FBfrac(complnd), 'lfrac', &
         !      is_local%wrap%FBfrac(compglc), 'lfrac', &
         !      is_local%wrap%RH(complnd,compglc,mapconsf), rc=rc)
         ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    endif

    !---------------------------------------
    !--- update time varying fractions
    !---------------------------------------

    call med_fraction_set(gcomp, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- check fractions
    !---------------------------------------

    do n1 = 1,ncomps
      if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n1),rc=rc)) then
        call med_fraction_check(is_local%wrap%FBfrac(n1), compname(n1)//'frac init', rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_fraction_init

  !-----------------------------------------------------------------------------

  subroutine med_fraction_set(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! This subroutine initializes the fractions

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Field)            :: field
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: dataPtr(:)
    real(ESMF_KIND_R8), pointer :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                     :: i,j,n
    character(len=*),parameter  :: subname='(med_fraction_set)'

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

    !---------------------------------------
    !--- update ice fraction
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', dataPtr1, rc=rc)
      ! If 'ifrac' and 'ofrac' exist then set 'ifrac' to 'ofrac' and 'ofrac' to 1.-'ifrac'
      if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
        call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', dataPtr3, rc=rc)
        if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_ifrac' , dataPtr2, rc=rc)
          if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
            dataPtr1 = dataPtr2
            dataPtr3 = 1.0_ESMF_KIND_R8 - dataPtr1
          endif
        endif
      endif

      if (is_local%wrap%comp_present(compocn)) then
         call shr_nuopc_methods_FB_FieldRegrid(&
              is_local%wrap%FBfrac(compice), 'ifrac', &
              is_local%wrap%FBfrac(compocn), 'ifrac', &
              is_local%wrap%RH(compice,compocn,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_FieldRegrid(&
             is_local%wrap%FBfrac(compice), 'ofrac', &
             is_local%wrap%FBfrac(compocn), 'ofrac', &
             is_local%wrap%RH(compice,compocn,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      if (is_local%wrap%comp_present(compatm)) then
         call shr_nuopc_methods_FB_FieldRegrid(&
              is_local%wrap%FBfrac(compice), 'ifrac', &
              is_local%wrap%FBfrac(compatm), 'ifrac', &
              is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_FieldRegrid(&
             is_local%wrap%FBfrac(compice), 'ofrac', &
             is_local%wrap%FBfrac(compatm), 'ofrac', &
             is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (atm_frac_correct) then
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ifrac', dataPtr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', dataPtr2, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', dataPtr3, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          where (dataPtr1 + dataPtr2 > 0.0_ESMF_KIND_R8)
            dataPtr1 = dataPtr1 * ((1.0_ESMF_KIND_R8 - dataPtr3)/(dataPtr2+dataPtr1))
            dataPtr2 = dataPtr2 * ((1.0_ESMF_KIND_R8 - dataPtr3)/(dataPtr2+dataPtr1))
          elsewhere
            dataPtr1 = 0.0_ESMF_KIND_R8
            dataPtr2 = 0.0_ESMF_KIND_R8
          end where
        endif
      endif
    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_fraction_set

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  subroutine med_fraction_check(FB, string, rc)
    type(ESMF_FieldBundle)  :: FB
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! This subroutine initializes the fractions

    character(len=*),parameter  :: subname='(med_fraction_check)'

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    !--- write fraction info
    !---------------------------------------

    call shr_nuopc_methods_FB_diagnose(FB, subname//trim(string), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_fraction_check

  !-----------------------------------------------------------------------------

end module med_fraction_mod
