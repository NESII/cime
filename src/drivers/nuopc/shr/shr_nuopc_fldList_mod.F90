!================================================================================
module shr_nuopc_fldList_mod

  use shr_kind_mod,only : r8 => shr_kind_r8
  use shr_kind_mod,only : CL => SHR_KIND_CL, CX => SHR_KIND_CX, CS => SHR_KIND_CS
  use shr_sys_mod, only : shr_sys_abort
  use shr_log_mod, only : loglev  => shr_log_Level
  use shr_log_mod, only : logunit => shr_log_Unit
  use shr_mpi_mod, only : shr_mpi_bcast
  use shr_string_mod, only : shr_string_listGetNum, shr_string_listGetName
  use seq_flds_mod, only : seq_flds_lookup, seq_flds_get_num_entries, seq_flds_get_entry
  use seq_flds_mod, only : seq_flds_scalar_name, seq_flds_scalar_num
  use ESMF
  use NUOPC

  implicit none
  save
  private

  public :: shr_nuopc_fldList_setDict_fromseqflds
  public :: shr_nuopc_fldList_fromseqflds
  public :: shr_nuopc_fldList_Type
  public :: shr_nuopc_fldList_Zero
  public :: shr_nuopc_fldList_Add
  public :: shr_nuopc_fldList_Advertise
  public :: shr_nuopc_fldList_Realize
  public :: shr_nuopc_fldList_SetScalarField
!  public :: shr_nuopc_fldList_CopyScalarToState
!  public :: shr_nuopc_fldList_CopyStateToScalar

  type shr_nuopc_fldList_Type
    integer :: num
    character(CS),pointer :: stdname(:)       => null()
    character(CS),pointer :: shortname(:)     => null()
    character(CS),pointer :: transferOffer(:) => null()  ! provide or accept
    character(CS),pointer :: mapping(:)       => null()
  end type shr_nuopc_fldList_Type

  integer :: dbrc
  character(len=CL) :: infostr

!================================================================================
contains
!================================================================================

  subroutine shr_nuopc_fldList_setDict_fromseqflds(rc)
    ! ----------------------------------------------
    ! Build NUOPC dictionary from seq_flds data
    ! ----------------------------------------------
    integer,          intent(inout)  :: rc

    ! local variables
    integer :: n, num
    character(CS) :: stdname
    character(CS) :: units
    character(len=*), parameter :: subname='(shr_nuopc_fldList_setDict_fromseqflds)'

    rc = ESMF_SUCCESS

    call seq_flds_get_num_entries(num)

    do n = 1,num
      call seq_flds_get_entry(n,shortname=stdname,units=units)
      if (.not.NUOPC_FieldDictionaryHasEntry(stdname)) then
        call ESMF_LogWrite(subname//': Add:'//trim(stdname), ESMF_LOGMSG_INFO, rc=rc)
        call NUOPC_FieldDictionaryAddEntry(standardName=stdname, canonicalUnits=units, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo

  end subroutine shr_nuopc_fldList_setDict_fromseqflds

!================================================================================

  subroutine shr_nuopc_fldList_fromseqflds(fldlist, seq_flds_list, transferOffer, tag, mapping, rc)
    ! ----------------------------------------------
    ! Build fldlist from seq_flds list
    ! ----------------------------------------------
    type(shr_nuopc_fldList_Type), intent(inout)  :: fldlist
    character(len=*)            , intent(in)     :: seq_flds_list
    character(len=*), intent(in)     :: transferOffer
    character(len=*), intent(in)     :: tag
    integer,          intent(inout)  :: rc
    character(len=*), intent(in), optional :: mapping

    ! local variables
    integer :: n, num
    character(len=CS) :: stdname
    character(len=CS) :: lmapping
    character(len=*), parameter :: subname='(shr_nuopc_fldList_fromseqflds)'

    rc = ESMF_SUCCESS

    if (present(mapping)) then
      lmapping = mapping
    else
      lmapping = "undefined"
    endif

    num = shr_string_listGetNum(seq_flds_list)
    do n = 1,num
      call shr_string_listGetName(seq_flds_list, n, stdname)
      call shr_nuopc_fldList_Add(fldlist, stdname, transferOffer, tag, mapping=lmapping, rc=rc)
    enddo

  end subroutine shr_nuopc_fldList_fromseqflds

!================================================================================

  subroutine shr_nuopc_fldList_Zero(fldlist, rc)
    ! ----------------------------------------------
    ! Zero out list of field information
    ! ----------------------------------------------
    type(shr_nuopc_fldList_Type), intent(inout)  :: fldlist
    integer,          intent(inout)  :: rc

    ! local variables
    integer :: num
    character(len=*), parameter :: subname='(shr_nuopc_fldList_Zero)'

    rc = ESMF_SUCCESS

    fldlist%num = 0
    if (associated(fldlist%stdname))       deallocate(fldlist%stdname)
    if (associated(fldlist%shortname))     deallocate(fldlist%shortname)
    if (associated(fldlist%transferOffer)) deallocate(fldlist%transferOffer)
    if (associated(fldlist%mapping))       deallocate(fldlist%mapping)
    allocate(fldlist%stdname(10))
    allocate(fldlist%shortname(10))
    allocate(fldlist%transferOffer(10))
    allocate(fldlist%mapping(10))

  end subroutine shr_nuopc_fldList_Zero

!================================================================================

  subroutine shr_nuopc_fldList_Add(fldlist, stdname, transferOffer, tag, shortname, mapping, rc)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    type(shr_nuopc_fldList_Type), intent(inout)  :: fldlist
    character(len=*), intent(in)     :: stdname
    character(len=*), intent(in)     :: transferOffer
    character(len=*), intent(in)     :: tag
    integer,          intent(inout)  :: rc
    character(len=*), intent(in)   , optional :: shortname
    character(len=*), intent(in)   , optional :: mapping

    ! local variables
    integer :: num,n
    character(len=CS),pointer :: copystr(:)
    character(len=*), parameter :: subname='(shr_nuopc_fldList_Add)'

    rc = ESMF_SUCCESS

    ! make sure fldlist has been zero-ed out first
    if (.not.associated(fldlist%stdname) .or. .not.associated(fldlist%shortname) .or. &
        .not.associated(fldlist%transferOffer) .or. .not.associated(fldlist%mapping)) then
      call ESMF_LogWrite(trim(subname)//":"//trim(tag)//" ERROR in fldlist, call shr_nuopc_fldList_Zero first "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      rc = ESMF_FAILURE
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! make sure fldlist arrays are consistent
    if ((size(fldlist%stdname) /= size(fldlist%shortname)) .or. &
        (size(fldlist%stdname) /= size(fldlist%transferOffer)) .or. &
        (size(fldlist%stdname) /= size(fldlist%mapping))) then
      call ESMF_LogWrite(trim(subname)//":"//trim(tag)//" ERROR in fldlist, size of arrays out of sync "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      rc = ESMF_FAILURE
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! make sure fldlist num has not gotten ahead of fldlist arrays
    if (fldlist%num > size(fldlist%stdname)) then
      call ESMF_LogWrite(trim(subname)//":"//trim(tag)//" ERROR fldlist%num gt size of fldlist%stdname "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      rc = ESMF_FAILURE
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! increase size of fldlist arrays
    if (fldlist%num == size(fldlist%stdname)) then
      num = fldlist%num
      allocate(copystr(num))

      copystr(1:num) = fldlist%stdname(1:num)
      deallocate(fldlist%stdname)
      allocate(fldlist%stdname(num+10))
      fldlist%stdname(1:num) = copystr(1:num)

      copystr(1:num) = fldlist%shortname(1:num)
      deallocate(fldlist%shortname)
      allocate(fldlist%shortname(num+10))
      fldlist%shortname(1:num) = copystr(1:num)

      copystr(1:num) = fldlist%transferOffer(1:num)
      deallocate(fldlist%transferOffer)
      allocate(fldlist%transferOffer(num+10))
      fldlist%transferOffer(1:num) = copystr(1:num)

      copystr(1:num) = fldlist%mapping(1:num)
      deallocate(fldlist%mapping)
      allocate(fldlist%mapping(num+10))
      fldlist%mapping(1:num) = copystr(1:num)

      deallocate(copystr)
    endif

    ! fill in the new entry
    fldlist%num = fldlist%num + 1
    num = fldlist%num

    fldlist%stdname(num)      = trim(stdname)
    if (present(shortname)) then
       fldlist%shortname(num) = trim(shortname)
    else
       fldlist%shortname(num) = trim(stdname)
    endif
    fldlist%transferOffer(num)= trim(transferOffer)
    if (present(mapping)) then
       fldlist%mapping(num)   = trim(mapping)
    else
       fldlist%mapping(num)   = 'undefined'
    endif

  end subroutine shr_nuopc_fldList_Add

!================================================================================

  subroutine shr_nuopc_fldList_Advertise(state, fldlist, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(shr_nuopc_fldList_Type), intent(inout) :: fldlist
    character(len=*), intent(in)                :: tag
    integer,          intent(inout)             :: rc

    ! local variables
    integer       :: n
    character(CS) :: units
    character(ESMF_MAXSTR), pointer :: StandardNameList(:)
    character(ESMF_MAXSTR), pointer :: ConnectedList(:)
    character(ESMF_MAXSTR), pointer :: NameSpaceList(:)
    character(ESMF_MAXSTR), pointer :: itemNameList(:)
    character(len=*),parameter  :: subname='(shr_nuopc_fldList_Advertise)'

    rc = ESMF_SUCCESS

    do n = 1, fldlist%num

      call ESMF_LogWrite(subname//':'//trim(tag)//':'//trim(fldlist%stdname(n)), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! adding to dictionary here not a good idea because will only be added on this PET.  needs
! to be added to all PETs at the driver level.
!      if (.not.NUOPC_FieldDictionaryHasEntry(fldlist%stdname(n))) then
!        call ESMF_LogWrite(subname//':DictionaryAddEntry:'//trim(tag)//':'//trim(fldlist%stdname(n))//':'//trim(units), ESMF_LOGMSG_INFO, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!        call seq_flds_lookup(fldlist%stdname(n), units=units)
!        call NUOPC_FieldDictionaryAddEntry(standardName=fldlist%stdname(n), canonicalUnits=units, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!        if (.not.NUOPC_FieldDictionaryHasEntry(fldlist%shortname(n))) then
!          call ESMF_LogWrite(subname//':DictionaryAddEntry:'//trim(tag)//':'//trim(fldlist%shortname(n))//':'//trim(units), ESMF_LOGMSG_INFO, rc=rc)
!          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!          call NUOPC_FieldDictionaryAddEntry(standardName=fldlist%shortname(n), canonicalUnits=units, rc=rc)
!          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!!          if (fldlist%stdname(n) /= fldlist%shortname(n)) then
!!            call NUOPC_FieldDictionarySetSyno( (/fldlist%stdname(n), fldlist%shortname(n)/), rc=rc)
!!            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!!          endif
!        endif
!      endif

      call NUOPC_Advertise(state, &
        standardName=trim(fldlist%stdname(n)), &
        shortname=trim(fldlist%shortname(n)), &
        name=trim(fldlist%shortname(n)), &
        TransferOfferGeomObject=fldlist%transferOffer(n), &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    enddo

#if (1 == 0)
    if (fldlist%num > 0) then
      nullify(StandardNameList)
      nullify(ConnectedList)
      nullify(NameSpaceList)
      nullify(ItemNameList)

      call NUOPC_GetStateMemberLists(state, StandardNameList=StandardNameList, ConnectedList=ConnectedList, &
        NamespaceList=NamespaceList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      do n = 1,size(StandardNameList)
        call ESMF_LogWrite(trim(subname)//trim(tag)//" StandardNameList = "//trim(StandardNameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo
      do n = 1,size(ConnectedList)
        call ESMF_LogWrite(trim(subname)//trim(tag)//" ConnectedList = "//trim(ConnectedList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo
      do n = 1,size(NamespaceList)
        call ESMF_LogWrite(trim(subname)//trim(tag)//" NamespaceList = "//trim(NamespaceList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo
      do n = 1,size(ItemnameList)
        call ESMF_LogWrite(trim(subname)//trim(tag)//" ItemnameList = "//trim(ItemnameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo
    endif
#endif

  end subroutine shr_nuopc_fldList_Advertise

!================================================================================

  subroutine shr_nuopc_fldList_Realize(state, grid, mesh, fldlist, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in),optional        :: grid
    type(ESMF_Mesh), intent(in),optional        :: mesh
    type(shr_nuopc_fldList_Type), intent(inout) :: fldlist
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    ! local variables
    integer           :: n
    integer           :: itemCount
    type(ESMF_Field)  :: field
    character(ESMF_MAXSTR)      :: transferAction
    character(ESMF_MAXSTR), pointer :: StandardNameList(:)
    character(ESMF_MAXSTR), pointer :: ConnectedList(:)
    character(ESMF_MAXSTR), pointer :: NameSpaceList(:)
    character(ESMF_MAXSTR), pointer :: itemNameList(:)
    character(len=*),parameter  :: subname='(shr_nuopc_fldList_Realize)'

    rc = ESMF_SUCCESS

    if (present(grid) .and. present(mesh)) then
      call ESMF_LogWrite(trim(subname)//trim(tag)//": ERROR both grid and mesh not allowed", ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    nullify(StandardNameList)
    nullify(ConnectedList)
    nullify(NameSpaceList)
    nullify(ItemNameList)

    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(infostr,'(i6)') itemCount
    call ESMF_LogWrite(trim(subname)//trim(tag)//" count = "//trim(infostr), ESMF_LOGMSG_INFO, rc=dbrc)
    if (itemCount > 0) then
      allocate(itemNameList(itemCount))
      call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      do n = 1,itemCount
        call ESMF_LogWrite(trim(subname)//trim(tag)//" itemNameList = "//trim(itemNameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo
      deallocate(itemNameList)
    endif

#if (1 == 0)
    call NUOPC_GetStateMemberLists(state, StandardNameList=StandardNameList, ConnectedList=ConnectedList, &
      NamespaceList=NamespaceList, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(infostr,'(i6)') size(StandardNameList)
    call ESMF_LogWrite(trim(subname)//trim(tag)//" size = "//trim(infostr), ESMF_LOGMSG_INFO, rc=dbrc)

    do n = 1,size(StandardNameList)
      call ESMF_LogWrite(trim(subname)//trim(tag)//" StandardNameList = "//trim(StandardNameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    do n = 1,size(ConnectedList)
      call ESMF_LogWrite(trim(subname)//trim(tag)//" ConnectedList = "//trim(ConnectedList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    do n = 1,size(NamespaceList)
      call ESMF_LogWrite(trim(subname)//trim(tag)//" NamespaceList = "//trim(NamespaceList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    do n = 1,size(ItemnameList)
      call ESMF_LogWrite(trim(subname)//trim(tag)//" ItemnameList = "//trim(ItemnameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
#endif

    do n = 1, fldlist%num

!      call ESMF_LogWrite(subname//' fld = '//trim(fldlist%shortname(n)), ESMF_LOGMSG_INFO, rc=dbrc)

      if (NUOPC_IsConnected(state, fieldName=fldlist%shortname(n))) then

        call ESMF_StateGet(state, field=field, itemName=fldlist%shortname(n), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call NUOPC_GetAttribute(field, name="TransferActionGeomObject", value=transferAction, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (trim(transferAction) == "accept") then
          call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(fldlist%shortname(n))//" is connected, grid/mesh TBD", ESMF_LOGMSG_INFO, rc=dbrc)

        else   ! provide

          if (fldlist%shortname(n) == trim(seq_flds_scalar_name)) then
            call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(fldlist%shortname(n))//" is connected on root pe", &
              ESMF_LOGMSG_INFO, rc=dbrc)
            call shr_nuopc_fldList_SetScalarField(field, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          elseif (present(grid)) then
            call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(fldlist%shortname(n))//" is connected using grid", &
              ESMF_LOGMSG_INFO, rc=dbrc)
            field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=fldlist%shortname(n),rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          elseif (present(mesh)) then
            call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(fldlist%shortname(n))//" is connected using mesh", &
              ESMF_LOGMSG_INFO, rc=dbrc)
            field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=fldlist%shortname(n),meshloc=ESMF_MESHLOC_ELEMENT,rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          else
            call ESMF_LogWrite(trim(subname)//trim(tag)//": ERROR grid or mesh expected", ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif

          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!          call ESMF_FieldPrint(field=field, rc=rc)
!          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif

      else

        call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(fldlist%shortname(n)) // " is not connected.", &
          ESMF_LOGMSG_INFO, rc=dbrc)
        call ESMF_StateRemove(state, (/fldlist%shortname(n)/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      endif

    enddo

    call ESMF_LogWrite(subname//' done ', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine shr_nuopc_fldList_Realize

!================================================================================

  subroutine shr_nuopc_fldList_SetScalarField(field, rc)
    ! ----------------------------------------------
    ! create a field with scalar data on the root pe
    ! ----------------------------------------------
    type(ESMF_Field), intent(inout)  :: field
    integer,          intent(inout)  :: rc

    ! local variables
    type(ESMF_Distgrid) :: distgrid
    type(ESMF_Grid)     :: grid
    character(len=*), parameter :: subname='(shr_nuopc_fldList_SetScalarField)'

    rc = ESMF_SUCCESS

    ! create a DistGrid with a single index space element, which gets mapped
    ! onto DE 0.

    distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    grid = ESMF_GridCreate(distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    field = ESMF_FieldCreate(name=trim(seq_flds_scalar_name), &
         grid=grid, &
         typekind=ESMF_TYPEKIND_R8, &
         ungriddedLBound=(/1/), ungriddedUBound=(/seq_flds_scalar_num/), & ! num of scalar values
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine shr_nuopc_fldList_SetScalarField

!================================================================================

  subroutine shr_nuopc_fldList_CopyStateToScalar(State, data, mpicom, rc)
    ! ----------------------------------------------
    ! Copy scalar data from State to local data on root then broadcast
    ! ----------------------------------------------
    type(ESMF_State), intent(in)     :: State
    real(r8),         intent(inout)  :: data(:)
    integer,          intent(in)     :: mpicom
    integer,          intent(inout)  :: rc

    ! local variables
    integer :: mytask
    type(ESMF_Field) :: field
    real(ESMF_KIND_R8), pointer :: farrayptr(:)
    character(len=*), parameter :: subname='(shr_nuopc_fldList_CopyStateToScalar)'

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, mytask, rc)
    call ESMF_StateGet(State, itemName=trim(seq_flds_scalar_name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (size(data) < seq_flds_scalar_num .or. size(farrayptr) < seq_flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR on data size", ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      data(1:seq_flds_scalar_num) = farrayptr(1:seq_flds_scalar_num)
    endif
    call shr_mpi_bcast(data, mpicom)

  end subroutine shr_nuopc_fldList_CopyStateToScalar

!================================================================================

  subroutine shr_nuopc_fldList_CopyScalarToState(data, State, mpicom, rc)
    ! ----------------------------------------------
    ! Copy local scalar data into State, root only
    ! ----------------------------------------------
    real(r8),         intent(in)     :: data(:)
    type(ESMF_State), intent(inout)  :: State
    integer,          intent(in)     :: mpicom
    integer,          intent(inout)  :: rc

    ! local variables
    integer :: mytask
    type(ESMF_Field) :: field
    real(ESMF_KIND_R8), pointer :: farrayptr(:)
    character(len=*), parameter :: subname='(shr_nuopc_fldList_CopyScalarToState)'

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, mytask, rc)
    call ESMF_StateGet(State, itemName=trim(seq_flds_scalar_name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (size(data) < seq_flds_scalar_num .or. size(farrayptr) < seq_flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR on data size", ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      farrayptr(1:seq_flds_scalar_num) = data(1:seq_flds_scalar_num)
    endif

  end subroutine shr_nuopc_fldList_CopyScalarToState

!================================================================================

end module shr_nuopc_fldList_mod
