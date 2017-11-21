module med_internalstate_mod

  !-----------------------------------------------------------------------------
  ! Mediator Internal State Datatype.
  !-----------------------------------------------------------------------------

  use ESMF
  use shr_nuopc_fldList_mod, only: shr_nuopc_fldList_type

  implicit none

  public

  !--- Component arrays ---

  ! This defines the components and med_mapping_allowed is a starting point for what is
  ! allowed in this coupled system.  It will be revised further after the system
  ! starts, but any coupling set to false will never be allowed.  As new connections
  ! are allowed, just update the table below.  The rows are the destination of
  ! coupling, the columns are the source of coupling.  So, the second column
  ! indicates which models the atm is coupled to.  And the second row indicates
  ! which models are coupled to the atm.
  ! The mediator is not connected to any components in CESM because the mediator
  ! doesn't have it's own grid and only acts as a hub.

  integer, parameter :: ncomps=8
  integer, parameter :: compmed=1
  integer, parameter :: compatm=2
  integer, parameter :: complnd=3
  integer, parameter :: compocn=4
  integer, parameter :: compice=5
  integer, parameter :: comprof=6
  integer, parameter :: compwav=7
  integer, parameter :: compglc=8

  character(len=*),parameter :: compname(ncomps) = &
   (/ 'med', &
      'atm', &
      'lnd', &
      'ocn', &
      'ice', &
      'rof', &
      'wav', &
      'glc' /)

  ! tcraig, turned off glc2ocn and glc2ice for time being
  logical, parameter :: med_coupling_allowed(ncomps,ncomps) = &
   (/ .false., .false., .false., .false., .false., .false., .false., .false., &  ! med
      .false., .false., .true. , .true. , .true. , .false., .false., .false., &  ! atm
      .false., .true. , .false., .false., .false., .true. , .false., .true. , &  ! lnd
      .false., .true. , .false., .false., .true. , .true. , .true. , .false., &  ! ocn
      .false., .true. , .false., .true. , .false., .true. , .false., .false., &  ! ice
      .false., .false., .true. , .false., .false., .false., .false., .false., &  ! rof
      .false., .true. , .false., .true. , .true. , .false., .false., .false., &  ! wav
      .false., .false., .true. , .false., .false., .false., .false., .false.  /) ! glc
!       med      atm      lnd      ocn      ice      rof      wav      glc

  !--- RH arrays ---
  integer, parameter :: nmappers=5
  integer, parameter :: mapbilnr=1
  integer, parameter :: mapconsf=2
  integer, parameter :: mapconsd=3
  integer, parameter :: mappatch=4
  integer, parameter :: mapfcopy=5

  !--- fld lists ---
  type (shr_nuopc_fldList_Type) :: fldsTo(ncomps)
  type (shr_nuopc_fldList_Type) :: fldsFr(ncomps)
  type (shr_nuopc_fldList_Type) :: fldsAtmOcn

  ! private internal state to keep instance data
  type InternalStateStruct

    ! NState_Imp and NState_Exp are the standard NUOPC coupling datatypes
    ! FBImp and FBExp are the internal mediator datatypes
    ! NState_Exp(n) = FBExp(n), copied in the connector prep phase
    ! FBImp(n,n) = NState_Imp(n), copied in connector post phase
    ! FBImp(n,k) is the FBImp(n,n) interpolated to grid k
    ! RH(n,k,m) is a RH from grid n to grid k, map type m

    logical               :: comp_present(ncomps)               ! comp present flag
    logical               :: med_coupling_active(ncomps,ncomps) ! computes the active coupling
    type(ESMF_State)      :: NStateImp(ncomps)                  ! Import data from various component, on their grid
    type(ESMF_State)      :: NStateExp(ncomps)                  ! Export data to various component, on their grid
    type(ESMF_FieldBundle):: FBImp(ncomps,ncomps)               ! Import data from various components interpolated to various grids
    type(ESMF_FieldBundle):: FBImpAccum(ncomps)                 ! Accumulator for various components import on their grid
    integer               :: FBImpAccumcnt(ncomps)              ! Accumulator counter for each FBImpAccum
    type(ESMF_FieldBundle):: FBExp(ncomps)                      ! Export data for various components, on their grid
    type(ESMF_FieldBundle):: FBExpAccum(ncomps)                 ! Accumulator for various components export on their grid
    integer               :: FBExpAccumcnt(ncomps)              ! Accumulator counter for each FBExpAccum
    integer               :: conn_prep_cnt(ncomps)              ! Connector prep count
    integer               :: conn_post_cnt(ncomps)              ! Connector post count
    type(ESMF_FieldBundle):: FBfrac(ncomps)                     ! Fraction data for various components, on their grid
    type(ESMF_RouteHandle):: RH(ncomps,ncomps,nmappers)         ! Routehandles for pairs of components and different mappers

    type(ESMF_RouteHandle):: RH_r2ol_consf                      ! rof to ocn liquid
    type(ESMF_RouteHandle):: RH_r2oi_consf                      ! rof to ocn frozen
    type(ESMF_FieldBundle):: FBAtmOcn_o                         ! Atm/Ocn flux fields on ocn grid
    type(ESMF_FieldBundle):: FBAtmOcn_a                         ! Atm/Ocn flux fields on atm grid
    type(ESMF_FieldBundle):: FBaccumAtmOcn                      ! accumulator of atm export data
    integer               :: mpicom
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type

  !-----------------------------------------------------------------------------

end module med_internalstate_mod
