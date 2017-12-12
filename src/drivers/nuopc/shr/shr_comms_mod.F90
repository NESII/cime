module shr_comms_mod

  !---------------------------------------------------------------------
  !
  ! Purpose: Set up necessary communications
  !          Note that if no MPI, will call MCTs fake version
  !          (including mpif.h) will be utilized
  !
  !---------------------------------------------------------------------


  ! NOTE: If all atmospheres are identical in number of processes,
  ! number of threads, and grid layout, we should check that the
  ! user-provided number of processes and threads are consistent
  ! (or else, only accept one entry for these quantities when reading
  ! the namelist).  ARE OTHER PROTECTIONS/CHECKS NEEDED???

  use ESMF
  use mct_mod     , only : mct_world_init
  use shr_sys_mod , only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod , only : shr_mpi_chkerr, shr_mpi_bcast, shr_mpi_max
  use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
  use esmf        , only : ESMF_LogKind_Flag, ESMF_LOGKIND_NONE
  use esmf        , only : ESMF_LOGKIND_SINGLE, ESMF_LOGKIND_MULTI

  implicit none

  private
#include <mpif.h>

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public shr_comms_init
  public shr_comms_getinfo
  public shr_comms_setnthreads
  public shr_comms_getnthreads

  private comp_init
  private pelayout_init
  private mkname
  private printcomms
  private set_shr_comms

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  integer, public, parameter :: default_logunit = 6
  integer, public :: logunit  = default_logunit ! log unit number
  integer, public :: loglevel = 1               ! log level
  integer, public :: global_mype = -1           ! To be initialized

  ! Note - NUM_COMP_INST_XXX are cpp variables set in buildlib.csm_share

  integer, parameter :: ncomptypes = 8  ! total number of component types
  integer, parameter :: ncouplers  = 1  ! number of couplers

  integer, parameter, public :: num_inst_atm = NUM_COMP_INST_ATM
  integer, parameter, public :: num_inst_lnd = NUM_COMP_INST_LND
  integer, parameter, public :: num_inst_ocn = NUM_COMP_INST_OCN
  integer, parameter, public :: num_inst_ice = NUM_COMP_INST_ICE
  integer, parameter, public :: num_inst_glc = NUM_COMP_INST_GLC
  integer, parameter, public :: num_inst_wav = NUM_COMP_INST_WAV
  integer, parameter, public :: num_inst_rof = NUM_COMP_INST_ROF
  integer, parameter, public :: num_inst_esp = NUM_COMP_INST_ESP

  integer, parameter, public :: num_inst_total= num_inst_atm + &
                                num_inst_lnd + &
                                num_inst_ocn + &
                                num_inst_ice + &
                                num_inst_glc + &
                                num_inst_wav + &
                                num_inst_rof + &
                                num_inst_esp + 1

  integer, public :: num_inst_min, num_inst_max
  integer, public :: num_inst_xao    ! for xao flux
  integer, public :: num_inst_frc    ! for fractions
  integer, public :: num_inst_driver = 1

  integer, parameter, public :: &
       num_inst_phys = num_inst_atm + num_inst_lnd + &
                       num_inst_ocn + num_inst_ice + &
                       num_inst_glc + num_inst_rof + &
                       num_inst_wav + num_inst_esp

  integer, parameter, public :: &
       num_cpl_phys  = num_inst_atm + num_inst_lnd + &
                       num_inst_ocn + num_inst_ice + &
                       num_inst_glc + num_inst_rof + &
                       num_inst_wav + num_inst_esp

 !integer, parameter :: ncomps = (1 + ncouplers + 2*ncomptypes + num_inst_phys + num_cpl_phys)
  integer, parameter :: ncomps = 1 + ncouplers + num_inst_phys

  integer, public :: GLOID
  integer, public :: MEDID
  integer, public :: ATMID(num_inst_atm)
  integer, public :: LNDID(num_inst_lnd)
  integer, public :: OCNID(num_inst_ocn)
  integer, public :: ICEID(num_inst_ice)
  integer, public :: GLCID(num_inst_glc)
  integer, public :: ROFID(num_inst_rof)
  integer, public :: WAVID(num_inst_wav)
  integer, public :: ESPID(num_inst_esp)

  type(ESMF_LogKind_Flag), public :: esmf_logfile_kind

  integer, parameter, public :: shr_comms_namelen=16

  ! suffix for log and timing files if multi coupler driver
  character(len=shr_comms_namelen), public  :: cpl_inst_tag

  type shr_comms_type
     character(len=shr_comms_namelen) :: name     ! my name
     character(len=shr_comms_namelen) :: suffix   ! recommended suffix
     integer :: inst      ! my inst index
     integer :: ID        ! my id number
     integer :: mpicom    ! mpicom
     integer :: npes      ! number of mpi tasks in comm
     integer :: nthreads  ! number of omp threads per task
     integer :: iam       ! my task number in mpicom
     logical :: iamroot   ! am i the root task in mpicom
     integer :: gloiam    ! my task number in global_comm
     integer :: gloroot   ! the global task number of each comps root on all pes
     integer :: pethreads ! max number of threads on my task
     logical :: set       ! has this datatype been set
    integer, pointer :: petlist(:) ! esmf pet list
  end type shr_comms_type

  type(shr_comms_type) :: shr_comms(ncomps)

  character(*), parameter :: layout_concurrent = 'concurrent'
  character(*), parameter :: layout_sequential = 'sequential'

  character(*), parameter :: F11 = "(a,a,'(',i3,' ',a,')',a,   3i6,' (',a,i6,')',' (',a,i3,')','(',a,a,')')"
  character(*), parameter :: F12 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')','(',a,2i6,')')"
  character(*), parameter :: F13 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')')"
  character(*), parameter :: F14 = "(a,a,'(',i3,' ',a,')',a,    6x,' (',a,i6,')',' (',a,i3,')')"

  ! Exposed for use in the esp component, please don't use this elsewhere
  integer, public :: global_comm
  integer         :: driver_comm

  character(len=32), public :: &
       atm_layout, lnd_layout, ice_layout, glc_layout, rof_layout, &
       ocn_layout, wav_layout, esp_layout

  logical :: shr_comms_initialized = .false.  ! whether this module has been initialized

!=======================================================================
contains
!======================================================================

  !---------------------------------------------------------
  subroutine shr_comms_init(global_comm_in, driver_comm_in, nmlfile, &
       comp_id, comp_comm, comp_comm_iam, comp_iamin, comp_name, maxthreads, &
       drv_comm_id)

    ! Arguments
    integer                         , intent(in)           :: global_comm_in
    integer                         , intent(in)           :: driver_comm_in
    character(len=*)                , intent(in)           :: nmlfile
    integer                         , intent(out)          :: comp_id(num_inst_total)
    integer                         , intent(out)          :: comp_comm(num_inst_total)
    integer                         , intent(out)          :: comp_comm_iam(num_inst_total)
    logical                         , intent(out)          :: comp_iamin(num_inst_total)
    character(len=shr_comms_namelen) , intent(out)          :: comp_name(num_inst_total)
    integer                         , intent(out)          :: maxthreads
    integer                         , intent(in), optional :: drv_comm_id
    !
    ! Local variables
    logical                 :: error_state
    logical                 :: iamroot
    integer                 :: nu
    integer                 :: ierr, n, count, it
    integer                 :: mpicom_GLOID
    integer                 :: mype,numpes,myncomps,max_threads,gloroot, global_numpes
    integer                 :: pelist(3,1) ! start, stop, stride for group
    integer, pointer        :: comps(:)    ! array with component ids
    integer, pointer        :: comms(:)    ! array with mpicoms
    logical                 :: output_perf = .false.  ! require timing data output for this pe
    character(len=(shr_comms_namelen+1)*(num_inst_phys+1)) :: complist

    integer :: &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, &
         esp_ntasks, esp_rootpe, esp_pestride, esp_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads

    namelist /ccsm_pes/  &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, atm_layout, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, lnd_layout, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, ice_layout, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, glc_layout, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, wav_layout, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, rof_layout, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, ocn_layout, &
         esp_ntasks, esp_rootpe, esp_pestride, esp_nthreads, esp_layout, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads

    character(*), parameter :: u_FILE_u = __FILE__
    character(*), parameter :: subName =   '(shr_comms_init) '
    !----------------------------------------------------------

    ! make sure this is first pass and set comms unset
    if (shr_comms_initialized) then
       write(logunit,*) trim(subname),' ERROR shr_comms_init already called '
       call shr_sys_abort()
    endif
    shr_comms_initialized = .true.

    global_comm = global_comm_in
    driver_comm = driver_comm_in

    !-----------------------------------------
    ! Initialize shr_comms elements
    !-----------------------------------------

    do n = 1,ncomps
       shr_comms(n)%name = 'unknown'
       shr_comms(n)%suffix = ' '
       shr_comms(n)%inst = 0
       shr_comms(n)%mpicom = MPI_COMM_NULL    ! do some initialization here
       shr_comms(n)%iam = -1
       shr_comms(n)%iamroot = .false.
       shr_comms(n)%npes = -1
       shr_comms(n)%nthreads = -1
       shr_comms(n)%gloiam = -1
       shr_comms(n)%gloroot = -1
       shr_comms(n)%pethreads = -1
    enddo

    !-----------------------------------------
    ! Initialize MPI,  Note that if no MPI, will call ESMF serial version
    !-----------------------------------------

    call mpi_comm_size(GLOBAL_COMM_IN, global_numpes , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    call mpi_comm_rank(DRIVER_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank driver')

    call mpi_comm_size(DRIVER_COMM, numpes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size driver')

    if (mod(global_numpes, numpes) .ne. 0) then
       write(logunit,*) trim(subname),' ERROR: numpes driver: ', numpes, ' should divide global_numpes: ',global_numpes
       call shr_sys_abort(trim(subname)//' ERROR decomposition error ')
    endif

    !-----------------------------------------
    ! Initialize gloiam on all IDs
    !-----------------------------------------

    global_mype = mype

    do n = 1,ncomps
       shr_comms(n)%gloiam = mype
    enddo

    !-----------------------------------------
    ! Set ntasks, rootpe, pestride, nthreads for all components
    !-----------------------------------------

    if (mype == 0) then

       ! Set up default component process parameters
       call pelayout_init(numpes, atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, atm_layout)
       call pelayout_init(numpes, lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, lnd_layout)
       call pelayout_init(numpes, ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, ice_layout)
       call pelayout_init(numpes, ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, ocn_layout)
       call pelayout_init(numpes, rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, rof_layout)
       call pelayout_init(numpes, wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, wav_layout)
       call pelayout_init(numpes, glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, glc_layout)
       call pelayout_init(numpes, esp_ntasks, esp_rootpe, esp_pestride, esp_nthreads, esp_layout)
       call pelayout_init(numpes, cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads)

       ! Read namelist if it exists
       nu = shr_file_getUnit()
       open(nu, file=trim(nmlfile), status='old', iostat=ierr)
       if (ierr == 0) then
          ierr = 1
          do while( ierr > 0 )
             read(nu, nml=ccsm_pes, iostat=ierr)
          end do
          close(nu)
       end if
       call shr_file_freeUnit(nu)

    end if

    call shr_mpi_bcast(atm_nthreads, DRIVER_COMM, 'atm_nthreads')
    call shr_mpi_bcast(lnd_nthreads, DRIVER_COMM, 'lnd_nthreads')
    call shr_mpi_bcast(ocn_nthreads, DRIVER_COMM, 'ocn_nthreads')
    call shr_mpi_bcast(ice_nthreads, DRIVER_COMM, 'ice_nthreads')
    call shr_mpi_bcast(glc_nthreads, DRIVER_COMM, 'glc_nthreads')
    call shr_mpi_bcast(wav_nthreads, DRIVER_COMM, 'wav_nthreads')
    call shr_mpi_bcast(rof_nthreads, DRIVER_COMM, 'rof_nthreads')
    call shr_mpi_bcast(esp_nthreads, DRIVER_COMM, 'esp_nthreads')
    call shr_mpi_bcast(cpl_nthreads, DRIVER_COMM, 'cpl_nthreads')

    call shr_mpi_bcast(atm_layout, DRIVER_COMM, 'atm_layout')
    call shr_mpi_bcast(lnd_layout, DRIVER_COMM, 'lnd_layout')
    call shr_mpi_bcast(ocn_layout, DRIVER_COMM, 'ocn_layout')
    call shr_mpi_bcast(ice_layout, DRIVER_COMM, 'ice_layout')
    call shr_mpi_bcast(glc_layout, DRIVER_COMM, 'glc_layout')
    call shr_mpi_bcast(wav_layout, DRIVER_COMM, 'wav_layout')
    call shr_mpi_bcast(rof_layout, DRIVER_COMM, 'rof_layout')
    call shr_mpi_bcast(esp_layout, DRIVER_COMM, 'esp_layout')

    !-----------------------------------------
    ! compute some other num_inst values
    !-----------------------------------------

    num_inst_xao = max(num_inst_atm, num_inst_ocn)
    num_inst_frc = num_inst_ice

    !-----------------------------------------
    ! compute num_inst_min, num_inst_max
    ! instances must be either 1 or a constant across components
    ! checks for prognostic/present consistency in the driver
    !-----------------------------------------

    error_state = .false.

    num_inst_min = min(num_inst_atm, num_inst_lnd, num_inst_ocn,&
         num_inst_ice, num_inst_glc, num_inst_wav, num_inst_rof,&
         num_inst_esp)

    num_inst_max = max(num_inst_atm, num_inst_lnd, num_inst_ocn,&
         num_inst_ice, num_inst_glc, num_inst_wav, num_inst_rof,&
         num_inst_esp)

    if (num_inst_min /= num_inst_max .and. num_inst_min /= 1) error_state = .true.
    if (num_inst_atm /= num_inst_min .and. num_inst_atm /= num_inst_max) error_state = .true.
    if (num_inst_lnd /= num_inst_min .and. num_inst_lnd /= num_inst_max) error_state = .true.
    if (num_inst_ocn /= num_inst_min .and. num_inst_ocn /= num_inst_max) error_state = .true.
    if (num_inst_ice /= num_inst_min .and. num_inst_ice /= num_inst_max) error_state = .true.
    if (num_inst_glc /= num_inst_min .and. num_inst_glc /= num_inst_max) error_state = .true.
    if (num_inst_wav /= num_inst_min .and. num_inst_wav /= num_inst_max) error_state = .true.
    if (num_inst_rof /= num_inst_min .and. num_inst_rof /= num_inst_max) error_state = .true.
    if (num_inst_esp /= num_inst_min .and. num_inst_esp /= num_inst_max) error_state = .true.

    if (error_state) then
       write(logunit,*) trim(subname),' ERROR: num_inst inconsistent'
       write(logunit,*) num_inst_atm, num_inst_lnd, num_inst_ocn,&
                        num_inst_ice, num_inst_glc, num_inst_wav, num_inst_rof,&
                        num_inst_esp, num_inst_min, num_inst_max
       call shr_sys_abort(trim(subname)//' ERROR: num_inst inconsistent')
    endif

    !-----------------------------------------
    ! Initialize IDs
    !-----------------------------------------

    count = 0

    count = count + 1
    GLOID = count
    if (mype == 0) then
       pelist(1,1) = 0
       pelist(2,1) = numpes-1
       pelist(3,1) = 1
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
    call set_shr_comms(GLOID, pelist=pelist, iname='GLOBAL')

    count = count + 1
    MEDID = count
    if (mype == 0) then
       pelist(1,1) = cpl_rootpe
       pelist(2,1) = cpl_rootpe + (cpl_ntasks -1) * cpl_pestride
       pelist(3,1) = cpl_pestride
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
    call set_shr_comms(MEDID, pelist=pelist, nthreads=cpl_nthreads, iname='MED')

    call comp_init(driver_comm, atm_rootpe, atm_nthreads, atm_layout, atm_ntasks, atm_pestride, num_inst_atm, &
         ATMID, 'ATM', count, drv_comm_id)
    call comp_init(driver_comm, lnd_rootpe, lnd_nthreads, lnd_layout, lnd_ntasks, lnd_pestride, num_inst_lnd, &
         LNDID, 'LND', count, drv_comm_id)
    call comp_init(driver_comm, ice_rootpe, ice_nthreads, ice_layout, ice_ntasks, ice_pestride, num_inst_ice, &
         ICEID, 'ICE', count, drv_comm_id)
    call comp_init(driver_comm, ocn_rootpe, ocn_nthreads, ocn_layout, ocn_ntasks, ocn_pestride, num_inst_ocn, &
         OCNID, 'OCN', count, drv_comm_id)
    call comp_init(driver_comm, rof_rootpe, rof_nthreads, rof_layout, rof_ntasks, rof_pestride, num_inst_rof, &
         ROFID, 'ROF', count, drv_comm_id)
    call comp_init(driver_comm, glc_rootpe, glc_nthreads, glc_layout, glc_ntasks, glc_pestride, num_inst_glc, &
         GLCID, 'GLC', count, drv_comm_id)
    call comp_init(driver_comm, wav_rootpe, wav_nthreads, wav_layout, wav_ntasks, wav_pestride, num_inst_wav, &
         WAVID, 'WAV', count, drv_comm_id)
    call comp_init(driver_comm, esp_rootpe, esp_nthreads, esp_layout, esp_ntasks, esp_pestride, num_inst_esp, &
         ESPID, 'ESP', count, drv_comm_id)

    if (count /= ncomps) then
       write(logunit,*) trim(subname),' ERROR in ID count ',count,ncomps
       call shr_sys_abort(trim(subname)//' ERROR in ID count')
    endif

    !-----------------------------------------
    ! Count the total number of threads
    !-----------------------------------------

    max_threads = -1
    do n = 1,ncomps
       max_threads = max(max_threads,shr_comms(n)%nthreads)
    enddo
    do n = 1,ncomps
       shr_comms(n)%pethreads = max_threads
    enddo

    !-----------------------------------------
    ! Compute each components root pe global id and broadcast so all pes have info
    !-----------------------------------------

    do n = 1,ncomps
       gloroot = -999
       if (shr_comms(n)%iamroot) gloroot = shr_comms(n)%gloiam
       call shr_mpi_max(gloroot,shr_comms(n)%gloroot,DRIVER_COMM, trim(subname)//' gloroot',all=.true.)
    enddo

    !----------------------------------------------
    ! Initialize MCT
    !----------------------------------------------

    ! add up valid comps on local pe
    myncomps = 0
    do n = 1,ncomps
       if (shr_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
       endif
    enddo

    ! set comps and comms
    allocate(comps(myncomps), comms(myncomps), stat=ierr)
    if( ierr/=0) call shr_sys_abort(subName // 'allocate comps comms')

    myncomps = 0
    do n = 1,ncomps
       if (shr_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
          if (myncomps > size(comps)) then
             write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
             call shr_sys_abort()
          endif
          comps(myncomps) = shr_comms(n)%ID
          comms(myncomps) = shr_comms(n)%mpicom
       endif
    enddo

    if (myncomps /= size(comps)) then
       write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
       call shr_sys_abort()
    endif

    call mct_world_init(ncomps, DRIVER_COMM, comms, comps)
    deallocate(comps,comms)

    call printcomms()

    !--------------------------------------------------
    ! Set arrays comp_id, comp_iamin, comp_name
    ! Set output_perf and complist
    !--------------------------------------------------

    call shr_comms_getinfo(ID=GLOID, mpicom=mpicom_GLOID, iamroot=iamroot)
    if (iamroot) output_perf = .true.

    it=1
    call shr_comms_getinfo(ID=MEDID, mpicom=comp_comm(it), iamroot=iamroot, &
         iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
    if (comp_comm_iam(it)) complist = trim(complist)//' cpl'
    if (iamroot) output_perf = .true.

    do n = 1,num_inst_atm
       it=it+1
       comp_id(it) = ATMID(n)
       call shr_comms_getinfo(ID=ATMID(n), mpicom=comp_comm(it), iamroot=iamroot, &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(ATMID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_lnd
       it=it+1
       comp_id(it) = LNDID(n)
       call shr_comms_getinfo(ID=LNDID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(LNDID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_ocn
       it=it+1
       comp_id(it) = OCNID(n)
       call shr_comms_getinfo(ID=OCNID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(OCNID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_ice
       it=it+1
       comp_id(it) = ICEID(n)
       call shr_comms_getinfo(ID=ICEID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(ICEID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_glc
       it=it+1
       comp_id(it) = GLCID(n)
       call shr_comms_getinfo(ID=GLCID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(GLCID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_rof
       it=it+1
       comp_id(it) = ROFID(n)
       call shr_comms_getinfo(ID=ROFID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(ROFID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_wav
       it=it+1
       comp_id(it) = WAVID(n)
       call shr_comms_getinfo(ID=WAVID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(WAVID(n)))
       if (iamroot) output_perf = .true.
    enddo

    do n = 1,num_inst_esp
       it=it+1
       comp_id(it) = ESPID(n)
       call shr_comms_getinfo(ID=ESPID(n), mpicom=comp_comm(it), iamroot=iamroot,  &
            iam=comp_comm_iam(it), iamin=comp_iamin(it), name=comp_name(it))
       if (comp_comm_iam(it)) complist = trim(complist)//' '//trim(comp_name(ESPID(n)))
       if (iamroot) output_perf = .true.
    enddo

    ! ESP components do not use the coupler (they are 'external')

  end subroutine shr_comms_init

  !---------------------------------------------------------
  subroutine comp_init(driver_comm, comp_rootpe, comp_nthreads, comp_layout, &
       comp_ntasks, comp_pestride, num_inst_comp, COMPID, name, count, drv_comm_id)

    ! Arguments
    integer          , intent(in)           :: driver_comm
    integer          , intent(in)           :: comp_rootpe
    integer          , intent(in)           :: comp_nthreads
    character(len=*) , intent(in)           :: comp_layout
    integer          , intent(in)           :: comp_ntasks
    integer          , intent(in)           :: comp_pestride
    integer          , intent(in)           :: num_inst_comp
    integer          , intent(out)          :: COMPID(num_inst_comp)
    integer          , intent(inout)        :: count
    integer          , intent(in), optional :: drv_comm_id
    character(len=*) , intent(in)           :: name

    ! Local variables
    integer :: comp_inst_tasks
    integer :: ntasks, ntask, cnt
    integer :: droot
    integer :: current_task_rootpe
    integer :: cmin(num_inst_comp), cmax(num_inst_comp), cstr(num_inst_comp)
    integer :: n
    integer :: pelist (3,1)
    integer :: ierr
    integer :: mype
    character(len=*), parameter :: subname = "comp_init"
    !----------------------------------------------------------

    call mpi_comm_rank(driver_comm, mype, ierr)

    do n = 1, num_inst_comp
       count = count + 1
       COMPID(n) = count
    enddo

    if (mype == 0) then
       !--- validation of inputs ---
       ! rootpes >= 0
       ! Determine the process layout
       !
       ! We will assign comp_ntasks / num_inst_comp tasks to each component
       ! instance.  (This may lead to unallocated tasks if comp_ntasks is
       ! not an integer multiple of num_inst_comp.)

       if (comp_rootpe < 0) then
          call shr_sys_abort(trim(subname)//' ERROR: rootpes must be >= 0 for component '//trim(name))
       endif

       if (trim(comp_layout) == trim(layout_concurrent)) then
          comp_inst_tasks = comp_ntasks / num_inst_comp
          droot = (comp_inst_tasks * comp_pestride)
       elseif (trim(comp_layout) == trim(layout_sequential)) then
          comp_inst_tasks = comp_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid comp_layout for component '//trim(name))
       endif

       current_task_rootpe = comp_rootpe
       do n = 1, num_inst_comp
          cmin(n) = current_task_rootpe
          cmax(n) = current_task_rootpe + ((comp_inst_tasks - 1) * comp_pestride)
          cstr(n) = comp_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do
    endif

    do n = 1, num_inst_comp
       if (mype==0) then
          pelist(1,1) = cmin(n)
          pelist(2,1) = cmax(n)
          pelist(3,1) = cstr(n)
       endif
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, DRIVER_COMM, ierr)

       ntasks = ((pelist(2,1) - pelist(1,1)) / pelist(3,1)) + 1
       allocate(shr_comms(COMPID(n))%petlist(ntasks))
       cnt = 0
       do ntask = pelist(1,1),pelist(2,1),pelist(3,1)
          cnt = cnt + 1
          if (cnt > ntasks) then
             write(logunit,*) subname,' ERROR in petlist init ',ntasks,pelist(1:3,1),ntask,cnt
             call shr_sys_abort(subname//' ERROR in petlist init')
          endif
          shr_comms(COMPID(n))%petlist(cnt) = ntask
       enddo

       if (present(drv_comm_id)) then
          print *,__FILE__,__LINE__,drv_comm_id
          call set_shr_comms(COMPID(n), pelist, comp_nthreads, name, drv_comm_id)
       else
          call set_shr_comms(COMPID(n), pelist, comp_nthreads, name, n, num_inst_comp)
       endif
    enddo

  end subroutine comp_init

  !---------------------------------------------------------
  subroutine pelayout_init(numpes, ntasks, rootpe, pestride, nthreads, layout)

    ! Arguments
    integer          , intent(in)            :: numpes
    integer          , intent(out)           :: ntasks
    integer          , intent(out)           :: rootpe
    integer          , intent(out)           :: pestride
    integer          , intent(out)           :: nthreads
    character(len=*) , intent(out), optional :: layout
    !----------------------------------------------------------

    ntasks = numpes
    rootpe = 0
    pestride = 1
    nthreads = 1
    if (present(layout)) then
       layout = trim(layout_concurrent)
    endif
  end subroutine pelayout_init

  !---------------------------------------------------------
  subroutine set_shr_comms(ID, pelist, nthreads, iname, inst, tinst)

    ! Arguments
    integer          ,intent(IN)          :: ID
    integer          ,intent(IN)          :: pelist(:,:)
    integer          ,intent(IN),optional :: nthreads
    character(len=*) ,intent(IN),optional :: iname ! name of component
    integer          ,intent(IN),optional :: inst  ! instance of component
    integer          ,intent(IN),optional :: tinst ! total number of instances for this component

    ! Local variables
    integer :: mpigrp_world
    integer :: mpigrp
    integer :: mpicom
    integer :: ierr
    logical :: set_suffix
    character(len=shr_comms_namelen) :: cname
    character(*),parameter :: subName =   '(set_shr_comms) '
    !---------------------------------------------

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif

    ! TODO: Can the following be removed in place of ESMF VM queries?
    ! The arrays generated by shr_comms_init are needed for initializing PIO
    ! But do we need the PIO initialization 2 with ESMF?

    call mpi_comm_group(driver_comm, mpigrp_world, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')

    call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')

    call mpi_comm_create(driver_comm, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    shr_comms(ID)%ID = ID

    if (present(inst)) then
       shr_comms(ID)%inst = inst
       set_suffix = .true.
    else
       shr_comms(ID)%inst = 1
       set_suffix = .false.
    endif

    if (present(tinst)) then
       if (tinst == 1) set_suffix = .false.
    endif

    if (present(iname)) then
       shr_comms(ID)%name = trim(iname)
       if (set_suffix) then
          call mkname(cname,iname,shr_comms(ID)%inst)
          shr_comms(ID)%name = trim(cname)
       endif
    endif
    write(6,*)'DEBUG: shr_comms(ID)%name = ',shr_comms(ID)%name

    if (set_suffix) then
       call mkname(cname, '_', shr_comms(ID)%inst)
       shr_comms(ID)%suffix = trim(cname)
    else
       shr_comms(ID)%suffix = ' '
    endif

    shr_comms(ID)%mpicom = mpicom
    if (present(nthreads)) then
       shr_comms(ID)%nthreads = nthreads
    else
       shr_comms(ID)%nthreads = 1
    endif

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom, shr_comms(ID)%npes, ierr)
       call shr_mpi_chkerr(ierr, subname//' mpi_comm_size')

       call mpi_comm_rank(mpicom, shr_comms(ID)%iam, ierr)
       call shr_mpi_chkerr(ierr, subname//' mpi_comm_rank')

       if (shr_comms(ID)%iam == 0) then
          shr_comms(ID)%iamroot = .true.
       else
          shr_comms(ID)%iamroot = .false.
       endif
    else
       shr_comms(ID)%npes = -1
       shr_comms(ID)%iam = -1
       shr_comms(ID)%nthreads = 1
       shr_comms(ID)%iamroot = .false.
    endif

    if (shr_comms(ID)%iamroot) then
       write(logunit, F11) trim(subname), '  initialize ID ', ID, shr_comms(ID)%name, &
            ' pelist   =', pelist, ' npes =', shr_comms(ID)%npes, ' nthreads =', shr_comms(ID)%nthreads, &
            ' suffix =', trim(shr_comms(ID)%suffix)
    endif

  end subroutine set_shr_comms

  !---------------------------------------------------------
  subroutine printcomms()

    ! Local variables
    integer :: n, mype, npes, ierr
    character(*), parameter :: subName =   '(printcomms) '
    !----------------------------------------------------------

    call mpi_comm_size(DRIVER_COMM, npes  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    call mpi_comm_rank(DRIVER_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    call shr_sys_flush(logunit)
    call mpi_barrier(DRIVER_COMM,ierr)
    if (mype == 0) then
       do n = 1,ncomps
          write(logunit,'(a,4i6,2x,3a)') trim(subName),n, &
               shr_comms(n)%gloroot,shr_comms(n)%npes,shr_comms(n)%nthreads, &
               trim(shr_comms(n)%name),':',trim(shr_comms(n)%suffix)
       enddo
       call shr_sys_flush(logunit)
    endif
  end subroutine printcomms

  !---------------------------------------------------------
  subroutine shr_comms_getinfo(ID, mpicom, nthreads, iam, iamroot, iamin, &
       gloiam, gloroot, pethreads, inst, name, suffix, petlist)

    ! Arguments
    integer ,intent(in)            :: ID
    integer ,intent(out), optional :: mpicom
    integer ,intent(out), optional :: nthreads
    integer ,intent(out), optional :: iam
    logical ,intent(out), optional :: iamroot
    logical ,intent(out), optional :: iamin
    integer ,intent(out), optional :: gloiam
    integer ,intent(out), optional :: gloroot
    integer ,intent(out), optional :: pethreads
    integer ,intent(out), optional :: inst
    character(len=shr_comms_namelen) , intent(out), optional :: name
    character(len=shr_comms_namelen) , intent(out), optional :: suffix
    integer , pointer   , optional :: petlist(:)

    ! Local variables
    character(*),parameter :: subName =   '(shr_comms_getinfo) '
    !---------------------------------------------

    ! Negative ID means there is no comm, return default or inactive values
    if ((ID == 0) .or. (ID > ncomps)) then
       write(logunit,*) subname,' ID out of range, return ',ID
       return
    endif
    if (present(mpicom)) then
       if (ID > 0) then
          mpicom = shr_comms(ID)%mpicom
       else
          mpicom = MPI_COMM_NULL
       end if
    endif
    if (present(nthreads)) then
       if (ID > 0) then
          nthreads = shr_comms(ID)%nthreads
       else
          nthreads = 1
       end if
    endif
    if (present(iam)) then
       if (ID > 0) then
          iam = shr_comms(ID)%iam
       else
          iam = -1
       end if
    endif
    if (present(iamroot)) then
       if (ID > 0) then
          iamroot = shr_comms(ID)%iamroot
       else
          iamroot = .false.
       end if
    endif
    if (present(iamin)) then
       if ((ID < 1)) then
          iamin = .false.
       else if (shr_comms(ID)%iam >= 0) then
          iamin = .true.
       else
          iamin = .false.
       endif
    end if
    if (present(gloiam)) then
       if (ID > 0) then
          gloiam = shr_comms(ID)%gloiam
       else
          gloiam = -1
       end if
    endif
    if (present(gloroot)) then
       if (ID > 0) then
          gloroot = shr_comms(ID)%gloroot
       else
          gloroot = -1
       end if
    endif
    if (present(pethreads)) then
       if (ID > 0) then
          pethreads = shr_comms(ID)%pethreads
       else
          pethreads = 1
       end if
    endif
    if(present(name)) then
       if (ID > 0) then
          name = trim(shr_comms(ID)%name)
       else
          name = ''
       end if
    end if
    if(present(inst)) then
       if (ID > 0) then
          inst = shr_comms(ID)%inst
       else
          inst = 0
       end if
    end if
    if (present(suffix)) then
       if (ID > 0) then
          suffix = trim(shr_comms(ID)%suffix)
       else
          suffix = ''
       end if
    end if
    if (present(petlist)) then
       if (ID > 0) then
          petlist => shr_comms(ID)%petlist
       else
          nullify(petlist)
       end if
    end if

  end subroutine shr_comms_getinfo

  !---------------------------------------------------------
  subroutine shr_comms_setnthreads(nthreads)

    integer,intent(in) :: nthreads
    character(*),parameter :: subName = '(shr_comms_setnthreads) '
    !---------------------------------------------

#ifdef _OPENMP
    if (nthreads < 1) then
       call shr_sys_abort(subname//' ERROR: nthreads less than one')
    endif
    call omp_set_num_threads(nthreads)
#endif
  end subroutine shr_comms_setnthreads

  !---------------------------------------------------------
  integer function shr_comms_getnthreads()

    integer :: omp_get_num_threads
    character(*),parameter :: subName =   '(shr_comms_getnthreads) '
    !---------------------------------------------

    shr_comms_getnthreads = -1
#ifdef _OPENMP
    !$OMP PARALLEL
    shr_comms_getnthreads = omp_get_num_threads()
    !$OMP END PARALLEL
#endif
  end function shr_comms_getnthreads

  !---------------------------------------------------------
  subroutine mkname(oname,str1,num)

    ! Arguments
    character(len=*),intent(out) :: oname
    character(len=*),intent(in)  :: str1
    integer,intent(in)           :: num

    ! Local variables
    character(len=8) :: cnum
    character(*),parameter :: subName = '(mkname) '
    !---------------------------------------------

    write(cnum,'(i4.4)') num
    if (len_trim(str1) + len_trim(cnum) > len(oname)) then
       write(logunit,*) trim(subname),' ERROR in str lens ',len(oname),trim(str1),trim(cnum)
       call shr_sys_abort(trim(subname))
    endif
    oname = trim(str1)//trim(cnum)
  end subroutine mkname

end module shr_comms_mod
