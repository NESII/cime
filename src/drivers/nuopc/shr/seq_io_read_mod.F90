!===============================================================================
! !MODULE: seq_io_read_mod -- reads integer, real arrays and chacter of driver files
!
! Current Problems
!  - the original use of seq_io will now ONLY work with the cpl because
!    of hardwiring cpl_io_type and cpl_io_iosystem.  want the original
!    io capabilities to be usable by any component
!  - the init1 method depends on shr_comms for name consistency but shr_comms_init
!    wants to be called after init1 so the global_comm can be modified for
!    async IO.  this needs to be reconciled.
!  - this routine stores information for all components but most methods are
!    hardwired to work only for the coupler.  should all the components info
!    be stored here or should this be more a general set of methods that are
!    reusable as it's original intent.

module seq_io_read_mod

  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in => shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl, cs => shr_kind_cs
  use shr_pio_mod,  only: shr_pio_getiosys, shr_pio_getiotype
  use shr_sys_mod       ! system calls
  use shr_comms_mod
  use mct_mod           ! mct wrappers
  use pio

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public seq_io_read

  interface seq_io_read
     module procedure seq_io_read_int
     module procedure seq_io_read_int1d
     module procedure seq_io_read_r8
     module procedure seq_io_read_r81d
     module procedure seq_io_read_char
  end interface

  !-------------------------------------------------------------------------------
  ! Local data
  !-------------------------------------------------------------------------------
   character(*) , parameter :: prefix = "seq_io_"
   character(*) , parameter :: version ='cpl7v10'
   character(*) , parameter :: version0='cpl7v00'
   character(CL)            :: charvar   ! buffer for string read/write

!=================================================================================
contains
!=================================================================================

  !===============================================================================
  subroutine seq_io_read_int(filename,pioid,idata,dname)

    ! !DESCRIPTION: Read scalar integer from netcdf file

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    integer         ,intent(inout) :: idata    ! integer data
    character(len=*),intent(in)    :: dname    ! name of data

    integer :: i1d(1)
    character(*),parameter :: subName = '(seq_io_read_int) '
    !-------------------------------------------------------------------------------

    call seq_io_read_int1d(filename,pioid,i1d,dname)
    idata = i1d(1)
  end subroutine seq_io_read_int

  !===============================================================================
  subroutine seq_io_read_int1d(filename,pioid,idata,dname)

    ! !IROUTINE: seq_io_read_int1d - read 1d integer from netcdf file

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*),intent(in)   :: filename ! file
    type(file_desc_t)             :: pioid
    integer(in)     ,intent(inout):: idata(:)  ! integer data
    character(len=*),intent(in)   :: dname    ! name of data

    integer(in)                     :: rcode
    type(var_desc_t)                :: varid
    logical                         :: exists
    character(CL)                   :: name1
    character(*),parameter          :: subName = '(seq_io_read_int1d) '
    logical :: addprefix
    !-------------------------------------------------------------------------------
    call seq_io_read_openfile(filename,pioid,addprefix)

    if (addprefix) then
       name1 = trim(prefix)//trim(dname)
    else
       name1 = trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,idata)
  end subroutine seq_io_read_int1d

  !===============================================================================
  subroutine seq_io_read_openfile(filename,pioid,addprefix)

    character(len=*) , intent(in)  :: filename
    type(file_desc_t)              :: pioid
    logical          , intent(out) :: addprefix

    logical                          :: exists
    integer(in)                      :: iam,mpicom
    type(iosystem_desc_t) , pointer  :: cpl_io_subsystem
    character(len=shr_comms_namelen) :: cpl_name
    integer(in)                      :: cpl_pio_iotype
    logical, save                    :: laddprefix
    integer                          :: rcode
    character(CL)                    :: lversion
    character(*),parameter           :: subName = '(seq_io_read_openfile) '

    if(.not. pio_file_is_open(pioid)) then
       cpl_io_subsystem => shr_pio_getiosys(cpl_name)
       cpl_pio_iotype   =  shr_pio_getiotype(cpl_name)

       call shr_comms_getinfo(MEDID, iam=iam, mpicom=mpicom, name=cpl_name)
       if (iam==0) inquire(file=trim(filename),exist=exists)

       call shr_mpi_bcast(exists,mpicom,'seq_io_read_openfile')
       if (exists) then
          rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
          call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
          rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
          call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
          if (trim(lversion) == trim(version)) then
             laddprefix=.false.
          else
             laddprefix=.true.
          endif
       else
          if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename)
          call shr_sys_abort()
       endif
    endif
    addprefix = laddprefix

  end subroutine seq_io_read_openfile

  !===============================================================================
  subroutine seq_io_read_r8(filename,pioid,rdata,dname)

    ! !DESCRIPTION: Read scalar double from netcdf file

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    real(r8)        ,intent(inout) :: rdata    ! real data
    character(len=*),intent(in)    :: dname    ! name of data

    real(r8) :: r1d(1)
    character(*),parameter :: subName = '(seq_io_read_r8) '
    !-------------------------------------------------------------------------------

    call seq_io_read_r81d(filename,pioid,r1d,dname)
    rdata = r1d(1)
  end subroutine seq_io_read_r8

  !===============================================================================
  subroutine seq_io_read_r81d(filename,pioid,rdata,dname)

    ! !DESCRIPTION: Read 1d double array from netcdf file

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    real(r8)        ,intent(inout) :: rdata(:) ! real data
    character(len=*),intent(in)    :: dname    ! name of data

    type(var_desc_t)                :: varid
    character(CL)                   :: name1
    character(*),parameter          :: subName = '(seq_io_read_r81d) '
    logical :: addprefix
    integer :: rcode
    !-------------------------------------------------------------------------------

    call seq_io_read_openfile(filename,pioid,addprefix)

    if (addprefix) then
       name1 = trim(prefix)//trim(dname)
    else
       name1 = trim(dname)
    endif

    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,rdata)

  end subroutine seq_io_read_r81d

  !===============================================================================
  subroutine seq_io_read_char(filename,pioid,rdata,dname)

    ! !DESCRIPTION: Read char string from netcdf file

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*),intent(in)    :: filename ! file
    type(file_desc_t) :: pioid
    character(len=*),intent(inout) :: rdata    ! character data
    character(len=*),intent(in)    :: dname    ! name of data

    type(var_desc_t)                :: varid
    character(CL)                   :: name1
    character(*),parameter          :: subName = '(seq_io_read_char) '
    logical :: addprefix
    integer :: rcode
    !-------------------------------------------------------------------------------

    call seq_io_read_openfile(filename,pioid,addprefix)

    if (addprefix) then
       name1 = trim(prefix)//trim(dname)
    else
       name1 = trim(dname)
    endif

    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,charvar)
    rdata = trim(charvar)

  end subroutine seq_io_read_char

end module seq_io_read_mod
