module startup_log

  implicit none
  character(len=*), parameter :: STARTUP_SEPARATOR = &
       '**********************************************************************'

contains

  logical function is_stage_child()
    implicit none
    character(len=16) :: child_flag
    integer :: stat

    child_flag=''
    call get_environment_variable('STRATA_STAGE_CHILD',child_flag,status=stat)
    is_stage_child = stat.eq.0 .and. len_trim(child_flag).gt.0 .and. &
         trim(adjustl(child_flag)).ne.'0'
  end function is_stage_child


  subroutine print_title_card
    implicit none

    if (is_stage_child()) return

    write(*,'(A)') ''
    write(*,'(A)') '****  ***** ****   ***  *****  ***    *   *  ***  *   *  *** '
    write(*,'(A)') '*       *   *   * *   *   *   *   *   **  * *   * *   * *   *'
    write(*,'(A)') '****    *   ****  *****   *   *****   * * * *   * *   * *****'
    write(*,'(A)') '    *   *   * *   *   *   *   *   *   *  ** *   *  * *  *   *'
    write(*,'(A)') '****    *   *  ** *   *   *   *   *   *   *  ***    *   *   *'
    write(*,'(A)') STARTUP_SEPARATOR
    write(*,'(A)') ''
  end subroutine print_title_card


  subroutine print_startup_separator
    implicit none
    write(*,'(A)') STARTUP_SEPARATOR
  end subroutine print_startup_separator


  subroutine print_stage_banner(label,start_day,end_day)
    implicit none
    character(len=*), intent(in) :: label
    real*8, intent(in) :: start_day,end_day

    write(*,'(A)') ''
    write(*,'(A)') STARTUP_SEPARATOR
    write(*,'(A,A)') ' STARTING NEW SIMULATION: ',trim(label)
    write(*,'(A,F10.5,A,F10.5,A)') ' Absolute time: ',start_day,' d  -> ',end_day,' d'
    write(*,'(A)') STARTUP_SEPARATOR
  end subroutine print_stage_banner

end module startup_log
