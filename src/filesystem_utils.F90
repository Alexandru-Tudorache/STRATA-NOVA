module filesystem_utils
  implicit none
  private
  public :: shell_quote, reset_output_directory, check_output_directory
contains
  ! quote filesystem paths before passing them to the shell
  function shell_quote(text) result(quoted)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: quoted
    integer :: i
    quoted = "'"
    do i=1,len_trim(text)
       if (text(i:i).eq."'") then
          quoted = quoted // "'" // achar(34) // "'" // achar(34) // "'"
       else
          quoted = quoted // text(i:i)
       endif
    enddo
    quoted = quoted // "'"
  end function shell_quote

  subroutine check_output_directory(path, protected_paths)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: protected_paths(:)
    character(len=:), allocatable :: command
    character(len=1024) :: executable
    integer :: stat, i

    if (len_trim(path).eq.0) error stop 'Empty output path'
    call get_command_argument(0,executable,status=stat)
    if (stat.ne.0) error stop 'Cannot resolve executable path'
    command='p=$(realpath -m -- '//shell_quote(trim(path))//') || exit 2; '// &
         'cwd=$(realpath -e -- .) || exit 2; '// &
         'test "$p" != "/" || exit 4; '// &
         'for tree in tables src .git; do q=$(realpath -m -- "$tree") || exit 2; '// &
         'case "$p/" in "$q/"*) exit 5;; esac; done; '// &
         'test ! -L '//shell_quote(trim(path))//' || exit 3; '// &
         'case "$cwd/" in "$p/"*) exit 4;; esac; '// &
         'for input in '//shell_quote(trim(executable))//' tables src .git'
    do i=1,size(protected_paths)
       if (len_trim(protected_paths(i)).gt.0) command=command//' '//shell_quote(trim(protected_paths(i)))
    enddo
    command=command//'; do q=$(realpath -m -- "$input") || exit 2; '// &
         'case "$q/" in "$p/"*) exit 5;; esac; done'
    call execute_command_line('bash -c '//shell_quote(command),wait=.true.,exitstat=stat)
    if (stat.ne.0) then
       write(*,*) 'ERROR: refusing unsafe output directory: ',trim(path)
       stop 1
    endif
  end subroutine check_output_directory

  subroutine reset_output_directory(path, protected_paths)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: protected_paths(:)
    integer :: stat

    call check_output_directory(path,protected_paths)
    call execute_command_line('find '//shell_quote(trim(path))// &
         ' -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +',wait=.true.,exitstat=stat)
    if (stat.ne.0) error stop 'Could not reset output directory'
  end subroutine reset_output_directory
end module filesystem_utils
