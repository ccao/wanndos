MODULE input
  !
  use constants
  !
  implicit none
  !
  real(dp) sigma
  real(dp) netot
  character(len=80) seed
  logical lpdos
  real(dp) emin, emax, de
  integer nkx, nky, nkz
  !
CONTAINS
  !
 SUBROUTINE read_input
  !
  !************ INPUT FILE *************
  !** file name: wanndos.inp
  !line 1: seed name
  !line 2: nqx nqy nqz
  !line 3: sigma
  !line 4: emin emax de
  !line 5: netot
  !line 6: pdos (?)
  !*************************************
  !
  use constants, only : fin
  use para
  !
  implicit none
  !
  integer mode
  !
  if (inode.eq.0) then
    open(unit=fin, file="wanndos.inp")
    !
    read(fin, *) seed
    !
    read(fin, *) nkx, nky, nkz
    read(fin, *) sigma
    read(fin, *) emin, emax, de
    read(fin, *) netot
    read(fin, *) mode
    !
  endif
  !
  CALL para_sync(nkx)
  CALL para_sync(nky)
  CALL para_sync(nkz)
  CALL para_sync(sigma)
  CALL para_sync(emin)
  CALL para_sync(emax)
  CALL para_sync(de)
  CALL para_sync(netot)
  CALL para_sync(mode)
  !
  if (mode.eq.0) then 
    lpdos = .false.
  else
    lpdos = .true.
  endif
  !
 END SUBROUTINE
  !
END MODULE
