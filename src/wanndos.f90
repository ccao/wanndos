PROGRAM wanndos
  !
  use para,     only : init_para, inode, distribute_k, finalize_para, first_k, last_k, para_merge
  use wanndata, only : read_ham, norb, finalize_wann
  use dosdata
  use input,    only : read_input, seed, lpdos
  !
  implicit none
  !
  integer ik
  !
  CALL init_para
  CALL read_input
  CALL read_ham(seed)
  !
  nproj=norb
  !
  CALL init_dos
  CALL distribute_k(nkpt)
  !
  do ik=first_k, last_k
    !
    CALL interpolate_bands(kvec(:, ik))
    !
    CALL calc_dos_contrib
    !
  enddo
  !
  CALL para_merge(dos, nedos, 2)
  !
  dos(:,:)=dos(:,:)/nkpt
  !
  if (inode.eq.0) CALL search_fermi
  !
  if (lpdos) then
    CALL para_merge(pdos, nedos, nproj)
    pdos(:,:)=pdos(:,:)/nkpt
  endif
  !
  if (inode.eq.0) CALL output_dos
  !
  CALL finalize_wann
  !
  CALL finalize_dos
  !
  CALL finalize_para
  !
END PROGRAM
