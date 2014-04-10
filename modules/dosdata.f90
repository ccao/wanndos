MODULE dosdata
  !
  use constants
  !
  implicit none
  !
  real(dp) ef
  real(dp), allocatable :: emesh(:)
  real(dp), allocatable :: dos(:, :)
  real(dp), allocatable :: pdos(:, :)
  integer nedos, nproj
  !
  integer nkpt
  real(dp), allocatable :: kvec(:, :)
  !
  real(dp), allocatable :: eig(:)
  complex(dp), allocatable :: work(:, :)
  !
CONTAINS

SUBROUTINE search_fermi
  !
  use constants,  only : dp, stdout
  use input,      only : netot
  !
  implicit none
  !
  integer ien
  real(dp) x
  !
  if ((dos(nedos, 2).lt.netot).or.(dos(1, 2).gt.netot)) then
    write(stdout, *) " # Fermi level cannot be determined !"
  else
    do ien=1,nedos-1
      if((dos(ien, 2).le.netot).and.(dos(ien+1, 2).ge.netot)) then
        x=(netot-dos(ien,2))/(dos(ien+1,2)-dos(ien,2))
        ef=emesh(ien)+x*(emesh(ien+1)-emesh(ien))
        write(stdout, '(A,1F22.16)') " # Fermi level :", ef
        exit
      endif
    enddo
  endif
  !
END SUBROUTINE

SUBROUTINE calc_dos_contrib
  !
  use constants,  only : dp
  use input,      only : sigma, lpdos
  !
  implicit none
  !
  real(dp) x
  integer ien, ii, iproj
  real(dp) gaussian, erf_2
  !
  do ien=1, nedos
    do ii=1, nproj
      x=(emesh(ien)-eig(ii))/sigma
      if (abs(x).lt.40) dos(ien, 1) = dos(ien, 1)+2.0*gaussian(x)/sigma
      dos(ien, 2) = dos(ien, 2)+erf(x)+1.0
      if(lpdos) then
        if (abs(x).lt.40) then
          do iproj=1, nproj
            pdos(ien, iproj) = pdos(ien, iproj) + 2.0*gaussian(x)*real(work(iproj, ii)*conjg(work(iproj,ii)))/sigma
          enddo
        endif
      endif
    enddo
  enddo
  !
END SUBROUTINE

SUBROUTINE init_dos
  !
  use input,  only : lpdos, emin, emax, de, nkx, nky, nkz
  !
  implicit none
  !
  integer ikx, iky, ikz, ii
  !
  nedos=int((emax-emin)/de)+1
  de=(emax-emin)/(nedos-1)
  !
  allocate(emesh(1:nedos))
  !
  do ii=1, nedos
    emesh(ii)=emin+(ii-1)*de
  enddo
  !
  allocate(dos(1:nedos, 2))
  dos(:,:)=0.d0
  if (lpdos) then
    allocate(pdos(1:nedos, 1:nproj))
    pdos(:,:)=0.d0
  endif
  !
  nkpt=nkx*nky*nkz
  !
  allocate(kvec(1:3, 1:nkpt))
  !
  do ikx=1, nkx
    do iky=1, nky
      do ikz=1, nkz
        ii=(ikz-1)*nkx*nky+(iky-1)*nkx+ikx
        kvec(1, ii)=(ikx-1.d0)/nkx
        kvec(2, ii)=(iky-1.d0)/nky
        kvec(3, ii)=(ikz-1.d0)/nkz
      enddo
    enddo
  enddo
  !
  allocate(eig(1:nproj))
  allocate(work(1:nproj, 1:nproj))
  !
END SUBROUTINE

SUBROUTINE output_dos
  !
  use constants, only : stdout, fout
  use input,     only : lpdos
  !
  implicit none
  !
  integer ii, ip
  !
  if(lpdos) open(unit=fout, file="pdos.dat")
  !
  do ii=1, nedos
    write(stdout, '(1F9.3,2F18.6)') emesh(ii), dos(ii, :)
    if (lpdos) then
      write(fout, '(1F9.3,50F18.6)') emesh(ii), pdos(ii, :)
    endif
  enddo
  !
  if(lpdos) close(unit=fout)
  !
END SUBROUTINE

SUBROUTINE finalize_dos
  !
  implicit none
  !
  if(allocated(emesh)) deallocate(emesh)
  if(allocated(dos)) deallocate(dos)
  if(allocated(pdos)) deallocate(pdos)
  if(allocated(kvec)) deallocate(kvec)
  if(allocated(eig)) deallocate(eig)
  if(allocated(work)) deallocate(work)
  !
END SUBROUTINE

END MODULE

