!Code to generate 'G' vectors for the supercell
!For now, we have assumed that the reciprocal lattice vectors and the real space lattice vectors
!are known for both the primitive cell and the suoercell
program generate_Gvectors
  implicit none
  !
  real*4, allocatable :: pclatvec(:,:),sclatvec(:,:)
  integer :: m,ratio
  real*4 :: pcvol,scvol
  real*4,dimension(3) :: pc_prod, sc_prod
  !
  allocate(pclatvec(3,3),sclatvec(3,3))
  !
  open(201,file='./data/input.in')
  do m = 1,7
    read(201,*)
  end do
  do m = 1,3
    read(201,*)pclatvec(:,m)
  end do
  read(201,*)
  do m = 1,3
    read(201,*)sclatvec(:,m)
  end do
  close(201)
  !
  !Calculating volume of the primitive cell
  pc_prod = 0
  sc_prod = 0
  !
  pc_prod(1) = pclatvec(2,2)*pclatvec(3,3)-pclatvec(3,2)*pclatvec(2,3)
  pc_prod(2) = pclatvec(3,2)*pclatvec(1,3)-pclatvec(1,2)*pclatvec(3,3)
  pc_prod(3) = pclatvec(1,2)*pclatvec(2,3)-pclatvec(2,2)*pclatvec(1,3)
  pcvol = dot_product(pclatvec(:,1),pc_prod)
  !
  !Calculating volume of the supercell
  sc_prod(1) = sclatvec(2,2)*sclatvec(3,3)-sclatvec(3,2)*sclatvec(2,3)
  sc_prod(2) = sclatvec(3,2)*sclatvec(1,3)-sclatvec(1,2)*sclatvec(3,3)
  sc_prod(3) = sclatvec(1,2)*sclatvec(2,3)-sclatvec(2,2)*sclatvec(1,3)
  scvol = dot_product(sclatvec(:,1),sc_prod)
  !
  ratio = nint(scvol/pcvol)
  call ggen(ratio)
end program generate_Gvectors
!
!
subroutine ggen(ratio)
  implicit none
  !
  integer :: ratio
  real*4,allocatable :: pcrecip(:,:),screcip(:,:)
  integer :: n,m1,m2,m3,val
  real*4,dimension(3) :: recip
  !
  !
  allocate(pcrecip(3,3),screcip(3,3))
  !
  open(202,file='./data/input.in')
  read(202,*)
  read(202,*)
  read(202,*)
  do n = 1,3
    read(202,*)pcrecip(:,n)
  end do
  read(202,*)
  do n = 1,3
    read(202,*)screcip(:,n)
  end do
  close(202)
  open(203,file='./sc_Gvectors.in')
  if (ratio%2 .eq. 0) then
    val = ratio/2
  else
    val = (ratio-1)/2
  end if
  !
  do m1 = -val,val
    do m2 = -val,val
      do m3 = -val,val
        recip(1) = m1*screcip(1,1) + m2*screcip(1,2) + m3*screcip(1,3)
        recip(2) = m1*screcip(2,1) + m2*screcip(2,2) + m3*screcip(2,3)
        recip(3) = m1*screcip(3,1) + m2*screcip(3,2) + m3*screcip(3,3)
        write(203,*)recip
      end do
    end do
  end do
  close(203)
  deallocate(pcrecip)
  deallocate(screcip)
end subroutine ggen
