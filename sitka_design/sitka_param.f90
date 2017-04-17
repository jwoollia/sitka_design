module sitka_param
implicit none
integer, parameter, private :: dp=kind(0.0d0)
integer, parameter :: tcro=110, tpar=65
integer, dimension(tpar) :: idp, nfs
integer, dimension(tpar) :: nn, mm, pp, qq ! nn, modified total count; mm, maternal excess; pp, site 2 excess; qq has progeny
integer, dimension(tcro,3) :: icr
integer, dimension(:), allocatable :: keep ! list of families with first nfam selected, dimension ncro
integer :: npar,ncro,nfam ! npar parents producing ncro crosses, with nfam crosses selected
real(kind=dp), dimension(3) :: c ! coefficients
real(kind=dp) :: fval

save
end module sitka_param
