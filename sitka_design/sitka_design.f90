program sitka_design
use sitka_param
use jaw_seeds
use sort_kit
use an_function
use anneal
implicit none
integer, parameter :: dp=kind(0.0d0)
integer, dimension(3) :: icdum
integer, dimension(2) :: iflag
integer :: ii,ierr,idum,jj,tstep,nover,nlim
real(kind=dp), dimension(3) :: ddum
real(kind=dp) :: tfact
character(len=10) :: adate,atime
!
call date_and_time(date=adate,time=atime)
open(16,file='log'//atime//'.txt')
open(15,file='in_param.txt')
read(15,*) nfam
read(15,*) c(:)
read(15,*) tfact,tstep,nover,nlim
write(16,'(a,i4)') ' ... families to be selected ... ',nfam
write(16,'(a,3f5.0)') ' ... coefficients for utility function ... ',c(:)
close(15)
!
open(11,file='parents.txt') 
ii=0
parents: do 
    read(11,*,iostat=ierr) idum
    if(ierr<0) exit
    ii=ii+1
    if(ii>tpar) then; write(16,*) 'increase tpar'; end if
    idp(ii)=idum
end do parents
npar=ii
write(16,*) ' ... parents ... ',npar
close(11)
open(12,file='crosses.txt') 
nfs(:)=0
jj=0
crosses: do 
    read(12,*,iostat=ierr) icdum(:)
    if(ierr<0) exit
    jj=jj+1
    if(jj>tcro) then; write(16,*) 'increase tcro'; end if
    icr(jj,:)=icdum(:)
    if(icr(jj,2)==icr(jj,3)) write(16,*) ' ... selfing ... ',icr(jj,2)
    iflag(:)=1
    parent_check: do ii=1,npar
        if(iflag(1)==1) then
            if(idp(ii)==icr(jj,2)) then
                nfs(ii)=nfs(ii)+1
                icr(jj,2)=ii ! code as in parent list
                iflag(1)=0
                if(sum(iflag(:))==0) exit parent_check
            end if
        end if
        if(iflag(2)==1) then
            if(idp(ii)==icr(jj,3)) then
                nfs(ii)=nfs(ii)+1
                icr(jj,3)=ii ! code as in parent list
                iflag(2)=0
                if(sum(iflag(:))==0) exit parent_check
            end if
        end if
    end do parent_check
    if(sum(iflag(:))/=0) write(16,*) ' ... parent search failure ... ',icr(jj,3)
end do crosses
ncro=jj
write(16,*) ' ... crosses ... ', ncro
close(12)
    !write(16,'(2i6)') ((idp(ii),nfs(ii)), ii=1,npar)
if(ncro<nfam) then
    write(16,*) ' ... more families requested than crosses ... ', nfam, ncro
    stop
elseif(ncro==nfam) then
    write(16,*) ' ... all families required ... ', nfam, ncro
    stop
else ! anneal
    allocate(keep(ncro))
end if
! only here if require selection of families
call seed_set(20,16)
! initialise the solution
call ranperm(keep)
    !write(16,*) keep(1:nfam)
! characterise the solution
fval=an_fun(keep)
    !f_initial: do ii=1,npar
    !    if(qq(ii)==0) cycle
    !    write(16,'(4i5)') idp(ii),nn(ii),mm(ii),pp(ii)
    !end do f_initial
    !write(16,'(a,f6.0)') ' ... initial function value ... ',fval
    !write(16,'(a,4f6.0)') ' ... dot products ... ',dot_product(nn,nn),dot_product(mm,mm),dot_product(pp,pp)
! go
call an_engine(keep,fval,t_factr=tfact,t_steps=tstep,n_over=nover,n_limit=nlim)
! final
write(16,'(a,f6.0)') ' ... final function value ... ',fval
ddum(1)=dot_product(nn,nn); ddum(2)=dot_product(mm,mm); ddum(3)=dot_product(pp,pp);
write(16,'(a,4f6.0)') ' ... dot products ... ',ddum(:)
f_final_par: do ii=1,npar
    if(qq(ii)==0) cycle
    write(16,'(4i5)') idp(ii),nn(ii),mm(ii),pp(ii)
end do f_final_par
f_final_cro: do ii=1,nfam
    write(16,'(3i5)') icr(keep(ii),1),idp(icr(keep(ii),2)),idp(icr(keep(ii),3))
end do f_final_cro
end program sitka_design