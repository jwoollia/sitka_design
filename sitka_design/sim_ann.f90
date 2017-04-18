!     Last change:  JAW  21 Feb 2001   12:36 pm
!<<<<<<<<<<<<<<<<<
module an_function
!<<<<<<<<<<<<<<<<<
USE sitka_param
IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)
save
contains
!
!========================
subroutine transition(tk)
!========================
integer, dimension(:), intent(inout) :: tk
integer, dimension(2) :: iswap
integer :: idum
real(kind=dp), dimension(2) :: u
call random_number(u(1:2))
iswap(1)=ceiling(real(nfam)*u(1))
iswap(2)=ceiling(real(ncro-nfam)*u(2))
idum=tk(iswap(1))
tk(iswap(1))=tk(iswap(2)+nfam)
tk(iswap(2)+nfam)=idum
end subroutine transition
!
!==================
function an_fun(tk)
!==================
REAL(KIND=dp) :: an_fun
integer, dimension(:), intent(inout) :: tk
integer :: i
nn(:)=0; mm(:)=0; pp(:)=0; qq(:)=0
do i=1,nfam
    nn(icr(tk(i),2))=nn(icr(tk(i),2))+1 ! increment the parent's total count from the cross
    nn(icr(tk(i),3))=nn(icr(tk(i),3))+1
    mm(icr(tk(i),2))=mm(icr(tk(i),2))+1 ! increment the parent's maternal excess from the cross
    mm(icr(tk(i),3))=mm(icr(tk(i),3))-1
    pp(icr(tk(i),2))=pp(icr(tk(i),2))+2*icr(tk(i),1)-3 ! increment the parent's site 2 excess from the cross
    pp(icr(tk(i),3))=pp(icr(tk(i),3))+2*icr(tk(i),1)-3
end do
where(nn>0)
    qq=1
    nn=nn-2
end where
an_fun=c(1)*dot_product(nn,nn)+c(2)*dot_product(mm,mm)+c(3)*dot_product(pp,pp)
end function an_fun
!
!>>>>>>>>>>>>>>>>>>>>>
END module an_function
!>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<
module anneal
!<<<<<<<<<<<<
USE sitka_param
USE an_function
INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)
save
contains
!
!===========================================================================
subroutine an_engine(best,fbest,t_ini,t_factr,t_steps,n_over,n_limit,r_unit)
!===========================================================================
!
! b villanueva 09/00 ... modified jaw
! annealing parameters are
! tini   = initial temperature
! temp   = temperature
! tsteps = number of temperature steps
! tfactr = factor by which temp is reduced on each step
! nover  = maximum number of modifications tried at any temperature
! nlimit = maximum number of successful changes before continuing
! rfile  = report file
!
implicit none
INTEGER, DIMENSION(:), INTENT(INOUT) :: best
REAL(KIND=dp), INTENT(OUT) :: fbest
INTEGER, INTENT(IN), OPTIONAL :: t_steps,n_over,n_limit,r_unit
REAL(KIND=dp), INTENT(IN), OPTIONAL :: t_ini,t_factr
INTEGER :: tsteps,nover,nlimit,nsucc,ic,jc,runit,rflag,msucc
REAL(KIND=dp) :: tini,tfactr,del,score,temp,u
INTEGER, DIMENSION(SIZE(best)) :: trial,rxval
REAL(KIND=dp) :: rfval
!
IF(PRESENT(t_ini)) THEN; tini=t_ini; ELSE; tini=1_dp; END if
IF(PRESENT(t_factr)) THEN; tfactr=t_factr; ELSE; tfactr=0.9_dp; END if
IF(PRESENT(t_steps)) THEN; tsteps=t_steps; ELSE; tsteps=200; END if
IF(PRESENT(n_over)) THEN; nover=n_over; ELSE; nover=100; END if
IF(PRESENT(n_limit)) THEN; nlimit=n_limit; ELSE; nlimit=20; END if
IF(PRESENT(r_unit)) THEN; runit=r_unit; END if
!
if(present(r_unit)) then
    write(runit,'(a,f12.6)') ' ... ANNEAL initial temperature        ... ',tini
    write(runit,'(a,f12.6)') ' ... ANNEAL temperature decay rate     ... ',tfactr
    write(runit,'(a,i12)')   ' ... ANNEAL decay steps                ... ',tsteps
    write(runit,'(a,i12)')   ' ... ANNEAL maximum trials per step    ... ',nover
    write(runit,'(a,i12)')   ' ... ANNEAL maximum successes per step ... ',nover
end if
!
temp=tini
fbest=an_fun(best)
rflag=nlimit
msucc=nlimit
across_t: do ic=1,tsteps      ! try up to tsteps temperature steps
    nsucc=0
    within_t: do jc=1,nover
        trial=best
        call transition(trial)
        score=an_fun(trial)
        del=score-fbest
        call random_number(u)
        if((del<=0_dp).or.(u<exp(-del/temp)))then
            nsucc=nsucc+1
            if(fbest>score) then
                rfval=fbest
                rxval=best
            end if
            fbest=score
            best=trial
        end if
        if(nsucc>=nlimit) exit within_t
    end do within_t
    if(present(r_unit)) then
        if(nsucc<msucc) msucc=nsucc
        if(nsucc<rflag) then
            write(runit,'(a,i8,a,i8)') ' ... ANNEAL successful modifications first < ',rflag,' after step ',ic
            rflag=rflag/2
        end if
    end if
    if(nsucc==0) exit across_t
    temp=temp*tfactr
end do across_t
if(fbest>rfval) then
    fbest=rfval
    best=rxval
end if
if(present(r_unit)) then
    write(runit,'(a,i8)') ' ... ANNEAL completed at step ',min(ic,tsteps)  
    write(runit,'(a,i8)') ' ... ANNEAL minimum successful modifications per step ',msucc
end if
end subroutine an_engine
!
!>>>>>>>>>>>>>>>>
END module anneal
!>>>>>>>>>>>>>>>>
