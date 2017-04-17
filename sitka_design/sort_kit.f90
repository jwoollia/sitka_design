!     Last change:  JAW  18 Jan 2009    6:11 pm
!<<<<<<<<<<<<<<
module sort_kit
!<<<<<<<<<<<<<<
implicit none
save
!
INTEGER, PARAMETER, private :: dp=KIND(1.0d0)
TYPE sort_tree
  INTEGER :: node
  REAL(KIND=dp), POINTER :: item
  TYPE(sort_tree), POINTER :: left
  TYPE(sort_tree), POINTER :: right
END TYPE sort_tree
!
contains
!
!=======================
subroutine ranperm(perm)
!=======================
INTEGER, DIMENSION(:), INTENT(OUT) :: perm
REAL(KIND=dp), DIMENSION(SIZE(perm)) :: xdum
call random_number(xdum)
call rank_dp(xdum,perm)
end subroutine ranperm
!
!=====================================
subroutine rank_dp(array,rank_to_item)
!=====================================
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: array
INTEGER, DIMENSION(:), INTENT(OUT), TARGET :: rank_to_item
INTEGER, DIMENSION(:), pointer :: prank
TYPE(sort_tree), pointer :: ranking
INTEGER :: i
IF(SIZE(array)/=SIZE(rank_to_item)) then
  PRINT *, 'ranking error'
END if
do i=1,SIZE(array)
  call grow(array(i),ranking,i)
END do
i=0
prank=>rank_to_item
call flower(ranking,prank,i)
IF(i/=SIZE(array)) then
  PRINT *, 'not enough flowers in the forest'
END if
END subroutine rank_dp
!
!======================================
RECURSIVE subroutine grow(xi,ranking,i)
!======================================
REAL(KIND=dp), target :: xi
TYPE (sort_tree), POINTER :: ranking
INTEGER :: i
IF(.not.ASSOCIATED(ranking)) THEN ! need to allocate
  ALLOCATE (ranking)
  ranking%item=>xi
  ranking%node=i
  NULLIFY (ranking%left)
  NULLIFY (ranking%right)
elseif(xi<=ranking%item) then
  call grow(xi,ranking%left,i)
else
  call grow(xi,ranking%right,i)
END if
end subroutine grow
!
!==============================================
RECURSIVE subroutine flower(ranking,rank_vec,j)
!==============================================
TYPE(sort_tree), POINTER :: ranking
INTEGER, DIMENSION(:), POINTER :: rank_vec
INTEGER :: j
IF(associated(ranking)) then
  call flower(ranking%right,rank_vec,j) ! going to the right first means highest has rank 1
  IF(associated(ranking%right)) DEALLOCATE(ranking%right)
  j=j+1
  rank_vec(j)=ranking%node
  call flower(ranking%left,rank_vec,j)
  IF(associated(ranking%left)) DEALLOCATE(ranking%left)
  DEALLOCATE(ranking)
END if
end subroutine flower
!
!=====================
subroutine test_rank() ! because I always forget how it works
!=====================
REAL(KIND=dp), DIMENSION(5) :: xx
INTEGER, DIMENSION(5) :: jxx
INTEGER :: j
xx=(/1.3d0,2.6d0,0.3d0,-1.2d0,1.5d0/)
CALL rank_dp(xx,jxx)
PRINT '(f8.4,2I4,f8.4)', (xx(j),jxx(j),j,xx(jxx(j)),j=1,5)
end subroutine test_rank
!
!=====================
subroutine test_perm() ! because I always forget how it works
!=====================
INTEGER, DIMENSION(5) :: jxx
CALL ranperm(jxx)
PRINT '(5i4)', jxx
end subroutine test_perm
!
!>>>>>>>>>>>>>>>>>>
END module sort_kit
!>>>>>>>>>>>>>>>>>>
