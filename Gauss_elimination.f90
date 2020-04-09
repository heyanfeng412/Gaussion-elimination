Program Gaussian_elimination
implicit none

Integer::n=4 !the order of system of linear equations
Integer i,j,k !i is the row number, j is the colume number, k is the step number
Real(8),Allocatable :: a(:,:),b(:),x(:),c(:)
Real(8) d

Allocate(a(1:n,1:n),b(1:n),x(1:n),c(1:n))

!create two files to save the value of a and b
Open(unit=10,FILE='a_value',FORM='formatted',status='old')
Open(unit=11,FILE='b_value',FORM='formatted',status='old')

do i=1,n,1
	do j=1,n,1
	read(10,*) a(i,j)
	end do
end do

do i=1,n,1
	read(11,*) b(i)
end do

close(10)
close(11)

do k=1,n-1,1  !n=4 when k=1

	do i=k+1,n,1 !assume i=2
		if(a(k,k) /= 0) then
			b(i)=b(i)-a(i,k)/a(k,k)*b(k)
		else
			goto 100
		endif !a(i,k) a(2,1) will be changed when j=1
		d=a(i,k)
		do j=1,n,1
			a(i,j)=a(i,j)-a(k,j)*(d/a(k,k))
		enddo
	enddo
enddo

!until there, we have finished step 1-3 if n=4
!calculate the value of x now
do i=n,1,-1
	do j=1,n,1
		if(j /= i ) then
			c(i)=c(i)+a(i,j)*x(j)
		else
			cycle
		endif
	enddo
	x(i)=(b(i)-c(i))/a(i,i)
enddo

!write the result of x(i)
do i=1,n,1
	write(*,*) "x(",i,"):",x(i)
enddo
		   
100 stop

end


