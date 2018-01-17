	program readwrite

	implicit none
	integer m,n,i,j,k,numit
	parameter (m = 600)
	parameter (n = 450)
	parameter (numit = 1)

	
	real im(0:m+1,0:n+1),old(0:m+1,0:n+1),new(0:m+1,0:n+1)
	real buf(m,n)
	
	do j=0,n+1
	  do i=0,m+1
	    im(i,j)=0.0
	    old(i,j)=0.0
	  end do
	end do

	call datread('input_600x450.dat',buf,m,n)
	do j=1,n
	  do i=1,m
	    im(i,j)=buf(i,j)
	    old(i,j)=buf(i,j)
	  end do
	end do

	do k=1,numit
	  do j=1,n
	    do i=1,m
              new(i,j)=-(old(i-1,j)+old(i+1,j)+old(i,j-1)+old(i,j+1)-4*im(i,j))
	    end do
	  end do

	  do j=1,n
	    do i=1,m
	      old(i,j)=new(i,j)
	    end do
	  end do
	end do
	
	do j=1,n
	  do i=1,m
	    buf(i,j)=old(i,j)
	  end do
	end do
	
	call dump('image.pgm',buf,m,n)
	
	end

