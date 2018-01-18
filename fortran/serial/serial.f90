	program readwrite

	implicit none
	integer m,n,i,j,k,numit,nn,mm
	parameter (m = 600)
	parameter (n = 450)
	parameter (numit = 10000)
        character (len=32) :: filename
	
	real im(0:m+1,0:n+1),old(0:m+1,0:n+1),new(0:m+1,0:n+1)
	real buf(m,n)

        im(:,:)  = 0.e0
        old(:,:) = 0.e0	

        write(filename,'(a,I3,a,I3,a)'),"edge",m,"x",n,".dat"
	call datread(filename,buf,m,n)

        im(1:m,1:n) = buf(1:m,1:n)
        old(1:m,1:n) = buf(1:m,1:n)


	do k = 1, numit
          new(1:m,1:n) = 0.25e0*(old(0:m-1,1:n) + old(2:m+1,1:n) + old(1:m,0:n-1) + old(1:m,2:n+1) - im(1:m,1:n))
	  old(:,:)=new(:,:)
          if ( mod(k,1000) == 0 ) then
          write(*,*) 'Iteration # ',k
          end if
	end do
	
	buf(1:m,1:n)=old(1:m,1:n)
	
	call pgmwrite('image.pgm',buf,m,n)
	
	end

