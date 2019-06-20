! Angelica de Oliveira-Costa & Max Tegmark 2007
! - Please reference the following paper : http://adsabs.harvard.edu/abs/2008MNRAS.388..247D


	program gsm
	implicit  none
	integer ncomp, nside, npix
   	parameter(ncomp=3,nside=512,npix=12*512**2)
   	! number of pixels = 12*512^2 = 3145728 = number of lines in in/out files
   	! 12 - is base HealPIX N_phi=4 , N_theta=3 times 512^2 ( N_side=512 - every Quadrilateral side is divided to 512 parts )
	integer ncomploaded, i, j, lnblnk, idx, n_spaces
	real*8 nu, f(ncomp+1), norm, A(npix,ncomp), t, nu_min, nu_max,d_nu
	character*60 infile, outfile, outfile_basename, nu_str
	character*60 outdir,comline
	
	! default place where output maps will be saved :
	outdir = 'skymaps/'	        

	print *,'Input frequency range at which to make a map [in MHz]'
	print *,'Lower limit :'
	read *,nu_min
	print *,'Upper limit :'
	read *,nu_max
	print *,'Enter the step [in MHz] :'
	read *,d_nu
	print *,'Base file name of file in which to save the map ?'
	read *,outfile_basename
	print *,'Enter name of output directory [default = ',outdir(1:lnblnk(outdir)),']'
	read *,outdir

        print *,'-------------------------------------------------------'
        print *,'PARAMETERS :'
        print *,'-------------------------------------------------------'                        
	print *,'Making sky map at frequeny range    ',nu_min,'-',nu_max,' [MHz]';
	print *,'Frequency step                      ',d_nu,' [MHz]';
	        
	idx = lnblnk(outfile_basename)
	print *,'Outfile basename____________________',outfile_basename(1:idx)
	print *,'Outdir______________________________',outdir(1:lnblnk(outdir))
	print *,'-------------------------------------------------------'

        ! creating output directory :
	print *,'Creating output directory ',outdir(1:lnblnk(outdir))
	comline = 'mkdir -p ' // outdir
	if (system(comline).ne.0)then
	   stop 'DEATH ERROR : Could not create output directory !'
	endif
	

	call LoadComponents(ncomploaded)
	print *,'Loaded ',ncomploaded,' components (expected ',ncomp,')'
	if (ncomploaded.ne.ncomp) then 
		stop 'DEATH ERROR: WRONG NUMBER OF COMPONENTS LOADED'
	endif

! loading components in every sky pixel :
	infile = 'component_maps_408locked.dat'
	print *,'Loading ',infile(1:lnblnk(infile))
	open(2,file=infile,status='old')
	do i=1,npix
	  		read(2,*) (A(i,j),j=1,ncomp)
	end do	
	close(2)
	print *,'All components read from all ',npix,' pixels'

	
	! main loop over frequencies :
	nu=nu_min
!	while ( nu .le. nu_max) do		   
!	do 1001 nu=nu_min, nu_max , d_nu
	do nu=nu_min, nu_max , d_nu
		print *,'########################### nu = ',nu,' ###########################'

		call ComputeComponents(nu,ncomp,f)
		norm = f(ncomp+1)

! Moved before loop :
! loading components in every sky pixel :
!		infile = 'component_maps_408locked.dat'
!		print *,'Loading ',infile(1:lnblnk(infile))
!		open(2,file=infile,status='old')
!		do i=1,npix
! 	  		read(2,*) (A(i,j),j=1,ncomp)
!		end do	
!		close(2)
!		print *,'All components read from all ',npix,' pixels'

		write(nu_str,'(F6.2)') nu
		print *,'nu_str=|',trim(adjustl(nu_str)),'|'
		
		n_spaces=lnblnk(nu_str)-len_trim(ADJUSTL(nu_str))
		nu_str = repeat( '0', n_spaces ) // adjustl(nu_str)
		print *,'final nu_str=',nu_str
		      
		
!		return;
      outfile=outdir(1:lnblnk(outdir)) // '/' // outfile_basename(1:lnblnk(outfile_basename)) // '_ni' // trim(adjustl(nu_str)) // '.out'
! calculating intensity at every pixel and save to output file :
! outfile containes temperatures for each HealPix pixel :
		print *,'Calculating temperatures and saving to file ',outfile(1:lnblnk(outfile))
		open(3,file=outfile)	 	 
		do i=1,npix
		  t = 0
		  do j=1,ncomp
		    t = t + f(j)*A(i,j) 
		  end do
		  t = norm*t
		  write(3,*) t
		end do
		close(3)				

!		nu = nu + d_nu
! 1001	continue
	end do
	
        return
 777    call usage
       end
	 
	subroutine LoadComponents(ncomp) ! ncomp = Number of components to load
	! Load the principal components from file and spline them for later use.
	! The "extra" component (ncomp+1) is the overall scaling - we spline its logarithm.
	implicit none
	integer  nmax, ncompmax, n, ncomp, idx
	parameter(nmax=1000,ncompmax=11)
	real*8 x(nmax), y(nmax,ncompmax+1), ypp(nmax,ncompmax+1)
	common/PCA/x, y, ypp, n
	integer  i, lnblnk
	real*8   xn, scaling, tmp(nmax), yp0, yp1
	character*80 infile, comline, comline2, outline
	!
	infile = 'components.dat'
	! Count number of columns in the infile:'
	comline = 'head -1 '//infile(1:lnblnk(infile))//'|wc|cut -c9-16 >qaz_cols.dat'
!	comline = comline2 // '|wc|cut -c9-16 >qaz_cols.dat'
!	comline = comline2 // '|wc'
	print *,'DEBUG = ',comline
!	comline = 'head -1 '//infile(1:lnblnk(infile))//' | wc | cut -c9-16 >qaz_cols.dat'
	!print *,'###'//comline(1:lnblnk(comline))//'###'
	if (system(comline).ne.0) stop 'DEATH ERROR COUNTING COLUMNS'
	open (2,file='qaz_cols.dat',status='old',err=777)
	read (2,*,end=777,err=777) n
	close(2)
	ncomp = n - 2
	if (ncomp.lt.0       ) stop 'DEATH ERROR: TOO FEW  COMPONENTS.'
	if (ncomp.gt.ncompmax) stop 'DEATH ERROR: TOO MANY COMPONENTS.'
	n = 0
	open (2,file=infile,status='old')
555	read (2,*,end=666) xn, scaling, (tmp(i),i=1,ncomp)
	n = n + 1
!	print *,'DEBUG : n = ',n
	IF (n.gt.nmax) THEN 
		pause 'n>nmax in LoadVector'
	ENDIF
	x(n) = log(xn) ! We'll spline against lg(nu)
	do i=1,ncomp
	   y(n,i) = tmp(i)
	end do
	y(n,ncomp+1) = log(scaling)
	goto 555
666	close(2)
	idx = lnblnk(infile)
!	outline = ' components read from '//infile(1:idx)//' with '//n//' spline points'
!	print *,ncomp,outline
	print *,ncomp,' components read from ',infile(1:idx),' with ',n,' spline points'
!	print *,ncomp,' components read from ',infile(1:idx),' with ',n,' spline points'
	yp0 = 1.d30 ! Imposes y''=0 at starting point
	yp1 = 1.d30 ! Imposes y''=0 at endpoint
	do i=1,ncomp+1
	   call myspline_r8(x,y(1,i),n,yp0,yp1,ypp(1,i))
	end do 
	return
777	stop 'DEATH ERROR 2 COUNTING COLUMNS'	
	end 
	   
	subroutine ComputeComponents(nu,ncomp,a) ! Computes the principal components at frequency nu
	implicit none
	integer nmax, ncompmax, n, ncomp
	parameter(nmax=1000,ncompmax=11)
	real*8 x(nmax), y(nmax,ncompmax+1), ypp(nmax,ncompmax+1)
	common/PCA/x, y, ypp, n
	integer i
	real*8 a(ncompmax+1), nu, lnnu, scaling
	lnnu = log(nu)
	do i=1,ncomp+1
	  call mysplint_r8(x,y(1,i),ypp(1,i),n,lnnu,a(i))
	end do
	a(ncomp+1) = exp(a(ncomp+1)) ! The overall scaling factor
	return
	end
 	
      SUBROUTINE myspline_r8 (x,y,n,yp1,ypn,y2)
      ! From numerical recipes.
      INTEGER   n, NMAX, i, k
      REAL*8    yp1,ypn,x(n),y(n),y2(n)
      PARAMETER(NMAX=10000)	! Increased from 500 by Max
      REAL*8    p,qn,sig,un,u(NMAX)
      if (N.gt.NMAX)    pause 'SPLINE NMAX DEATH ERROR'          ! Added by Max
      if (x(1).gt.x(n)) pause 'SPLINE WARNING: x NOT INCREASING' ! Added by Max
      if (yp1.gt..99e30) then
        y2(1)=0.d0
        u (1)=0.d0
      else
        y2(1)=-0.5d0
        u (1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig  =(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p    = sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u (i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE MYSPLINT_r8 (XA,YA,Y2A,N,X,Y)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! From numerical recipes.
      ! Modified to be more robust when
      ! extrapolating - It is linear if X lies outside 
      ! the bounds.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer  N, KLO, KHI, K
      REAL*8   XA(N),YA(N),Y2A(N), X, Y, H, A, B
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.d0) PAUSE 'Bad XA input.'
      if ((x-xa(1))*(x-xa(N)).gt.0.d0) then
        ! Outside bounds; do LINEAR extrapolation rather than cubic
        A = (YA(KHI)-YA(KLO))/H
        Y = ya(KLO) + A * (x-xa(KLO))
      else
        ! Within bounds; do cubic interpolation
        A=(XA(KHI)-X)/H
        B=(X-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      end if
      RETURN
      END
