      subroutine virt_part(itimestep, ntotal,nvirt,hsml,volume,x,vx,
     &           density,u,p,itype,Tem) 

c----------------------------------------------------------------------
c     Subroutine to determine the information of virtual particles
c     Here only the Monaghan type virtual particles 

c     itimestep : Current time step                                 [in]
c     ntotal : Number of particles                                  [in]
c     nvirt  : Number of virtual particles                         [out]
c     hsml   : Smoothing Length                                 [in/out]
c     volume : Particle volume                                  [in/out]
c     x      : Coordinates of all particles                     [in/out]
c     vx     : Velocities of all particles                      [in/out]
c     density: water density                                    [in/out]
c     u      : internal energy                                  [in/out]
c     itype  : type of particles                                [in/out]

      implicit none
      include 'param.inc.txt'
      integer itimestep, ntotal, nvirt, itype(maxn)
      double precision hsml(maxn),volume(maxn),x(dim,maxn),vx(dim,maxn),
     &                 density(maxn), u(maxn), p(maxn), Tem(maxn)
      integer i, j, d, im, mp
      double precision xl
	double precision space_x  ! l_j  
       
	nvirt = 0
	space_x = 100.0 

	do i = 1,53
	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = - space_x 
	  vx(1, ntotal + nvirt) = 0.
	  density (ntotal + nvirt) =  i * space_x * 0.001
	  volume(ntotal + nvirt) = space_x * density(i)
	  p(ntotal + nvirt) = 0.
 	  Tem(ntotal + nvirt) = 0.
	  u(ntotal + nvirt) = 0.
	  itype(ntotal + nvirt) = -2
	  hsml(ntotal + nvirt) = space_x * 1.2
	enddo

	do i = 1,53
	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = - 2*space_x   
	  vx(1, ntotal + nvirt) = 0.
	  density (ntotal + nvirt) =  -0.0025 + i*space_x * 0.001
	  volume(ntotal + nvirt) = space_x * density(i)
	  p(ntotal + nvirt) = 0.
	  Tem(ntotal + nvirt) = 0.
	  u(ntotal + nvirt) = 0.
	  itype(ntotal + nvirt) = -2
	  hsml(ntotal + nvirt) = space_x * 1.2
	enddo


      if (mod(itimestep,save_step).eq.0) then
        open(1,file="../xv_vp.dat")
        open(2,file="../state_vp.dat")
        open(3,file="../other_vp.dat")            
       
	  write(1,*) nvirt
        
	  do i = ntotal + 1, ntotal + nvirt         
          write(1,1001) i, (x(d,i), d=1,dim), (vx(d,i), d=1,dim)              
          write(2,1002) i, volume(i), density(i), p(i), u(i), Tem(i)
          write(3,1003) i, itype(i), hsml(i)                               
        enddo 
	        
1001    format(1x, I6, 6(2x, e14.8))
1002    format(1x, I6, 8(2x, e14.8)) 
1003    format(1x, I6, 2x, I4, 2x, e14.8)
  
        close(1)
        close(2) 
        close(3) 
      endif 

      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
         print *,' >> Statistics: Virtual boundary particles:'
         print *,'          Number of virtual particles:', nvirt
        endif     
      endif

      end