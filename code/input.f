      subroutine input(x, vx, volume, density, p, u, itype, hsml, 
     &	         ntotal,Tem)
      
c----------------------------------------------------------------------
c     Subroutine for loading or generating initial particle information

c     x       :coordinates of particles                            [out]
c     vx      :velocities of particles                             [out]
c     volume  :volume of particles                                 [out]
c     Tem    : Temperature                                         [out]
c     density :water density of particles                          [out]
c     p       :pressure  of particles                              [out]
c     u       :internal energy of particles                        [out]
c     itype   :types of particles                                  [out]
c     hsml    :smoothing lengths of particles                      [out]
c     ntotal  :total particle number                               [out]

      implicit none     
      include 'param.inc.txt'

      integer itype(maxn), ntotal , itimestep       
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn), 
     &                 p(maxn), u(maxn), hsml(maxn), density(maxn),
     &                 Tem(maxn), dt  
      integer i, d, im        
          
      open(1,file="../ini_xv.dat")
	open(2,file="../ini_state.dat")
      open(3,file="../ini_other.dat") 
       
      if (Hydrau_Jump) call HydrauJump(x, vx, volume, density, p, u, 
     &                      itype, hsml, ntotal,Tem)                

      do i = 1, ntotal 
		write(1,1001) i, (x(d,i), d=1,dim), (vx(d,i),d=1,dim) 
          write(2,1002) i, volume(i), density(i), Tem(i)         
          write(3,1003) i, itype(i), hsml(i)    
      enddo   
1001    format(1x, I5, 6(2x, e14.8)) 
1002    format(1x, I5, 3(2x, e14.8)) 
1003    format(1x, I5, 2x, I2, 2x, e14.8)
 
      write(*,*)'  **************************************************'
      write(*,*)'      Initial particle configuration generated   '       
      write(*,*)'      Total number of particles   ', ntotal    	
      write(*,*)'  **************************************************' 

      close(1)
      close(2) 
      close(3) 

      end              
       
       
      subroutine HydrauJump(x,vx,volume,density,p,u,itype,hsml,
     &                    	ntotal,Tem) 
c----------------------------------------------------------------------     
c     x      :coordinates of particles                            [out]
c     vx     :velocities of particles                             [out]
c     volume :volume of particles                                 [out]
c     density:water density of particles                          [out]
c     p      :pressure  of particles                              [out]
c     u      :internal energy of particles                        [out]
c     itype  :types of particles                                  [out]
c             = 1   ideal gas
c             = 2   fluid
c     hsml   :smoothing lengths of particles                      [out]
c     ntotal :total particle number                               [out]

      implicit none     
      include 'param.inc.txt'
      
      integer itype(maxn), ntotal, itimestep 
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn),
     &     density(maxn), p(maxn), u(maxn), hsml(maxn), Tem(maxn)
      integer i, d
      double precision space_x  
      


	ntotal =100 
      space_x =100.0 
      
	    
      do i=1,ntotal
		hsml(i) = space_x * 1.2
		itype(i)=2
		do d = 1, dim
			x(d,i) = 0. 
			vx(d,i) = 0.
		enddo        
      enddo                
                
       do i=1,ntotal
		x(1,i) =  space_x * (i-1) 
	    vx(1,i) = 0.5                                    
          density(i) = 0.5 
          Tem(i) = 10.0
          p(i) = 0.5 * 9.81 * density(i) * 1000    
	    volume(i) = space_x * density(i)
	enddo

      open(7,file="../ini_all.txt")

      do i = 1, ntotal 
		write(7,1007) i, (x(d,i), d=1,dim), (vx(d,i),d=1,dim),
     &                    volume(i), density(i), Tem(i)           
      enddo   
1007  format(1x, I6, 5(2x, e14.8)) 

      close(7)
	  	               
      end   