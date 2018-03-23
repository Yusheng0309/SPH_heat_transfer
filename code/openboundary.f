      subroutine inflow(x, vx, volume, density, p, Tem,
     &              itype, hsml, ntotal, itimestep, storageip)

c----------------------------------------------------------------------     
c     This subroutine is used to generate boundary condition for the 
c     open flow field by using inflow condition 
c     x     :coordinates of particles                             [out]
c     vx    :velocities of particles                              [out]
c     mass  :mass of particles                                    [out]
c     p     :pressure  of particles                               [out]
c     Tem   :Temperature of particles                             [out]
c     u     :internal energy of particles                         [out]
c     itype :types of particles                                   [out]
c             =2   water
c     h     :smoothing lengths of particles                       [out]
c     ntotal:total particle number                                [out]
c     ninstep :how many times putting inflow particles            [out]
c     nstorage :the total number of storage particle               [in]
c     storageip :the number of the storage particles               [in]

      implicit none     
      include 'param.inc.txt'
      
      integer itype(maxn), ntotal, itimestep, ninstep
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn),
     &     density(maxn), p(maxn), hsml(maxn), Tem(maxn),
     &     storageip(maxn)
      integer i, j, d, k,inputn
      double precision dx,space_x, nstorage
    
c     setting inflow boundary condition and putting particles
c     to the flow domain

    
      dx = 100.
      space_x = 100.
      inputn = 1
c      if the storage particles are less than 40 (inflow particles
c      are put one time) , then putting new particles to the flow       
	  if(nstorage .le. inputn) then
	    ninstep = ninstep + 1 
	      do i = 1, inputn
             x(1, ntotal+i) = - (ninstep-1)*space_x
	       vx(1, ntotal+i) = 1.0
	       density(ntotal + i) = 0.5
	       volume(ntotal + i) = space_x * density(ntotal + i)
	       p(ntotal + i) = 0.
	       Tem(ntotal + i) = 1.0
	       itype(ntotal + i) = 2
	       hsml(ntotal + i) = space_x * 1.2 
           enddo
         ntotal =  ntotal + inputn

c      if the storage particles are more than 40 ,
c      then putting the particles to the flow by using
c      the storage particles 	   
	  else
	    ninstep = ninstep + 1 
	      do i = 1, inputn
             k = storageip(i)
             x(1, k) = - (ninstep-1)*space_x
	       vx(1, k) = 1.0 
	       p(k) = 0.
	       Tem(k) = 1.0
	       itype(k) = 2
	       hsml(k) = space_x * 1.2 
            enddo 
          nstorage = nstorage - inputn
	  endif       
      end	 

      subroutine outflow(x, vx, volume, density, p, Tem, 
     &                   itype, hsml, ntotal, ip, storageip)

c----------------------------------------------------------------------     
c     This subroutine is used to generate boundary condition for the 
c     open flow field by using inflow condition 
c     x     :coordinates of particles                             [out]
c     vx    :velocities of particles                              [out]
c     mass  :mass of particles                                    [out]
c     p     :pressure  of particles                               [out]
c     Tem   :Temperature of particles                             [out]
c     u     :internal energy of particles                         [out]
c     itype :types of particles                                   [out]
c          =2   water
c     h     :smoothing lengths of particles                       [out]
c     ntotal :total particle number                               [out]
c     ip     :tne number of outflow particles                      [in]
c     nstorage :the total number of storage particles             [out]
c     storageip :the number of the storage particles              [out]
      
	implicit none     
      include 'param.inc.txt'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn),
     &     density(maxn), p(maxn), hsml(maxn), Tem(maxn),
     &     storageip(maxn)
      integer i, j, d, k, ip, count,inputn,nstorage
      double precision xl, dx,space_x

c     setting outflow boundary condition and putting particles
c     out of the flow domain

      space_x = 100.
      inputn = 1
	nstorage = nstorage + 1
	count = count +1
c       if count > maxn ? 
        x(1, ip) = 15000 + (nstorage-1)*1 
	  vx(1, ip) = 0.  
	  p(ip) = 0.
	  Tem(ip) = 0.
	  itype(ip) = 3
        storageip(count) = ip
      open(10,file="../data/nstorage.txt")
      write(10,*) 'nstorage  count'
	write(10,1010) nstorage, count

1010  format(1x, I6, 5x, I6)	  	         
      end	 
	
	