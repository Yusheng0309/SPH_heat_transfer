      subroutine output(x,vx,volume,density,p,u,c,itype,hsml,ntotal,nf, 
     &             Tem) 

c----------------------------------------------------------------------           
c     Subroutine for saving particle information to external disk file

c     x      :coordinates of particles                             [in]
c     vx     :velocities of particles                              [in]
c     volume :volume of particles                                  [in]
c     density:water density of particles                           [in]
c     p      :pressure  of particles                               [in]
c     u      :internal energy of particles                         [in]
c     Tem    : Temperature                                         [in]
c     c      :sound velocity of particles                          [in]
c     itype  :types of particles                                   [in]
c     hsml   :smoothing lengths of particles                       [in]
c     ntotal :total particle number                                [in]

      implicit none     
      include 'param.inc.txt'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn), 
     &       density(maxn),p(maxn),u(maxn),Tem(maxn),c(maxn),hsml(maxn),
     &       nx(dim,maxn)
      integer i, d, nf      
c ------- normalize      
	do i=1,ntotal ! U*t    !ntotal*space_x
	nx(1,i)=(x(1,i)-1.0*5000)/0.99000000E+04
	enddo 
c--------
      open(1,file="../hydraulic jump/f_xv.dat")

      write(1,*) ntotal
      do i=1,ntotal         
		write(nf+1,1001) i, (x(d,i),d=1,dim), (vx(d,i),d=1,dim), 
     &	                 volume(i), density(i), Tem(i)
      enddo
      
1001  format(1x,I6,5(2x,e14.8))   
                                        
      close(1)
	   
      end           