      subroutine viscosity(ntotal,itype,x,density,eta)

c----------------------------------------------------------------------
c     Subroutine to define the fluid particle viscosity
 
c     ntotal  : Number of particles                                 [in]
c     itype   : Type of particle                                    [in]
c     x       : Coordinates of all particles                        [in]
c     density : water density                                       [in]
c     eta     : Dynamic viscosity                                  [out]

      implicit none
      include 'param.inc.txt'
      
      integer ntotal,i,itype(maxn)
      double precision x(dim,maxn),density(maxn),eta(maxn)

      do i=1,ntotal
        if (abs(itype(i)).eq.1) then
          eta(i)=0.
        else if (abs(itype(i)).eq.2) then
          eta(i)=1.0e-3
        endif  
      enddo  
 
      end