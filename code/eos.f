      subroutine p_gas(density, u, p, c)
      
c----------------------------------------------------------------------
c     Gamma law EOS: subroutine to calculate the pressure and sound  
 
c     density: water density                                        [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
          
      implicit none
      double precision density, u, p, c   
      double precision gamma 
          
c      For air (idea gas)

      gamma=1.4
      p = (gamma-1) * density * u     
      c = sqrt((gamma-1) * u) 
     
      end 

c****************************************************************************
	             
      subroutine p_art_water(density, p, c)
      
c----------------------------------------------------------------------
c     Artificial equation of state for the artificial compressibility 

c     density: water density                                        [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
c     Equation of state for artificial compressibility   

      implicit none
      double precision density, u, p, c
      double precision gamma, density0

c     Artificial EOS, Form 1 (Monaghan, 1994)
c      gamma=7.
c      density0=1000.       
c      b = 1.013e5
c      p = b*((density/density0)**gamma-1)      
c      c = 1480.

c     Artificial EOS, Form 2 (Morris, 1997)
    
	c = sqrt(9.81*density)
      p = 0.5*9.81*density*density

      end 