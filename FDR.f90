         subroutine FDR(I,J,FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ)
         USE VAR
         IMPLICIT NONE
         integer I,J
         REAL*8 FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ,RIJINV
         REAL*8 OMMEGA,FRXIJ,FRYIJ,FRZIJ,FDXIJ,FDYIJ,FDZIJ
         REAL*8 FDFAC,FRFAC
         REAL*8 EXIJ,EYIJ,EZIJ,VXIJ,VYIJ,VZIJ
         real R2S
        
!        CM1 = MASS(ITYPE(I))
!        CM2 = MASS(ITYPE(J))
!        MW = CM1*CM2/(CM1+CM2)
         RIJINV = 1.D0/RIJ
         OMMEGA = 1.D0-RIJ/RCUT    
     
         EXIJ = RXIJ*RIJINV
         EYIJ = RYIJ*RIJINV
         EZIJ = RZIJ*RIJINV
        
         VXIJ = VOX(I)-VOX(J)
         VYIJ = VOY(I)-VOY(J)
         VZIJ = VOZ(I)-VOZ(J)
        
         FDFAC = -gamma*OMMEGA**2*(EXIJ*VXIJ+EYIJ*VYIJ+EZIJ*VZIJ)
         FDXIJ = FDFAC*EXIJ
         FDYIJ = FDFAC*EYIJ
         FDZIJ = FDFAC*EZIJ
         
         FRFAC = sigma*OMMEGA*rfac*(2.*R2S()-1.)/dsqrt(DT)
         FRXIJ = FRFAC*EXIJ
         FRYIJ = FRFAC*EYIJ
         FRZIJ = FRFAC*EZIJ
         
         FXIJ = FXIJ + FDXIJ + FRXIJ
         FYIJ = FYIJ + FDYIJ + FRYIJ
         FZIJ = FZIJ + FDZIJ + FRZIJ
        
         RETURN
         END
         
        
         
         

           
