!This used to be in RDCOOR but I don't want it there, here temporarily

     !--------------------------------------------

     ! Divide beads and atoms in two different part to be processed separately

     NUM_BEAD = 0
     NUM_VS   = 0
     NUM_BA  = 0

     if(MTS_CHECK .eq. 0)then
        if(virtsite .eq. 0)then

           ! From 1 to NUM_BEAD we have Bead
           ! From NUM_BEAD+1 to NUM_BA we have beads connected to atoms (not considered in MTS) 
           ! From NUM_BA + 1 to NATOMS we have atoms 

           DO I = 1,NATOMS
              NOBOUND = 0
              IF(TYPE_LABEL(I) .EQ. 2) THEN
                 DO JP = 1,NBONDS(I)
                    IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                       NOBOUND = 1                              
                    END IF
                 END DO
                 IF(NOBOUND .EQ. 0) THEN
                    NUM_BEAD = NUM_BEAD+1
                    INDEX_AB(NUM_BEAD) = I
                 END IF
              END IF
           END DO

           PP = NUM_BEAD

           DO I = 1,NATOMS
              NOBOUND = 0
              IF (TYPE_LABEL(I) .EQ. 2) THEN
                 DO JP = 1,NBONDS(I)
                    IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                       NOBOUND = 1                              
                    END IF
                 END DO
                 IF(NOBOUND .EQ. 1) THEN
                    PP = PP+1
                    INDEX_AB(PP) = I
                 END IF
              END IF
           END DO

           NUM_BA = PP

           DO I = 1,NATOMS
              IF (TYPE_LABEL(I) .EQ. 1) THEN
                 PP = PP + 1
                 INDEX_AB(PP) = I
              END IF
           END DO
        else

           ! From 1 to NUM_BEAD we have Bead
           ! From NUM_BEAD+1 to NUM_BA we have beads connected to atoms (not considered in MTS) 
           ! From NUM_BA + 1 to NUM_VS we have Virtual Site
           ! From NUM_VS+1 to NATOMS we have atoms (no VS) 

           DO I = 1,NATOMS
              NOBOUND = 0
              IF(TYPE_LABEL(I) .EQ. 2) THEN
                 DO JP = 1,NBONDS(I)
                    IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                       NOBOUND = 1                              
                    END IF
                 END DO
                 IF(NOBOUND .EQ. 0) THEN
                    NUM_BEAD = NUM_BEAD+1
                    INDEX_AB(NUM_BEAD) = I
                 END IF
              END IF
           END DO

           NUM_BA = NUM_BEAD

           DO I = 1,NATOMS
              NOBOUND = 0
              IF (TYPE_LABEL(I) .EQ. 2) THEN
                 DO JP = 1,NBONDS(I)
                    IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                       NOBOUND = 1                              
                    END IF
                 END DO
                 IF(NOBOUND .EQ. 1) THEN
                    NUM_BA = NUM_BA+1
                    INDEX_AB(NUM_BA) = I
                 END IF
              END IF
           END DO

           NUM_VS = NUM_BA

           DO I = 1,NATOMS
              IF (TYPE_LABEL(I) .EQ. 1 .and.  VIRT_VS_IND(I) .ne. 0) THEN
                 NUM_VS = NUM_VS + 1
                 INDEX_AB(NUM_VS) = I
              END IF
           END DO

           PP = NUM_VS

           DO I = 1,NATOMS
              IF (TYPE_LABEL(I) .EQ. 1  .and.  VIRT_VS_IND(I) .eq. 0) THEN
                 PP = PP + 1
                 INDEX_AB(PP) = I
              END IF
           END DO

        end if !if(virtsite .eq. 0)
     else

        if(virtsite .eq. 0)then

           ! From 1 to NUM_BEAD we have Bead
           ! From NUM_BEAD+1 to NATOMS we have atoms 

           DO I = 1,NATOMS
              IF(TYPE_LABEL(I) .EQ. 2) THEN
                 NUM_BEAD = NUM_BEAD+1
                 INDEX_AB(NUM_BEAD) = I
              END IF
           END DO

           PP = NUM_BEAD

           DO I = 1,NATOMS
              IF (TYPE_LABEL(I) .EQ. 1) THEN
                 PP = PP + 1
                 INDEX_AB(PP) = I
              END IF
           END DO

        else

           ! From 1 to NUM_BEAD we have Bead
           ! From NUM_BEAD+1 to NUM_VS we have virtual site 
           ! From NUM_VS+1 to NATOMS we have atoms 

           DO I = 1,NATOMS
              IF(TYPE_LABEL(I) .EQ. 2) THEN
                 NUM_BEAD = NUM_BEAD+1
                 INDEX_AB(NUM_BEAD) = I
              END IF
           END DO
           NUM_VS = NUM_BEAD
           DO I = 1,NATOMS
              IF (TYPE_LABEL(I) .EQ. 1 .and.  VIRT_VS_IND(I) .ne. 0) THEN
                 NUM_VS = NUM_VS + 1
                 INDEX_AB(NUM_VS) = I
              END IF
           END DO
           PP = NUM_VS
           DO I = 1,NATOMS
              IF (TYPE_LABEL(I) .EQ. 1  .and.  VIRT_VS_IND(I) .eq. 0) THEN
                 PP = PP + 1
                 INDEX_AB(PP) = I
              END IF
           END DO

        end if ! VIRTSITE
     end if ! MTS