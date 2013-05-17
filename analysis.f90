subroutine analysis ()

use var

implicit none

integer :: i,j,l,k,inA,inF,q,vsite,endVS,begVS

OPEN ( unit=213 , FILE = 'analysis.conf',status='replace', form='UNFORMATTED', access='SEQUENTIAL')

write(213)natoms

write(213)ntype
do i=1,ntype
    write(213)label(i),name_label(i)
end do

if(ibrdescr .EQ. 0)then
    write(213)0
    if(VIRTSITE .EQ. 0)then
        write(213)0
    else
        write(213)1  
    end if  
else
    write(213)1  
end if

write(213)nmol


!name_mol='PS_AT'

l=0
vsite=1
begVS=1
endVS=0

if(ibrdescr .eq. 100) then !FIX THIS
    write(213)nvirta
!    write(213)name_label
    if(VIRTSITE .EQ. 0)then
        do i=1,nmol
            write(213)natm(i),name_mol(i)
            do j=1,natm(i)
                l=l+1
                write(213)itype(l),mass(itype(l)),nbonds(l),type_label(l)
                write(213) (JBOND(L,K),K=1,NBONDS(L))
            end do
            write(213)virtNmol(i)
            endVS = endVS+virtNmol(i)
            do k=begVS,endVS
                write(213)init_numbcomp(k),VITYPE(k),TOTBMASS(k)
                do q=1,init_numbcomp(k)
                    write(213)COMPCOM(k,q)
                end do
            end do
            begVS=begVS+virtNmol(i)
        end do
    else
        do i=1,nmol
            write(213)natm(i),name_mol(i)
            do j=1,natm(i)
                l=l+1
                write(213)itype(l),mass(itype(l)),nbonds(l),type_label(l)
                write(213) (JBOND(L,K),K=1,NBONDS(L))
            end do
            write(213)virtNmol(i)
            endVS = endVS+virtNmol(i)
            do k=begVS,endVS
                write(213)INDEX_VSITE(k),init_numbcomp(k),VITYPE(k),TOTBMASS(k)
                do q=1,init_numbcomp(k)
                    write(213)INDX_ATM(k,q)
                end do
            end do
            begVS=begVS+virtNmol(i)
        end do
    end if
else
    do i=1,nmol
        write(213)natm(i),name_mol
        do j=1,natm(i)

        end do
    end do
end if

 close(213) 
end subroutine
